#!/usr/bin/env python3
"""Batch GNINA docking (CPU or GPU); writes standardized dock_score CSV."""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "docking_open_source.yaml"


def load_config(config_path: Path, target_key: str) -> tuple[dict, dict]:
    cfg = yaml.safe_load(config_path.read_text())
    return cfg["targets"][target_key], cfg.get("gnina", {})


def parse_gnina_log(log_text: str, score_mode: str) -> tuple[float | None, str]:
    """Parse GNINA log table. score_mode: affinity | cnnaff."""
    affinities: list[float] = []
    cnnaffs: list[float] = []
    for line in log_text.splitlines():
        # mode | affinity | intramol | CNN | CNNaff
        m = re.match(
            r"\s*\d+\s+([-0-9.]+)\s+[-0-9.]+\s+[-0-9.]+\s+([-0-9.]+)",
            line,
        )
        if m:
            affinities.append(float(m.group(1)))
            cnnaffs.append(float(m.group(2)))
        else:
            m2 = re.match(r"\s*\d+\s+([-0-9.]+)\s+", line)
            if m2:
                affinities.append(float(m2.group(1)))
    if score_mode == "cnnaff" and cnnaffs:
        # Higher CNNaffinity (pK) = better → store negative for lower-is-better convention
        return -max(cnnaffs), "docked"
    if affinities:
        return min(affinities), "docked"
    return None, "no_pose"


def dock_one(
    gnina_bin: str,
    receptor: str,
    ligand_pdbqt: str,
    center: list[float],
    size: list[float],
    out_sdf: str,
    exhaustiveness: int,
    num_modes: int,
    cpu: int,
    no_gpu: bool,
    cnn_scoring: str,
    score_mode: str,
    extra_args: list[str],
) -> tuple[float | None, str, str]:
    out_dir = Path(out_sdf).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(out_sdf).with_suffix(".log")
    cmd = [
        gnina_bin,
        "-r",
        receptor,
        "-l",
        ligand_pdbqt,
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--cpu",
        str(cpu),
        "--cnn_scoring",
        cnn_scoring,
        "-o",
        out_sdf,
        "--log",
        str(log_path),
    ]
    if no_gpu:
        cmd.append("--no_gpu")
    cmd.extend(extra_args)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        log_text = ""
        if log_path.exists():
            log_text = log_path.read_text()
        log_text += "\n" + (proc.stdout or "") + "\n" + (proc.stderr or "")
        if proc.returncode != 0 and not log_path.exists():
            return None, "gnina_error", (proc.stderr or proc.stdout or "")[:500]
        score, status = parse_gnina_log(log_text, score_mode)
        if score is None:
            return None, "no_pose", str(log_path)
        return score, status, str(log_path)
    except subprocess.TimeoutExpired:
        return None, "timeout", ""
    except FileNotFoundError:
        return None, "gnina_not_found", ""


def _worker(args_tuple):
    return dock_one(*args_tuple)


def run_batch(
    manifest_csv: Path,
    receptor_pdbqt: Path,
    center: list[float],
    size: list[float],
    output_dir: Path,
    pdb_id: str,
    gnina_cfg: dict,
    jobs: int = 1,
    limit: int | None = None,
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_csv)
    manifest = manifest[manifest["status"] == "prepared"].copy()
    if limit:
        manifest = manifest.head(limit)

    gnina_bin = gnina_cfg.get("binary", "gnina")
    bin_path = PROJECT_ROOT / "tools" / "gnina"
    if gnina_bin == "gnina" and bin_path.exists():
        gnina_bin = str(bin_path)

    score_mode = gnina_cfg.get("score_mode", "affinity")
    tasks = []
    for _, row in manifest.iterrows():
        rid = row["repurposing_id"]
        lig = row["pdbqt"]
        out_sdf = output_dir / "poses" / f"{rid}_out.sdf"
        tasks.append(
            (
                gnina_bin,
                str(receptor_pdbqt),
                lig,
                center,
                size,
                str(out_sdf),
                int(gnina_cfg.get("exhaustiveness", 16)),
                int(gnina_cfg.get("num_modes", 9)),
                int(gnina_cfg.get("cpu", 1)),
                bool(gnina_cfg.get("no_gpu", True)),
                str(gnina_cfg.get("cnn_scoring", "rescore")),
                score_mode,
                list(gnina_cfg.get("extra_args", [])),
            )
        )

    results: list[tuple] = []
    if jobs <= 1:
        for t in tasks:
            results.append(_worker(t))
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(_worker, t): i for i, t in enumerate(tasks)}
            buf = [None] * len(tasks)
            for fut in as_completed(futs):
                buf[futs[fut]] = fut.result()
            results = buf

    records = []
    for (_, _, _, _, _, out_sdf, *_), (score, status, log_path) in zip(tasks, results):
        rid = Path(out_sdf).stem.replace("_out", "")
        mrow = manifest[manifest["repurposing_id"] == rid].iloc[0]
        records.append(
            {
                "repurposing_id": rid,
                "canonical_smiles": mrow["canonical_smiles"],
                "dock_score": score,
                "glide_score_xp": score,
                "docking_status": status,
                "pdb_id": pdb_id,
                "docking_engine": "gnina_cpu" if gnina_cfg.get("no_gpu", True) else "gnina",
                "pose_file": out_sdf if status == "docked" else None,
                "log_file": log_path,
            }
        )

    out_df = pd.DataFrame(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"docking_{pdb_id.lower()}_gnina.csv"
    out_df.to_csv(csv_path, index=False)
    summary = {
        "pdb_id": pdb_id,
        "n_docked": int((out_df["docking_status"] == "docked").sum()),
        "n_total": int(len(out_df)),
        "output_csv": str(csv_path),
        "engine": "gnina",
        "no_gpu": gnina_cfg.get("no_gpu", True),
        "exhaustiveness": gnina_cfg.get("exhaustiveness", 16),
    }
    (output_dir / "docking_summary.json").write_text(json.dumps(summary, indent=2))
    return out_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch GNINA docking (WSL CPU-friendly)")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--target", type=str, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1, help="Parallel workers (CPU: 1-4 recommended)")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    target, gnina_cfg = load_config(args.config, args.target)
    receptor = PROJECT_ROOT / target["prepared_receptor"]
    if not receptor.exists():
        raise FileNotFoundError(f"Receptor missing: {receptor}. Run prepare_receptor_vina.py first.")

    run_batch(
        args.manifest,
        receptor,
        target["center"],
        target["size"],
        args.output_dir,
        target["pdb_id"],
        gnina_cfg,
        jobs=args.jobs,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
