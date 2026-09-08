#!/usr/bin/env python3
"""Batch AutoDock Vina docking; writes standardized dock_score CSV."""
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
    return cfg["targets"][target_key], cfg.get("vina", {})


def parse_vina_log(log_text: str) -> tuple[float | None, str]:
    scores = []
    for line in log_text.splitlines():
        m = re.match(r"\s*(\d+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", line)
        if m:
            scores.append(float(m.group(2)))
    if not scores:
        return None, "no_pose"
    return min(scores), "docked"


def dock_one(
    vina_bin: str,
    receptor: str,
    ligand_pdbqt: str,
    center: list[float],
    size: list[float],
    out_pdbqt: str,
    exhaustiveness: int,
    num_modes: int,
    energy_range: float,
    cpu: int,
) -> tuple[float | None, str, str]:
    out_dir = Path(out_pdbqt).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(out_pdbqt).with_suffix(".log")
    cmd = [
        vina_bin,
        "--receptor",
        receptor,
        "--ligand",
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
        "--energy_range",
        str(energy_range),
        "--cpu",
        str(cpu),
        "--out",
        out_pdbqt,
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        log_path.write_text(log_text)
        if proc.returncode != 0:
            return None, "vina_error", log_text[:500]
        score, status = parse_vina_log(log_text)
        return score, status, str(log_path)
    except subprocess.TimeoutExpired:
        return None, "timeout", ""
    except FileNotFoundError:
        return None, "vina_not_found", ""


def _worker(args_tuple):
    return dock_one(*args_tuple)


def run_batch(
    manifest_csv: Path,
    receptor_pdbqt: Path,
    center: list[float],
    size: list[float],
    output_dir: Path,
    pdb_id: str,
    vina_cfg: dict,
    jobs: int = 4,
    limit: int | None = None,
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_csv)
    manifest = manifest[manifest["status"] == "prepared"].copy()
    if limit:
        manifest = manifest.head(limit)

    vina_bin = vina_cfg.get("binary", "vina")
    if vina_bin == "vina" and Path("/tmp/vina").exists():
        vina_bin = "/tmp/vina"

    tasks = []
    rows = []
    for _, row in manifest.iterrows():
        rid = row["repurposing_id"]
        lig = row["pdbqt"]
        out_pose = output_dir / "poses" / f"{rid}_out.pdbqt"
        tasks.append(
            (
                vina_bin,
                str(receptor_pdbqt),
                lig,
                center,
                size,
                str(out_pose),
                int(vina_cfg.get("exhaustiveness", 32)),
                int(vina_cfg.get("num_modes", 9)),
                float(vina_cfg.get("energy_range", 3.0)),
                int(vina_cfg.get("cpu", 1)),
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
    for (_, _, _, _, _, out_pose, *_), (score, status, log_path) in zip(tasks, results):
        rid = Path(out_pose).stem.replace("_out", "")
        mrow = manifest[manifest["repurposing_id"] == rid].iloc[0]
        records.append(
            {
                "repurposing_id": rid,
                "canonical_smiles": mrow["canonical_smiles"],
                "dock_score": score,
                "docking_status": status,
                "pdb_id": pdb_id,
                "docking_engine": "vina",
                "pose_file": out_pose if status == "docked" else None,
                "log_file": log_path,
            }
        )

    out_df = pd.DataFrame(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"docking_{pdb_id.lower()}_vina.csv"
    out_df.to_csv(csv_path, index=False)
    summary = {
        "pdb_id": pdb_id,
        "n_docked": int((out_df["docking_status"] == "docked").sum()),
        "n_total": int(len(out_df)),
        "output_csv": str(csv_path),
        "engine": "vina",
        "exhaustiveness": vina_cfg.get("exhaustiveness", 32),
    }
    (output_dir / "docking_summary.json").write_text(json.dumps(summary, indent=2))
    return out_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch Vina docking")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--target", type=str, required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="ligand_manifest.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    parser.add_argument("--limit", type=int, default=None, help="Dock first N ligands (smoke test)")
    args = parser.parse_args()

    target, vina_cfg = load_config(args.config, args.target)
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
        vina_cfg,
        jobs=args.jobs,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
