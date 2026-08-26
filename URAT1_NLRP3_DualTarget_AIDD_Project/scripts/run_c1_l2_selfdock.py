#!/usr/bin/env python3
"""C1 L2 driver: self-dock lesinurad@9DKB and NP3-146@7ALV across seeds 42/43/44.

Uses config/docking_c1_cpu.yaml (no_gpu). Does NOT open L3.
Readouts parsed from SDF (C1_P2star), not batch max-CNNaffinity CSV.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from parse_c1_sdf_readouts import evaluate_selfdock  # noqa: E402


def run_gnina(
    gnina: Path,
    receptor: Path,
    ligand: Path,
    center: list[float],
    size: list[float],
    out_sdf: Path,
    exhaustiveness: int,
    num_modes: int,
    cpu: int,
    seed: int,
    no_gpu: bool,
) -> None:
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina),
        "-r",
        str(receptor),
        "-l",
        str(ligand),
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
        "rescore",
        "--seed",
        str(seed),
        "-o",
        str(out_sdf),
        "--log",
        str(log),
    ]
    if no_gpu:
        cmd.append("--no_gpu")
    print("RUN:", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    (out_sdf.parent / (out_sdf.stem + "_stdout.txt")).write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    if proc.returncode != 0 and not out_sdf.exists():
        raise RuntimeError(f"gnina failed rc={proc.returncode}: {(proc.stderr or '')[:500]}")


def main() -> None:
    cfg = yaml.safe_load((PROJECT_ROOT / "config/docking_c1_cpu.yaml").read_text())
    gnina_cfg = cfg["gnina"]
    gnina = PROJECT_ROOT / "tools" / "gnina"
    if not gnina.exists():
        gnina = Path(gnina_cfg.get("binary", "gnina"))

    seeds = [42, 43, 44]
    out_root = PROJECT_ROOT / "data/campaigns/c1/02_selfdock"
    refs = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
    arg_json = refs / "arg477_coords.json"
    les_ref = refs / "lesinurad_crystal_ref.sdf"
    rm5_ref = refs / "NP3-146_RM5_crystal_ref.sdf"

    # ligand pdbqts
    fr_manifest = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/forced_recovery/ligand_manifest.csv"
    np3_manifest = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/np3_146/ligand_manifest.csv"
    fr = pd.read_csv(fr_manifest)
    les_row = fr[fr["repurposing_id"] == "lesinurad"].iloc[0]
    np3 = pd.read_csv(np3_manifest)
    np3_row = np3.iloc[0]

    jobs = [
        {
            "ligand_id": "lesinurad",
            "target": "urat1_9dkb",
            "pdbqt": Path(les_row["pdbqt"]),
            "ref": les_ref,
            "use_arg": True,
        },
        {
            "ligand_id": "NP3-146",
            "target": "nlrp3_7alv",
            "pdbqt": Path(np3_row["pdbqt"]),
            "ref": rm5_ref,
            "use_arg": False,
        },
    ]

    results = []
    for job in jobs:
        tcfg = cfg["targets"][job["target"]]
        receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
        for seed in seeds:
            out_sdf = out_root / job["target"] / f"seed{seed}" / f"{job['ligand_id']}_out.sdf"
            if not (out_sdf.exists() and out_sdf.stat().st_size > 0):
                run_gnina(
                    gnina,
                    receptor,
                    job["pdbqt"],
                    tcfg["center"],
                    tcfg["size"],
                    out_sdf,
                    int(gnina_cfg["exhaustiveness"]),
                    int(gnina_cfg["num_modes"]),
                    int(gnina_cfg.get("cpu", 4)),
                    seed,
                    bool(gnina_cfg.get("no_gpu", True)),
                )
            else:
                print(f"SKIP existing {out_sdf}", flush=True)
            metrics = evaluate_selfdock(
                out_sdf,
                job["ref"],
                arg_json if job["use_arg"] else None,
                job["ligand_id"],
                seed,
                job["target"],
            )
            metrics.pop("_poses", None)
            (out_sdf.parent / f"{job['ligand_id']}_metrics.json").write_text(
                json.dumps(metrics, indent=2)
            )
            results.append(metrics)
            print(json.dumps(metrics, indent=2), flush=True)

    df = pd.DataFrame(results)
    df.to_csv(out_root / "l2_selfdock_metrics.csv", index=False)

    # Aggregate pass/fail
    def seed_pass(ligand, target):
        sub = [r for r in results if r["ligand_id"] == ligand and r["target"] == target]
        return {
            "n_seeds": len(sub),
            "n_pass": sum(1 for r in sub if r.get("pass") is True),
            "per_seed": {r["seed"]: r.get("pass") for r in sub},
            "rmsd": {r["seed"]: r.get("rmsd_cnnscore_selected") for r in sub},
            "acid_arg477": {r["seed"]: r.get("acid_arg477_min_A") for r in sub},
        }

    summary = {
        "lesinurad_9DKB": seed_pass("lesinurad", "urat1_9dkb"),
        "NP3_146_7ALV": seed_pass("NP3-146", "nlrp3_7alv"),
        "criteria": {
            "lesinurad": "CNNscore-selected RMSD<=2.0 A AND acid-Arg477<=4.0 A",
            "NP3-146": "CNNscore-selected RMSD<=2.0 A",
            "seeds": seeds,
        },
        "l3_blocked_until": "lesinurad passes on primary seed 42 (campaign: L2 gate)",
    }
    # Campaign: L2 pass for opening L3 — require seed 42 pass for lesinurad; NP3 failure demotes NLRP3 arm
    les42 = next(r for r in results if r["ligand_id"] == "lesinurad" and r["seed"] == 42)
    np342 = next(r for r in results if r["ligand_id"] == "NP3-146" and r["seed"] == 42)
    summary["gate_seed42"] = {
        "lesinurad_pass": bool(les42.get("pass")),
        "NP3_146_pass": bool(np342.get("pass")),
        "allow_L3_rank_track": bool(les42.get("pass")),
        "nlrp3_structure_exploratory_only": not bool(np342.get("pass")),
    }
    (out_root / "l2_gate_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
