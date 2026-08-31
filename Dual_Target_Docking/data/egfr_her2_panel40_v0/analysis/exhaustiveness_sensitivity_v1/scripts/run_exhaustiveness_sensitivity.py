#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0")
OUT = ROOT / "analysis" / "exhaustiveness_sensitivity_v1"
TABLES = OUT / "tables"
POSES = OUT / "poses"
LOGS = OUT / "logs"
SCRIPTS = OUT / "scripts"
VINA = Path("/home/gwj/miniconda3/bin/vina")
LIGAND_DIR = Path("/mnt/d/CADD paper exercise/dual target docking/Maestro doc/vina_docking/ligands_pdbqt")
RECEPTOR_DIR = ROOT / "receptors"
BOX_JSON = ROOT / "boxes" / "all_boxes.json"

SEED_FIXED_PRIMARY = 20260727
SEED_REPLICATES = [20260727, 7, 42]
EXHAUST_LIST = [8, 16, 32]
N_MODES = 9
ENERGY_RANGE = 3
TARGETS = ["3POZ", "3RCD"]
KEY_NOISE_LIGANDS = ["EH40_01", "EH40_18", "EH40_23"]


def read_subset() -> list[str]:
    with (TABLES / "subset_ligands.csv").open() as fh:
        return [row["ligand_id"] for row in csv.DictReader(fh)]


def split_models(all_pose_path: Path, pose_dir: Path) -> int:
    text = all_pose_path.read_text(errors="ignore")
    blocks = []
    current = []
    in_model = False
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if current:
                blocks.append(current)
            current = [line]
            in_model = True
        elif line.startswith("ENDMDL"):
            current.append(line)
            blocks.append(current)
            current = []
            in_model = False
        elif in_model:
            current.append(line)
    if current:
        blocks.append(current)
    if not blocks and text.strip():
        blocks = [[text]]
    pose_dir.mkdir(parents=True, exist_ok=True)
    for idx, block in enumerate(blocks, start=1):
        (pose_dir / f"mode_{idx:02d}.pdbqt").write_text("".join(block))
    return len(blocks)


def build_jobs() -> list[dict]:
    subset = read_subset()
    jobs = []
    for lig in subset:
        for target in TARGETS:
            for exhaustiveness in EXHAUST_LIST:
                jobs.append(
                    {
                        "experiment": "A",
                        "ligand_id": lig,
                        "target": target,
                        "exhaustiveness": exhaustiveness,
                        "seed": SEED_FIXED_PRIMARY,
                    }
                )
    for lig in KEY_NOISE_LIGANDS:
        for target in TARGETS:
            for exhaustiveness in [8, 16]:
                for seed in SEED_REPLICATES:
                    jobs.append(
                        {
                            "experiment": "B",
                            "ligand_id": lig,
                            "target": target,
                            "exhaustiveness": exhaustiveness,
                            "seed": seed,
                        }
                    )
    return jobs


def run_one(job: dict) -> dict:
    target = job["target"]
    lig = job["ligand_id"]
    exhaustiveness = job["exhaustiveness"]
    seed = job["seed"]
    receptor = RECEPTOR_DIR / f"{target}_receptor.pdbqt"
    ligand = LIGAND_DIR / f"{lig}.pdbqt"
    all_pose_dir = POSES / f"E{exhaustiveness}_seed{seed}" / target / lig
    pose_dir = all_pose_dir
    all_pose_dir.mkdir(parents=True, exist_ok=True)
    all_pose_path = all_pose_dir / "all_modes.pdbqt"
    log_path = LOGS / f"E{exhaustiveness}_seed{seed}" / f"{target}_{lig}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with BOX_JSON.open() as fh:
        boxes = json.load(fh)
    box = boxes[target]
    cmd = [
        str(VINA),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--out",
        str(all_pose_path),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(N_MODES),
        "--energy_range",
        str(ENERGY_RANGE),
        "--seed",
        str(seed),
        "--cpu",
        "1",
    ]
    if not VINA.exists():
        return {**job, "status": "fail", "reason": f"vina_missing:{VINA}"}
    if not receptor.exists() or not ligand.exists():
        return {
            **job,
            "status": "fail",
            "reason": f"missing_input receptor={receptor.exists()} ligand={ligand.exists()}",
        }
    with log_path.open("w") as log_fh:
        proc = subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not all_pose_path.exists():
        return {
            **job,
            "status": "fail",
            "reason": f"vina_returncode={proc.returncode}",
            "log_path": str(log_path),
            "pose_dir": str(pose_dir),
        }
    n_modes = split_models(all_pose_path, pose_dir)
    return {
        **job,
        "status": "success",
        "reason": "",
        "log_path": str(log_path),
        "pose_dir": str(pose_dir),
        "all_pose_path": str(all_pose_path),
        "n_modes_written": n_modes,
    }


def main() -> int:
    LOGS.mkdir(parents=True, exist_ok=True)
    POSES.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs()
    workers = min(8, max(1, os.cpu_count() or 1))
    results = []
    print(f"Running {len(jobs)} jobs with {workers} workers", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        fut_map = {ex.submit(run_one, job): job for job in jobs}
        for idx, fut in enumerate(as_completed(fut_map), start=1):
            result = fut.result()
            results.append(result)
            print(
                f"[{idx}/{len(jobs)}] {result['experiment']} "
                f"{result['target']} {result['ligand_id']} "
                f"E{result['exhaustiveness']} seed{result['seed']} -> {result['status']}",
                flush=True,
            )
    results.sort(key=lambda r: (r["experiment"], r["exhaustiveness"], r["seed"], r["target"], r["ligand_id"]))
    out_csv = TABLES / "job_status.csv"
    with out_csv.open("w", newline="") as fh:
        fieldnames = [
            "experiment",
            "target",
            "ligand_id",
            "exhaustiveness",
            "seed",
            "status",
            "reason",
            "log_path",
            "pose_dir",
            "all_pose_path",
            "n_modes_written",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    summary = {
        "n_jobs": len(results),
        "n_success": sum(r["status"] == "success" for r in results),
        "n_fail": sum(r["status"] != "success" for r in results),
        "workers": workers,
    }
    (TABLES / "job_status_summary.json").write_text(json.dumps(summary, indent=2))
    print("SUMMARY", json.dumps(summary, ensure_ascii=False), flush=True)
    return 0 if summary["n_fail"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
