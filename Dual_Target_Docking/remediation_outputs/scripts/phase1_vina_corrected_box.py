#!/usr/bin/env python3
"""Rerun EGFR/HER2 primary-panel Vina with corrected heavy-atom boxes.

Box is the only intended protocol change. Ligand PDBQT, receptor PDBQT,
Vina 1.2.7, exhaustiveness 8, num_modes 9, energy_range 3, and as-run seeds
are held fixed. Outputs stay under remediation_outputs/.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs"
P1 = OUT / "phase1_boxes"
DOCK = OUT / "phase1_vina"
LIGDIR = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0/ligands_pdbqt"
)
RECDIR = ROOT / "data/egfr_her2_panel40_v0/receptors"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
SEEDS40 = ROOT / "data/egfr_her2_panel40_v0/protocol/seeds_as_run.csv"
VINA = "/home/gwj/miniconda3/bin/vina"
EXHAUST = 8
N_MODES = 9
ENERGY_RANGE = 3
WORKERS = 6
EH120_SEED = 20260727


def load_boxes() -> dict:
    boxes = {}
    for pdb in ("3POZ", "3RCD"):
        boxes[pdb] = json.loads((P1 / f"{pdb}_box_corrected.json").read_text())
    return boxes


def load_seeds() -> dict[tuple[str, str], int]:
    seeds = {}
    with SEEDS40.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            seeds[(row["target"], row["ligand_id"])] = int(row["seed"])
    return seeds


def seed_for(target: str, lig: str, as_run: dict[tuple[str, str], int]) -> int:
    if lig.startswith("EH40_"):
        return as_run[(target, lig)]
    if lig.startswith("EH120_"):
        return EH120_SEED
    raise KeyError(lig)


def parse_mode1(out_pdbqt: Path) -> float | None:
    text = out_pdbqt.read_text(errors="replace")
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", text)
    if not m:
        return None
    return float(m.group(1))


def write_conf(target: str, lig: str, box: dict, seed: int, ligand: Path) -> tuple[Path, Path]:
    conf_dir = DOCK / "confs"
    out_dir = DOCK / "out"
    log_dir = DOCK / "logs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    rec = RECDIR / f"{target}_receptor.pdbqt"
    out = out_dir / f"{target}_{lig}_out.pdbqt"
    log = log_dir / f"{target}_{lig}.log"
    conf = conf_dir / f"{target}_{lig}.txt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {ligand}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                f"energy_range = {ENERGY_RANGE}",
                "cpu = 1",
                f"seed = {seed}",
                f"out = {out}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return conf, out, log, rec


def run_one(job: dict) -> dict:
    conf, out, log, rec = (
        Path(job["conf"]),
        Path(job["out"]),
        Path(job["log"]),
        Path(job["receptor"]),
    )
    if out.exists() and out.stat().st_size > 0:
        aff = parse_mode1(out)
        if aff is not None:
            job = dict(job)
            job["affinity"] = aff
            job["status"] = "cached"
            return job
    cmd = [VINA, "--config", str(conf)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write_text((proc.stdout or "") + (proc.stderr or ""), encoding="utf-8")
    job = dict(job)
    if proc.returncode != 0 or not out.exists():
        job["status"] = f"FAIL rc={proc.returncode}"
        job["affinity"] = ""
        return job
    aff = parse_mode1(out)
    job["affinity"] = aff if aff is not None else ""
    job["status"] = "ok" if aff is not None else "NO_MODE1"
    return job


def build_jobs() -> list[dict]:
    boxes = load_boxes()
    as_run = load_seeds()
    ligs = []
    with PANEL.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ligs.append(row["panel_id"])
    jobs = []
    missing = []
    for lig in ligs:
        ligand = LIGDIR / f"{lig}.pdbqt"
        if not ligand.exists():
            missing.append(str(ligand))
            continue
        for target in ("3POZ", "3RCD"):
            seed = seed_for(target, lig, as_run)
            rec = RECDIR / f"{target}_receptor.pdbqt"
            conf, out, log, rec = write_conf(target, lig, boxes[target], seed, ligand)
            jobs.append(
                {
                    "target": target,
                    "ligand": lig,
                    "seed": seed,
                    "exhaustiveness": EXHAUST,
                    "num_modes": N_MODES,
                    "energy_range": ENERGY_RANGE,
                    "ligand_pdbqt": str(ligand),
                    "receptor": str(rec),
                    "conf": str(conf),
                    "out": str(out),
                    "log": str(log),
                    "vina": VINA,
                }
            )
    if missing:
        raise SystemExit(f"missing {len(missing)} ligand pdbqt files, e.g. {missing[:3]}")
    man = DOCK / "job_manifest.csv"
    DOCK.mkdir(parents=True, exist_ok=True)
    with man.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(jobs[0].keys()))
        w.writeheader()
        w.writerows(jobs)
    print(f"jobs={len(jobs)} ligands={len(ligs)} manifest={man}")
    return jobs


def main() -> int:
    if not Path(VINA).exists():
        print("vina missing", VINA, file=sys.stderr)
        return 2
    jobs = build_jobs()
    if os.environ.get("PREP_ONLY") == "1":
        return 0
    results = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(run_one, job) for job in jobs]
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            results.append(rec)
            done += 1
            if done % 10 == 0 or rec["status"] != "ok" and rec["status"] != "cached":
                print(f"[{done}/{len(jobs)}] {rec['target']} {rec['ligand']} {rec['status']} {rec.get('affinity','')}", flush=True)
    results.sort(key=lambda r: (r["ligand"], r["target"]))
    out_csv = DOCK / "scores_vina_mode1_corrected_box.csv"
    fields = list(results[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    n_ok = sum(1 for r in results if r["status"] in {"ok", "cached"})
    print(f"done {n_ok}/{len(results)} -> {out_csv}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
