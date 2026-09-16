#!/usr/bin/env python3
"""EGFR/HER2 five-seed Vina with canonical heavy-atom boxes.

Seeds 20260811–20260814 are redocked. Seed 20260727 reuses the corrected
production mode-1 scores (as-run EH40 seeds + EH120 seed 20260727).
Outputs stay under remediation_outputs/ until analysis CSVs are written.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs/phase_fiveseed"
LIGDIR = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0/ligands_pdbqt"
)
RECDIR = ROOT / "data/egfr_her2_panel40_v0/receptors"
BOXDIR = ROOT / "data/egfr_her2_panel40_v0/boxes"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
PROD = ROOT / "remediation_outputs/phase1_vina/scores_vina_mode1_corrected_box.csv"
VINA = "/home/gwj/miniconda3/bin/vina"
PRIMARY = 20260727
NEW_SEEDS = [20260811, 20260812, 20260813, 20260814]
WORKERS = int(os.environ.get("FIVESEED_WORKERS", "6"))


def parse_mode1(path: Path) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def dock(target: str, lig: str, seed: int, box: dict) -> dict:
    rec = RECDIR / f"{target}_receptor.pdbqt"
    ligand = LIGDIR / f"{lig}.pdbqt"
    out_dir = OUT / "out" / f"seed_{seed}"
    log_dir = OUT / "logs" / f"seed_{seed}"
    conf_dir = OUT / "confs" / f"seed_{seed}"
    for d in (out_dir, log_dir, conf_dir):
        d.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{target}_{lig}_out.pdbqt"
    log = log_dir / f"{target}_{lig}.log"
    conf = conf_dir / f"{target}_{lig}.txt"
    rec_row = {
        "pair": "EGFR/HER2",
        "seed": seed,
        "ligand": lig,
        "pdb": target,
        "status": "",
        "vina_mode1": "",
    }
    if out.exists() and parse_mode1(out) is not None:
        rec_row["status"] = "cached"
        rec_row["vina_mode1"] = parse_mode1(out)
        return rec_row
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
                "exhaustiveness = 8",
                "num_modes = 9",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {seed}",
                f"out = {out}",
            ]
        )
        + "\n"
    )
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text((proc.stdout or "") + (proc.stderr or ""))
    aff = parse_mode1(out) if out.exists() else None
    rec_row["status"] = "ok" if aff is not None else f"FAIL rc={proc.returncode}"
    rec_row["vina_mode1"] = aff if aff is not None else ""
    return rec_row


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    boxes = {pdb: json.loads((BOXDIR / f"{pdb}_box.json").read_text()) for pdb in ("3POZ", "3RCD")}
    ligs = [r["panel_id"] for r in csv.DictReader(PANEL.open())]
    cls = {r["panel_id"]: r["class"] for r in csv.DictReader(PANEL.open())}
    rows = []
    prod = {}
    with PROD.open() as fh:
        for r in csv.DictReader(fh):
            if r["status"] in {"ok", "cached"}:
                prod.setdefault(r["ligand"], {})[r["target"]] = float(r["affinity"])
    for lig in ligs:
        for pdb in ("3POZ", "3RCD"):
            rows.append(
                {
                    "pair": "EGFR/HER2",
                    "seed": PRIMARY,
                    "ligand": lig,
                    "class": cls[lig],
                    "pdb": pdb,
                    "status": "production_corrected_box",
                    "vina_mode1": prod[lig][pdb],
                }
            )
    jobs = [(pdb, lig, seed) for seed in NEW_SEEDS for lig in ligs for pdb in ("3POZ", "3RCD")]
    print(f"fiveseed new jobs={len(jobs)} workers={WORKERS}", flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(dock, pdb, lig, seed, boxes[pdb]) for pdb, lig, seed in jobs]
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            rec["class"] = cls[rec["ligand"]]
            rows.append(rec)
            done += 1
            if done % 20 == 0 or rec["status"] not in {"ok", "cached"}:
                print(f"[{done}/{len(jobs)}] {rec['pdb']} {rec['ligand']} seed={rec['seed']} {rec['status']} {rec.get('vina_mode1','')}", flush=True)
    dest = OUT / "scores_vina_mode1_fiveseed.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair", "seed", "ligand", "class", "pdb", "status", "vina_mode1"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["seed"], r["ligand"], r["pdb"])))
    n_ok = sum(1 for r in rows if r["status"] in {"ok", "cached", "production_corrected_box"})
    print(f"done {n_ok}/{len(rows)} -> {dest}")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
