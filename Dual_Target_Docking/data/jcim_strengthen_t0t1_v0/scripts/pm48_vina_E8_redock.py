#!/usr/bin/env python3
"""T1.1 B1a: Re-dock PM48 RDKit ligands at E=8 (4L23/4JT6). Reuses existing meeko ligands."""
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0")
REPO = Path(__file__).resolve().parents[3]
PANEL = REPO / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv"
SEED = 20260727
EXHAUST = 8
N_MODES = 9
TARGETS = ["4L23", "4JT6"]
VINA = "/home/gwj/miniconda3/bin/vina"
MAX_WORKERS = 6


def load_box(target: str) -> dict:
    return json.loads((ROOT / "boxes" / f"{target}_box.json").read_text())


def parse_affinities(log_text: str) -> list[float]:
    affs = []
    for line in log_text.splitlines():
        m = re.search(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            affs.append(float(m.group(1)))
    return affs


def run_one(target: str, lig: str, ligand_pdbqt: Path) -> dict:
    box = load_box(target)
    rec = ROOT / "receptors" / f"{target}_receptor.pdbqt"
    out = ROOT / "logs" / "vina_E8" / f"{target}_{lig}_out.pdbqt"
    log = ROOT / "logs" / "vina_E8" / f"{target}_{lig}.log"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        VINA,
        "--receptor", str(rec),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(EXHAUST),
        "--num_modes", str(N_MODES),
        "--seed", str(SEED),
        "--cpu", "1",
        "--out", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write_text(proc.stdout + "\n" + proc.stderr)
    affs = parse_affinities(proc.stdout + proc.stderr)
    best = min(affs) if affs else None
    return {
        "target": target,
        "ligand": lig,
        "status": "ok" if proc.returncode == 0 and best is not None else "fail",
        "best_affinity": best,
        "n_modes": len(affs),
    }


def main():
    rows = list(csv.DictReader(PANEL.open()))
    jobs = []
    for r in rows:
        lig = r["panel_id"]
        pdbqt = ROOT / "ligands_pdbqt" / f"{lig}.pdbqt"
        if not pdbqt.exists():
            print(f"MISSING pdbqt {lig}", file=sys.stderr)
            continue
        for t in TARGETS:
            jobs.append((t, lig, pdbqt))

    results = []
    print(f"E8 redock jobs: {len(jobs)}", flush=True)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(run_one, t, l, p): (t, l) for t, l, p in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            results.append(res)
            print(f"[{i}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']} {res.get('best_affinity')}", flush=True)

    # pivot to wide scores
    by_lig = {}
    for r in results:
        if r["status"] != "ok":
            continue
        by_lig.setdefault(r["ligand"], {"ligand": r["ligand"]})
        by_lig[r["ligand"]][f"{r['target']}_affinity_E8"] = r["best_affinity"]

    out_dir = ROOT / "tables"
    out_dir.mkdir(exist_ok=True)
    fields = ["ligand", "4L23_affinity_E8", "4JT6_affinity_E8"]
    with (out_dir / "scores_vina_E8_best.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for lig in sorted(by_lig):
            w.writerow(by_lig[lig])

    # also copy to repo
    repo_tab = REPO / "data/pik3ca_mtor_panel48_rdkit_v0/tables"
    repo_tab.mkdir(parents=True, exist_ok=True)
    (repo_tab / "scores_vina_E8_best.csv").write_text((out_dir / "scores_vina_E8_best.csv").read_text())

    ok = sum(1 for r in results if r["status"] == "ok")
    print(f"done ok={ok}/{len(results)}")
    return 0 if ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
