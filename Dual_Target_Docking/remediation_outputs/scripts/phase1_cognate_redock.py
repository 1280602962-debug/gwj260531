#!/usr/bin/env python3
"""Cognate redock of 3POZ/3RCD ligands into corrected heavy-atom boxes.

Same Vina protocol as primary (E=8, n_modes=9, energy_range=3). Seed 20260727
to match other cognate QC jobs. Outputs stay under remediation_outputs.
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
P1 = ROOT / "remediation_outputs/phase1_boxes"
OUT = ROOT / "remediation_outputs/phase1_cognate_redock"
VINA = "/home/gwj/miniconda3/bin/vina"
SEED = 20260727

# Prefer already-prepared cognate ligand pdbqt if present in the panel pack.
CANDIDATES = {
    "3POZ": [
        ROOT / "data/egfr_her2_panel40_v0/cognate/3POZ_ligand.pdbqt",
        ROOT / "data/egfr_her2_panel40_v0/ligands_pdbqt/03_ligand.pdbqt",
    ],
    "3RCD": [
        ROOT / "data/egfr_her2_panel40_v0/cognate/3RCD_ligand.pdbqt",
        ROOT / "data/egfr_her2_panel40_v0/ligands_pdbqt/03_ligand.pdbqt",
    ],
}


def parse_mode1(path: Path) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def find_ligand(pdb: str) -> Path | None:
    for p in CANDIDATES.get(pdb, []):
        if p.exists():
            return p
    # search common cognate locations
    hits = list((ROOT / "data/egfr_her2_panel40_v0").rglob(f"*{pdb}*ligand*.pdbqt"))
    return hits[0] if hits else None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    recdir = ROOT / "data/egfr_her2_panel40_v0/receptors"
    rows = []
    for pdb in ("3POZ", "3RCD"):
        box = json.loads((P1 / f"{pdb}_box_corrected.json").read_text())
        lig = find_ligand(pdb)
        rec = recdir / f"{pdb}_receptor.pdbqt"
        if lig is None or not rec.exists():
            rows.append({"pdb": pdb, "status": "missing_input", "ligand": str(lig), "affinity": ""})
            continue
        out = OUT / f"{pdb}_cognate_corrected_out.pdbqt"
        log = OUT / f"{pdb}_cognate_corrected.log"
        conf = OUT / f"{pdb}_cognate_corrected.txt"
        conf.write_text(
            "\n".join(
                [
                    f"receptor = {rec}",
                    f"ligand = {lig}",
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
                    f"seed = {SEED}",
                    f"out = {out}",
                ]
            )
            + "\n"
        )
        proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
        log.write_text((proc.stdout or "") + (proc.stderr or ""))
        aff = parse_mode1(out) if out.exists() else None
        rows.append(
            {
                "pdb": pdb,
                "status": "ok" if aff is not None else f"FAIL rc={proc.returncode}",
                "ligand": str(lig),
                "affinity": aff if aff is not None else "",
                "seed": SEED,
            }
        )
        print(rows[-1])
    dest = OUT / "cognate_redock_corrected_box.csv"
    import csv

    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
