#!/usr/bin/env python3
"""Assemble one primary chemically mapped RMSD table for all 14 receptors.

Does not redock. Reads locked CalcRMS or production-QC tables already in git.
The earlier Hungarian assignment for the eight added receptors is not used here.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_novelty_v0/tables/primary_cognate_rmsd_calcrrms_v1.csv"

ORDER = [
    ("PIK3CA", "4L23"),
    ("mTOR", "4JT6"),
    ("AChE", "4EY7"),
    ("BChE", "4BDS"),
    ("EGFR", "3POZ"),
    ("HER2", "3RCD"),
    ("F2", "4UDW"),
    ("F10", "2JKH"),
    ("JAK1", "6N7A"),
    ("TYK2", "3LXP"),
    ("JAK2", "8BXH"),
    ("PPARG", "9V8H"),
    ("PPARA", "6LXA"),
    ("PPARD", "5U3Q"),
]


def read_csv(rel: str) -> list[dict]:
    path = ROOT / rel
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def r3(value: float) -> str:
    return f"{float(value):.3f}"


def main() -> int:
    rank = read_csv("data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv")
    added = read_csv(
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/layer3_cognate_rmsd_calcrrms_v1.csv"
    )
    pm = read_csv("data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16.csv")

    by_pdb: dict[str, dict] = {}
    seen = set()
    for row in rank:
        pdb = row["pdb"]
        if pdb in seen:
            continue
        seen.add(pdb)
        by_pdb[pdb] = {
            "protein": row["target"],
            "pdb": pdb,
            "exhaustiveness": int(row["exhaustiveness"]),
            "n_modes": int(row["n_modes_deposited"]),
            "top1_A": r3(row["best_top1_A"]),
            "top3_A": r3(row["best_top3_A"]),
            "best_A": r3(row["best_all_deposited_A"]),
            "mapping_method": row["mapping_method"],
            "rmsd_method": row["rmsd_method"],
            "source_table": "cognate_rank_rmsd_reaudit_v1.csv",
            "pose_in_git": 1,
            "note": "EGFR 3POZ and HER2 3RCD poses are reconstructed QC, not recovered production files."
            if pdb in {"3POZ", "3RCD"}
            else "",
        }

    for row in added:
        by_pdb[row["pdb"]] = {
            "protein": row["protein"],
            "pdb": row["pdb"],
            "exhaustiveness": 8,
            "n_modes": int(row["n_modes"]),
            "top1_A": r3(row["calcrrms_top1_A"]),
            "top3_A": r3(row["calcrrms_top3_A"]),
            "best_A": r3(row["calcrrms_best_A"]),
            "mapping_method": row["mapping_method"],
            "rmsd_method": "RDKit symmetry-aware CalcRMS; no superposition",
            "source_table": "layer3_cognate_rmsd_calcrrms_v1.csv",
            "pose_in_git": 1,
            "note": "2JKH/BI7 used CCD SMILES because the OpenBabel SDF had invalid nitrogen valence."
            if row["pdb"] == "2JKH"
            else "",
        }

    for row in pm:
        if row["seed"] != "20260727":
            continue
        pdb = row["target"]
        best = r3(row["rmsd_best_of_9"])
        top1 = r3(row["rmsd_mode1"])
        best_mode = int(row["best_of_9_mode"])
        top3 = best if best_mode <= 3 else ""
        by_pdb[pdb] = {
            "protein": "PIK3CA" if pdb == "4L23" else "mTOR",
            "pdb": pdb,
            "exhaustiveness": 16,
            "n_modes": 9,
            "top1_A": top1,
            "top3_A": top3,
            "best_A": best,
            "mapping_method": "Meeko SMILES index map + automorphism min CalcRMS",
            "rmsd_method": "RDKit symmetry-aware CalcRMS; no superposition",
            "source_table": "pm48_01_rmsd_E16.csv",
            "pose_in_git": 0,
            "note": "Historical production QC; cognate SDF/PDBQT not in git. Not a new docking run.",
        }

    rows = []
    for protein, pdb in ORDER:
        rec = dict(by_pdb[pdb])
        rec["protein"] = protein
        rec["pass_best_lt2"] = int(float(rec["best_A"]) < 2.0)
        rec["pass_top1_lt2"] = int(float(rec["top1_A"]) < 2.0)
        rows.append(rec)

    fields = [
        "protein",
        "pdb",
        "exhaustiveness",
        "n_modes",
        "top1_A",
        "top3_A",
        "best_A",
        "pass_best_lt2",
        "pass_top1_lt2",
        "mapping_method",
        "rmsd_method",
        "source_table",
        "pose_in_git",
        "note",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 14:
        raise SystemExit(f"expected 14 receptors, got {len(rows)}")
    if not all(r["pass_best_lt2"] == 1 for r in rows):
        raise SystemExit("search-coverage gate failed for at least one receptor")
    print("wrote", OUT, "n=", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
