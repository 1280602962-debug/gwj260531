#!/usr/bin/env python3
"""Class-wise chemistry audit on the current eligible score master.

Descriptive only. Same descriptors as the earlier four-pair confounding audit.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import ECFP_NBITS, ECFP_RADIUS, PRIMARY_PAIRS, add_io_args, io_paths, is_primary_row  # noqa: E402

CANON = ROOT / "results/canonical"
MASTER = CANON / "current_score_master.csv"
OUT = CANON / "class_chemistry_summary.csv"
CLASSES = ("dual", "A_only", "B_only", "neither")


def qsummary(values):
    x = np.asarray(values, dtype=float)
    if len(x) == 0:
        return "", "", ""
    return (
        round(float(np.median(x)), 4),
        round(float(np.percentile(x, 25)), 4),
        round(float(np.percentile(x, 75)), 4),
    )


def main() -> int:
    import argparse

    global MASTER, OUT
    parser = argparse.ArgumentParser(description="Class-wise chemistry audit on the eligible score master.")
    add_io_args(parser)
    args = parser.parse_args()
    outdir, MASTER = io_paths(args.outdir, args.master)
    outdir.mkdir(parents=True, exist_ok=True)
    OUT = outdir / "class_chemistry_summary.csv"

    rows = list(csv.DictReader(MASTER.open(encoding="utf-8-sig", newline="")))
    out = []
    for pair in PRIMARY_PAIRS:
        recs = []
        for r in rows:
            if r.get("pair") != pair or not is_primary_row(r):
                continue
            mol = Chem.MolFromSmiles(r.get("smiles") or "")
            if mol is None:
                continue
            fp = Chem.RDKFingerprint(mol) if False else None
            from rdkit.Chem import AllChem

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, ECFP_RADIUS, nBits=ECFP_NBITS)
            recs.append(
                {
                    "ligand": r["ligand_id"],
                    "cls": r["primary_class_theta6"],
                    "mol": mol,
                    "fp": fp,
                    "heavy": Descriptors.HeavyAtomCount(mol),
                    "mw": Descriptors.MolWt(mol),
                    "clogp": Descriptors.MolLogP(mol),
                    "tpsa": Descriptors.TPSA(mol),
                    "charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
                    "rot": Lipinski.NumRotatableBonds(mol),
                    "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "__acyclic__",
                }
            )
        dual = [r for r in recs if r["cls"] == "dual"]
        dual_fps = [r["fp"] for r in dual]
        for cls in CLASSES:
            group = [r for r in recs if r["cls"] == cls]
            if not group:
                continue
            nn = []
            for r in group:
                comparison = dual_fps if cls != "dual" else [d["fp"] for d in dual if d["ligand"] != r["ligand"]]
                if not comparison:
                    nn.append(float("nan"))
                else:
                    nn.append(max(DataStructs.TanimotoSimilarity(r["fp"], fp) for fp in comparison))
            sc = Counter(r["scaffold"] for r in group)
            hm, hq1, hq3 = qsummary([r["heavy"] for r in group])
            mw, mw1, mw3 = qsummary([r["mw"] for r in group])
            lp, lp1, lp3 = qsummary([r["clogp"] for r in group])
            tp, tp1, tp3 = qsummary([r["tpsa"] for r in group])
            ch, ch1, ch3 = qsummary([r["charge"] for r in group])
            rt, rt1, rt3 = qsummary([r["rot"] for r in group])
            nnm, nn1, nn3 = qsummary([x for x in nn if np.isfinite(x)])
            out.append(
                {
                    "pair": pair,
                    "class": cls,
                    "n_ligands": len(group),
                    "n_murcko_scaffolds": len(sc),
                    "singleton_scaffold_ligand_fraction": round(sum(1 for n in sc.values() if n == 1) / len(group), 4),
                    "heavy_median": hm,
                    "heavy_q1": hq1,
                    "heavy_q3": hq3,
                    "mw_median": mw,
                    "mw_q1": mw1,
                    "mw_q3": mw3,
                    "clogp_median": lp,
                    "clogp_q1": lp1,
                    "clogp_q3": lp3,
                    "tpsa_median": tp,
                    "tpsa_q1": tp1,
                    "tpsa_q3": tp3,
                    "formal_charge_median": ch,
                    "formal_charge_q1": ch1,
                    "formal_charge_q3": ch3,
                    "rotatable_bonds_median": rt,
                    "rotatable_bonds_q1": rt1,
                    "rotatable_bonds_q3": rt3,
                    "nearest_dual_ecfp4_median": nnm,
                    "nearest_dual_ecfp4_q1": nn1,
                    "nearest_dual_ecfp4_q3": nn3,
                    "nearest_dual_note": "nearest other dual" if cls == "dual" else "nearest member of dual class",
                    "interpretation": "descriptive post-hoc confounding audit on current eligible scores; no multiplicity-adjusted tests",
                }
            )
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(out[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("wrote", OUT, "n=", len(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
