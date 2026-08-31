#!/usr/bin/env python3
"""Describe class-wise chemical-space differences in the frozen four-pair panels.

This is a confounding audit, not a model-selection exercise.  It reports robust
descriptor summaries, Bemis--Murcko scaffold diversity, and nearest-neighbour
ECFP4 similarity to the dual class.  No hypothesis tests are performed because
the audit is post hoc and spans correlated descriptors and multiple classes.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

from benchmark_formulation_v1 import ROOT, SPEC, assemble

OUT = ROOT / "data" / "jcim_novelty_v0" / "tables"
CLASSES = ("dual", "A_only", "B_only", "neither")


def qsummary(values):
    x = np.asarray(values, dtype=float)
    return (
        round(float(np.median(x)), 4),
        round(float(np.percentile(x, 25)), 4),
        round(float(np.percentile(x, 75)), 4),
    )


def write_csv(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scaffold_key(mol, ligand):
    """Keep acyclic molecules together instead of treating an empty Murcko string as a parse failure."""
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return f"__parse_failure_{ligand}"
    return scaffold or "__acyclic__"


def main():
    rows = []
    for pair, cfg in SPEC.items():
        records = assemble(pair, cfg)
        dual = [r for r in records if r["cls"] == "dual"]
        dual_fps = [r["fp"] for r in dual]

        for cls in CLASSES:
            group = [r for r in records if r["cls"] == cls]
            if not group:
                continue
            formal_charge = []
            rotatable = []
            nn_dual = []
            scaffolds = []
            for r in group:
                mol = Chem.MolFromSmiles(r["smiles"])
                formal_charge.append(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
                rotatable.append(Lipinski.NumRotatableBonds(mol))
                scaffolds.append(scaffold_key(mol, r["ligand"]))
                comparison = dual_fps
                if cls == "dual":
                    comparison = [d["fp"] for d in dual if d["ligand"] != r["ligand"]]
                nn_dual.append(
                    max(DataStructs.TanimotoSimilarity(r["fp"], fp) for fp in comparison)
                    if comparison else float("nan")
                )

            scaffold_counts = Counter(scaffolds)
            descriptor_values = {
                "heavy": [r["heavy"] for r in group],
                "mw": [r["mw"] for r in group],
                "clogp": [r["clogp"] for r in group],
                "tpsa": [r["tpsa"] for r in group],
                "formal_charge": formal_charge,
                "rotatable_bonds": rotatable,
                "nearest_dual_ecfp4": nn_dual,
            }
            result = {
                "pair": pair,
                "class": cls,
                "n_ligands": len(group),
                "n_murcko_scaffolds": len(scaffold_counts),
                "singleton_scaffold_ligand_fraction": round(
                    sum(n for n in scaffold_counts.values() if n == 1) / len(group), 4
                ),
            }
            for name, values in descriptor_values.items():
                med, q1, q3 = qsummary(values)
                result[f"{name}_median"] = med
                result[f"{name}_q1"] = q1
                result[f"{name}_q3"] = q3
            result["nearest_dual_note"] = (
                "nearest other dual" if cls == "dual" else "nearest member of dual class"
            )
            result["interpretation"] = (
                "descriptive post-hoc confounding audit; correlated descriptors; no multiplicity-adjusted tests"
            )
            rows.append(result)

    path = OUT / "class_chemistry_summary_v1.csv"
    write_csv(path, rows)
    print(path)
    for row in rows:
        print(
            row["pair"], row["class"], "n=", row["n_ligands"],
            "scaffolds=", row["n_murcko_scaffolds"],
            "MW=", row["mw_median"],
            "cLogP=", row["clogp_median"],
            "TPSA=", row["tpsa_median"],
            "NNdual=", row["nearest_dual_ecfp4_median"],
        )


if __name__ == "__main__":
    main()
