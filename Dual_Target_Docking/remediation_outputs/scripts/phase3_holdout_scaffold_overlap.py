#!/usr/bin/env python3
"""Main vs unused-pool holdout Bemis–Murcko overlap (post-fix wording support).

Does not delete holdout. Reports shared scaffolds. Also writes an optional
scaffold-disjoint holdout *subset* as exploratory sensitivity only.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs"
CANON = ROOT / "audit_outputs/canonical_ligand_table.csv"


def murcko(smiles: str) -> str:
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return ""


def main() -> int:
    rows = list(csv.DictReader(CANON.open(newline="", encoding="utf-8")))
    by_pair: dict[str, dict[str, list]] = defaultdict(lambda: {"main": [], "holdout": []})
    for r in rows:
        pair = r["pair"]
        status = (r.get("main_or_holdout") or r.get("main_holdout_status") or r.get("analysis_set") or "").lower()
        if "holdout" in status:
            bucket = "holdout"
        elif "main" in status:
            bucket = "main"
        else:
            continue
        scaf = r.get("scaffold_id") or r.get("murcko_scaffold") or murcko(r.get("canonical_smiles") or r.get("smiles") or "")
        by_pair[pair][bucket].append(
            {
                "ligand_id": r.get("ligand_id"),
                "chembl": r.get("molecule_chembl_id"),
                "smiles": r.get("canonical_smiles") or r.get("smiles") or "",
                "scaffold": scaf,
                "cls": r.get("primary_class_theta6") or r.get("primary_state_theta6") or r.get("panel_class"),
            }
        )

    out_rows = []
    expl = []
    for pair, d in by_pair.items():
        main_scaf = {x["scaffold"] for x in d["main"] if x["scaffold"]}
        hold_scaf = {x["scaffold"] for x in d["holdout"] if x["scaffold"]}
        shared = sorted(main_scaf & hold_scaf)
        n_hold = len(hold_scaf)
        frac = (len(shared) / n_hold) if n_hold else ""
        out_rows.append(
            {
                "pair": pair,
                "n_main_ligands": len(d["main"]),
                "n_holdout_ligands": len(d["holdout"]),
                "n_main_scaffolds": len(main_scaf),
                "n_holdout_scaffolds": n_hold,
                "n_shared_scaffolds": len(shared),
                "shared_scaffold_fraction_of_holdout": frac,
                "shared_scaffolds": ";".join(shared),
                "definition": "unused-pool internal holdout; not scaffold-disjoint",
            }
        )
        # exploratory: drop holdout ligands whose Murcko is in main
        kept = [x for x in d["holdout"] if x["scaffold"] and x["scaffold"] not in main_scaf]
        counts: dict[str, int] = defaultdict(int)
        for x in kept:
            counts[x["cls"] or "unknown"] += 1
        expl.append(
            {
                "pair": pair,
                "n_holdout_original": len(d["holdout"]),
                "n_holdout_scaffold_disjoint": len(kept),
                "n_removed_shared_scaffold": len(d["holdout"]) - len(kept),
                "D": counts.get("dual", 0),
                "A": counts.get("A_only", 0),
                "B": counts.get("B_only", 0),
                "N": counts.get("neither", 0),
                "note": "exploratory sensitivity only; does not replace unused-pool holdout",
            }
        )

    path = OUT / "holdout_scaffold_overlap_post_fix.csv"
    fields = [
        "pair",
        "n_main_ligands",
        "n_holdout_ligands",
        "n_main_scaffolds",
        "n_holdout_scaffolds",
        "n_shared_scaffolds",
        "shared_scaffold_fraction_of_holdout",
        "shared_scaffolds",
        "definition",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    epath = OUT / "phase3_holdout" / "scaffold_disjoint_holdout_subset_counts.csv"
    epath.parent.mkdir(parents=True, exist_ok=True)
    with epath.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(expl[0].keys()))
        w.writeheader()
        w.writerows(expl)
    print(f"wrote {path}")
    print(f"wrote {epath}")
    for r in out_rows:
        print(r["pair"], "shared", r["n_shared_scaffolds"], "frac", r["shared_scaffold_fraction_of_holdout"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
