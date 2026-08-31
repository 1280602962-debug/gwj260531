#!/usr/bin/env python3
"""Audit class/property patterns and rank-extreme sensitivity of docking failures."""
from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"

SPECS = {
    "AChE/BChE": {
        "panel": "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "scores": "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "status": "data/ache_bche_panel_v0/tables/job_status.csv",
        "id": "panel_id",
        "score_id": "ligand",
        "A": "vina_ACHE",
        "B": "vina_BCHE",
    },
    "PIK3CA/PIK3CB": {
        "panel": "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "scores": "data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        "status": "data/pik3ca_pik3cb_panel_v0/tables/job_status.csv",
        "id": "panel_id",
        "score_id": "ligand",
        "A": "vina_PIK3CA",
        "B": "vina_PIK3CB",
    },
}


def read_csv(path):
    with (ROOT / path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def score(value):
    return -float(value) if value not in ("", None) else None


def auroc(pos, neg):
    wins = sum(p > n for p in pos for n in neg)
    ties = sum(p == n for p in pos for n in neg)
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def largest_fragment(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(fragments, key=lambda m: m.GetNumHeavyAtoms())


def main():
    failure_rows, sensitivity_rows = [], []
    for pair, cfg in SPECS.items():
        panel = read_csv(cfg["panel"])
        scores = {r[cfg["score_id"]]: r for r in read_csv(cfg["scores"])}
        statuses = read_csv(cfg["status"])
        status_by_ligand = {}
        for row in statuses:
            if row["status"] not in {"success", "exists"}:
                status_by_ligand.setdefault(row["ligand"], []).append(
                    f"{row['target']}:{row['status']}:{row.get('reason', '')}".rstrip(":")
                )

        records = []
        for ligand in panel:
            ligand_id = ligand[cfg["id"]]
            scored = scores.get(ligand_id, {})
            a, b = score(scored.get(cfg["A"])), score(scored.get(cfg["B"]))
            records.append({"ligand": ligand_id, "class": ligand["class"], "A": a, "B": b})
            if a is None or b is None:
                mol = largest_fragment(ligand["smiles"])
                failure_rows.append(
                    {
                        "pair": pair,
                        "ligand": ligand_id,
                        "class": ligand["class"],
                        "score_A_available": int(a is not None),
                        "score_B_available": int(b is not None),
                        "heavy_atoms": mol.GetNumHeavyAtoms(),
                        "mw": round(Descriptors.MolWt(mol), 3),
                        "clogp": round(Descriptors.MolLogP(mol), 3),
                        "tpsa": round(Descriptors.TPSA(mol), 3),
                        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
                        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
                        "status_reason": ";".join(status_by_ligand.get(ligand_id, [])),
                        "note": "largest organic/heavy-atom fragment; descriptive coverage audit",
                    }
                )

        complete = [r for r in records if r["A"] is not None and r["B"] is not None]
        for contrast, negative_class, end in (
            ("D_vs_A_pocketB", "A_only", "B"),
            ("D_vs_B_pocketA", "B_only", "A"),
        ):
            complete_pos = [r[end] for r in complete if r["class"] == "dual"]
            complete_neg = [r[end] for r in complete if r["class"] == negative_class]
            arm_pos = [r[end] for r in records if r["class"] == "dual" and r[end] is not None]
            arm_neg = [r[end] for r in records if r["class"] == negative_class and r[end] is not None]
            total_pos = sum(r["class"] == "dual" for r in records)
            total_neg = sum(r["class"] == negative_class for r in records)
            observed_pairs = len(arm_pos) * len(arm_neg)
            total_pairs = total_pos * total_neg
            observed_wins = auroc(arm_pos, arm_neg) * observed_pairs
            sensitivity_rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "complete_case_n_pos": len(complete_pos),
                    "complete_case_n_neg": len(complete_neg),
                    "complete_case_auroc": round(auroc(complete_pos, complete_neg), 4),
                    "arm_available_n_pos": len(arm_pos),
                    "arm_available_n_neg": len(arm_neg),
                    "arm_available_auroc": round(auroc(arm_pos, arm_neg), 4),
                    "full_panel_n_pos": total_pos,
                    "full_panel_n_neg": total_neg,
                    "rank_extreme_lower_bound": round(observed_wins / total_pairs, 4),
                    "rank_extreme_upper_bound": round(
                        (observed_wins + total_pairs - observed_pairs) / total_pairs, 4
                    ),
                    "note": "bounds assign every missing positive/negative comparison against the claimed direction; not an imputation model",
                }
            )

    write_csv(TAB / "docking_failed_ligand_properties_v1.csv", failure_rows)
    write_csv(TAB / "docking_failure_rank_extreme_v1.csv", sensitivity_rows)
    for row in sensitivity_rows:
        print(row)


if __name__ == "__main__":
    main()
