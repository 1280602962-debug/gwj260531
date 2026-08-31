#!/usr/bin/env python3
"""Audit usable-pChEMBL overlap and source-document concentration.

Absence from a cached molecule->pChEMBL map means no usable cached pChEMBL,
not experimental inactivity. Document concentration uses retained records from
the dated high-confidence current-ChEMBL view.
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PUBLIC = ROOT / "data" / "public_pair_selection"
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"

PAIR_MAPS = {
    "EGFR/HER2": ("mols_EGFR.json", "mols_HER2.json", "CHEMBL203", "CHEMBL1824"),
    "AChE/BChE": ("mols_ACHE.json", "mols_BCHE.json", "CHEMBL220", "CHEMBL1914"),
    "PIK3CA/PIK3CB": ("mols_PIK3CA.json", "mols_PIK3CB.json", "CHEMBL4005", "CHEMBL3145"),
    "PIK3CA/mTOR": ("mols_PIK3CA.json", "mols_MTOR.json", "CHEMBL4005", "CHEMBL2842"),
}


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    overlap_rows = []
    for pair, (file_a, file_b, _, _) in PAIR_MAPS.items():
        map_a = json.loads((PUBLIC / file_a).read_text(encoding="utf-8"))
        map_b = json.loads((PUBLIC / file_b).read_text(encoding="utf-8"))
        ids_a, ids_b = set(map_a), set(map_b)
        both = ids_a & ids_b
        union = ids_a | ids_b
        overlap_rows.append(
            {
                "pair": pair,
                "n_usable_pchembl_A": len(ids_a),
                "n_usable_pchembl_B": len(ids_b),
                "n_both": len(both),
                "n_A_only_measured": len(ids_a - ids_b),
                "n_B_only_measured": len(ids_b - ids_a),
                "fraction_union_measured_both": round(len(both) / len(union), 6),
                "jaccard_usable_pchembl": round(len(both) / len(union), 6),
                "note": "map overlap; absence means no usable cached pChEMBL, not inactivity",
            }
        )

    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    activities = [r for r in read_csv(TAB / "high_confidence_activity_audit_v1.csv") if r["keep"] == "1"]
    activity_index: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in activities:
        activity_index[(row["molecule_chembl_id"], row["target_chembl_id"])].append(row)

    document_rows = []
    for pair, (_, _, target_a, target_b) in PAIR_MAPS.items():
        for cls in ("dual", "A_only", "B_only", "neither"):
            kept = [r for r in labels if r["pair"] == pair and r["high_conf_class_theta6"] == cls]
            documents = []
            ligand_docs: dict[str, set[str]] = defaultdict(set)
            for ligand in kept:
                molecule = ligand["molecule_chembl_id"]
                for target in (target_a, target_b):
                    for activity in activity_index[(molecule, target)]:
                        doc = activity.get("document_chembl_id") or "__missing_document__"
                        documents.append(doc)
                        ligand_docs[molecule].add(doc)
            counts = Counter(documents)
            top_doc, top_records = counts.most_common(1)[0] if counts else ("", 0)
            top_ligands = sum(top_doc in docs for docs in ligand_docs.values()) if top_doc else 0
            document_rows.append(
                {
                    "pair": pair,
                    "class": cls,
                    "n_ligands": len(kept),
                    "n_retained_activity_records": len(documents),
                    "n_unique_documents": len(counts),
                    "top_document": top_doc,
                    "top_document_record_fraction": round(top_records / len(documents), 4) if documents else "",
                    "top_document_ligand_fraction": round(top_ligands / len(kept), 4) if kept else "",
                    "note": "current high-confidence retained records; descriptive source concentration",
                }
            )

    write_csv(TAB / "complete_case_usable_pchembl_overlap_v1.csv", overlap_rows)
    write_csv(TAB / "source_document_concentration_v1.csv", document_rows)
    print(json.dumps(overlap_rows, indent=2))
    print(json.dumps(document_rows, indent=2))


if __name__ == "__main__":
    main()
