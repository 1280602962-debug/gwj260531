#!/usr/bin/env python3
"""Machine assay-context audit for priority DualFourClass ligands.

This is a risk-stratified extraction, not a completed human paper review.
Human include/exclude decisions remain empty until local reading of source
documents. The script does not change frozen labels.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "EGFR/HER2": ("CHEMBL203", "CHEMBL1824", "EGFR", "HER2"),
    "AChE/BChE": ("CHEMBL220", "CHEMBL1914", "AChE", "BChE"),
    "PIK3CA/PIK3CB": ("CHEMBL4005", "CHEMBL3145", "PIK3CA", "PIK3CB"),
    "PIK3CA/mTOR": ("CHEMBL4005", "CHEMBL2842", "PIK3CA", "mTOR"),
}
THETA = 6.0
BORDER_DELTA = 0.30


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        if value in ("", None):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    if not pos or not neg:
        return float("nan")
    p = np.asarray(pos, float)
    n = np.asarray(neg, float)
    delta = p[:, None] - n[None, :]
    return float(((delta > 0).sum() + 0.5 * (delta == 0).sum()) / (len(p) * len(n)))


def truthy(value) -> bool:
    return str(value) in {"1", "True", "true"}


def loo_influence(recs, score_key, pos_cls, neg_cls):
    pos = [rec for rec in recs if rec["cls"] == pos_cls]
    neg = [rec for rec in recs if rec["cls"] == neg_cls]
    base = auroc([r[score_key] for r in pos], [r[score_key] for r in neg])
    out = {}
    for rec in pos + neg:
        p = [r[score_key] for r in pos if r["ligand"] != rec["ligand"]]
        n = [r[score_key] for r in neg if r["ligand"] != rec["ligand"]]
        value = auroc(p, n)
        out[rec["ligand"]] = None if value != value or base != base else float(base - value)
    return base, out


def mutation_hint(text: str) -> str:
    lowered = (text or "").lower()
    hits = []
    for token in (
        "wild-type",
        "wild type",
        "wt ",
        "mutant",
        "mutation",
        "l858r",
        "t790m",
        "exon 19",
        "egfrviii",
        "h1047r",
        "e545k",
        "kinase-dead",
        "truncated",
        "delta",
    ):
        if token in lowered:
            hits.append(token.strip())
    return ";".join(dict.fromkeys(hits))


def main() -> None:
    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    activities = read_csv(TAB / "high_confidence_activity_audit_v1.csv")
    kept = [row for row in activities if truthy(row.get("keep"))]
    by_mol_target: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in kept:
        by_mol_target[(row["molecule_chembl_id"], row["target_chembl_id"])].append(row)

    top_docs = set()
    concentration = read_csv(TAB / "source_document_concentration_v1.csv")
    for row in concentration:
        if row.get("top_document"):
            top_docs.add((row["pair"], row["class"], row["top_document"]))

    recs_by_pair = defaultdict(list)
    for row in labels:
        recs_by_pair[row["pair"]].append(
            {
                "ligand": row["ligand"],
                "cls": row["frozen_class"],
                "vina_A": fnum(row["vina_A"]),
                "vina_B": fnum(row["vina_B"]),
            }
        )
    influence = {}
    for pair, recs in recs_by_pair.items():
        for contrast, pos_cls, neg_cls, score_key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            _, deltas = loo_influence(recs, score_key, pos_cls, neg_cls)
            for ligand, delta in deltas.items():
                prev = influence.get((pair, ligand), 0.0)
                if delta is not None and abs(delta) > abs(prev):
                    influence[(pair, ligand)] = delta

    abs_inf = sorted((abs(v), k) for k, v in influence.items() if v is not None)
    top_influence = {k for _, k in abs_inf[-20:]}

    audit_rows = []
    ligand_rows = []
    for row in labels:
        pair = row["pair"]
        target_a, target_b, name_a, name_b = TARGETS[pair]
        mol = row["molecule_chembl_id"]
        cls = row["frozen_class"]
        recs_a = by_mol_target[(mol, target_a)]
        recs_b = by_mol_target[(mol, target_b)]
        p_a = fnum(row["high_conf_max_A"])
        p_b = fnum(row["high_conf_max_B"])
        types_a = sorted({r["standard_type"] for r in recs_a if r.get("standard_type")})
        types_b = sorted({r["standard_type"] for r in recs_b if r.get("standard_type")})
        assays_a = sorted({r["assay_type"] for r in recs_a if r.get("assay_type")})
        assays_b = sorted({r["assay_type"] for r in recs_b if r.get("assay_type")})
        orgs_a = sorted({r["assay_organism"] for r in recs_a if r.get("assay_organism")})
        orgs_b = sorted({r["assay_organism"] for r in recs_b if r.get("assay_organism")})
        docs = sorted(
            {r["document_chembl_id"] for r in recs_a + recs_b if r.get("document_chembl_id")}
        )
        flags = []
        if pair == "EGFR/HER2" and cls in {"dual", "A_only", "B_only"}:
            flags.append("egfr_directional_priority")
        if pair == "PIK3CA/mTOR" and cls == "neither":
            flags.append("pm_neither_single_document")
        if p_a is not None and abs(p_a - THETA) <= BORDER_DELTA:
            flags.append("borderline_pA")
        if p_b is not None and abs(p_b - THETA) <= BORDER_DELTA:
            flags.append("borderline_pB")
        if len(types_a) > 1 or len(types_b) > 1:
            flags.append("mixed_endpoint")
        if ("B" in assays_a and "F" in assays_a) or ("B" in assays_b and "F" in assays_b):
            flags.append("biochem_and_functional")
        if any(org and org != "Homo sapiens" for org in orgs_a + orgs_b):
            flags.append("non_human_assay_organism")
        if any((pair, cls, doc) in top_docs for doc in docs):
            flags.append("top_document_series")
        if (pair, row["ligand"]) in top_influence:
            flags.append("high_auroc_influence")
        if not flags:
            continue

        inf = influence.get((pair, row["ligand"]))
        ligand_rows.append(
            {
                "pair": pair,
                "ligand": row["ligand"],
                "molecule_chembl_id": mol,
                "frozen_class": cls,
                "high_conf_class": row["high_conf_class_theta6"],
                "pchembl_A": row["high_conf_max_A"],
                "pchembl_B": row["high_conf_max_B"],
                "n_records_A": len(recs_a),
                "n_records_B": len(recs_b),
                "standard_types_A": ";".join(types_a),
                "standard_types_B": ";".join(types_b),
                "assay_types_A": ";".join(assays_a),
                "assay_types_B": ";".join(assays_b),
                "organisms_A": ";".join(orgs_a),
                "organisms_B": ";".join(orgs_b),
                "n_documents": len(docs),
                "documents": ";".join(docs),
                "loo_auroc_delta": "" if inf is None else round(inf, 4),
                "priority_flags": ";".join(flags),
                "machine_incomparability_flag": int(
                    "mixed_endpoint" in flags
                    or "biochem_and_functional" in flags
                    or "non_human_assay_organism" in flags
                ),
                "human_include_exclude": "",
                "human_reviewed_class": "",
                "human_rationale": "",
                "reviewed_by": "",
                "review_date": "",
                "note": "machine extraction only; frozen label unchanged until human review",
            }
        )
        for target_id, target_name, recs in (
            (target_a, name_a, recs_a),
            (target_b, name_b, recs_b),
        ):
            for rec in recs:
                audit_rows.append(
                    {
                        "pair": pair,
                        "ligand": row["ligand"],
                        "molecule_chembl_id": mol,
                        "frozen_class": cls,
                        "target_name": target_name,
                        "target_chembl_id": target_id,
                        "assay_chembl_id": rec.get("assay_chembl_id", ""),
                        "document_chembl_id": rec.get("document_chembl_id", ""),
                        "assay_type": rec.get("assay_type", ""),
                        "standard_type": rec.get("standard_type", ""),
                        "standard_relation": rec.get("standard_relation", ""),
                        "pchembl_value": rec.get("pchembl_value", ""),
                        "assay_organism": rec.get("assay_organism", ""),
                        "confidence_score": rec.get("confidence_score", ""),
                        "protein_construct": "",
                        "wildtype_or_mutant": mutation_hint(
                            " ".join(
                                [
                                    rec.get("assay_type", ""),
                                    rec.get("standard_type", ""),
                                    rec.get("assay_organism", ""),
                                ]
                            )
                        ),
                        "priority_flags": ";".join(flags),
                        "original_label": cls,
                        "human_include_exclude": "",
                        "human_reviewed_label": "",
                        "human_rationale": "",
                        "incomparable_record": "",
                        "note": "construct/mutation require source-document reading; API fields deposited here",
                    }
                )

    write_csv(TAB / "assay_context_priority_ligands_v1.csv", ligand_rows)
    write_csv(TAB / "assay_context_audit.csv", audit_rows)

    flag_counts = Counter()
    for row in ligand_rows:
        for flag in row["priority_flags"].split(";"):
            flag_counts[flag] += 1
    lines = [
        "# Assay-context machine audit",
        "",
        f"Priority ligands: {len(ligand_rows)} / {len(labels)} scored panel compounds.",
        f"Activity rows extracted: {len(audit_rows)}.",
        "",
        "This file is not a completed human audit. Include/exclude columns are empty.",
        "Frozen labels must not be overwritten from this machine extraction.",
        "",
        "| flag | n_ligands |",
        "|------|----------:|",
    ]
    for flag, count in flag_counts.most_common():
        lines.append(f"| {flag} | {count} |")
    lines += [
        "",
        "Human review SOP: `docs/ASSAY_CONTEXT_HUMAN_REVIEW_SOP.md`.",
        "If any reviewed label changes, recompute Table 2 before claiming robustness.",
        "",
    ]
    (ANALYSIS / "ASSAY_CONTEXT_MACHINE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"priority ligands {len(ligand_rows)}; activity rows {len(audit_rows)}")


if __name__ == "__main__":
    main()
