#!/usr/bin/env python3
"""Regression checks for the reviewer-facing revision analyses."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"


def rows(name):
    with (TAB / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def rows_at(relative_path):
    with (ROOT / relative_path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def one(data, **keys):
    found = [row for row in data if all(row[key] == value for key, value in keys.items())]
    assert len(found) == 1, (keys, len(found))
    return found[0]


def near(value, expected, tolerance=5e-4):
    assert abs(float(value) - expected) <= tolerance, (value, expected)


def main():
    high = rows("high_confidence_summary_v1.csv")
    assert sum(int(r["n_frozen_scored"]) for r in high) == 352
    assert sum(int(r["n_class_matches_frozen"]) for r in high) == 352
    for pair, expected in {
        "EGFR/HER2": 0.4297,
        "AChE/BChE": 0.6058,
        "PIK3CA/PIK3CB": 0.5000,
        "PIK3CA/mTOR": 0.6921,
    }.items():
        near(one(high, pair=pair)["summary_min"], expected)

    primary = rows_at(
        "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    )
    expected_primary = {
        "EGFR/HER2": (0.4297, 0.2818, 0.5775),
        "AChE/BChE": (0.6058, 0.4370, 0.7303),
        "PIK3CA/PIK3CB": (0.5000, 0.3502, 0.6495),
        "PIK3CA/mTOR": (0.6921, 0.4702, 0.8133),
    }
    for pair, (point, lo, hi) in expected_primary.items():
        record = one(primary, pair=pair, label_rule="theta_6.0")
        near(record["pocket_matched_summary_min"], point)
        near(record["ci_lo"], lo)
        near(record["ci_hi"], hi)

    ml_compare = rows_at(
        "data/jcim_strengthen_t0t1_v0/tables/ligand_ml_scaffold_vs_random_v1.csv"
    )
    mean_delta = sum(float(r["delta_random_minus_scaffold"]) for r in ml_compare) / len(
        ml_compare
    )
    near(mean_delta, 0.0257875, 1e-7)
    near(one(ml_compare, pair="EGFR/HER2", contrast="D_vs_B")["auroc_scaffold_GroupKFold"], 0.8895)

    equal = rows("formulation_equal_score_negative_v1.csv")
    near(
        one(equal, pair="EGFR/HER2", contrast="D_vs_B_or_neither_pocketA")[
            "delta_neither_minus_selective"
        ],
        0.3783,
    )

    overlap = rows("complete_case_usable_pchembl_overlap_v1.csv")
    expected_overlap = {
        "EGFR/HER2": 0.145119,
        "AChE/BChE": 0.340172,
        "PIK3CA/PIK3CB": 0.233349,
        "PIK3CA/mTOR": 0.265252,
    }
    for pair, expected in expected_overlap.items():
        near(one(overlap, pair=pair)["fraction_union_measured_both"], expected, 1e-6)

    documents = rows("source_document_concentration_v1.csv")
    pm_neither = one(documents, pair="PIK3CA/mTOR", **{"class": "neither"})
    assert (pm_neither["n_ligands"], pm_neither["n_unique_documents"]) == ("4", "1")
    near(pm_neither["top_document_record_fraction"], 1.0)

    failures = rows("docking_failure_rank_extreme_v1.csv")
    near(one(failures, pair="AChE/BChE", contrast="D_vs_A_pocketB")["rank_extreme_lower_bound"], 0.5599)
    near(one(failures, pair="PIK3CA/PIK3CB", contrast="D_vs_A_pocketB")["arm_available_auroc"], 0.6952)

    cognate = rows("cognate_rank_rmsd_reaudit_v1.csv")
    near(one(cognate, pdb="4BDS", pose_rank="1")["best_top1_A"], 4.7941)
    near(one(cognate, pdb="4BDS", pose_rank="1")["best_top3_A"], 0.3856)
    near(one(cognate, pdb="2WXF", pose_rank="1")["best_top1_A"], 0.4048)

    chemistry = rows("class_chemistry_summary_v1.csv")
    near(one(chemistry, pair="AChE/BChE", **{"class": "dual"})["tpsa_median"], 76.02)
    near(one(chemistry, pair="AChE/BChE", **{"class": "A_only"})["nearest_dual_ecfp4_median"], 0.236)

    manuscript = (ROOT / "docs" / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    required_phrases = (
        "Four-Pair Formulation Audit",
        "14.5%–34.0%",
        "rank-extreme lower bounds",
        "not a top-ranked-pose validation",
        "do not define a target-general reliability boundary",
        "is not claimed as external validation",
        "not stably estimable",
    )
    for phrase in required_phrases:
        assert phrase in manuscript, phrase

    blocked = rows("document_blocked_cv_summary_v1.csv")
    near(one(blocked, pair="EGFR/HER2", contrast="D_vs_B")["rank_auroc_full"], 0.4297)
    near(one(blocked, pair="EGFR/HER2", contrast="D_vs_B")["ecfp4_auroc_oof"], 0.6228)
    assert one(blocked, pair="PIK3CA/mTOR", contrast="D_vs_B")["status"] == "cannot_stably_estimate"
    assert sum(r["status"] == "ok" for r in blocked) == 7

    time_split = rows("time_split_auroc_v1.csv")
    assert all(
        r["packaged_as_external_validation"] == "0"
        for r in time_split
        if r["cutoff_year"] == "2018"
    )
    counts = rows("time_split_class_counts_v1.csv")
    primary_egfr = one(
        counts, cutoff_year="2018", pair="EGFR/HER2", split="test_on_or_after"
    )
    assert (primary_egfr["n_dual"], primary_egfr["n_A_only"], primary_egfr["n_B_only"]) == (
        "6",
        "3",
        "14",
    )

    assay = rows("assay_context_priority_ligands_v1.csv")
    assert len(assay) == 186
    # After local human review, include/exclude must be filled for every priority ligand.
    assert all(r["human_include_exclude"] in {"include", "exclude", "uncertain"} for r in assay)
    assert all(r.get("reviewed_by", "") != "" for r in assay)

    zh = (ROOT / "docs" / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    assert "confidence≥8 与 Homo sapiens 过滤未重建" not in zh
    assert "不作为外部验证" in zh
    print("revision validation: PASS")


if __name__ == "__main__":
    main()
