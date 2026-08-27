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
    near(one(cognate, pdb="3POZ", pose_rank="1")["best_top1_A"], 9.5054)
    near(one(cognate, pdb="3POZ", pose_rank="1")["best_top3_A"], 6.227)
    near(one(cognate, pdb="3POZ", pose_rank="1")["best_all_deposited_A"], 0.7599)
    assert one(cognate, pdb="3POZ", pose_rank="1")["pass_top1_lt2"] == "0"
    assert one(cognate, pdb="3POZ", pose_rank="1")["pass_top3_lt2"] == "0"
    near(one(cognate, pdb="3RCD", pose_rank="1")["best_top1_A"], 1.8546)
    assert one(cognate, pdb="3RCD", pose_rank="1")["pass_top1_lt2"] == "1"

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
        "reconstructed QC",
        "179 include / 7 uncertain / 0 exclude",
        "17 unique pairs",
        "AND-like dual filter",
        "does not replace Table 2",
        "zero pairs meeting the pre-frozen primary external gate",
        "applicability stress-test",
        "not formal pose-gold validation",
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
    assert all(r["human_include_exclude"] in {"include", "exclude", "uncertain"} for r in assay)
    assert all(r.get("reviewed_by", "") != "" for r in assay)
    assert sum(r["human_include_exclude"] == "include" for r in assay) == 179
    assert sum(r["human_include_exclude"] == "uncertain" for r in assay) == 7
    assert sum(r["human_include_exclude"] == "exclude" for r in assay) == 0
    assert all(r["human_reviewed_class"] == r["frozen_class"] for r in assay)

    bdb = TAB / "bindingdb_independence_summary_v1.csv"
    if bdb.exists() and bdb.stat().st_size > 0:
        summary = rows("bindingdb_independence_summary_v1.csv")
        assert all(r["packaged_as_external_validation"] == "0" for r in summary)

    native = rows("external_slice_summary_v1.csv")
    assert all(r["packaged_as_external_evaluation"] == "0" for r in native)
    assert all(r["gate"] == "insufficient" for r in native)
    egfr_native = one(native, pair="EGFR/HER2")
    assert (egfr_native["n_dual"], egfr_native["n_A_only"], egfr_native["n_B_only"]) == (
        "180",
        "10",
        "20",
    )
    ache_native = one(native, pair="AChE/BChE")
    assert (ache_native["n_dual"], ache_native["n_A_only"], ache_native["n_B_only"]) == (
        "4",
        "8",
        "14",
    )
    flow = rows("external_candidate_flow.csv")
    native_egfr = one(flow, pair="EGFR/HER2", layer="native_paired_theta6")
    assert native_egfr["n_dual"] == "371"
    mcl1_panel = one(rows("mcl1_bclxl_panel_freeze_v1.csv"), pair="MCL1/Bcl-xL")
    assert (
        mcl1_panel["panel_dual"],
        mcl1_panel["panel_A_only"],
        mcl1_panel["panel_B_only"],
        mcl1_panel["panel_neither"],
    ) == ("24", "24", "24", "24")
    assert mcl1_panel["docked"] == "1"
    assert "applicability_stress_test" in mcl1_panel["pose_gold_gate"]
    assert "fail_3WIZ" in mcl1_panel["pose_gold_gate"]
    rec = one(rows("mcl1_bclxl_receptor_freeze_v1.csv"), pdb_id="3WIY")
    near(rec["resolution_A"], 2.15)
    assert rec["primary_chain"] == "A"
    assert rec["pose_gold_gate_run"] == "1"
    mcl1_auroc = rows_at("data/mcl1_bclxl_panel_v0/tables/formulation_auroc_MBX_v1.csv")
    assert all(r["panel_role"] == "applicability_stress_test" for r in mcl1_auroc)
    near(one(mcl1_auroc, contrast="summary_min")["auroc"], 0.6087)
    near(one(mcl1_auroc, contrast="D_vs_A_pocketB")["auroc"], 0.7935)
    near(one(mcl1_auroc, contrast="D_vs_B_pocketA")["auroc"], 0.6087)
    near(one(mcl1_auroc, contrast="D_vs_neither_mean")["auroc"], 0.6277)
    lc6 = rows_at("data/mcl1_bclxl_panel_v0/tables/cognate_qc_lc6_v1.csv")
    near(one(lc6, pdb="3WIY")["rmsd_best_of_top3"], 1.689)
    assert one(lc6, pdb="3WIY")["best_top3_pass_lt_2"] == "1"
    near(one(lc6, pdb="3WIZ")["rmsd_best_of_top3"], 2.011)
    assert one(lc6, pdb="3WIZ")["best_top3_pass_lt_2"] == "0"
    assert "not_topology_aware" in one(lc6, pdb="3WIY")["rmsd_method"]
    assert "not a formal pose-gold result" in one(lc6, pdb="3WIY")["note"]

    census = rows("theta6_pair_census_v1.csv")
    assert len(census) == 49
    assert sum(int(r["directional_n10"]) for r in census) == 17
    assert sum(int(r["formulation_n10"]) for r in census) == 17
    assert sum(int(r["docked_in_this_paper"]) for r in census) == 4

    caliper = rows("property_caliper_match_v1.csv")
    near(one(caliper, pair="EGFR/HER2", contrast="D_vs_B_pocketA", caliper_sd="1.0")["auroc_matched"], 0.5664)
    near(one(caliper, pair="AChE/BChE", contrast="D_vs_B_pocketA", caliper_sd="1.0")["auroc_matched"], 0.4615)

    and_rows = rows("and_filter_operating_point_v1.csv")
    near(
        one(and_rows, pair="EGFR/HER2", score="vina_worst", dual_percentile="50")["precision_dual"],
        0.2979,
    )
    near(
        one(and_rows, pair="EGFR/HER2", score="vina_worst", dual_percentile="90")["hardneg_fraction_pass"],
        0.8696,
    )

    ligand = rows("ligand_only_fullmap_auroc_v1.csv")
    near(one(ligand, pair="EGFR/HER2", contrast="D_vs_neither")["ecfp4_groupkfold_auroc"], 0.9214)
    near(one(ligand, pair="EGFR/HER2", contrast="D_vs_B")["ecfp4_groupkfold_auroc"], 0.8636)
    near(one(ligand, pair="EGFR/HER2", contrast="summary_min_ecfp4")["ecfp4_groupkfold_auroc"], 0.8013)

    zh = (ROOT / "docs" / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    assert "confidence≥8 与 Homo sapiens 过滤未重建" not in zh
    assert "不作为外部验证" in zh
    assert "纳入/排除与构建体/突变核查栏仍为空" not in zh
    assert "人工纳入/排除仍待本地阅读原文" not in zh
    print("revision validation: PASS")


if __name__ == "__main__":
    main()
