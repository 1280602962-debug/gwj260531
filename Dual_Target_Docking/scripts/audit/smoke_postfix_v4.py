#!/usr/bin/env python3
"""Sub-minute post-fix smoke checks. Ground truth is canonical CSVs/JSON, not manuscript text."""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from canonical_paths_v4 import (  # noqa: E402
    ALLOWED_CLASSES,
    DT_ROOT,
    OFFICIAL_BOXES,
    PRIMARY_PAIRS,
    THETA_PRIMARY,
    WITHDRAWN_PAIRS,
    fail,
    is_legacy_path,
    load_csv,
    load_json,
    near,
    one,
    require_exists,
    require_finite,
    required_paths,
)

PASS: list[str] = []


def ok(name: str, detail: str = "") -> None:
    extra = f" ({detail})" if detail else ""
    print(f"PASS  {name}{extra}")
    PASS.append(name)


def main() -> int:
    t0 = time.monotonic()

    style_text = require_exists("figures/jcim_article/scripts/jcim_figure_style.py").read_text(encoding="utf-8")
    if "PRIMARY_PAIRS" not in style_text:
        fail("jcim_figure_style.py missing PRIMARY_PAIRS")
    if any(p not in style_text for p in PRIMARY_PAIRS):
        fail("PRIMARY_PAIRS text missing a locked pair")
    if "PIK3CA/PIK3CB" in style_text and 'PRIMARY_PAIRS = [' in style_text:
        block = style_text.split("PRIMARY_PAIRS = [", 1)[1].split("]", 1)[0]
        if "PIK3CA/PIK3CB" in block:
            fail("withdrawn PIK3CA/PIK3CB is inside PRIMARY_PAIRS")
    ok("primary_pairs_locked", f"n={len(PRIMARY_PAIRS)}")

    unified = load_csv(
        "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
        ["pair", "label_rule", "n_dual", "n_A_only", "n_B_only", "pocket_matched_summary_min", "auroc_D_vs_A", "auroc_D_vs_B"],
        ["pair", "label_rule"],
    )
    theta_pairs = {r["pair"] for r in unified if r["label_rule"] == THETA_PRIMARY}
    if "EGFR/HER2" not in theta_pairs or "AChE/BChE" not in theta_pairs:
        fail(f"theta=6 primary labels missing EGFR/AChE; have {sorted(theta_pairs)}")
    ok("theta6_primary_labels", THETA_PRIMARY)

    membership = load_csv(
        "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv",
        ["pair", "ligand", "cls", "pA", "pB", "score_A", "score_B"],
        ["pair", "ligand"],
    )
    pairs = tuple(dict.fromkeys(r["pair"] for r in membership))
    if pairs != PRIMARY_PAIRS and set(pairs) != set(PRIMARY_PAIRS):
        fail(f"membership pair set {set(pairs)} != PRIMARY_PAIRS {set(PRIMARY_PAIRS)}")
    if set(pairs) != set(PRIMARY_PAIRS) or len(pairs) != 8:
        fail(f"expected exactly 8 primary pairs, got {pairs}")
    if any(p in pairs for p in WITHDRAWN_PAIRS):
        fail("withdrawn PIK3CA/PIK3CB entered current primary pair set")
    ok("exactly_8_primary_pairs", ",".join(PRIMARY_PAIRS))

    classes = {r["cls"] for r in membership}
    if classes != ALLOWED_CLASSES:
        fail(f"allowed classes limited to D/A/B/N; got {classes}")
    ok("allowed_classes_DABN")

    for row in membership:
        require_finite("review_scored_membership_v1.csv", row, ["pA", "pB", "score_A", "score_B"])
        if float(row["pA"]) >= 6.0 and float(row["pB"]) >= 6.0:
            if row["cls"] != "dual":
                fail(f"dual-positive rule broken for {row['pair']} {row['ligand']}")
    ok("dual_positive_theta6")

    egfr_mem = [r for r in membership if r["pair"] == "EGFR/HER2"]
    if {r["cls"] for r in egfr_mem} != ALLOWED_CLASSES:
        fail("EGFR membership classes are not D/A/B/N")
    n_egfr = (
        sum(1 for r in egfr_mem if r["cls"] == "dual"),
        sum(1 for r in egfr_mem if r["cls"] == "A_only"),
        sum(1 for r in egfr_mem if r["cls"] == "B_only"),
    )
    if n_egfr != (28, 38, 32):
        fail(f"EGFR membership n {n_egfr} != (28, 38, 32)")
    ok("EGFR_membership_n_28_38_32")

    master = load_csv(
        "remediation_outputs/canonical_tables/post_fix_master_metrics.csv",
        ["pair", "analysis_set", "estimand", "score_definition", "independent"],
        None,
    )
    mode1_rows = [
        r
        for r in master
        if r["pair"] == "EGFR/HER2"
        and r["analysis_set"] == "main"
        and r["estimand"] in {"AUROC_D_vs_A_pocketB", "AUROC_D_vs_B_pocketA"}
    ]
    da_row = [r for r in mode1_rows if r["estimand"] == "AUROC_D_vs_A_pocketB"]
    db_row = [r for r in mode1_rows if r["estimand"] == "AUROC_D_vs_B_pocketA"]
    if len(da_row) != 1 or not da_row[0]["estimand"].endswith("pocketB"):
        fail("D vs A-only must be the pocket-B estimand")
    if len(db_row) != 1 or not db_row[0]["estimand"].endswith("pocketA"):
        fail("D vs B-only must be the pocket-A estimand")
    proto = [
        r
        for r in master
        if r["pair"] == "JAK1/JAK2"
        and r["estimand"] == "AUROC_D_vs_A_pocketB"
        and r["analysis_set"] == "main"
    ]
    if len(proto) != 1 or "mode1" not in proto[0]["score_definition"]:
        fail("locked protocol score_definition is not Vina mode1")
    ok("mode1_primary_score")
    ok("DA_uses_pocketB_DB_uses_pocketA")

    equal = load_csv(
        "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv",
        [
            "pair",
            "contrast",
            "score",
            "auroc_dual_vs_selective",
            "auroc_dual_vs_neither",
            "delta_neither_minus_selective",
            "delta_ci_lo",
            "delta_ci_hi",
        ],
        ["pair", "contrast"],
    )
    pocket_a = one(equal, "formulation_equal_score_negative_v1.csv", pair="EGFR/HER2", contrast="D_vs_B_or_neither_pocketA")
    require_finite("formulation_equal_score_negative_v1.csv", pocket_a, ["delta_neither_minus_selective", "auroc_dual_vs_selective", "auroc_dual_vs_neither"])
    if pocket_a["score"] != "vina_A":
        fail(f"fixed-score channel must be identical vina_A, got {pocket_a['score']}")
    pocket_b = one(equal, "formulation_equal_score_negative_v1.csv", pair="EGFR/HER2", contrast="D_vs_A_or_neither_pocketB")
    if pocket_b["score"] != "vina_B":
        fail(f"fixed-score channel must be identical vina_B, got {pocket_b['score']}")
    # Same score for selective and neither arms of each contrast (fixed-score channel).
    near(
        float(pocket_a["delta_neither_minus_selective"]),
        float(pocket_a["auroc_dual_vs_neither"]) - float(pocket_a["auroc_dual_vs_selective"]),
        5e-4,
        "EGFR pocket-A delta identity",
    )
    near(float(pocket_a["delta_neither_minus_selective"]), 0.4621, 5e-4, "EGFR fixed delta")
    ok("fixed_score_channel_identical")
    ok("EGFR_fixed_delta_0.4621")

    egfr_u = one(unified, "unified_threshold_sensitivity_v2.csv", pair="EGFR/HER2", label_rule=THETA_PRIMARY)
    require_finite("unified_threshold_sensitivity_v2.csv", egfr_u, ["pocket_matched_summary_min"])
    near(float(egfr_u["pocket_matched_summary_min"]), 0.3237, 5e-4, "EGFR summary_min")
    smin = min(float(egfr_u["auroc_D_vs_A"]), float(egfr_u["auroc_D_vs_B"]))
    near(float(egfr_u["pocket_matched_summary_min"]), smin, 5e-4, "EGFR summary_min identity")
    ok("EGFR_summary_min_0.3237")
    ok("summary_min_is_min_of_two_arms", f"{smin:.4f}")

    ache_u = one(unified, "unified_threshold_sensitivity_v2.csv", pair="AChE/BChE", label_rule=THETA_PRIMARY)
    n_scored = (int(ache_u["n_dual"]), int(ache_u["n_A_only"]), int(ache_u["n_B_only"]))
    if n_scored != (27, 26, 28):
        fail(f"AChE n_scored {n_scored} != (27, 26, 28)")
    ok("AChE_n_scored_27_26_28")

    ranking = load_csv(
        "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv",
        [
            "pair",
            "n_ranked",
            "n_dual",
            "top_k",
            "top_dual",
            "top_A_only",
            "top_B_only",
            "top_neither",
            "ef_dual_10pct",
            "ranking_rule",
            "tie_rule",
        ],
        ["pair"],
    )
    rank_pairs = [r["pair"] for r in ranking]
    if rank_pairs != list(PRIMARY_PAIRS):
        fail(f"ranking pair order {rank_pairs} != {list(PRIMARY_PAIRS)}")
    for row in ranking:
        n = int(row["n_ranked"])
        k = int(row["top_k"])
        if k != math.ceil(0.1 * n):
            fail(f"{row['pair']}: top_k={k} != ceil(0.1*{n})")
        top_sum = int(row["top_dual"]) + int(row["top_A_only"]) + int(row["top_B_only"]) + int(row["top_neither"])
        if top_sum != k:
            fail(f"{row['pair']}: top D+A+B+N={top_sum} != k={k}")
        ef = (int(row["top_dual"]) / k) / (int(row["n_dual"]) / n)
        if abs(ef - float(row["ef_dual_10pct"])) > 1e-9:
            fail(f"{row['pair']}: EF formula mismatch {ef} vs {row['ef_dual_10pct']}")
        if "k=ceil(0.10 n)" not in row["ranking_rule"] or "EF=(top_dual/k)/(n_dual/n)" not in row["ranking_rule"]:
            fail(f"{row['pair']}: ranking_rule does not record EF/k formula")
        if "ligand ID" not in row["tie_rule"]:
            fail(f"{row['pair']}: tie-breaking is not deterministic by ligand ID")
    ok("top10_k_ceil_0.1n")
    ok("top_DABN_equals_k")
    ok("EF_formula")

    leak = load_csv(
        "audit_outputs/cv_leakage_audit.csv",
        ["pair", "contrast", "fold", "overlap_count"],
        ["pair", "contrast", "fold"],
    )
    leak_pairs = {r["pair"] for r in leak}
    if leak_pairs != set(PRIMARY_PAIRS):
        fail(f"ECFP leakage table pairs {leak_pairs} != primary set")
    if any(int(r["overlap_count"]) != 0 for r in leak):
        fail("ECFP scaffold leakage overlap_count != 0")
    ok("ECFP_scaffold_leakage_0")

    holdout = load_csv(
        "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv",
        ["pair", "ligand", "chembl", "cls"],
        ["pair", "ligand"],
    )
    main_ids = {(r["pair"], r["ligand"]) for r in membership}
    hold_ids = {(r["pair"], r["ligand"]) for r in holdout}
    if main_ids & hold_ids:
        fail(f"main/holdout exact ligand-id overlap: {sorted(main_ids & hold_ids)[:5]}")
    ache_main = load_csv(
        "data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        ["panel_id", "molecule_chembl_id"],
        ["panel_id"],
    )
    hoab = load_csv(
        "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv",
        ["holdout_id", "molecule_chembl_id"],
        ["holdout_id"],
    )
    chembl_overlap = {r["molecule_chembl_id"] for r in ache_main} & {r["molecule_chembl_id"] for r in hoab}
    if chembl_overlap:
        fail(f"AChE main/holdout ChEMBL overlap: {sorted(chembl_overlap)[:5]}")
    ok("main_holdout_exact_ligand_overlap_0")

    wrong = load_csv(
        "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv",
        ["pair", "set", "delta_ci_lo", "delta_ci_hi", "ci_excludes_zero", "n_dual", "n_A_only", "n_B_only"],
        ["pair", "set"],
    )
    egfr_w = one(wrong, "wrong_pocket_paired_delta_bootstrap_v1.csv", pair="EGFR/HER2", set="main_panel")
    ache_w = one(wrong, "wrong_pocket_paired_delta_bootstrap_v1.csv", pair="AChE/BChE", set="main_panel")
    require_finite("wrong_pocket_paired_delta_bootstrap_v1.csv", egfr_w, ["delta_ci_lo", "delta_ci_hi"])
    require_finite("wrong_pocket_paired_delta_bootstrap_v1.csv", ache_w, ["delta_ci_lo", "delta_ci_hi"])
    egfr_lo, egfr_hi = float(egfr_w["delta_ci_lo"]), float(egfr_w["delta_ci_hi"])
    if egfr_lo > 0 or egfr_hi < 0:
        fail(f"EGFR matched CI should include 0, got [{egfr_lo}, {egfr_hi}]")
    if egfr_w["ci_excludes_zero"].lower() == "true":
        fail("EGFR ci_excludes_zero must be False")
    ache_lo, ache_hi = float(ache_w["delta_ci_lo"]), float(ache_w["delta_ci_hi"])
    if ache_lo <= 0 <= ache_hi:
        fail(f"AChE matched CI should exclude 0, got [{ache_lo}, {ache_hi}]")
    if ache_w["ci_excludes_zero"].lower() != "true":
        fail("AChE ci_excludes_zero must be True")
    if (int(ache_w["n_dual"]), int(ache_w["n_A_only"]), int(ache_w["n_B_only"])) != (27, 26, 28):
        fail("AChE wrong-pocket n_scored != 27/26/28")
    ok("EGFR_matched_CI_includes_0")
    ok("AChE_matched_CI_excludes_0")

    for pdb, expect in OFFICIAL_BOXES.items():
        box = load_json(f"remediation_outputs/phase1_boxes/{pdb}_box_corrected.json")
        for key, value in expect.items():
            got = box.get(key)
            if isinstance(value, float):
                if abs(float(got) - value) > 1e-3:
                    fail(f"{pdb} {key}={got} != official {value}")
            elif got != value:
                fail(f"{pdb} {key}={got!r} != official {value!r}")
    ok("corrected_box_JSON_matches_official_3POZ_3RCD")

    for rel in required_paths():
        path = DT_ROOT / rel
        if path.exists() and is_legacy_path(path):
            fail(f"canonical required path is legacy: {rel}")
    analysis_inputs = [
        "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
        "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv",
        "remediation_outputs/phase1_boxes/3POZ_box_corrected.json",
        "remediation_outputs/canonical_tables/post_fix_master_metrics.csv",
        "figures/jcim_article/plotted_values_postfix.json",
    ]
    for rel in analysis_inputs:
        path = require_exists(rel)
        if is_legacy_path(path):
            fail(f"active analysis input resolves under data/_legacy_archive: {rel}")
    ok("no_active_input_under_legacy_archive")
    ok("no_withdrawn_PIK3CA_PIK3CB")

    elapsed = time.monotonic() - t0
    if elapsed >= 60:
        fail(f"smoke exceeded 1 minute ({elapsed:.1f}s)")
    print(f"smoke_postfix_v4: {len(PASS)} PASS in {elapsed:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
