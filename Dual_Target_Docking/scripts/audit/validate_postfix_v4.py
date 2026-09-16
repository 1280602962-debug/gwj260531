#!/usr/bin/env python3
"""Post-fix canonical validation. Ground truth is machine-readable files, not manuscript text."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from canonical_paths_v4 import (  # noqa: E402
    ALLOWED_CLASSES,
    DT_ROOT,
    FIGURE_STEMS,
    OFFICIAL_BOXES,
    PRIMARY_PAIRS,
    THETA_PRIMARY,
    WITHDRAWN_PAIRS,
    fail,
    load_csv,
    load_json,
    near,
    one,
    require_exists,
    require_finite,
    sha256_raw,
)

PASS: list[str] = []


def ok(name: str, detail: str = "") -> None:
    extra = f" ({detail})" if detail else ""
    print(f"PASS  {name}{extra}")
    PASS.append(name)


def main() -> int:
    lock = require_exists("docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md")
    lock_text = lock.read_text(encoding="utf-8")
    if "FIGURE_TABLE_LOCK_POSTFIX_V4.md" not in lock_text:
        fail("V4 lock file does not identify itself")
    ok("v4_lock_present")

    audit_md = require_exists("remediation_outputs/POST_FIX_AUDIT_REPORT.md").read_text(encoding="utf-8")
    if "98 PASS" not in audit_md or "0 WARNING" not in audit_md or "0 FAIL" not in audit_md:
        fail("canonical post-fix scientific audit is not 98 PASS / 0 WARNING / 0 FAIL")
    ok("postfix_scientific_audit_98_0_0")

    membership = load_csv(
        "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv",
        ["pair", "ligand", "cls", "score_A", "score_B"],
        ["pair", "ligand"],
    )
    if set(r["pair"] for r in membership) != set(PRIMARY_PAIRS):
        fail("membership pair set is not an exact match to PRIMARY_PAIRS")
    if set(r["cls"] for r in membership) != ALLOWED_CLASSES:
        fail("membership classes are not limited to dual/A_only/B_only/neither")
    if any(p in {r["pair"] for r in membership} for p in WITHDRAWN_PAIRS):
        fail("withdrawn pair in membership")
    ok("pair_set_exact_match")

    unified = load_csv(
        "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
        ["pair", "label_rule", "pocket_matched_summary_min", "auroc_D_vs_A", "auroc_D_vs_B"],
        ["pair", "label_rule"],
    )
    master = load_csv(
        "remediation_outputs/canonical_tables/post_fix_master_metrics.csv",
        ["pair", "estimand", "independent"],
        None,
    )
    for pair in ("EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"):
        u = one(unified, "unified_threshold_sensitivity_v2.csv", pair=pair, label_rule=THETA_PRIMARY)
        require_finite("unified_threshold_sensitivity_v2.csv", u, ["auroc_D_vs_A", "auroc_D_vs_B", "pocket_matched_summary_min"])
        smin = min(float(u["auroc_D_vs_A"]), float(u["auroc_D_vs_B"]))
        near(float(u["pocket_matched_summary_min"]), smin, 5e-4, f"{pair} summary_min identity")
        m = [
            r
            for r in master
            if r["pair"] == pair and r["estimand"] == "summary_min" and r["analysis_set"] == "main"
        ]
        if len(m) != 1:
            fail(f"{pair}: expected 1 master summary_min row, got {len(m)}")
        require_finite("post_fix_master_metrics.csv", m[0], ["independent"])
        near(float(m[0]["independent"]), smin, 5e-4, f"{pair} master vs unified summary_min")
    ok("summary_min_matches_min_of_arms_in_canonical_tables")

    ranking = load_csv(
        "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv",
        ["pair", "n_ranked", "n_dual", "top_k", "top_dual", "top_A_only", "top_B_only", "top_neither", "ef_dual_10pct"],
        ["pair"],
    )
    if [r["pair"] for r in ranking] != list(PRIMARY_PAIRS):
        fail("ranking pair order is not the locked eight-pair order")
    for row in ranking:
        n, k = int(row["n_ranked"]), int(row["top_k"])
        if k != math.ceil(0.1 * n):
            fail(f"{row['pair']} k != ceil(0.1 n)")
        if int(row["top_dual"]) + int(row["top_A_only"]) + int(row["top_B_only"]) + int(row["top_neither"]) != k:
            fail(f"{row['pair']} top composition != k")
        ef = (int(row["top_dual"]) / k) / (int(row["n_dual"]) / n)
        near(float(row["ef_dual_10pct"]), ef, 1e-9, f"{row['pair']} EF")
    ok("ranking_operating_point_k_and_EF")

    rmsd = load_csv(
        "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv",
        ["pdb", "calcrrms_top1_A", "pose_status"],
        None,
    )
    p3 = one(rmsd, "all14_cognate_rmsd_calcrrms_v1.csv", pdb="3POZ")
    r3 = one(rmsd, "all14_cognate_rmsd_calcrrms_v1.csv", pdb="3RCD")
    near(float(p3["calcrrms_top1_A"]), 1.019, 5e-4, "3POZ top1")
    near(float(r3["calcrrms_top1_A"]), 1.947, 5e-4, "3RCD top1")
    if p3["pose_status"] != "canonical_heavy_atom_redock":
        fail("3POZ pose_status is not canonical_heavy_atom_redock")
    ok("cognate_RMSD_3POZ_1.019_3RCD_1.947")

    for pdb, expect in OFFICIAL_BOXES.items():
        box = load_json(f"remediation_outputs/phase1_boxes/{pdb}_box_corrected.json")
        for key, value in expect.items():
            got = box.get(key)
            if isinstance(value, float):
                near(float(got), value, 1e-3, f"{pdb}.{key}")
            elif got != value:
                fail(f"{pdb}.{key}={got!r}")
    ok("box_json_official")

    gnina = load_csv(
        "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv",
        ["pair", "engine", "contrast", "auroc"],
        None,
    )
    g_smin = one(
        gnina,
        "independent_dock_formulation_v1.csv",
        pair="EGFR/HER2",
        engine="gnina_dock_mode1",
        contrast="summary_min",
    )
    g_neither = one(
        gnina,
        "independent_dock_formulation_v1.csv",
        pair="EGFR/HER2",
        engine="gnina_dock_mode1",
        contrast="D_vs_neither_mean",
    )
    near(float(g_smin["auroc"]), 0.2645, 5e-4, "GNINA EGFR summary_min")
    near(float(g_neither["auroc"]), 0.7370, 5e-4, "GNINA EGFR D-vs-neither")
    ok("GNINA_EGFR_frozen_values")

    plotted = load_json("figures/jcim_article/plotted_values_postfix.json")
    recs = plotted.get("records")
    if not isinstance(recs, list) or not recs:
        fail("plotted_values_postfix.json missing records")
    def rec(fig, panel, pair, metric):
        hits = [r for r in recs if r.get("figure")==fig and r.get("panel")==panel and r.get("pair")==pair and r.get("metric")==metric]
        if len(hits) != 1:
            fail(f"plotted record missing/duplicate {fig}/{panel}/{pair}/{metric}")
        return hits[0]
    checks = [
        ("Figure 2", "C", "EGFR/HER2", "summary_min", 0.324, 0.002),
        ("Figure 2", "A", "EGFR/HER2", "fixed_delta_pocketA", 0.462, 0.01),
        ("Figure S3", "A", "EGFR 3POZ", "top1_rmsd_A", 1.019, 0.001),
    ]
    for fig, panel, pair, metric, expect, tol in checks:
        value = float(rec(fig, panel, pair, metric)["raw_value"])
        near(value, expect, tol, f"plotted {fig}{panel} {metric}")
    ok("plotted_values_postfix_matches_canonical_tables")

    for name, stem in FIGURE_STEMS.items():
        exts = ("png", "tif") if name == "TOC" else ("pdf", "png", "tif")
        for ext in exts:
            require_exists(f"figures/jcim_article/{stem}.{ext}")
    ok("current_figures_present")

    for rel, _role in [
        ("docs/MANUSCRIPT_JCIM_EN.md", "manuscript"),
        ("docs/MANUSCRIPT_JCIM_ZH.md", "manuscript"),
        ("docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md", "si"),
        ("docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md", "si"),
    ]:
        require_exists(rel)
    ok("current_manuscript_SI_present")

    # Publication copies in the pack must match sources when the pack exists.
    pack_en = DT_ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md"
    src_en = DT_ROOT / "docs/MANUSCRIPT_JCIM_EN.md"
    if pack_en.is_file():
        if sha256_raw(pack_en) != sha256_raw(src_en):
            fail("submission_pack EN manuscript hash != docs source")
        ok("pack_EN_manuscript_hash_parity")
    else:
        fail("submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md missing")

    print(f"validate_postfix_v4: {len(PASS)} PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
