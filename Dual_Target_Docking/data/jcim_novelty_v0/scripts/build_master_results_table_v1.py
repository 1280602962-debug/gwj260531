#!/usr/bin/env python3
"""Build MASTER_RESULTS_TABLE.csv from frozen DualFourClass sources.

No new docking. Numbers are copied from deposited CSVs and from the A4
summary/AUROC tables. EGFR frozen Table 2 (cached max) is kept distinct
from the A4 API-refetched max (one cache/API mismatch: EH120_060).
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
SR = ROOT / "data" / "jcim_structure_robust_v0" / "tables"
TH = ROOT / "data" / "jcim_strengthen_t0t1_v0" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

FIELDS = [
    "block",
    "manuscript_table",
    "pair",
    "setting",
    "metric",
    "value",
    "ci_lo",
    "ci_hi",
    "n_scored",
    "n_dual",
    "n_A_only",
    "n_B_only",
    "n_flip",
    "label_agreement",
    "note",
    "source_file",
]


def _read(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def _f(v):
    if v in ("", None):
        return ""
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def row(**kwargs) -> dict:
    out = {k: "" for k in FIELDS}
    out.update(kwargs)
    for k in ("value", "ci_lo", "ci_hi", "label_agreement"):
        if out[k] not in ("", None):
            out[k] = f"{float(out[k]):.6g}"
    return out


def main() -> None:
    theta = _read(TH / "unified_threshold_sensitivity_v2.csv")
    theta6 = {r["pair"]: r for r in theta if r["label_rule"] == "theta_6.0"}
    form = _read(TAB / "formulation_conventional_vs_directional_v1.csv")
    agg = _read(TAB / "aggregation_min_mean_geometric_harmonic_v1.csv")
    a4s = _read(TAB / "assay_max_vs_median_summary_v1.csv")
    a4a = _read(TAB / "assay_max_vs_median_auroc_v1.csv")
    census = _read(TAB / "docking_failure_census_v1.csv")
    desc = _read(TAB / "descriptor_all_four_directional_v1.csv")
    for r in desc:
        expected = max(float(r[f"{name}_summary_min"]) for name in ("heavy", "mw", "clogp", "tpsa"))
        actual = float(r["best_single_descriptor_summary_min"])
        if abs(actual - expected) > 1e-12:
            raise ValueError(
                f"{r['pair']}: best descriptor value {actual} does not equal max {expected}"
            )

    pm_jps = _read(SR / "pocket_matched_PM48_alt4JPS_v1.csv")[0]
    pm_dxt = _read(SR / "pocket_matched_PM48_alt5DXT_v1.csv")[0]
    pm_jsx = _read(SR / "pocket_matched_PM48_alt4JSX_v1.csv")[0]
    pab_jps = _read(SR / "pocket_matched_PAB_alt4JPS_v1.csv")[0]
    pab_dxt = _read(SR / "pocket_matched_PAB_alt5DXT_v1.csv")[0]

    rows: list[dict] = []

    # --- Table 2 frozen directional ---
    for pair, r in theta6.items():
        rows.append(
            row(
                block="primary_directional",
                manuscript_table="Table 2",
                pair=pair,
                setting="frozen_cached_max_theta6",
                metric="summary_min",
                value=r["pocket_matched_summary_min"],
                ci_lo=r["ci_lo"],
                ci_hi=r["ci_hi"],
                n_dual=r["n_dual"],
                n_A_only=r["n_A_only"],
                n_B_only=r["n_B_only"],
                n_scored=int(r["n_dual"]) + int(r["n_A_only"]) + int(r["n_B_only"]),
                note="Primary endpoint. Cached max pChEMBL labels. Neither excluded.",
                source_file="data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
            )
        )
        rows.append(
            row(
                block="primary_directional",
                manuscript_table="Table 2",
                pair=pair,
                setting="frozen_cached_max_theta6",
                metric="auroc_D_vs_A_pocketB",
                value=r["auroc_D_vs_A"],
                n_dual=r["n_dual"],
                n_A_only=r["n_A_only"],
                source_file="data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
            )
        )
        rows.append(
            row(
                block="primary_directional",
                manuscript_table="Table 2",
                pair=pair,
                setting="frozen_cached_max_theta6",
                metric="auroc_D_vs_B_pocketA",
                value=r["auroc_D_vs_B"],
                n_dual=r["n_dual"],
                n_B_only=r["n_B_only"],
                source_file="data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
            )
        )

    # --- Table 3 formulation ---
    for r in form:
        if r["contrast"] in ("D_vs_neither_mean", "D_vs_all_nondual_mean", "summary_min"):
            rows.append(
                row(
                    block="formulation",
                    manuscript_table="Table 3",
                    pair=r["pair"],
                    setting=r["formulation"],
                    metric=r["contrast"],
                    value=r["auroc"],
                    ci_lo=r["ci_lo"],
                    ci_hi=r["ci_hi"],
                    n_dual=r["n_pos"],
                    note=r.get("note", ""),
                    source_file="data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv",
                )
            )

    # --- aggregation sensitivity ---
    for r in agg:
        for metric in ("summary_min", "summary_mean", "summary_geometric", "summary_harmonic"):
            rows.append(
                row(
                    block="aggregation_sensitivity",
                    manuscript_table="Table S26",
                    pair=r["pair"],
                    setting="frozen_cached_max_theta6",
                    metric=metric,
                    value=r[metric],
                    n_dual=r["n_dual"],
                    n_A_only=r["n_A_only"],
                    n_B_only=r["n_B_only"],
                    note="Pair ranking unchanged under min/arithmetic/geometric/harmonic.",
                    source_file="data/jcim_novelty_v0/tables/aggregation_min_mean_geometric_harmonic_v1.csv",
                )
            )

    # --- A4 max vs median ---
    a4s_by = {r["pair"]: r for r in a4s}
    a4_agree_rows = []
    for r in a4a:
        pair = r["pair"]
        s = a4s_by[pair]
        n = int(s["n_ligands_scored"])
        n_flip = int(s["n_class_flip_theta6"])
        agree = 1.0 - n_flip / n
        frozen = theta6[pair]
        rows.append(
            row(
                block="assay_aggregation_A4",
                manuscript_table="Table S29",
                pair=pair,
                setting=f"api_{r['aggregation']}_theta6",
                metric="summary_min",
                value=r["summary_min"],
                ci_lo=r["ci_lo_D_vs_B"] if float(r["summary_min"]) == float(r["auroc_D_vs_B"]) else r["ci_lo_D_vs_A"],
                ci_hi=r["ci_hi_D_vs_B"] if float(r["summary_min"]) == float(r["auroc_D_vs_B"]) else r["ci_hi_D_vs_A"],
                n_scored=n,
                n_dual=r["n_dual"],
                n_A_only=r["n_A_only"],
                n_B_only=r["n_B_only"],
                n_flip=n_flip if r["aggregation"] == "median" else "",
                label_agreement=agree if r["aggregation"] == "median" else "",
                note=(
                    "API-refetched pChEMBL. Do not mix EGFR API-max 0.417 with frozen Table 2 0.430 "
                    "(EH120_060 cache/API mismatch on HER2)."
                    if pair == "EGFR/HER2"
                    else "API-refetched pChEMBL; cache max matches API max on this pair."
                ),
                source_file="data/jcim_novelty_v0/tables/assay_max_vs_median_auroc_v1.csv",
            )
        )
        if r["aggregation"] == "median":
            a4_agree_rows.append(
                {
                    "pair": pair,
                    "n_scored": n,
                    "n_cache_matches_both_max": s["n_cache_matches_both_max"],
                    "n_any_end_max_ne_median": s["n_any_end_max_ne_median"],
                    "n_class_flip_theta6": n_flip,
                    "label_agreement": f"{agree:.6g}",
                    "pct_class_flip": f"{100 * n_flip / n:.2f}",
                    "frozen_summary_min": frozen["pocket_matched_summary_min"],
                    "api_max_summary_min": next(
                        x["summary_min"] for x in a4a if x["pair"] == pair and x["aggregation"] == "max"
                    ),
                    "api_median_summary_min": r["summary_min"],
                    "api_max_n_dual_A_B": "/".join(
                        next(
                            f"{x['n_dual']}/{x['n_A_only']}/{x['n_B_only']}"
                            for x in a4a
                            if x["pair"] == pair and x["aggregation"] == "max"
                        ).split()
                    )
                    if False
                    else "",
                    "note": (
                        "EGFR frozen labels use cached max; API-max differs on EH120_060/CHEMBL24828 only. "
                        "Class-flip rate uses scored n, not construction n. Numeric max≠median is more common than class flips."
                    ),
                }
            )

    # rewrite agreement rows cleanly
    a4_agree_rows = []
    a4a_max = {r["pair"]: r for r in a4a if r["aggregation"] == "max"}
    a4a_med = {r["pair"]: r for r in a4a if r["aggregation"] == "median"}
    for pair, s in a4s_by.items():
        n = int(s["n_ligands_scored"])
        n_flip = int(s["n_class_flip_theta6"])
        agree = 1.0 - n_flip / n
        mx = a4a_max[pair]
        md = a4a_med[pair]
        a4_agree_rows.append(
            {
                "pair": pair,
                "n_scored": n,
                "n_cache_matches_both_max": s["n_cache_matches_both_max"],
                "n_any_end_max_ne_median": s["n_any_end_max_ne_median"],
                "n_class_flip_theta6": n_flip,
                "label_agreement": f"{agree:.6g}",
                "pct_class_flip": f"{100.0 * n_flip / n:.2f}",
                "frozen_summary_min": theta6[pair]["pocket_matched_summary_min"],
                "api_max_summary_min": mx["summary_min"],
                "api_median_summary_min": md["summary_min"],
                "frozen_n_dual_A_B": f"{theta6[pair]['n_dual']}/{theta6[pair]['n_A_only']}/{theta6[pair]['n_B_only']}",
                "api_max_n_dual_A_B": f"{mx['n_dual']}/{mx['n_A_only']}/{mx['n_B_only']}",
                "api_median_n_dual_A_B": f"{md['n_dual']}/{md['n_A_only']}/{md['n_B_only']}",
                "delta_api_median_minus_api_max": f"{float(md['summary_min']) - float(mx['summary_min']):.6g}",
                "delta_api_median_minus_frozen": f"{float(md['summary_min']) - float(theta6[pair]['pocket_matched_summary_min']):.6g}",
                "note": (
                    "Do not mix EGFR frozen 0.4297 with API-max 0.4170. "
                    "One cache/API mismatch: EH120_060 CHEMBL24828 (frozen A_only, API-max dual). "
                    if pair == "EGFR/HER2"
                    else "Frozen cached max equals API-max on this pair. "
                )
                + "Report label agreement = 1 - n_flip/n_scored, not flip count alone. "
                "Numeric max≠median does not imply a class flip at θ=6.0.",
            }
        )

    agree_path = TAB / "assay_max_vs_median_agreement_v1.csv"
    with agree_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=list(a4_agree_rows[0].keys()), lineterminator="\n"
        )
        w.writeheader()
        w.writerows(a4_agree_rows)

    # --- B5 / PM receptor realization ---
    two_pair = [
        {
            "pair": "PIK3CA/mTOR",
            "pik3ca_receptor": "4L23",
            "kept_pocket_B": "4JT6",
            "n_attempted": 48,
            "n_successful": 48,
            "n_failed": 0,
            "n_dual": 18,
            "n_A_only": 14,
            "n_B_only": 12,
            "auroc_D_vs_A": 0.7143,
            "auroc_D_vs_B": 0.6921,
            "summary_min": 0.6921,
            "ci_lo": 0.4638,
            "ci_hi": 0.8015,
            "delta_vs_original": 0.0,
            "weak_arm": "D_vs_B",
            "failed_ligand": "",
            "note": "Original frozen panel. Exhaustiveness 16. Source: unified_threshold_sensitivity_v2.csv",
        },
        {
            "pair": "PIK3CA/mTOR",
            "pik3ca_receptor": "4JPS",
            "kept_pocket_B": "4JT6",
            "n_attempted": 48,
            "n_successful": 48,
            "n_failed": 0,
            "n_dual": pm_jps["n_dual"],
            "n_A_only": pm_jps["n_A_only"],
            "n_B_only": pm_jps["n_B_only"],
            "auroc_D_vs_A": pm_jps["auroc_D_vs_A"],
            "auroc_D_vs_B": pm_jps["auroc_D_vs_B"],
            "summary_min": pm_jps["summary_min"],
            "ci_lo": pm_jps["summary_min_ci_lo"],
            "ci_hi": pm_jps["summary_min_ci_hi"],
            "delta_vs_original": float(pm_jps["summary_min"]) - 0.6921,
            "weak_arm": "D_vs_B",
            "failed_ligand": "",
            "note": "Pocket A replaced; pocket B frozen 4JT6. Exhaustiveness 16. Use deposited CSV CIs.",
        },
        {
            "pair": "PIK3CA/mTOR",
            "pik3ca_receptor": "5DXT",
            "kept_pocket_B": "4JT6",
            "n_attempted": 48,
            "n_successful": 48,
            "n_failed": 0,
            "n_dual": pm_dxt["n_dual"],
            "n_A_only": pm_dxt["n_A_only"],
            "n_B_only": pm_dxt["n_B_only"],
            "auroc_D_vs_A": pm_dxt["auroc_D_vs_A"],
            "auroc_D_vs_B": pm_dxt["auroc_D_vs_B"],
            "summary_min": pm_dxt["summary_min"],
            "ci_lo": pm_dxt["summary_min_ci_lo"],
            "ci_hi": pm_dxt["summary_min_ci_hi"],
            "delta_vs_original": float(pm_dxt["summary_min"]) - 0.6921,
            "weak_arm": "D_vs_B",
            "failed_ligand": "",
            "note": "Pocket A replaced; pocket B frozen 4JT6. Exhaustiveness 16. Use deposited CSV CIs.",
        },
        {
            "pair": "PIK3CA/mTOR",
            "pik3ca_receptor": "4L23",
            "kept_pocket_B": "4JSX",
            "n_attempted": 48,
            "n_successful": 48,
            "n_failed": 0,
            "n_dual": pm_jsx["n_dual"],
            "n_A_only": pm_jsx["n_A_only"],
            "n_B_only": pm_jsx["n_B_only"],
            "auroc_D_vs_A": pm_jsx["auroc_D_vs_A"],
            "auroc_D_vs_B": pm_jsx["auroc_D_vs_B"],
            "summary_min": pm_jsx["summary_min"],
            "ci_lo": pm_jsx["summary_min_ci_lo"],
            "ci_hi": pm_jsx["summary_min_ci_hi"],
            "delta_vs_original": float(pm_jsx["summary_min"]) - 0.6921,
            "weak_arm": "D_vs_A",
            "failed_ligand": "",
            "note": "Pocket B replaced (mTOR 4JSX); pocket A frozen 4L23. Not a PIK3CA swap.",
        },
        {
            "pair": "PIK3CA/PIK3CB",
            "pik3ca_receptor": "4L23",
            "kept_pocket_B": "2WXF",
            "n_attempted": 100,
            "n_successful": 99,
            "n_failed": 1,
            "n_dual": 28,
            "n_A_only": 27,
            "n_B_only": 28,
            "auroc_D_vs_A": 0.6905,
            "auroc_D_vs_B": 0.5,
            "summary_min": 0.5,
            "ci_lo": 0.3468,
            "ci_hi": 0.648,
            "delta_vs_original": 0.0,
            "weak_arm": "D_vs_B",
            "failed_ligand": "PAB_034 A_only timeout_900s_torsdof=23 on 4L23; 2WXF success. Not a label filter.",
            "note": "Original frozen panel. Exhaustiveness 8. Same 99-ligand AUROC set as Table 2.",
        },
        {
            "pair": "PIK3CA/PIK3CB",
            "pik3ca_receptor": "4JPS",
            "kept_pocket_B": "2WXF",
            "n_attempted": 100,
            "n_successful": 99,
            "n_failed": 1,
            "n_dual": pab_jps["n_dual"],
            "n_A_only": pab_jps["n_A_only"],
            "n_B_only": pab_jps["n_B_only"],
            "auroc_D_vs_A": pab_jps["auroc_D_vs_A"],
            "auroc_D_vs_B": pab_jps["auroc_D_vs_B"],
            "summary_min": pab_jps["summary_min"],
            "ci_lo": pab_jps["summary_min_ci_lo"],
            "ci_hi": pab_jps["summary_min_ci_hi"],
            "delta_vs_original": float(pab_jps["summary_min"]) - 0.5,
            "weak_arm": "D_vs_A",
            "failed_ligand": "PAB_034 A_only timeout_600s (668.6 s) on 4JPS; 2WXF frozen success. Same ligand as original 4L23 timeout.",
            "note": "Pocket A replaced; pocket B frozen 2WXF. Exhaustiveness 8. Weak arm switches from D/B to D/A. Use deposited CSV CIs.",
        },
        {
            "pair": "PIK3CA/PIK3CB",
            "pik3ca_receptor": "5DXT",
            "kept_pocket_B": "2WXF",
            "n_attempted": 100,
            "n_successful": 99,
            "n_failed": 1,
            "n_dual": pab_dxt["n_dual"],
            "n_A_only": pab_dxt["n_A_only"],
            "n_B_only": pab_dxt["n_B_only"],
            "auroc_D_vs_A": pab_dxt["auroc_D_vs_A"],
            "auroc_D_vs_B": pab_dxt["auroc_D_vs_B"],
            "summary_min": pab_dxt["summary_min"],
            "ci_lo": pab_dxt["summary_min_ci_lo"],
            "ci_hi": pab_dxt["summary_min_ci_hi"],
            "delta_vs_original": float(pab_dxt["summary_min"]) - 0.5,
            "weak_arm": "D_vs_B",
            "failed_ligand": "PAB_034 A_only timeout_600s (665.0 s) on 5DXT; 2WXF frozen success. Same ligand as original 4L23 timeout.",
            "note": "Pocket A replaced; pocket B frozen 2WXF. Exhaustiveness 8. D/B 0.6849 ≈ D/A 0.6905. Use deposited CSV CIs.",
        },
    ]
    two_path = SR / "receptor_realization_two_pair_v1.csv"
    with two_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(two_pair[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(two_pair)

    for r in two_pair:
        rows.append(
            row(
                block="receptor_realization_B5",
                manuscript_table="Table S9 / Table S30",
                pair=r["pair"],
                setting=f"{r['pik3ca_receptor']}/{r['kept_pocket_B']}",
                metric="summary_min",
                value=r["summary_min"],
                ci_lo=r["ci_lo"],
                ci_hi=r["ci_hi"],
                n_scored=r["n_successful"],
                n_dual=r["n_dual"],
                n_A_only=r["n_A_only"],
                n_B_only=r["n_B_only"],
                note=r["note"] + " " + r["failed_ligand"],
                source_file="data/jcim_structure_robust_v0/tables/receptor_realization_two_pair_v1.csv",
            )
        )

    # --- docking census ---
    for r in census:
        rows.append(
            row(
                block="docking_census",
                manuscript_table="Table S27",
                pair=r["pair"],
                setting=r["set"],
                metric="n_success_both_ends",
                value=r["n_success_both_ends"],
                n_scored=r["n_attempted"],
                note=r["note"] + (
                    " failed=" + r["failed_ligands"] if r.get("failed_ligands") else ""
                ),
                source_file="data/jcim_novelty_v0/tables/docking_failure_census_v1.csv",
            )
        )

    # --- descriptors ---
    for r in desc:
        rows.append(
            row(
                block="descriptor_reference",
                manuscript_table="Table 2 / Table S28",
                pair=r["pair"],
                setting="prespecified_descriptor",
                metric="best_single_descriptor_summary_min",
                value=r["best_single_descriptor_summary_min"],
                note="See Table S28 for all four descriptors; highest is a reference, not a confirmatory test.",
                source_file="data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv",
            )
        )

    det = _read(TAB / "detectable_effect_simulation_v1.csv")
    for r in det:
        if r["contrast"] != "summary_min" or r["true_auroc"] not in ("0.60", "0.70", "0.75"):
            continue
        rows.append(
            row(
                block="detectable_effect",
                manuscript_table="Table S31",
                pair=r["pair"],
                setting=f"true_auroc_{r['true_auroc']}",
                metric="p_summary_min_ci_excludes_0.5",
                value=r["p_ci_excludes_0p5"],
                n_dual=r["n_pos"],
                note="Binormal simulation; N_MC=1000; ligand bootstrap B=2000; not observed power. n_neg column in source CSV is n_A/n_B.",
                source_file="data/jcim_novelty_v0/tables/detectable_effect_simulation_v1.csv",
            )
        )

    out_path = TAB / "MASTER_RESULTS_TABLE.csv"
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {out_path} ({len(rows)} rows)")
    print(f"wrote {agree_path}")
    print(f"wrote {two_path}")


if __name__ == "__main__":
    main()
