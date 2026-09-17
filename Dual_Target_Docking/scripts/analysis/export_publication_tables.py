#!/usr/bin/env python3
"""Write derived publication CSVs from results/canonical (no new estimands)."""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    direc = read_csv(CANON / "primary_directional_auroc.csv")
    da = {(r["pair"], r["estimand"]): r for r in direc}
    labels = read_csv(CANON / "label_aggregation_sensitivity.csv")
    ranking = read_csv(CANON / "top10_operating_points.csv")
    matched = read_csv(CANON / "matched_mismatched_pocket.csv")
    fixed = read_csv(CANON / "fixed_score_negative_class_delta.csv")
    inc = read_csv(CANON / "ecfp4_incremental_information.csv")
    desc = {r["pair"]: r for r in read_csv(CANON / "descriptor_baselines.csv")}
    hold = {r["pair"]: r for r in read_csv(CANON / "holdout_metrics.csv")}
    cluster = read_csv(CANON / "cluster_bootstrap_sensitivity.csv")
    rmsd = read_csv(CANON / "cognate_rmsd.csv")

    # Unified threshold: all eight pairs, all label rules from canonical sensitivity.
    unified = []
    for r in labels:
        unified.append(
            {
                "pair": r["pair"],
                "label_rule": r["label_rule"],
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "n_neither_excluded": "",
                "underpowered": int(min(int(r["n_dual"]), int(r["n_A_only"]), int(r["n_B_only"])) < 8),
                "pocket_matched_summary_min": r["summary_min"],
                "auroc_D_vs_A": r["auroc_D_vs_A"],
                "auroc_D_vs_B": r["auroc_D_vs_B"],
                "ci_lo": r["ci_lo"],
                "ci_hi": r["ci_hi"],
                "bootstrap": r["bootstrap"],
                "seed": r["seed"],
            }
        )
    write_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv", unified)

    five = []
    for pair in ("F2/F10", "JAK1/TYK2", "JAK1/JAK2", "PPARG/PPARA", "PPARA/PPARD"):
        s = smin[pair]
        d = desc[pair]
        five.append(
            {
                "pair": pair,
                "label_rule": "theta_6.0",
                "bootstrap": "class_stratified_shared_dual",
                "n_dual": s["n_dual"],
                "n_A_only": s["n_A_only"],
                "n_B_only": s["n_B_only"],
                "summary_min": s["summary_min"],
                "ci_lo": s["ci_lo"],
                "ci_hi": s["ci_hi"],
                "auroc_D_vs_A_pocketB": s["auroc_D_vs_A_pocketB"],
                "auroc_D_vs_B_pocketA": s["auroc_D_vs_B_pocketA"],
                "best_single_descriptor": d["best_descriptor"],
                "best_single_descriptor_summary_min": d["best_descriptor_summary_min"],
                "note": "Derived from results/canonical; scheme-B primary",
            }
        )
    write_csv(
        ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv",
        five,
    )

    rank_out = []
    for r in ranking:
        rank_out.append(
            {
                **r,
                "ranking_rule": r["ranking_rule"],
                "tie_rule": r["tie_rule"],
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv", rank_out)

    wp = []
    for r in matched:
        wp.append(
            {
                "pair": r["pair"],
                "set": r["set"],
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "n_boot_ok": 2000,
                "matched_D_vs_A": r["matched_auroc_D_vs_A"],
                "matched_D_vs_B": r["matched_auroc_D_vs_B"],
                "matched_summary_min": r["matched_summary_min"],
                "wrong_D_vs_A": r["mismatched_auroc_D_vs_A"],
                "wrong_D_vs_B": r["mismatched_auroc_D_vs_B"],
                "wrong_summary_min": r["mismatched_summary_min"],
                "delta_matched_minus_wrong": r["delta"],
                "delta_ci_lo": r["delta_ci_lo"],
                "delta_ci_hi": r["delta_ci_hi"],
                "ci_excludes_zero": "True" if str(r["ci_excludes_zero"]) in {"1", "True"} else "False",
                "weaker_arm_switched": r["weaker_arm_switched"],
            }
        )
    write_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv", wp)

    write_csv(ROOT / "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv", fixed)
    write_csv(ROOT / "data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv", cluster)
    write_csv(ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv", rmsd)

    inc_out = []
    for r in inc:
        for model, col in (
            ("ECFP4", "cv_auroc_ECFP4"),
            ("docking", "cv_auroc_docking"),
            ("ECFP4+docking", "cv_auroc_ECFP4_docking"),
        ):
            inc_out.append(
                {
                    "pair": r["pair"],
                    "contrast": r["contrast"],
                    "model": model,
                    "n": r["n"],
                    "n_pos": r["n_pos"],
                    "n_neg": r["n_neg"],
                    "n_scaffolds": r["n_scaffolds"],
                    "cv_auroc": r[col],
                    "rank_auroc_docking": r["rank_auroc_docking"],
                    "delta_ECFP4_plus_docking_minus_ECFP4": r["delta_ECFP4_plus_docking_minus_ECFP4"],
                    "note": r["note"],
                }
            )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/incremental_information_v1.csv", inc_out)

    # Membership derived from master (corrected EGFR scores).
    master = read_csv(CANON / "current_score_master.csv")
    memb = []
    for r in master:
        if r["analysis_set"] != "main" or r["complete_case"] not in {"1", "True"}:
            continue
        memb.append(
            {
                "pair": r["pair"],
                "ligand": r["ligand_id"],
                "cls": r["primary_class_theta6"],
                "pA": r["pA"],
                "pB": r["pB"],
                "score_A": r["score_A"],
                "score_B": r["score_B"],
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv", memb)

    hold_out = []
    for pair, r in hold.items():
        hold_out.append(r)
    write_csv(ROOT / "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv", hold_out)

    print("exported publication CSVs from results/canonical")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
