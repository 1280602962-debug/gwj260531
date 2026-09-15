#!/usr/bin/env python3
"""Recompute EGFR/HER2 primary metrics from corrected-box Vina mode-1 scores.

Labels, membership, and ligand IDs come from the pre-fix canonical table.
Only affinities change. Does not peek at AUROC to choose ligands.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs"
CANON = ROOT / "audit_outputs/canonical_ligand_table.csv"
NEW_SCORES = OUT / "phase1_vina/scores_vina_mode1_corrected_box.csv"
N_BOOT = 2000
SEED = 20260729


def sklearn_auc(pos, neg) -> float:
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def boot_auc(pos, neg, rng) -> tuple[float, float]:
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    vals = []
    for _ in range(N_BOOT):
        p = rng.choice(pos, size=len(pos), replace=True)
        n = rng.choice(neg, size=len(neg), replace=True)
        try:
            vals.append(sklearn_auc(p, n))
        except ValueError:
            continue
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def boot_summary_min(dual_a, dual_b, a_only_b, b_only_a, rng) -> tuple[float, float, int]:
    dual_a = np.asarray(dual_a, float)
    dual_b = np.asarray(dual_b, float)
    a_only_b = np.asarray(a_only_b, float)
    b_only_a = np.asarray(b_only_a, float)
    n_d = len(dual_a)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.choice(n_d, size=n_d, replace=True)
        a_idx = rng.choice(len(a_only_b), size=len(a_only_b), replace=True)
        b_idx = rng.choice(len(b_only_a), size=len(b_only_a), replace=True)
        da = sklearn_auc(dual_b[idx], a_only_b[a_idx])
        db = sklearn_auc(dual_a[idx], b_only_a[b_idx])
        vals.append(min(da, db))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), len(vals)


def boot_delta_fixed(dual, sel, nei, rng) -> tuple[float, float]:
    dual = np.asarray(dual, float)
    sel = np.asarray(sel, float)
    nei = np.asarray(nei, float)
    vals = []
    for _ in range(N_BOOT):
        d = rng.choice(dual, size=len(dual), replace=True)
        s = rng.choice(sel, size=len(sel), replace=True)
        n = rng.choice(nei, size=len(nei), replace=True)
        vals.append(sklearn_auc(d, n) - sklearn_auc(d, s))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def interpret_auroc(x: float, lo: float | None, hi: float | None) -> str:
    side = "above_0.5" if x > 0.5 else ("below_0.5" if x < 0.5 else "equal_0.5")
    if lo is None or hi is None or (isinstance(lo, str) and lo == ""):
        return side
    lo, hi = float(lo), float(hi)
    if lo > 0.5:
        ci = "CI_excludes_0.5_above"
    elif hi < 0.5:
        ci = "CI_excludes_0.5_below"
    else:
        ci = "CI_includes_0.5"
    return f"{side}; {ci}"


def main() -> int:
    if not NEW_SCORES.exists():
        print(f"waiting for {NEW_SCORES}")
        return 2
    scores = {}
    with NEW_SCORES.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["status"] not in {"ok", "cached"}:
                continue
            scores.setdefault(r["ligand"], {})[r["target"]] = float(r["affinity"])
    rows = []
    with CANON.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["pair"] != "EGFR/HER2":
                continue
            if r.get("main_or_holdout") != "main":
                continue
            if r.get("complete_case") not in ("1", "True", "true"):
                continue
            lig = r["ligand_id"]
            if lig not in scores or "3POZ" not in scores[lig] or "3RCD" not in scores[lig]:
                continue
            ea, eb = scores[lig]["3POZ"], scores[lig]["3RCD"]
            rec = dict(r)
            rec["vina_A_raw_new"] = ea
            rec["vina_B_raw_new"] = eb
            rec["score_A_new"] = -ea
            rec["score_B_new"] = -eb
            rec["score_mean_new"] = 0.5 * ((-ea) + (-eb))
            rows.append(rec)

    by = {}
    for r in rows:
        by.setdefault(r["primary_class_theta6"], []).append(r)
    dual, aonly, bonly, nei = by.get("dual", []), by.get("A_only", []), by.get("B_only", []), by.get("neither", [])
    rng = np.random.default_rng(SEED)

    da = sklearn_auc([r["score_B_new"] for r in dual], [r["score_B_new"] for r in aonly])
    db = sklearn_auc([r["score_A_new"] for r in dual], [r["score_A_new"] for r in bonly])
    da_lo, da_hi = boot_auc([r["score_B_new"] for r in dual], [r["score_B_new"] for r in aonly], rng)
    db_lo, db_hi = boot_auc([r["score_A_new"] for r in dual], [r["score_A_new"] for r in bonly], rng)
    smin = min(da, db)
    smin_lo, smin_hi, nvalid = boot_summary_min(
        [r["score_A_new"] for r in dual],
        [r["score_B_new"] for r in dual],
        [r["score_B_new"] for r in aonly],
        [r["score_A_new"] for r in bonly],
        rng,
    )
    d_n_a = sklearn_auc([r["score_A_new"] for r in dual], [r["score_A_new"] for r in nei])
    d_b_a = sklearn_auc([r["score_A_new"] for r in dual], [r["score_A_new"] for r in bonly])
    delta_a = d_n_a - d_b_a
    dlt_lo, dlt_hi = boot_delta_fixed(
        [r["score_A_new"] for r in dual],
        [r["score_A_new"] for r in bonly],
        [r["score_A_new"] for r in nei],
        rng,
    )
    d_n_mean = sklearn_auc([r["score_mean_new"] for r in dual], [r["score_mean_new"] for r in nei])
    d_ab_mean = sklearn_auc(
        [r["score_mean_new"] for r in dual],
        [r["score_mean_new"] for r in aonly + bonly],
    )
    # mismatched: swap S_A/S_B
    da_mm = sklearn_auc([r["score_A_new"] for r in dual], [r["score_A_new"] for r in aonly])
    db_mm = sklearn_auc([r["score_B_new"] for r in dual], [r["score_B_new"] for r in bonly])
    smin_mm = min(da_mm, db_mm)
    mm_delta = smin - smin_mm

    ranked = sorted(rows, key=lambda r: -r["score_mean_new"])
    k = math.ceil(0.10 * len(rows))
    top = ranked[:k]
    counts = {c: sum(1 for r in top if r["primary_class_theta6"] == c) for c in ("dual", "A_only", "B_only", "neither")}
    ef = (counts["dual"] / k) / (len(dual) / len(rows)) if dual and k else float("nan")

    pre = {
        "AUROC_D_vs_A_pocketB": (0.6663533834586466, 0.5237584218584672, 0.793095238095238),
        "AUROC_D_vs_B_pocketA": (0.4296875, 0.2817924317062248, 0.579011972770499),
        "summary_min": (0.4296875, 0.2817924317062248, 0.5775046404364821),
        "fixed_score_delta_pocketA": (0.3783482142857143, 0.20161830357142851, 0.5524553571428571),
        "AUROC_mean_D_vs_neither": (0.7559523809523809, 0.5625, 0.9197),
        "AUROC_mean_D_vs_AplusB": (0.5163265306122449, None, None),
        "delta_summary_min_matched_minus_mismatched": (0.16964285714285715, 0.059145422149122794, 0.2877330714894653),
        "EF_dual_10": (0.3571428571428572, None, None),
    }
    post = {
        "AUROC_D_vs_A_pocketB": (da, da_lo, da_hi),
        "AUROC_D_vs_B_pocketA": (db, db_lo, db_hi),
        "summary_min": (smin, smin_lo, smin_hi),
        "fixed_score_delta_pocketA": (delta_a, dlt_lo, dlt_hi),
        "AUROC_mean_D_vs_neither": (d_n_mean, None, None),
        "AUROC_mean_D_vs_AplusB": (d_ab_mean, None, None),
        "delta_summary_min_matched_minus_mismatched": (mm_delta, None, None),
        "EF_dual_10": (ef, None, None),
    }

    def qual_auroc(name, val, lo, hi):
        if "delta" in name.lower() or name.startswith("EF"):
            if "EF" in name:
                return "above_1" if val > 1 else ("below_1" if val < 1 else "equal_1")
            if lo is not None and hi is not None:
                if lo > 0:
                    return "positive; CI_excludes_0"
                if hi < 0:
                    return "negative; CI_excludes_0"
                return ("positive" if val > 0 else "negative") + "; CI_includes_0"
            return "positive" if val > 0 else "negative"
        return interpret_auroc(val, lo, hi)

    out_rows = []
    for metric, (nv, nlo, nhi) in post.items():
        ov, olo, ohi = pre[metric]
        q_old = qual_auroc(metric, ov, olo, ohi)
        q_new = qual_auroc(metric, nv, nlo, nhi)
        out_rows.append(
            {
                "metric": metric,
                "old": ov,
                "new": nv,
                "absolute_delta": nv - ov,
                "CI_old": "" if olo is None else f"{olo},{ohi}",
                "CI_new": "" if nlo is None else f"{nlo},{nhi}",
                "qualitative_interpretation_old": q_old,
                "qualitative_interpretation_new": q_new,
                "interpretation_changed": "yes" if q_old.split(";")[0].strip() != q_new.split(";")[0].strip() or ("CI_" in q_old and q_old != q_new) else "no",
                "n_dual": len(dual),
                "n_A": len(aonly),
                "n_B": len(bonly),
                "n_neither": len(nei),
                "n_complete": len(rows),
                "summary_min_boot_valid": nvalid,
            }
        )
    out_rows.append(
        {
            "metric": "Top10_D_A_B_N",
            "old": "1/5/5/0",
            "new": f"{counts['dual']}/{counts['A_only']}/{counts['B_only']}/{counts['neither']}",
            "absolute_delta": "",
            "CI_old": "",
            "CI_new": "",
            "qualitative_interpretation_old": "top10_majority_single_target",
            "qualitative_interpretation_new": (
                "top10_majority_single_target"
                if (counts["A_only"] + counts["B_only"]) > counts["dual"]
                else "top10_not_majority_single_target"
            ),
            "interpretation_changed": "no"
            if (counts["A_only"] + counts["B_only"]) > counts["dual"]
            else "yes",
            "n_dual": len(dual),
            "n_A": len(aonly),
            "n_B": len(bonly),
            "n_neither": len(nei),
            "n_complete": len(rows),
            "summary_min_boot_valid": nvalid,
        }
    )
    path = OUT / "egfr_her2_pre_vs_post_box_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"n={len(rows)} dual/A/B/N={len(dual)}/{len(aonly)}/{len(bonly)}/{len(nei)}")
    print(f"DA={da:.4f} DB={db:.4f} smin={smin:.4f} deltaA={delta_a:.4f}")
    print(f"DvsN={d_n_mean:.4f} Top10={counts} EF={ef:.3f} matched-mm={mm_delta:.4f}")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
