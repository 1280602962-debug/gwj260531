#!/usr/bin/env python3
"""Eight-pair post-fix master metrics.

Unchanged pairs are copied from the independent audit master table.
EGFR/HER2 uses canonical heavy-atom-box Vina scores.
AChE/BChE uses the no-ChEMBL-ID-prefix-cap panel (complete-case).
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs/canonical_tables"
AUDIT = ROOT / "audit_outputs/audit_master_metrics.csv"
EGFR_ABL = ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv"
EGFR_PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
ACHE_ABL = ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv"
ACHE_PANEL = ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict.csv"
N_BOOT = 2000
SEED = 20260729
UNCHANGED = (
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def sklearn_auc(pos, neg) -> float:
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def boot_auc(pos, neg, rng) -> tuple[float, float]:
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    vals = []
    for _ in range(N_BOOT):
        try:
            vals.append(
                sklearn_auc(rng.choice(pos, size=len(pos), replace=True), rng.choice(neg, size=len(neg), replace=True))
            )
        except ValueError:
            continue
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def boot_summary_min(dual_a, dual_b, a_only_b, b_only_a, rng):
    dual_a = np.asarray(dual_a, float)
    dual_b = np.asarray(dual_b, float)
    a_only_b = np.asarray(a_only_b, float)
    b_only_a = np.asarray(b_only_a, float)
    n_d = len(dual_a)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.choice(n_d, size=n_d, replace=True)
        ai = rng.choice(len(a_only_b), size=len(a_only_b), replace=True)
        bi = rng.choice(len(b_only_a), size=len(b_only_a), replace=True)
        vals.append(min(sklearn_auc(dual_b[idx], a_only_b[ai]), sklearn_auc(dual_a[idx], b_only_a[bi])))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def boot_delta_fixed(dual, sel, nei, rng):
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


def boot_matched(dual_a, dual_b, a_b, a_a, b_a, b_b, rng):
    dual_a = np.asarray(dual_a, float)
    dual_b = np.asarray(dual_b, float)
    a_b = np.asarray(a_b, float)
    a_a = np.asarray(a_a, float)
    b_a = np.asarray(b_a, float)
    b_b = np.asarray(b_b, float)
    n_d = len(dual_a)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.choice(n_d, size=n_d, replace=True)
        ai = rng.choice(len(a_b), size=len(a_b), replace=True)
        bi = rng.choice(len(b_a), size=len(b_a), replace=True)
        da = sklearn_auc(dual_b[idx], a_b[ai])
        db = sklearn_auc(dual_a[idx], b_a[bi])
        da_mm = sklearn_auc(dual_a[idx], a_a[ai])
        db_mm = sklearn_auc(dual_b[idx], b_b[bi])
        vals.append(min(da, db) - min(da_mm, db_mm))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def ci_vs_half(x, lo, hi) -> str:
    side = "above_0.5" if x > 0.5 else ("below_0.5" if x < 0.5 else "equal_0.5")
    if lo > 0.5:
        ci = "CI_excludes_0.5_above"
    elif hi < 0.5:
        ci = "CI_excludes_0.5_below"
    else:
        ci = "CI_includes_0.5"
    return f"{side}; {ci}"


def ci_vs_zero(x, lo, hi) -> str:
    sign = "positive" if x > 0 else ("negative" if x < 0 else "zero")
    if lo > 0:
        ci = "CI_excludes_0"
    elif hi < 0:
        ci = "CI_excludes_0"
    else:
        ci = "CI_includes_0"
    return f"{sign}; {ci}"


def load_pair(panel_path, abl_path, id_key, class_key, a_col, b_col, score_are_energy=True):
    panel = {r[id_key]: r for r in csv.DictReader(panel_path.open())}
    rows = []
    with abl_path.open() as fh:
        for r in csv.DictReader(fh):
            lig = r.get("ligand") or r.get("panel_id")
            try:
                ea = float(r[a_col])
                eb = float(r[b_col])
            except (TypeError, ValueError, KeyError):
                continue
            cls = r.get(class_key) or panel.get(lig, {}).get("class")
            if cls not in {"dual", "A_only", "B_only", "neither"}:
                continue
            sa = -ea if score_are_energy else ea
            sb = -eb if score_are_energy else eb
            rows.append(
                {
                    "ligand": lig,
                    "cls": cls,
                    "score_A": sa,
                    "score_B": sb,
                    "score_mean": 0.5 * (sa + sb),
                }
            )
    return rows


def metrics_for(pair, rows, rng):
    by = {}
    for r in rows:
        by.setdefault(r["cls"], []).append(r)
    dual, aonly, bonly, nei = by.get("dual", []), by.get("A_only", []), by.get("B_only", []), by.get("neither", [])
    da = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    da_lo, da_hi = boot_auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly], rng)
    db_lo, db_hi = boot_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly], rng)
    smin = min(da, db)
    smin_lo, smin_hi = boot_summary_min(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in bonly],
        rng,
    )
    d_n_a = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in nei])
    d_b_a = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    delta_a = d_n_a - d_b_a
    dlt_lo, dlt_hi = boot_delta_fixed(
        [r["score_A"] for r in dual],
        [r["score_A"] for r in bonly],
        [r["score_A"] for r in nei],
        rng,
    )
    d_n_b = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in nei])
    d_a_b = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    delta_b = d_n_b - d_a_b
    d_n_mean = sklearn_auc([r["score_mean"] for r in dual], [r["score_mean"] for r in nei])
    d_n_mean_lo, d_n_mean_hi = boot_auc([r["score_mean"] for r in dual], [r["score_mean"] for r in nei], rng)
    d_ab_mean = sklearn_auc([r["score_mean"] for r in dual], [r["score_mean"] for r in aonly + bonly])
    da_mm = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in aonly])
    db_mm = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in bonly])
    mm_delta = smin - min(da_mm, db_mm)
    mm_lo, mm_hi = boot_matched(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in aonly],
        [r["score_A"] for r in bonly],
        [r["score_B"] for r in bonly],
        rng,
    )
    ranked = sorted(rows, key=lambda r: (-r["score_mean"], r["ligand"]))
    k = math.ceil(0.10 * len(rows))
    top = ranked[:k]
    counts = {c: sum(1 for r in top if r["cls"] == c) for c in ("dual", "A_only", "B_only", "neither")}
    ef = (counts["dual"] / k) / (len(dual) / len(rows)) if dual and k else float("nan")
    recs = [
        ("AUROC_D_vs_A_pocketB", da, da_lo, da_hi, len(dual), len(aonly), "", ci_vs_half(da, da_lo, da_hi)),
        ("AUROC_D_vs_B_pocketA", db, db_lo, db_hi, len(dual), len(bonly), "", ci_vs_half(db, db_lo, db_hi)),
        ("summary_min", smin, smin_lo, smin_hi, len(dual), min(len(aonly), len(bonly)), "", ci_vs_half(smin, smin_lo, smin_hi)),
        ("fixed_score_D_vs_neither_pocketA", d_n_a, "", "", len(dual), len(nei), "", ""),
        ("fixed_score_D_vs_B_only_pocketA", d_b_a, "", "", len(dual), len(bonly), "", ""),
        ("fixed_score_delta_pocketA", delta_a, dlt_lo, dlt_hi, len(dual), len(bonly), len(nei), ci_vs_zero(delta_a, dlt_lo, dlt_hi)),
        ("fixed_score_delta_pocketB", delta_b, "", "", len(dual), len(aonly), len(nei), ""),
        ("AUROC_mean_D_vs_neither", d_n_mean, d_n_mean_lo, d_n_mean_hi, len(dual), len(nei), "", ci_vs_half(d_n_mean, d_n_mean_lo, d_n_mean_hi)),
        ("AUROC_mean_D_vs_AplusB", d_ab_mean, "", "", len(dual), len(aonly) + len(bonly), "", "above_0.5" if d_ab_mean > 0.5 else "below_0.5"),
        (
            "delta_summary_min_matched_minus_mismatched",
            mm_delta,
            mm_lo,
            mm_hi,
            len(dual),
            len(aonly) + len(bonly),
            "",
            ci_vs_zero(mm_delta, mm_lo, mm_hi),
        ),
        ("EF_dual_10", ef, "", "", len(dual), len(rows) - len(dual), "", "above_1" if ef > 1 else "below_1"),
        (
            "Top10_D_A_B_N",
            f"{counts['dual']}/{counts['A_only']}/{counts['B_only']}/{counts['neither']}",
            "",
            "",
            counts["dual"],
            k - counts["dual"],
            k,
            "",
        ),
    ]
    out = []
    for name, val, lo, hi, npos, nneg, n3, interp in recs:
        out.append(
            {
                "pair": pair,
                "analysis_set": "main",
                "estimand": name,
                "independent": val,
                "ci_lo": lo,
                "ci_hi": hi,
                "n_pos": npos,
                "n_neg": nneg,
                "n_neither": n3,
                "n_complete": len(rows),
                "interpretation": interp,
                "source": "post_fix recompute from canonical ablation",
            }
        )
    return out, {
        "pair": pair,
        "n_complete": len(rows),
        "n_dual": len(dual),
        "n_A": len(aonly),
        "n_B": len(bonly),
        "n_neither": len(nei),
        "AUROC_D_vs_A_pocketB": da,
        "AUROC_D_vs_B_pocketA": db,
        "summary_min": smin,
        "fixed_score_D_vs_neither_pocketA": d_n_a,
        "fixed_score_D_vs_B_only_pocketA": d_b_a,
        "fixed_score_delta_pocketA": delta_a,
        "fixed_score_delta_pocketA_ci": f"{dlt_lo},{dlt_hi}",
        "AUROC_mean_D_vs_neither": d_n_mean,
        "AUROC_mean_D_vs_AplusB": d_ab_mean,
        "matched_minus_mismatched": mm_delta,
        "matched_minus_mismatched_ci": f"{mm_lo},{mm_hi}",
        "matched_ci_includes_0": "yes" if mm_lo <= 0 <= mm_hi else "no",
        "Top10_D_A_B_N": f"{counts['dual']}/{counts['A_only']}/{counts['B_only']}/{counts['neither']}",
        "EF_dual_10": ef,
        "do_not_claim_matched_pocket_advantage": "yes" if mm_lo <= 0 <= mm_hi else "no",
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows = []
    with AUDIT.open() as fh:
        for r in csv.DictReader(fh):
            if r["pair"] in UNCHANGED:
                r["source"] = "pre_fix unchanged pair; copied from independent audit"
                rows.append(r)
    egfr_rows = load_pair(EGFR_PANEL, EGFR_ABL, "panel_id", "class", "3POZ_affinity", "3RCD_affinity")
    ache_rows = load_pair(ACHE_PANEL, ACHE_ABL, "panel_id", "class", "vina_ACHE", "vina_BCHE")
    egfr_m, egfr_sum = metrics_for("EGFR/HER2", egfr_rows, rng)
    ache_m, ache_sum = metrics_for("AChE/BChE", ache_rows, rng)
    rows.extend(egfr_m)
    rows.extend(ache_m)
    dest = OUT / "post_fix_master_metrics.csv"
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    cmp_path = OUT / "egfr_her2_corrected_comparison.csv"
    with cmp_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(egfr_sum.keys()))
        w.writeheader()
        w.writerow(egfr_sum)
        w.writerow(ache_sum)
    print("EGFR", egfr_sum)
    print("AChE", ache_sum)
    print("wrote", dest, "n_rows", len(rows))
    print("matched EGFR CI includes 0:", egfr_sum["matched_ci_includes_0"])
    print("matched AChE CI includes 0:", ache_sum["matched_ci_includes_0"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
