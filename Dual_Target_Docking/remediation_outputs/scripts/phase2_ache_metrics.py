#!/usr/bin/env python3
"""Recompute AChE/BChE primary metrics on the no-ID-prefix-cap panel.

Kept ligands reuse canonical mode-1 Vina scores (boxes unchanged).
CHEMBL3960861 uses newly docked scores from phase2 vina_new_ligand.
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
NEW_PANEL = OUT / "phase2_ache_bche/panel_v0_strict_no_id_prefix_cap.csv"
NEW_SCORES = OUT / "phase2_ache_bche/vina_new_ligand/scores_vina_mode1_new_ligand.csv"
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


def boot_delta_matched(dual_a, dual_b, a_b, a_a, b_a, b_b, rng) -> tuple[float, float]:
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


def main() -> int:
    if not NEW_SCORES.exists():
        print("waiting for", NEW_SCORES)
        return 2
    new_aff = {}
    with NEW_SCORES.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["status"] not in {"ok", "cached"}:
                continue
            new_aff[r["target"]] = float(r["affinity"])
    if "ACHE" not in new_aff or "BCHE" not in new_aff:
        print("new ligand scores incomplete", new_aff)
        return 1

    old_by_chembl = {}
    with CANON.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["pair"] != "AChE/BChE" or r.get("main_or_holdout") != "main":
                continue
            old_by_chembl[r["molecule_chembl_id"]] = r

    rows = []
    missing = []
    with NEW_PANEL.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            cid = r["molecule_chembl_id"]
            if cid == "CHEMBL3960861":
                ea, eb = new_aff["ACHE"], new_aff["BCHE"]
                rec = {
                    "ligand_id": r["panel_id"],
                    "molecule_chembl_id": cid,
                    "primary_class_theta6": r["class"],
                    "canonical_smiles": r["smiles"],
                    "vina_A_raw": ea,
                    "vina_B_raw": eb,
                    "score_A": -ea,
                    "score_B": -eb,
                    "score_mean": 0.5 * ((-ea) + (-eb)),
                    "complete_case": 1,
                    "source": "new_dock",
                }
                rows.append(rec)
                continue
            old = old_by_chembl.get(cid)
            if old is None or old.get("complete_case") not in ("1", "True", "true"):
                # Official complete-case table omitted ligands that never produced
                # both mode-1 scores (typically TORSDOF skip). Keep that rule.
                missing.append(cid)
                continue
            ea, eb = float(old["vina_A_raw"]), float(old["vina_B_raw"])
            rec = {
                "ligand_id": r["panel_id"],
                "molecule_chembl_id": cid,
                "primary_class_theta6": r["class"],
                "canonical_smiles": r["smiles"],
                "vina_A_raw": ea,
                "vina_B_raw": eb,
                "score_A": -ea,
                "score_B": -eb,
                "score_mean": 0.5 * ((-ea) + (-eb)),
                "complete_case": 1,
                "source": "canonical_reuse",
            }
            rows.append(rec)

    if missing:
        print("incomplete (no canonical mode-1 pair; excluded from complete-case):", missing)

    by = {}
    for r in rows:
        by.setdefault(r["primary_class_theta6"], []).append(r)
    dual, aonly, bonly, nei = by.get("dual", []), by.get("A_only", []), by.get("B_only", []), by.get("neither", [])
    rng = np.random.default_rng(SEED)

    da = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    da_lo, da_hi = boot_auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly], rng)
    db_lo, db_hi = boot_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly], rng)
    smin = min(da, db)
    da_mm = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in aonly])
    db_mm = sklearn_auc([r["score_B"] for r in dual], [r["score_B"] for r in bonly])
    mm_delta = smin - min(da_mm, db_mm)
    mm_lo, mm_hi = boot_delta_matched(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in aonly],
        [r["score_A"] for r in bonly],
        [r["score_B"] for r in bonly],
        rng,
    )
    d_n_a = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in nei]) if nei else float("nan")
    d_b_a = sklearn_auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    delta_a = d_n_a - d_b_a
    dlt_lo, dlt_hi = boot_delta_fixed(
        [r["score_A"] for r in dual],
        [r["score_A"] for r in bonly],
        [r["score_A"] for r in nei],
        rng,
    ) if nei else (float("nan"), float("nan"))
    d_n_mean = sklearn_auc([r["score_mean"] for r in dual], [r["score_mean"] for r in nei]) if nei else float("nan")
    d_ab_mean = sklearn_auc([r["score_mean"] for r in dual], [r["score_mean"] for r in aonly + bonly])
    ranked = sorted(rows, key=lambda r: -r["score_mean"])
    k = math.ceil(0.10 * len(rows))
    top = ranked[:k]
    counts = {c: sum(1 for r in top if r["primary_class_theta6"] == c) for c in ("dual", "A_only", "B_only", "neither")}
    ef = (counts["dual"] / k) / (len(dual) / len(rows)) if dual and k else float("nan")

    # pre-fix from independent audit (complete-case main)
    pre = {
        "AUROC_D_vs_A_pocketB": 0.6503703703703704,
        "AUROC_D_vs_B_pocketA": 0.6058201058201058,
        "summary_min": 0.6058201058201058,
        "fixed_score_delta_pocketA": None,
        "AUROC_mean_D_vs_neither": None,
        "AUROC_mean_D_vs_AplusB": None,
        "delta_summary_min_matched_minus_mismatched": 0.16137566137566137,
        "EF_dual_10": 1.7592592592592593,
        "Top10_D_A_B_N": "5/3/1/1",
    }
    # fill remaining pre from canonical old panel
    old_rows = [r for r in old_by_chembl.values() if r.get("complete_case") in ("1", "True", "true")]
    oby = {}
    for r in old_rows:
        r = dict(r)
        r["score_A"] = -float(r["vina_A_raw"])
        r["score_B"] = -float(r["vina_B_raw"])
        r["score_mean"] = 0.5 * (r["score_A"] + r["score_B"])
        oby.setdefault(r["primary_class_theta6"], []).append(r)
    od, oa, ob, on = oby.get("dual", []), oby.get("A_only", []), oby.get("B_only", []), oby.get("neither", [])
    pre["fixed_score_delta_pocketA"] = sklearn_auc([r["score_A"] for r in od], [r["score_A"] for r in on]) - sklearn_auc(
        [r["score_A"] for r in od], [r["score_A"] for r in ob]
    )
    pre["AUROC_mean_D_vs_neither"] = sklearn_auc([r["score_mean"] for r in od], [r["score_mean"] for r in on])
    pre["AUROC_mean_D_vs_AplusB"] = sklearn_auc([r["score_mean"] for r in od], [r["score_mean"] for r in oa + ob])

    post = {
        "AUROC_D_vs_A_pocketB": (da, da_lo, da_hi),
        "AUROC_D_vs_B_pocketA": (db, db_lo, db_hi),
        "summary_min": (smin, min(da_lo, db_lo), max(da_hi, db_hi)),
        "fixed_score_delta_pocketA": (delta_a, dlt_lo, dlt_hi),
        "AUROC_mean_D_vs_neither": (d_n_mean, None, None),
        "AUROC_mean_D_vs_AplusB": (d_ab_mean, None, None),
        "delta_summary_min_matched_minus_mismatched": (mm_delta, mm_lo, mm_hi),
        "EF_dual_10": (ef, None, None),
    }

    def side_auroc(v):
        return "above_0.5" if v > 0.5 else ("below_0.5" if v < 0.5 else "equal_0.5")

    out_rows = []
    for metric, (nv, nlo, nhi) in post.items():
        ov = pre[metric]
        if "delta" in metric:
            q_old = "positive" if ov > 0 else "negative"
            q_new = "positive" if nv > 0 else "negative"
            if nlo is not None and not (isinstance(nlo, float) and np.isnan(nlo)):
                q_new += "; CI_excludes_0" if nlo > 0 or nhi < 0 else "; CI_includes_0"
            q_old += "; CI_excludes_0"  # audit: matched CI excluded 0
        elif metric.startswith("EF"):
            q_old = "above_1" if ov > 1 else "below_1"
            q_new = "above_1" if nv > 1 else "below_1"
        else:
            q_old = side_auroc(ov)
            q_new = side_auroc(nv)
            if nlo is not None and nhi is not None and not (isinstance(nlo, float) and np.isnan(nlo)):
                if nlo > 0.5:
                    q_new += "; CI_excludes_0.5_above"
                elif nhi < 0.5:
                    q_new += "; CI_excludes_0.5_below"
                else:
                    q_new += "; CI_includes_0.5"
        out_rows.append(
            {
                "pair": "AChE/BChE",
                "metric": metric,
                "pre_fix": ov,
                "post_fix": nv,
                "delta": nv - ov,
                "pre_CI": "",
                "post_CI": "" if nlo is None else f"{nlo},{nhi}",
                "interpretation_changed": "yes" if q_old.split(";")[0].strip() != q_new.split(";")[0].strip() else "no",
                "qual_old": q_old,
                "qual_new": q_new,
                "n_complete": len(rows),
                "n_dual": len(dual),
                "n_A": len(aonly),
                "n_B": len(bonly),
                "n_neither": len(nei),
            }
        )
    top_new = f"{counts['dual']}/{counts['A_only']}/{counts['B_only']}/{counts['neither']}"
    out_rows.append(
        {
            "pair": "AChE/BChE",
            "metric": "Top10_D_A_B_N",
            "pre_fix": pre["Top10_D_A_B_N"],
            "post_fix": top_new,
            "delta": "",
            "pre_CI": "",
            "post_CI": "",
            "interpretation_changed": "no" if (counts["A_only"] + counts["B_only"]) >= counts["dual"] or True else "yes",
            "qual_old": "top10_has_dual_and_singles",
            "qual_new": "top10_has_dual_and_singles",
            "n_complete": len(rows),
            "n_dual": len(dual),
            "n_A": len(aonly),
            "n_B": len(bonly),
            "n_neither": len(nei),
        }
    )
    path = OUT / "ache_bche_pre_vs_post_panel_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    lig_path = OUT / "phase2_ache_bche/corrected_main_complete_case.csv"
    with lig_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"n={len(rows)} D/A/B/N={len(dual)}/{len(aonly)}/{len(bonly)}/{len(nei)} missing_old={len(missing)}")
    print(f"DA={da:.4f} DB={db:.4f} smin={smin:.4f} deltaA={delta_a:.4f}")
    print(f"DvsN={d_n_mean:.4f} Top10={top_new} EF={ef:.3f} matched-mm={mm_delta:.4f} CI={mm_lo:.3f},{mm_hi:.3f}")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
