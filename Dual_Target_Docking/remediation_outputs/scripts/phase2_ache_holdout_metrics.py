#!/usr/bin/env python3
"""Recompute AChE/BChE unused-pool holdout metrics after panel rebuild.

53 kept ligands reuse deposited HOAB mode-1 scores.
7 new ligands use phase2 vina_holdout_new scores.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
NEW_HO = ROOT / "remediation_outputs/phase2_ache_bche/holdout_panel_HOAB_corrected.csv"
CMP = ROOT / "remediation_outputs/phase2_ache_bche/holdout_old_vs_corrected.csv"
OLD_SCORES = ROOT / "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv"
NEW_SCORES = ROOT / "remediation_outputs/phase2_ache_bche/vina_holdout_new/scores_vina_mode1_new_holdout.csv"
OUT = ROOT / "remediation_outputs/phase2_ache_bche/ache_holdout_pre_vs_post.csv"
N_BOOT = 2000
SEED = 20260729


def auc(pos, neg):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def boot_matched(dual_a, dual_b, a_b, a_a, b_a, b_b, rng):
    dual_a, dual_b = np.asarray(dual_a, float), np.asarray(dual_b, float)
    a_b, a_a = np.asarray(a_b, float), np.asarray(a_a, float)
    b_a, b_b = np.asarray(b_a, float), np.asarray(b_b, float)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.choice(len(dual_a), size=len(dual_a), replace=True)
        ai = rng.choice(len(a_b), size=len(a_b), replace=True)
        bi = rng.choice(len(b_a), size=len(b_a), replace=True)
        da = auc(dual_b[idx], a_b[ai])
        db = auc(dual_a[idx], b_a[bi])
        da_mm = auc(dual_a[idx], a_a[ai])
        db_mm = auc(dual_b[idx], b_b[bi])
        vals.append(min(da, db) - min(da_mm, db_mm))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def metrics(rows):
    by = {}
    for r in rows:
        by.setdefault(r["class"], []).append(r)
    dual, aonly, bonly = by.get("dual", []), by.get("A_only", []), by.get("B_only", [])
    da = auc([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db = auc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    smin = min(da, db)
    da_mm = auc([r["score_A"] for r in dual], [r["score_A"] for r in aonly])
    db_mm = auc([r["score_B"] for r in dual], [r["score_B"] for r in bonly])
    mm = smin - min(da_mm, db_mm)
    rng = np.random.default_rng(SEED)
    lo, hi = boot_matched(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in aonly],
        [r["score_A"] for r in bonly],
        [r["score_B"] for r in bonly],
        rng,
    )
    return {
        "n": len(rows),
        "n_dual": len(dual),
        "n_A": len(aonly),
        "n_B": len(bonly),
        "AUROC_D_vs_A_pocketB": da,
        "AUROC_D_vs_B_pocketA": db,
        "summary_min": smin,
        "matched_minus_mismatched": mm,
        "mm_ci_lo": lo,
        "mm_ci_hi": hi,
        "mm_excludes_0": lo > 0 or hi < 0,
    }


def main() -> int:
    if not NEW_SCORES.exists():
        print("waiting for", NEW_SCORES)
        return 2
    new_aff = {}
    with NEW_SCORES.open() as fh:
        for r in csv.DictReader(fh):
            if r["status"] not in {"ok", "cached"}:
                continue
            new_aff.setdefault(r["holdout_id"], {})[r["target"]] = float(r["affinity"])
    old_by_chembl = {}
    with OLD_SCORES.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("pair") != "AChE/BChE":
                continue
            if r.get("vina_A_raw") in ("", None) or r.get("vina_B_raw") in ("", None):
                continue
            old_by_chembl[r["chembl"]] = r
    kept = {r["ligand_id"] for r in csv.DictReader(CMP.open()) if r["old_holdout"] == "yes" and r["new_holdout"] == "yes"}
    rows = []
    missing = []
    with NEW_HO.open() as fh:
        for r in csv.DictReader(fh):
            cid = r["molecule_chembl_id"]
            hid = r["holdout_id"]
            if cid in kept:
                old = old_by_chembl.get(cid)
                if old is None:
                    missing.append(cid)
                    continue
                ea, eb = float(old["vina_A_raw"]), float(old["vina_B_raw"])
            else:
                sc = new_aff.get(hid, {})
                if "4EY7" not in sc or "4BDS" not in sc:
                    missing.append(cid)
                    continue
                ea, eb = sc["4EY7"], sc["4BDS"]
            rows.append(
                {
                    "holdout_id": hid,
                    "molecule_chembl_id": cid,
                    "class": r["class"],
                    "score_A": -ea,
                    "score_B": -eb,
                }
            )
    if missing:
        print("incomplete holdout ligands", missing)
        return 1
    post = metrics(rows)
    # pre-fix from old complete holdout scores
    pre_rows = []
    with OLD_SCORES.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("pair") != "AChE/BChE":
                continue
            if r.get("vina_A_raw") in ("", None) or r.get("vina_B_raw") in ("", None):
                continue
            pre_rows.append(
                {
                    "class": r["cls"],
                    "score_A": -float(r["vina_A_raw"]),
                    "score_B": -float(r["vina_B_raw"]),
                }
            )
    pre = metrics(pre_rows)
    out_rows = []
    for k in (
        "AUROC_D_vs_A_pocketB",
        "AUROC_D_vs_B_pocketA",
        "summary_min",
        "matched_minus_mismatched",
    ):
        pv, nv = pre[k], post[k]
        if k == "matched_minus_mismatched":
            q_old = ("positive" if pv > 0 else "negative") + (
                "; CI_excludes_0" if pre["mm_excludes_0"] else "; CI_includes_0"
            )
            q_new = ("positive" if nv > 0 else "negative") + (
                "; CI_excludes_0" if post["mm_excludes_0"] else "; CI_includes_0"
            )
        else:
            q_old = "above_0.5" if pv > 0.5 else "below_0.5"
            q_new = "above_0.5" if nv > 0.5 else "below_0.5"
        out_rows.append(
            {
                "pair": "AChE/BChE",
                "set": "unused_pool_holdout",
                "metric": k,
                "pre_fix": pv,
                "post_fix": nv,
                "delta": nv - pv,
                "pre_CI": f"{pre['mm_ci_lo']},{pre['mm_ci_hi']}" if k == "matched_minus_mismatched" else "",
                "post_CI": f"{post['mm_ci_lo']},{post['mm_ci_hi']}" if k == "matched_minus_mismatched" else "",
                "interpretation_changed": "yes" if q_old != q_new else "no",
                "qual_old": q_old,
                "qual_new": q_new,
                "n_pre": pre["n"],
                "n_post": post["n"],
                "n_dual_post": post["n_dual"],
                "n_A_post": post["n_A"],
                "n_B_post": post["n_B"],
            }
        )
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print("pre", pre)
    print("post", post)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
