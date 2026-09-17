#!/usr/bin/env python3
"""Build current_score_master and independently recompute primary statistics.

Does not import production AUROC/bootstrap helpers. Score sources:
  EGFR/HER2  — current panel120 ablation (canonical heavy-atom box)
  AChE/BChE  — current ache ablation (no-ID-prefix panel, includes AB_056)
  remaining six pairs — review_scored_membership (unchanged Track B / PM48)
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
AUDIT = ROOT / "remediation_outputs" / "freeze_audit"
N_BOOT = 2000
SEED = 20260729
ORDER = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fnum(x):
    if x in (None, ""):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def auroc_mw(pos, neg) -> float:
    """Mann–Whitney AUROC; higher score is more favorable. Independent of sklearn."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (pos.size * neg.size))


def auroc_rankdata(pos, neg) -> float:
    """Second independent AUROC via midranks (not the pairwise implementation)."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    x = np.concatenate([pos, neg])
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, float)
    i = 0
    while i < x.size:
        j = i
        while j + 1 < x.size and x[order[j + 1]] == x[order[i]]:
            j += 1
        mid = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = mid
        i = j + 1
    r_pos = ranks[: pos.size].sum()
    return float((r_pos - pos.size * (pos.size + 1) / 2.0) / (pos.size * neg.size))


def theta6(pa, pb) -> str:
    a, b = pa >= 6.0, pb >= 6.0
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def strict6555(pa, pb):
    if pa >= 6.5 and pb >= 6.5:
        return "dual"
    if pa >= 6.5 and pb <= 5.5:
        return "A_only"
    if pb >= 6.5 and pa <= 5.5:
        return "B_only"
    if pa <= 5.5 and pb <= 5.5:
        return "neither"
    return "gray"


def pct(vals) -> tuple[float, float]:
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def boot_two_class(pos, neg, rng) -> tuple[float, float]:
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    vals = []
    for _ in range(N_BOOT):
        p = pos[rng.integers(0, pos.size, pos.size)]
        n = neg[rng.integers(0, neg.size, neg.size)]
        vals.append(auroc_mw(p, n))
    return pct(vals)


def boot_smin_stratified(dual_sa, dual_sb, a_sb, b_sa, rng) -> tuple[float, float]:
    dual_sa = np.asarray(dual_sa, float)
    dual_sb = np.asarray(dual_sb, float)
    a_sb = np.asarray(a_sb, float)
    b_sa = np.asarray(b_sa, float)
    nd, na, nb = dual_sa.size, a_sb.size, b_sa.size
    vals = []
    for _ in range(N_BOOT):
        idd = rng.integers(0, nd, nd)
        da = auroc_mw(dual_sb[idd], a_sb[rng.integers(0, na, na)])
        db = auroc_mw(dual_sa[idd], b_sa[rng.integers(0, nb, nb)])
        vals.append(min(da, db))
    return pct(vals)


def boot_delta_shared_dual(dual, sel, nei, rng) -> tuple[float, float]:
    dual = np.asarray(dual, float)
    sel = np.asarray(sel, float)
    nei = np.asarray(nei, float)
    vals = []
    for _ in range(N_BOOT):
        d = dual[rng.integers(0, dual.size, dual.size)]
        vals.append(auroc_mw(d, nei[rng.integers(0, nei.size, nei.size)]) - auroc_mw(d, sel[rng.integers(0, sel.size, sel.size)]))
    return pct(vals)


def boot_matched_mismatch(dual_sa, dual_sb, a_sa, a_sb, b_sa, b_sb, rng):
    dual_sa = np.asarray(dual_sa, float)
    dual_sb = np.asarray(dual_sb, float)
    a_sa = np.asarray(a_sa, float)
    a_sb = np.asarray(a_sb, float)
    b_sa = np.asarray(b_sa, float)
    b_sb = np.asarray(b_sb, float)
    nd, na, nb = dual_sa.size, a_sa.size, b_sa.size
    vals = []
    for _ in range(N_BOOT):
        idd = rng.integers(0, nd, nd)
        ida = rng.integers(0, na, na)
        idb = rng.integers(0, nb, nb)
        m_da = auroc_mw(dual_sb[idd], a_sb[ida])
        m_db = auroc_mw(dual_sa[idd], b_sa[idb])
        w_da = auroc_mw(dual_sa[idd], a_sa[ida])
        w_db = auroc_mw(dual_sb[idd], b_sb[idb])
        vals.append(min(m_da, m_db) - min(w_da, w_db))
    return pct(vals)


def build_master() -> list[dict]:
    memb = read_csv(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv")
    egfr_abl = {r["ligand"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")}
    egfr_panel = {r["panel_id"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")}
    ache_abl = {r["ligand"]: r for r in read_csv(ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv")}
    rows = []
    for r in memb:
        pair, lid, cls = r["pair"], r["ligand"], r["cls"]
        pA, pB = fnum(r["pA"]), fnum(r["pB"])
        sa, sb = fnum(r["score_A"]), fnum(r["score_B"])
        chembl = ""
        score_src = "review_scored_membership_v1"
        if pair == "EGFR/HER2":
            a = egfr_abl[lid]
            sa, sb = float(a["vina_3POZ_hb"]), float(a["vina_3RCD_hb"])
            pA, pB = fnum(egfr_panel[lid]["pchembl_EGFR"]), fnum(egfr_panel[lid]["pchembl_HER2"])
            chembl = egfr_panel[lid]["molecule_chembl_id"]
            cls = a["class"]
            score_src = "egfr_her2_panel120_v0/ablation_ligand_scores.csv"
        elif pair == "AChE/BChE":
            a = ache_abl[lid]
            sa, sb = float(a["vina_ACHE_hb"]), float(a["vina_BCHE_hb"])
            pA, pB = fnum(a["pchembl_ACHE"]), fnum(a["pchembl_BCHE"])
            chembl = a.get("molecule_chembl_id", "")
            cls = a["class"]
            score_src = "ache_bche_panel_v0/ablation_ligand_scores.csv"
        pred = strict6555(pA, pB) if pair == "AChE/BChE" else theta6(pA, pB)
        rows.append(
            {
                "pair": pair,
                "ligand_id": lid,
                "molecule_chembl_id": chembl,
                "analysis_set": "main",
                "pA": pA,
                "pB": pB,
                "class": cls,
                "class_from_pchembl": pred,
                "score_A": sa,
                "score_B": sb,
                "score_mean": (sa + sb) / 2.0,
                "score_worst": min(sa, sb),
                "complete_case": 1,
                "score_source": score_src,
            }
        )
    # AChE ligands present in ablation but missing from membership (should be none after AB_056)
    have = {(r["pair"], r["ligand_id"]) for r in rows}
    for lid, a in ache_abl.items():
        if ("AChE/BChE", lid) in have:
            continue
        sa, sb = fnum(a.get("vina_ACHE_hb")), fnum(a.get("vina_BCHE_hb"))
        if sa is None or sb is None:
            continue
        pA, pB = fnum(a.get("pchembl_ACHE")), fnum(a.get("pchembl_BCHE"))
        rows.append(
            {
                "pair": "AChE/BChE",
                "ligand_id": lid,
                "molecule_chembl_id": a.get("molecule_chembl_id", ""),
                "analysis_set": "main",
                "pA": pA,
                "pB": pB,
                "class": a["class"],
                "class_from_pchembl": strict6555(pA, pB),
                "score_A": sa,
                "score_B": sb,
                "score_mean": (sa + sb) / 2.0,
                "score_worst": min(sa, sb),
                "complete_case": 1,
                "score_source": "ache_bche_panel_v0/ablation_ligand_scores.csv",
            }
        )
    rows.sort(key=lambda r: (ORDER.index(r["pair"]) if r["pair"] in ORDER else 99, r["ligand_id"]))
    return rows


def pair_recs(master, pair):
    return [r for r in master if r["pair"] == pair]


def compute_pair(recs, pair):
    by = defaultdict(list)
    for r in recs:
        by[r["class"]].append(r)
    dual, aonly, bonly, neither = by["dual"], by["A_only"], by["B_only"], by["neither"]
    da = auroc_mw([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db = auroc_mw([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    da2 = auroc_rankdata([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db2 = auroc_rankdata([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    smin = min(da, db)
    d_n = auroc_mw([r["score_mean"] for r in dual], [r["score_mean"] for r in neither]) if neither else float("nan")
    d_n_a = auroc_mw([r["score_A"] for r in dual], [r["score_A"] for r in neither]) if neither else float("nan")
    d_n_b = auroc_mw([r["score_B"] for r in dual], [r["score_B"] for r in neither]) if neither else float("nan")
    rng = np.random.default_rng(SEED)
    da_lo, da_hi = boot_two_class([r["score_B"] for r in dual], [r["score_B"] for r in aonly], rng)
    rng = np.random.default_rng(SEED)
    db_lo, db_hi = boot_two_class([r["score_A"] for r in dual], [r["score_A"] for r in bonly], rng)
    rng = np.random.default_rng(SEED)
    sm_lo, sm_hi = boot_smin_stratified(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in bonly],
        rng,
    )
    rng = np.random.default_rng(SEED)
    dn_lo = dn_hi = float("nan")
    if neither:
        dn_lo, dn_hi = boot_two_class([r["score_mean"] for r in dual], [r["score_mean"] for r in neither], rng)
    delta_a = d_n_a - db
    delta_b = d_n_b - da
    rng = np.random.default_rng(SEED)
    dA_lo, dA_hi = boot_delta_shared_dual(
        [r["score_A"] for r in dual], [r["score_A"] for r in bonly], [r["score_A"] for r in neither], rng
    ) if neither else (float("nan"), float("nan"))
    rng = np.random.default_rng(SEED)
    dB_lo, dB_hi = boot_delta_shared_dual(
        [r["score_B"] for r in dual], [r["score_B"] for r in aonly], [r["score_B"] for r in neither], rng
    ) if neither else (float("nan"), float("nan"))
    m_da = da
    m_db = db
    w_da = auroc_mw([r["score_A"] for r in dual], [r["score_A"] for r in aonly])
    w_db = auroc_mw([r["score_B"] for r in dual], [r["score_B"] for r in bonly])
    w_smin = min(w_da, w_db)
    mm_delta = smin - w_smin
    rng = np.random.default_rng(SEED)
    mm_lo, mm_hi = boot_matched_mismatch(
        [r["score_A"] for r in dual],
        [r["score_B"] for r in dual],
        [r["score_A"] for r in aonly],
        [r["score_B"] for r in aonly],
        [r["score_A"] for r in bonly],
        [r["score_B"] for r in bonly],
        rng,
    )
    ranked = sorted(recs, key=lambda r: (-r["score_mean"], r["ligand_id"]))
    k = max(1, math.ceil(0.10 * len(ranked)))
    top = ranked[:k]
    ct = Counter(r["class"] for r in top)
    n_dual = len(dual)
    ef = (ct["dual"] / k) / (n_dual / len(recs))
    duals = [r for r in recs if r["class"] == "dual"]
    thr = float(np.median([r["score_worst"] for r in duals]))
    use = [r for r in recs if r["class"] in ("dual", "A_only", "B_only")]
    retained = [r for r in use if r["score_worst"] >= thr]
    cr = Counter(r["class"] for r in retained)
    class_mismatch = sum(r["class"] != r["class_from_pchembl"] and r["class_from_pchembl"] != "gray" for r in recs)
    return {
        "pair": pair,
        "n_dual": len(dual),
        "n_A_only": len(aonly),
        "n_B_only": len(bonly),
        "n_neither": len(neither),
        "n_complete_case": len(recs),
        "n_missing_score_A": 0,
        "n_missing_score_B": 0,
        "n_class_pchembl_mismatch": class_mismatch,
        "auroc_D_vs_A_pocketB": da,
        "auroc_D_vs_A_pocketB_rank": da2,
        "auroc_D_vs_B_pocketA": db,
        "auroc_D_vs_B_pocketA_rank": db2,
        "summary_min": smin,
        "summary_min_ci_lo": sm_lo,
        "summary_min_ci_hi": sm_hi,
        "ci_lo_D_vs_A": da_lo,
        "ci_hi_D_vs_A": da_hi,
        "ci_lo_D_vs_B": db_lo,
        "ci_hi_D_vs_B": db_hi,
        "auroc_mean_D_vs_neither": d_n,
        "ci_lo_D_vs_neither": dn_lo,
        "ci_hi_D_vs_neither": dn_hi,
        "fixed_score_D_vs_selective_pocketA": db,
        "fixed_score_D_vs_neither_pocketA": d_n_a,
        "fixed_score_delta_pocketA": delta_a,
        "fixed_score_delta_pocketA_ci_lo": dA_lo,
        "fixed_score_delta_pocketA_ci_hi": dA_hi,
        "fixed_score_D_vs_selective_pocketB": da,
        "fixed_score_D_vs_neither_pocketB": d_n_b,
        "fixed_score_delta_pocketB": delta_b,
        "fixed_score_delta_pocketB_ci_lo": dB_lo,
        "fixed_score_delta_pocketB_ci_hi": dB_hi,
        "matched_smin": smin,
        "mismatched_smin": w_smin,
        "matched_minus_mismatched": mm_delta,
        "matched_minus_mismatched_ci_lo": mm_lo,
        "matched_minus_mismatched_ci_hi": mm_hi,
        "matched_ci_excludes_0": bool(mm_lo > 0 or mm_hi < 0),
        "top_k": k,
        "top_dual": ct["dual"],
        "top_A_only": ct["A_only"],
        "top_B_only": ct["B_only"],
        "top_neither": ct["neither"],
        "ef_dual_10pct": ef,
        "and_n_input": len(use),
        "and_pass_dual": cr["dual"],
        "and_pass_A": cr["A_only"],
        "and_pass_B": cr["B_only"],
        "and_dual_precision": cr["dual"] / len(retained) if retained else float("nan"),
        "auroc_impl_agreement": abs(da - da2) < 1e-12 and abs(db - db2) < 1e-12,
        "n_boot": N_BOOT,
        "seed": SEED,
        "bootstrap": "class-stratified; summary_min and matched-mismatch share dual draws; B=2000 seed=20260729",
    }


def sync_membership(master: list[dict]) -> None:
    path = ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv"
    old = read_csv(path)
    by = {(r["pair"], r["ligand_id"]): r for r in master}
    out = []
    seen = set()
    for r in old:
        key = (r["pair"], r["ligand"])
        m = by.get(key)
        if m:
            r = dict(r)
            r["cls"] = m["class"]
            r["pA"] = m["pA"]
            r["pB"] = m["pB"]
            r["score_A"] = m["score_A"]
            r["score_B"] = m["score_B"]
        out.append(r)
        seen.add(key)
    for m in master:
        key = (m["pair"], m["ligand_id"])
        if key not in seen:
            out.append(
                {
                    "pair": m["pair"],
                    "ligand": m["ligand_id"],
                    "cls": m["class"],
                    "pA": m["pA"],
                    "pB": m["pB"],
                    "score_A": m["score_A"],
                    "score_B": m["score_B"],
                }
            )
    write_csv(path, out)


def patch_unified(stats: dict) -> None:
    path = ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    rows = read_csv(path)
    for r in rows:
        if r["pair"] in ("EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR") and r["label_rule"] == "theta_6.0":
            s = stats[r["pair"]]
            r["n_dual"] = s["n_dual"]
            r["n_A_only"] = s["n_A_only"]
            r["n_B_only"] = s["n_B_only"]
            r["pocket_matched_summary_min"] = f"{s['summary_min']:.4f}"
            r["auroc_D_vs_A"] = f"{s['auroc_D_vs_A_pocketB']:.4f}"
            r["auroc_D_vs_B"] = f"{s['auroc_D_vs_B_pocketA']:.4f}"
            r["ci_lo"] = f"{s['summary_min_ci_lo']:.4f}"
            r["ci_hi"] = f"{s['summary_min_ci_hi']:.4f}"
    write_csv(path, rows)


def patch_formulation(stats: dict) -> None:
    path = ROOT / "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv"
    rows = read_csv(path)
    for r in rows:
        if r["pair"] not in stats:
            continue
        s = stats[r["pair"]]
        c = r["contrast"]
        if c == "D_vs_A_pocketB":
            r["auroc"] = f"{s['auroc_D_vs_A_pocketB']:.4f}"
            r["ci_lo"] = f"{s['ci_lo_D_vs_A']:.4f}"
            r["ci_hi"] = f"{s['ci_hi_D_vs_A']:.4f}"
            r["n_pos"] = s["n_dual"]
            r["n_neg"] = s["n_A_only"]
        elif c == "D_vs_B_pocketA":
            r["auroc"] = f"{s['auroc_D_vs_B_pocketA']:.4f}"
            r["ci_lo"] = f"{s['ci_lo_D_vs_B']:.4f}"
            r["ci_hi"] = f"{s['ci_hi_D_vs_B']:.4f}"
            r["n_pos"] = s["n_dual"]
            r["n_neg"] = s["n_B_only"]
        elif c == "summary_min":
            r["auroc"] = f"{s['summary_min']:.4f}"
            r["n_pos"] = s["n_dual"]
            r["n_neg"] = s["n_B_only"] if s["auroc_D_vs_B_pocketA"] <= s["auroc_D_vs_A_pocketB"] else s["n_A_only"]
        elif c == "D_vs_neither_mean":
            r["auroc"] = f"{s['auroc_mean_D_vs_neither']:.4f}"
            r["ci_lo"] = f"{s['ci_lo_D_vs_neither']:.4f}"
            r["ci_hi"] = f"{s['ci_hi_D_vs_neither']:.4f}"
            r["n_pos"] = s["n_dual"]
            r["n_neg"] = s["n_neither"]
    write_csv(path, rows)

    path = ROOT / "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv"
    rows = read_csv(path)
    for r in rows:
        if r["pair"] not in stats:
            continue
        s = stats[r["pair"]]
        if r["contrast"] == "D_vs_B_or_neither_pocketA":
            r["n_dual"] = s["n_dual"]
            r["n_selective"] = s["n_B_only"]
            r["n_neither"] = s["n_neither"]
            r["auroc_dual_vs_selective"] = f"{s['fixed_score_D_vs_selective_pocketA']:.4f}"
            r["auroc_dual_vs_neither"] = f"{s['fixed_score_D_vs_neither_pocketA']:.4f}"
            r["delta_neither_minus_selective"] = f"{s['fixed_score_delta_pocketA']:.4f}"
            r["delta_ci_lo"] = f"{s['fixed_score_delta_pocketA_ci_lo']:.4f}"
            r["delta_ci_hi"] = f"{s['fixed_score_delta_pocketA_ci_hi']:.4f}"
        elif r["contrast"] == "D_vs_A_or_neither_pocketB":
            r["n_dual"] = s["n_dual"]
            r["n_selective"] = s["n_A_only"]
            r["n_neither"] = s["n_neither"]
            r["auroc_dual_vs_selective"] = f"{s['fixed_score_D_vs_selective_pocketB']:.4f}"
            r["auroc_dual_vs_neither"] = f"{s['fixed_score_D_vs_neither_pocketB']:.4f}"
            r["delta_neither_minus_selective"] = f"{s['fixed_score_delta_pocketB']:.4f}"
            r["delta_ci_lo"] = f"{s['fixed_score_delta_pocketB_ci_lo']:.4f}"
            r["delta_ci_hi"] = f"{s['fixed_score_delta_pocketB_ci_hi']:.4f}"
    write_csv(path, rows)


def patch_ranking_and_wp(stats: dict, master: list[dict]) -> None:
    path = ROOT / "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv"
    rows = read_csv(path)
    for r in rows:
        if r["pair"] not in stats:
            continue
        s = stats[r["pair"]]
        recs = pair_recs(master, r["pair"])
        ranked = sorted(recs, key=lambda x: (-x["score_mean"], x["ligand_id"]))
        k = s["top_k"]
        r["n_ranked"] = s["n_complete_case"]
        r["n_dual"] = s["n_dual"]
        r["n_A_only"] = s["n_A_only"]
        r["n_B_only"] = s["n_B_only"]
        r["n_neither"] = s["n_neither"]
        r["top_k"] = k
        r["top_dual"] = s["top_dual"]
        r["top_A_only"] = s["top_A_only"]
        r["top_B_only"] = s["top_B_only"]
        r["top_neither"] = s["top_neither"]
        r["top_dual_fraction"] = s["top_dual"] / k
        r["panel_dual_fraction"] = s["n_dual"] / s["n_complete_case"]
        r["ef_dual_10pct"] = s["ef_dual_10pct"]
        r["top_ligand_ids"] = ";".join(x["ligand_id"] for x in ranked[:k])
        r["n_filter_input"] = s["and_n_input"]
        r["retained_dual"] = s["and_pass_dual"]
        r["retained_A_only"] = s["and_pass_A"]
        r["retained_B_only"] = s["and_pass_B"]
        r["dual_precision"] = s["and_dual_precision"]
    write_csv(path, rows)

    path = ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv"
    rows = read_csv(path)
    for r in rows:
        if r["set"] != "main_panel" or r["pair"] not in stats:
            continue
        s = stats[r["pair"]]
        r["n_dual"] = s["n_dual"]
        r["n_A_only"] = s["n_A_only"]
        r["n_B_only"] = s["n_B_only"]
        r["matched_D_vs_A"] = f"{s['auroc_D_vs_A_pocketB']:.4f}"
        r["matched_D_vs_B"] = f"{s['auroc_D_vs_B_pocketA']:.4f}"
        r["matched_summary_min"] = f"{s['matched_smin']:.4f}"
        r["wrong_summary_min"] = f"{s['mismatched_smin']:.4f}"
        r["delta_matched_minus_wrong"] = f"{s['matched_minus_mismatched']:.4f}"
        r["delta_boot_mean"] = f"{s['matched_minus_mismatched']:.4f}"
        r["delta_ci_lo"] = f"{s['matched_minus_mismatched_ci_lo']:.4f}"
        r["delta_ci_hi"] = f"{s['matched_minus_mismatched_ci_hi']:.4f}"
        r["ci_excludes_zero"] = str(s["matched_ci_excludes_0"])
        r["point_matched_gt_wrong"] = str(s["matched_minus_mismatched"] > 0)
    write_csv(path, rows)


def patch_master_results(stats: dict) -> None:
    path = ROOT / "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv"
    rows = read_csv(path)
    for r in rows:
        if r.get("block") != "primary_directional" or r["pair"] not in stats:
            continue
        s = stats[r["pair"]]
        r["n_dual"] = s["n_dual"]
        r["n_A_only"] = s["n_A_only"]
        r["n_B_only"] = s["n_B_only"]
        r["n_scored"] = s["n_dual"] + s["n_A_only"] + s["n_B_only"]
        if r["metric"] == "summary_min":
            r["value"] = f"{s['summary_min']:.4f}"
            r["ci_lo"] = f"{s['summary_min_ci_lo']:.4f}"
            r["ci_hi"] = f"{s['summary_min_ci_hi']:.4f}"
        elif r["metric"] == "auroc_D_vs_A_pocketB":
            r["value"] = f"{s['auroc_D_vs_A_pocketB']:.4f}"
        elif r["metric"] == "auroc_D_vs_B_pocketA":
            r["value"] = f"{s['auroc_D_vs_B_pocketA']:.4f}"
    write_csv(path, rows)


def patch_canonical_ligand_level(master: list[dict]) -> None:
    path = ROOT / "data/jcim_independent_audit_v0/tables/canonical_ligand_level_v1.csv"
    if not path.exists():
        return
    rows = read_csv(path)
    by = {(r["pair"], r["ligand_id"]): r for r in master}
    out = []
    seen = set()
    for r in rows:
        key = (r["pair"], r["ligand_id"])
        m = by.get(key)
        if m and r.get("analysis_set") == "main":
            r["score_A"] = m["score_A"]
            r["score_B"] = m["score_B"]
            r["score_mean"] = m["score_mean"]
            r["panel_class"] = m["class"]
            r["membership_cls"] = m["class"]
            r["pChEMBL_A"] = m["pA"]
            r["pChEMBL_B"] = m["pB"]
            if m["molecule_chembl_id"]:
                r["molecule_chembl_id"] = m["molecule_chembl_id"]
        out.append(r)
        if r.get("analysis_set") == "main":
            seen.add(key)
    fields = list(out[0].keys()) if out else []
    for m in master:
        key = (m["pair"], m["ligand_id"])
        if key in seen:
            continue
        rec = {k: "" for k in fields}
        rec.update(
            {
                "analysis_set": "main",
                "pair": m["pair"],
                "ligand_id": m["ligand_id"],
                "molecule_chembl_id": m["molecule_chembl_id"],
                "pChEMBL_A": m["pA"],
                "pChEMBL_B": m["pB"],
                "panel_class": m["class"],
                "membership_cls": m["class"],
                "score_A": m["score_A"],
                "score_B": m["score_B"],
                "score_mean": m["score_mean"],
                "main_holdout_status": "main",
            }
        )
        out.append(rec)
    write_csv(path, out)


def write_canonical_stats(stats: dict, master: list[dict]) -> None:
    write_csv(CANON / "current_score_master.csv", master)
    dir_rows, smin_rows, fix_rows, rank_rows, mean_rows, cls_rows = [], [], [], [], [], []
    for pair in ORDER:
        s = stats[pair]
        cls_rows.append({k: s[k] for k in ("pair", "n_dual", "n_A_only", "n_B_only", "n_neither", "n_complete_case", "n_missing_score_A", "n_missing_score_B", "n_class_pchembl_mismatch")})
        dir_rows.append(
            {
                "pair": pair,
                "auroc_D_vs_A_pocketB": s["auroc_D_vs_A_pocketB"],
                "ci_lo": s["ci_lo_D_vs_A"],
                "ci_hi": s["ci_hi_D_vs_A"],
                "auroc_D_vs_B_pocketA": s["auroc_D_vs_B_pocketA"],
                "ci_lo_B": s["ci_lo_D_vs_B"],
                "ci_hi_B": s["ci_hi_D_vs_B"],
                "n_dual": s["n_dual"],
                "n_A_only": s["n_A_only"],
                "n_B_only": s["n_B_only"],
                "bootstrap": s["bootstrap"],
                "second_impl_agrees": s["auroc_impl_agreement"],
            }
        )
        smin_rows.append(
            {
                "pair": pair,
                "summary_min": s["summary_min"],
                "ci_lo": s["summary_min_ci_lo"],
                "ci_hi": s["summary_min_ci_hi"],
                "n_boot": N_BOOT,
                "seed": SEED,
                "bootstrap": s["bootstrap"],
            }
        )
        fix_rows.append(
            {
                "pair": pair,
                "pocket": "A",
                "auroc_dual_vs_selective": s["fixed_score_D_vs_selective_pocketA"],
                "auroc_dual_vs_neither": s["fixed_score_D_vs_neither_pocketA"],
                "delta": s["fixed_score_delta_pocketA"],
                "delta_ci_lo": s["fixed_score_delta_pocketA_ci_lo"],
                "delta_ci_hi": s["fixed_score_delta_pocketA_ci_hi"],
            }
        )
        fix_rows.append(
            {
                "pair": pair,
                "pocket": "B",
                "auroc_dual_vs_selective": s["fixed_score_D_vs_selective_pocketB"],
                "auroc_dual_vs_neither": s["fixed_score_D_vs_neither_pocketB"],
                "delta": s["fixed_score_delta_pocketB"],
                "delta_ci_lo": s["fixed_score_delta_pocketB_ci_lo"],
                "delta_ci_hi": s["fixed_score_delta_pocketB_ci_hi"],
            }
        )
        mean_rows.append(
            {
                "pair": pair,
                "auroc_mean_D_vs_neither": s["auroc_mean_D_vs_neither"],
                "ci_lo": s["ci_lo_D_vs_neither"],
                "ci_hi": s["ci_hi_D_vs_neither"],
                "n_neither": s["n_neither"],
            }
        )
        rank_rows.append(
            {
                "pair": pair,
                "n": s["n_complete_case"],
                "k": s["top_k"],
                "top_D_A_B_N": f"{s['top_dual']}/{s['top_A_only']}/{s['top_B_only']}/{s['top_neither']}",
                "ef_dual_10pct": s["ef_dual_10pct"],
            }
        )
    write_csv(CANON / "class_counts.csv", cls_rows)
    write_csv(CANON / "primary_directional_auroc.csv", dir_rows)
    write_csv(CANON / "primary_summary_min.csv", smin_rows)
    write_csv(CANON / "fixed_score_negative_class_delta.csv", fix_rows)
    write_csv(CANON / "two_pocket_mean_ranking.csv", mean_rows)
    write_csv(CANON / "top10_operating_points.csv", rank_rows)
    write_csv(
        CANON / "matched_minus_mismatched.csv",
        [
            {
                "pair": p,
                "delta": stats[p]["matched_minus_mismatched"],
                "ci_lo": stats[p]["matched_minus_mismatched_ci_lo"],
                "ci_hi": stats[p]["matched_minus_mismatched_ci_hi"],
                "ci_excludes_0": stats[p]["matched_ci_excludes_0"],
            }
            for p in ORDER
        ],
    )


def main() -> int:
    CANON.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    master = build_master()
    stats = {p: compute_pair(pair_recs(master, p), p) for p in ORDER}
    write_canonical_stats(stats, master)
    sync_membership(master)
    patch_unified(stats)
    patch_formulation(stats)
    patch_ranking_and_wp(stats, master)
    patch_master_results(stats)
    patch_canonical_ligand_level(master)
    summary = {
        p: {
            "counts": (stats[p]["n_dual"], stats[p]["n_A_only"], stats[p]["n_B_only"], stats[p]["n_neither"]),
            "da": stats[p]["auroc_D_vs_A_pocketB"],
            "db": stats[p]["auroc_D_vs_B_pocketA"],
            "smin": stats[p]["summary_min"],
            "smin_ci": (stats[p]["summary_min_ci_lo"], stats[p]["summary_min_ci_hi"]),
            "deltaA": stats[p]["fixed_score_delta_pocketA"],
            "mm": stats[p]["matched_minus_mismatched"],
            "mm_ci": (stats[p]["matched_minus_mismatched_ci_lo"], stats[p]["matched_minus_mismatched_ci_hi"]),
            "ef": stats[p]["ef_dual_10pct"],
            "impl_ok": stats[p]["auroc_impl_agreement"],
        }
        for p in ORDER
    }
    (AUDIT / "independent_primary_recompute.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    for p in ORDER:
        s = stats[p]
        print(
            f"{p}: {s['n_dual']}/{s['n_A_only']}/{s['n_B_only']}/{s['n_neither']} "
            f"D/A={s['auroc_D_vs_A_pocketB']:.4f} D/B={s['auroc_D_vs_B_pocketA']:.4f} "
            f"smin={s['summary_min']:.4f} [{s['summary_min_ci_lo']:.4f},{s['summary_min_ci_hi']:.4f}] "
            f"ΔA={s['fixed_score_delta_pocketA']:.3f} mm={s['matched_minus_mismatched']:.3f} "
            f"EF={s['ef_dual_10pct']:.3f} impl={s['auroc_impl_agreement']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
