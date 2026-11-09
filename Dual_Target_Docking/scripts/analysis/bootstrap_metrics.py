#!/usr/bin/env python3
"""Scheme-B class-stratified nonparametric percentile bootstrap.

Primary protocol (all eight pairs, all directional estimands):
  B = 2000
  seed = 20260729
  percentile 95% CI (2.5, 97.5)
  point estimate from the original complete sample (not the bootstrap mean)

Directional AUROC resamples each experimental-state class independently.
summary_min draws dual ligands once per replicate and applies that draw to
both pockets; A-only and B-only are resampled independently; min is taken
inside the replicate.

Fixed-score Δ shares the dual resample between dual-vs-selective and
dual-vs-neither.

Matched/mismatched Δ uses paired ligand resampling (each ligand keeps both
pocket scores).
"""
from __future__ import annotations

import math

import numpy as np

from analysis.analysis_config import N_BOOT, SEED, STRICT_HI, STRICT_LO, THETA_PRIMARY  # noqa: E402


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def percentile_ci(values, n_boot: int = N_BOOT) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < max(20, n_boot // 4):
        return float("nan"), float("nan")
    lo, hi = np.percentile(arr, [2.5, 97.5])
    return float(lo), float(hi)


def assign_fourclass(pA, pB, cut: float = THETA_PRIMARY) -> str | None:
    if pA is None or pB is None or (isinstance(pA, float) and math.isnan(pA)):
        return None
    try:
        a = float(pA) >= cut
        b = float(pB) >= cut
    except (TypeError, ValueError):
        return None
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def assign_strict(pA, pB) -> str | None:
    if pA is None or pB is None:
        return None
    try:
        pa, pb = float(pA), float(pB)
    except (TypeError, ValueError):
        return None
    if pa >= STRICT_HI and pb >= STRICT_HI:
        return "dual"
    if pa >= STRICT_HI and pb <= STRICT_LO:
        return "A_only"
    if pb >= STRICT_HI and pa <= STRICT_LO:
        return "B_only"
    if pa <= STRICT_LO and pb <= STRICT_LO:
        return "neither"
    return "gray"


def _choice(rng: np.random.Generator, idx: np.ndarray) -> np.ndarray:
    if idx.size == 0:
        return idx
    return rng.choice(idx, size=idx.size, replace=True)


def stratified_auroc_ci(pos, neg, n_boot: int = N_BOOT, seed: int = SEED):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    point = auroc(pos, neg)
    if pos.size == 0 or neg.size == 0:
        return point, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        vals[i] = auroc(rng.choice(pos, size=pos.size, replace=True), rng.choice(neg, size=neg.size, replace=True))
    lo, hi = percentile_ci(vals, n_boot)
    return float(point), lo, hi


def unstratified_auroc_ci(pos, neg, n_boot: int = N_BOOT, seed: int = SEED):
    """Ligand-pool bootstrap of the combined positive+negative set (sensitivity)."""
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    point = auroc(pos, neg)
    if pos.size == 0 or neg.size == 0:
        return point, float("nan"), float("nan")
    y = np.concatenate([np.ones(pos.size), np.zeros(neg.size)])
    s = np.concatenate([pos, neg])
    rng = np.random.default_rng(seed)
    vals = []
    n = s.size
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        p = s[idx][y[idx] == 1]
        nneg = s[idx][y[idx] == 0]
        if p.size == 0 or nneg.size == 0:
            continue
        vals.append(auroc(p, nneg))
    lo, hi = percentile_ci(vals, n_boot)
    return float(point), lo, hi


def summary_min_stratified(score_A, score_B, cls, n_boot: int = N_BOOT, seed: int = SEED):
    score_A = np.asarray(score_A, dtype=float)
    score_B = np.asarray(score_B, dtype=float)
    cls = np.asarray(cls)
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")
    da = auroc(score_B[d], score_B[a])
    db = auroc(score_A[d], score_A[b])
    point = min(da, db)
    weaker = "D_vs_A_pocketB" if da <= db else "D_vs_B_pocketA"
    rng = np.random.default_rng(seed)
    mins = np.empty(n_boot, dtype=float)
    das = np.empty(n_boot, dtype=float)
    dbs = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        di = _choice(rng, d)
        ai = _choice(rng, a)
        bi = _choice(rng, b)
        das[i] = auroc(score_B[di], score_B[ai])
        dbs[i] = auroc(score_A[di], score_A[bi])
        mins[i] = min(das[i], dbs[i])
    lo, hi = percentile_ci(mins, n_boot)
    da_lo, da_hi = percentile_ci(das, n_boot)
    db_lo, db_hi = percentile_ci(dbs, n_boot)
    return {
        "auroc_D_vs_A_pocketB": float(da),
        "auroc_D_vs_A_ci_lo": da_lo,
        "auroc_D_vs_A_ci_hi": da_hi,
        "auroc_D_vs_B_pocketA": float(db),
        "auroc_D_vs_B_ci_lo": db_lo,
        "auroc_D_vs_B_ci_hi": db_hi,
        "summary_min": float(point),
        "summary_min_ci_lo": lo,
        "summary_min_ci_hi": hi,
        "weaker_arm": weaker,
        "n_dual": int(d.size),
        "n_A_only": int(a.size),
        "n_B_only": int(b.size),
    }


def summary_min_unstratified(score_A, score_B, cls, n_boot: int = N_BOOT, seed: int = SEED):
    score_A = np.asarray(score_A, dtype=float)
    score_B = np.asarray(score_B, dtype=float)
    cls = np.asarray(cls)
    use = np.flatnonzero(np.isin(cls, ["dual", "A_only", "B_only"]))
    point_row = summary_min_stratified(score_A, score_B, cls, n_boot=1, seed=seed)
    rng = np.random.default_rng(seed)
    mins = []
    for _ in range(n_boot):
        idx = rng.choice(use, size=use.size, replace=True)
        lab = cls[idx]
        d = np.flatnonzero(lab == "dual")
        a = np.flatnonzero(lab == "A_only")
        b = np.flatnonzero(lab == "B_only")
        if min(d.size, a.size, b.size) == 0:
            continue
        da = auroc(score_B[idx][d], score_B[idx][a])
        db = auroc(score_A[idx][d], score_A[idx][b])
        mins.append(min(da, db))
    lo, hi = percentile_ci(mins, n_boot)
    return point_row["summary_min"], lo, hi, point_row


def fixed_score_delta_stratified(dual, selective, neither, n_boot: int = N_BOOT, seed: int = SEED):
    dual = np.asarray(dual, dtype=float)
    selective = np.asarray(selective, dtype=float)
    neither = np.asarray(neither, dtype=float)
    auc_s = auroc(dual, selective)
    auc_n = auroc(dual, neither)
    point = auc_n - auc_s
    if min(dual.size, selective.size, neither.size) == 0:
        return auc_s, auc_n, point, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        d = rng.choice(dual, size=dual.size, replace=True)
        s = rng.choice(selective, size=selective.size, replace=True)
        n = rng.choice(neither, size=neither.size, replace=True)
        deltas[i] = auroc(d, n) - auroc(d, s)
    lo, hi = percentile_ci(deltas, n_boot)
    return float(auc_s), float(auc_n), float(point), lo, hi


def fixed_score_delta_unstratified(dual, selective, neither, n_boot: int = N_BOOT, seed: int = SEED):
    dual = np.asarray(dual, dtype=float)
    selective = np.asarray(selective, dtype=float)
    neither = np.asarray(neither, dtype=float)
    auc_s = auroc(dual, selective)
    auc_n = auroc(dual, neither)
    point = auc_n - auc_s
    y = np.array(["dual"] * dual.size + ["sel"] * selective.size + ["nei"] * neither.size)
    s = np.concatenate([dual, selective, neither])
    rng = np.random.default_rng(seed)
    deltas = []
    n = s.size
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        lab = y[idx]
        d = s[idx][lab == "dual"]
        sel = s[idx][lab == "sel"]
        nei = s[idx][lab == "nei"]
        if min(d.size, sel.size, nei.size) == 0:
            continue
        deltas.append(auroc(d, nei) - auroc(d, sel))
    lo, hi = percentile_ci(deltas, n_boot)
    return float(auc_s), float(auc_n), float(point), lo, hi


def matched_mismatched_paired(score_A, score_B, cls, n_boot: int = N_BOOT, seed: int = SEED, stratified: bool = True):
    """Δ = matched summary_min − mismatched summary_min.

    Matched: D vs A uses pocket B; D vs B uses pocket A.
    Mismatched: D vs A uses pocket A; D vs B uses pocket B.
    Each ligand is drawn once per replicate (pairing). When stratified, D/A/B
    class sizes are preserved.
    """
    score_A = np.asarray(score_A, dtype=float)
    score_B = np.asarray(score_B, dtype=float)
    cls = np.asarray(cls)
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")

    def smin(idx_d, idx_a, idx_b, matched: bool):
        if matched:
            da = auroc(score_B[idx_d], score_B[idx_a])
            db = auroc(score_A[idx_d], score_A[idx_b])
        else:
            da = auroc(score_A[idx_d], score_A[idx_a])
            db = auroc(score_B[idx_d], score_B[idx_b])
        weaker = "D_vs_A" if da <= db else "D_vs_B"
        return min(da, db), da, db, weaker

    m_pt, m_da, m_db, m_weak = smin(d, a, b, True)
    u_pt, u_da, u_db, u_weak = smin(d, a, b, False)
    delta = m_pt - u_pt
    switched = m_weak != u_weak
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot, dtype=float)
    pool = np.concatenate([d, a, b])
    for i in range(n_boot):
        if stratified:
            di, ai, bi = _choice(rng, d), _choice(rng, a), _choice(rng, b)
        else:
            draw = rng.choice(pool, size=pool.size, replace=True)
            lab = cls[draw]
            di = draw[lab == "dual"]
            ai = draw[lab == "A_only"]
            bi = draw[lab == "B_only"]
            if min(di.size, ai.size, bi.size) == 0:
                deltas[i] = np.nan
                continue
        mv, *_ = smin(di, ai, bi, True)
        uv, *_ = smin(di, ai, bi, False)
        deltas[i] = mv - uv
    lo, hi = percentile_ci(deltas, n_boot)
    return {
        "matched_summary_min": float(m_pt),
        "matched_auroc_D_vs_A": float(m_da),
        "matched_auroc_D_vs_B": float(m_db),
        "matched_weaker_arm": m_weak,
        "mismatched_summary_min": float(u_pt),
        "mismatched_auroc_D_vs_A": float(u_da),
        "mismatched_auroc_D_vs_B": float(u_db),
        "mismatched_weaker_arm": u_weak,
        "delta": float(delta),
        "delta_ci_lo": lo,
        "delta_ci_hi": hi,
        "ci_excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)),
        "weaker_arm_switched": bool(switched),
        "n_dual": int(d.size),
        "n_A_only": int(a.size),
        "n_B_only": int(b.size),
    }


def cluster_delta(recs, group_key: str, n_boot: int = N_BOOT, seed: int = SEED):
    """Fixed-score Δ with group (scaffold or document) resampling.

    recs: dicts with cls in {dual, selective, neither}, score, and group_key.
    """
    from collections import defaultdict

    groups = defaultdict(list)
    for r in recs:
        groups[r[group_key]].append(r)
    gids = list(groups)
    dual = [r["score"] for r in recs if r["cls"] == "dual"]
    sel = [r["score"] for r in recs if r["cls"] == "selective"]
    nei = [r["score"] for r in recs if r["cls"] == "neither"]
    point = auroc(dual, nei) - auroc(dual, sel)
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(n_boot):
        draw = rng.choice(gids, size=len(gids), replace=True)
        bag = [r for gid in draw for r in groups[gid]]
        d = [r["score"] for r in bag if r["cls"] == "dual"]
        s = [r["score"] for r in bag if r["cls"] == "selective"]
        n = [r["score"] for r in bag if r["cls"] == "neither"]
        if min(len(d), len(s), len(n)) < 2:
            continue
        deltas.append(auroc(d, n) - auroc(d, s))
    lo, hi = percentile_ci(deltas, n_boot)
    return {
        "n_groups": len(gids),
        "n_dual": len(dual),
        "n_selective": len(sel),
        "n_neither": len(nei),
        "delta_point": float(point),
        "delta_ci_lo": lo,
        "delta_ci_hi": hi,
        "n_valid_boot": len(deltas),
        "excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)),
    }


def ci_crosses(lo, hi, value: float) -> bool:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return True
    return lo <= value <= hi
