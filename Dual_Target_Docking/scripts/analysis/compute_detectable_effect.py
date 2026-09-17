#!/usr/bin/env python3
"""Binormal detectable-effect simulation for the current eight pairs.

Zero docking. Class sizes come from current_score_master. Inner 95% CIs use
the same class-stratified shared-dual percentile bootstrap as Table 2
(B = 2000, seed 20260729). This is not observed/post-hoc power.

Contrasts:
  dual_vs_A_only / dual_vs_B_only / dual_vs_neither: class-stratified two-sample
  summary_min: dual ligands drawn once per replicate and applied to both pockets
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import N_BOOT, SEED  # noqa: E402

CANON = ROOT / "results" / "canonical"
MASTER = CANON / "current_score_master.csv"
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
TRUE_AUCS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)
N_MC_DEFAULT = 1000


def load_class_sizes() -> dict[str, dict[str, int]]:
    with MASTER.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    sizes: dict[str, dict[str, int]] = {}
    for pair in PAIRS:
        recs = [
            r
            for r in rows
            if r["pair"] == pair and r.get("analysis_set") == "main" and r.get("complete_case") == "1"
        ]
        counts = Counter(r.get("primary_class_theta6") or "" for r in recs)
        sizes[pair] = {
            "n_dual": int(counts["dual"]),
            "n_a": int(counts["A_only"]),
            "n_b": int(counts["B_only"]),
            "n_neither": int(counts["neither"]),
        }
        if min(sizes[pair]["n_dual"], sizes[pair]["n_a"], sizes[pair]["n_b"]) < 1:
            raise SystemExit(f"{pair} missing a required class: {sizes[pair]}")
    return sizes


def mu_from_auc(auc: float) -> float:
    a = float(auc)
    if a <= 0.0:
        return -np.inf
    if a >= 1.0:
        return np.inf
    return float(np.sqrt(2.0) * NormalDist().inv_cdf(a))


def batch_auroc(pos_b: np.ndarray, neg_b: np.ndarray) -> np.ndarray:
    diff = pos_b[:, :, None] - neg_b[:, None, :]
    n_p = pos_b.shape[1]
    n_n = neg_b.shape[1]
    return ((diff > 0).sum(axis=(1, 2)) + 0.5 * (diff == 0).sum(axis=(1, 2))) / (n_p * n_n)


def two_sample_boot_ci(pos: np.ndarray, neg: np.ndarray, rng: np.random.Generator, n_boot: int) -> tuple[float, float]:
    n_p, n_n = pos.shape[0], neg.shape[0]
    pb = pos[rng.integers(0, n_p, size=(n_boot, n_p))]
    nb = neg[rng.integers(0, n_n, size=(n_boot, n_n))]
    aucs = batch_auroc(pb, nb)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def summary_min_boot_ci(
    dual_sa: np.ndarray,
    dual_sb: np.ndarray,
    a_sb: np.ndarray,
    b_sa: np.ndarray,
    rng: np.random.Generator,
    n_boot: int,
) -> tuple[float, float]:
    nd, na, nb = dual_sa.shape[0], a_sb.shape[0], b_sa.shape[0]
    idd = rng.integers(0, nd, size=(n_boot, nd))
    ida = rng.integers(0, na, size=(n_boot, na))
    idb = rng.integers(0, nb, size=(n_boot, nb))
    auc_da = batch_auroc(dual_sb[idd], a_sb[ida])
    auc_db = batch_auroc(dual_sa[idd], b_sa[idb])
    smin = np.minimum(auc_da, auc_db)
    lo, hi = np.percentile(smin, [2.5, 97.5])
    return float(lo), float(hi)


def run_contrast(n_pos: int, n_neg: int, true_auc: float, rng: np.random.Generator, n_mc: int, n_boot: int) -> dict:
    mu = mu_from_auc(true_auc)
    excl = 0
    point = np.empty(n_mc, dtype=float)
    for i in range(n_mc):
        pos = rng.normal(mu, 1.0, size=n_pos)
        neg = rng.normal(0.0, 1.0, size=n_neg)
        lo, hi = two_sample_boot_ci(pos, neg, rng, n_boot)
        excl += int(lo > 0.5 or hi < 0.5)
        point[i] = float(batch_auroc(pos[None, :], neg[None, :])[0])
    p = excl / n_mc
    return dict(
        n_mc=n_mc,
        n_boot=n_boot,
        p_ci_excludes_0p5=p,
        se_binomial=float(np.sqrt(p * (1.0 - p) / n_mc)),
        mean_point_auroc=float(np.mean(point)),
    )


def run_summary_min(n_dual: int, n_a: int, n_b: int, true_auc: float, rng: np.random.Generator, n_mc: int, n_boot: int) -> dict:
    mu = mu_from_auc(true_auc)
    excl = 0
    point = np.empty(n_mc, dtype=float)
    for i in range(n_mc):
        dual_sa = rng.normal(mu, 1.0, size=n_dual)
        dual_sb = rng.normal(mu, 1.0, size=n_dual)
        a_sb = rng.normal(0.0, 1.0, size=n_a)
        b_sa = rng.normal(0.0, 1.0, size=n_b)
        lo, hi = summary_min_boot_ci(dual_sa, dual_sb, a_sb, b_sa, rng, n_boot)
        excl += int(lo > 0.5 or hi < 0.5)
        auc_da = float(batch_auroc(dual_sb[None, :], a_sb[None, :])[0])
        auc_db = float(batch_auroc(dual_sa[None, :], b_sa[None, :])[0])
        point[i] = min(auc_da, auc_db)
    p = excl / n_mc
    return dict(
        n_mc=n_mc,
        n_boot=n_boot,
        p_ci_excludes_0p5=p,
        se_binomial=float(np.sqrt(p * (1.0 - p) / n_mc)),
        mean_point_auroc=float(np.mean(point)),
    )


def fmt(stats: dict) -> dict:
    out = {}
    for key, val in stats.items():
        out[key] = f"{val:.6g}" if isinstance(val, float) else val
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mc", type=int, default=N_MC_DEFAULT)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    sizes = load_class_sizes()
    rng = np.random.default_rng(args.seed)
    rows = []
    for pair in PAIRS:
        n = sizes[pair]
        for true_auc in TRUE_AUCS:
            print(f"{pair} true={true_auc:.2f} n={n['n_dual']}/{n['n_a']}/{n['n_b']}/{n['n_neither']}", flush=True)
            for contrast, n_pos, n_neg in (
                ("dual_vs_A_only", n["n_dual"], n["n_a"]),
                ("dual_vs_B_only", n["n_dual"], n["n_b"]),
                ("dual_vs_neither", n["n_dual"], n["n_neither"]),
            ):
                stats = run_contrast(n_pos, n_neg, true_auc, rng, args.n_mc, args.n_boot)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "n_pos": n_pos,
                        "n_neg": n_neg,
                        "true_auroc": f"{true_auc:.2f}",
                        "bootstrap": "class_stratified_shared_dual",
                        "seed": args.seed,
                        **fmt(stats),
                    }
                )
            stats = run_summary_min(n["n_dual"], n["n_a"], n["n_b"], true_auc, rng, args.n_mc, args.n_boot)
            rows.append(
                {
                    "pair": pair,
                    "contrast": "summary_min",
                    "n_pos": n["n_dual"],
                    "n_neg": f"{n['n_a']}/{n['n_b']}",
                    "true_auroc": f"{true_auc:.2f}",
                    "bootstrap": "class_stratified_shared_dual",
                    "seed": args.seed,
                    **fmt(stats),
                }
            )

    out = CANON / "detectable_effect_simulation.csv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", out.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
