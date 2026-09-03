#!/usr/bin/env python3
"""Simulation-based detectable-effect analysis for DualFourClass-Bench.

No new docking. Uses the observed n_scored class sizes from Table 2 / Table 3
and the manuscript ligand-level bootstrap (B = 2000, percentile 95% CI).

This is not observed (post hoc) power. For each frozen (n_dual, n_neg) and a
grid of true AUROCs, it estimates the probability that the bootstrap CI
excludes 0.5 under a binormal score model.

summary_min uses class-preserving resampling of dual / A-only / B-only with
independent pocket-A and pocket-B score channels. Fixed class sizes are part
of the simulation design; this differs from the non-stratified empirical
bootstrap used for the canonical Table 2 interval.
"""
from __future__ import annotations

import csv
import json
from statistics import NormalDist
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_novelty_v0"
TAB = OUT / "tables"
AN = OUT / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
AN.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
N_MC = 1000
SEED = 20260729
TRUE_AUCS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)

# Table 2 n_scored (dual / A-only / B-only); Table 3 n_neither
PAIRS = {
    "EGFR/HER2": dict(n_dual=28, n_a=38, n_b=32, n_neither=12),
    "AChE/BChE": dict(n_dual=27, n_a=25, n_b=28, n_neither=15),
    "PIK3CA/PIK3CB": dict(n_dual=28, n_a=27, n_b=28, n_neither=16),
    "PIK3CA/mTOR": dict(n_dual=18, n_a=14, n_b=12, n_neither=4),
}


def mu_from_auc(auc: float) -> float:
    """Binormal equal-variance: pos~N(μ,1), neg~N(0,1) ⇒ AUROC = Φ(μ/√2)."""
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


def ci_excludes_half(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    return (lo > 0.5) | (hi < 0.5)


def two_sample_boot_ci(pos: np.ndarray, neg: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n_p, n_n = pos.shape[0], neg.shape[0]
    pb = pos[rng.integers(0, n_p, size=(N_BOOT, n_p))]
    nb = neg[rng.integers(0, n_n, size=(N_BOOT, n_n))]
    aucs = batch_auroc(pb, nb)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def summary_min_boot_ci(
    dual_sa: np.ndarray,
    dual_sb: np.ndarray,
    a_sb: np.ndarray,
    b_sa: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    nd, na, nb = dual_sa.shape[0], a_sb.shape[0], b_sa.shape[0]
    idd = rng.integers(0, nd, size=(N_BOOT, nd))
    ida = rng.integers(0, na, size=(N_BOOT, na))
    idb = rng.integers(0, nb, size=(N_BOOT, nb))
    auc_da = batch_auroc(dual_sb[idd], a_sb[ida])
    auc_db = batch_auroc(dual_sa[idd], b_sa[idb])
    smin = np.minimum(auc_da, auc_db)
    lo, hi = np.percentile(smin, [2.5, 97.5])
    return float(lo), float(hi)


def run_contrast(n_pos: int, n_neg: int, true_auc: float, rng: np.random.Generator) -> dict:
    mu = mu_from_auc(true_auc)
    excl = 0
    point = []
    for _ in range(N_MC):
        pos = rng.normal(mu, 1.0, size=n_pos)
        neg = rng.normal(0.0, 1.0, size=n_neg)
        lo, hi = two_sample_boot_ci(pos, neg, rng)
        excl += int(lo > 0.5 or hi < 0.5)
        point.append(float(batch_auroc(pos[None, :], neg[None, :])[0]))
    p = excl / N_MC
    se = float(np.sqrt(p * (1.0 - p) / N_MC))
    return dict(
        n_mc=N_MC,
        n_boot=N_BOOT,
        p_ci_excludes_0p5=p,
        se_binomial=se,
        mean_point_auroc=float(np.mean(point)),
    )


def run_summary_min(n_dual: int, n_a: int, n_b: int, true_auc: float, rng: np.random.Generator) -> dict:
    mu = mu_from_auc(true_auc)
    excl = 0
    point = []
    for _ in range(N_MC):
        dual_sa = rng.normal(mu, 1.0, size=n_dual)
        dual_sb = rng.normal(mu, 1.0, size=n_dual)
        a_sb = rng.normal(0.0, 1.0, size=n_a)
        b_sa = rng.normal(0.0, 1.0, size=n_b)
        lo, hi = summary_min_boot_ci(dual_sa, dual_sb, a_sb, b_sa, rng)
        excl += int(lo > 0.5 or hi < 0.5)
        auc_da = float(batch_auroc(dual_sb[None, :], a_sb[None, :])[0])
        auc_db = float(batch_auroc(dual_sa[None, :], b_sa[None, :])[0])
        point.append(min(auc_da, auc_db))
    p = excl / N_MC
    se = float(np.sqrt(p * (1.0 - p) / N_MC))
    return dict(
        n_mc=N_MC,
        n_boot=N_BOOT,
        p_ci_excludes_0p5=p,
        se_binomial=se,
        mean_point_auroc=float(np.mean(point)),
    )


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    for pair, n in PAIRS.items():
        for true_auc in TRUE_AUCS:
            print(f"{pair} true={true_auc:.2f}", flush=True)
            for contrast, n_pos, n_neg in (
                ("dual_vs_A_only", n["n_dual"], n["n_a"]),
                ("dual_vs_B_only", n["n_dual"], n["n_b"]),
                ("dual_vs_neither", n["n_dual"], n["n_neither"]),
            ):
                stats = run_contrast(n_pos, n_neg, true_auc, rng)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "n_pos": n_pos,
                        "n_neg": n_neg,
                        "true_auroc": f"{true_auc:.2f}",
                        **{k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in stats.items()},
                    }
                )
            stats = run_summary_min(n["n_dual"], n["n_a"], n["n_b"], true_auc, rng)
            rows.append(
                {
                    "pair": pair,
                    "contrast": "summary_min",
                    "n_pos": n["n_dual"],
                    "n_neg": f"{n['n_a']}/{n['n_b']}",
                    "true_auroc": f"{true_auc:.2f}",
                    **{k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in stats.items()},
                }
            )

    out_csv = TAB / "detectable_effect_simulation_v1.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    meta = {
        "n_mc": N_MC,
        "n_boot": N_BOOT,
        "seed": SEED,
        "score_model": "binormal_equal_variance",
        "bootstrap": "ligand-level class-preserving percentile 95% CI; fixed class sizes differ from the non-stratified empirical Table 2 bootstrap",
        "not": "observed/post-hoc power",
        "pairs": PAIRS,
        "source_csv": out_csv.relative_to(ROOT).as_posix(),
    }
    (AN / "DETECTABLE_EFFECT_SIMULATION_V1.md").write_text(
        _verdict_md(rows, meta), encoding="utf-8"
    )
    (TAB / "detectable_effect_simulation_v1.meta.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    print("wrote", out_csv)


def _verdict_md(rows: list[dict], meta: dict) -> str:
    def p(pair, contrast, auc):
        for r in rows:
            if r["pair"] == pair and r["contrast"] == contrast and r["true_auroc"] == f"{auc:.2f}":
                return float(r["p_ci_excludes_0p5"])
        raise KeyError((pair, contrast, auc))

    lines = [
        "# Detectable-effect simulation v1",
        "",
        "Zero docking. Binormal scores; ligand-level class-preserving bootstrap with fixed class sizes as part of the simulation design.",
        f"N_MC = {meta['n_mc']}; N_BOOT = {meta['n_boot']}; seed = {meta['seed']}.",
        "",
        "This is **not** observed power on the empirical AUROCs.",
        "",
        "## Probability that the 95% CI excludes 0.5 (`summary_min`)",
        "",
        "| Pair | n_scored | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for pair, n in PAIRS.items():
        ns = f"{n['n_dual']}/{n['n_a']}/{n['n_b']}"
        cells = " | ".join(f"{p(pair, 'summary_min', a):.3f}" for a in (0.55, 0.60, 0.65, 0.70, 0.75))
        lines.append(f"| {pair} | {ns} | {cells} |")
    lines += [
        "",
        "## Interpretation freeze",
        "",
        "- Current class sizes resolve **large** directional effects more readily than moderate ones.",
        "- Failure of an observed CI to exclude 0.5 does **not** establish equivalence to chance.",
        "- Dual versus neither uses a smaller negative set than the directional B-only/A-only arms on some pairs;",
        "  detectable-effect probabilities are therefore not interchangeable across formulations.",
        "",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
