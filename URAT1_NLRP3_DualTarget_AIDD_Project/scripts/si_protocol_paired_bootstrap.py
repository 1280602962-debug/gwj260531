#!/usr/bin/env python3
"""Paired bootstrap comparison of docking/rescoring protocols (P0-P5) on TrueDecoy.

The existing per-protocol bootstrap (data/si/protocol_enrichment_ci/) resamples
each protocol independently, so the resulting 95% CIs cannot be used to claim
that Pi*=P2 is statistically distinguishable from another protocol. This
script draws the SAME bootstrap resample of molecule rows for every protocol
in each iteration, and reports the resample-wise EF@1% difference against P2
(paired bootstrap), which is the correct way to test whether one protocol's
early-enrichment estimate is significantly different from another's on the
same benchmark. Does not change the locked Pi*, does not rescore anything.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MOL_SCORES = PROJECT_ROOT / "data" / "benchmarks" / "protocol_selection" / "mol_protocol_scores.csv"
OUT_DIR = PROJECT_ROOT / "data" / "si" / "protocol_paired_bootstrap"

PROTOCOLS = [
    ("P0", "P0_CNNscore", True),
    ("P1", "P1_vina_affinity", False),
    ("P2", "P2_CNNaffinity", True),
    ("P3", "P3_gnina_affinity", False),
    ("P4", "P4_RTM_vina", True),
    ("P5", "P5_RTM_gnina", True),
]
REF_PROTOCOL = "P2"
N_BOOT = 2000
SEED = 42
FRACTIONS = (0.01, 0.05)


def ef_at_frac(score: np.ndarray, y: np.ndarray, frac: float, higher_is_better: bool) -> float:
    mask = ~np.isnan(score)
    score, y = score[mask], y[mask]
    if len(y) == 0 or y.sum() == 0:
        return np.nan
    n_top = max(1, int(np.floor(frac * len(y))))
    order = np.argsort(-score if higher_is_better else score)
    hits = y[order[:n_top]].sum()
    return float((hits / n_top) / (y.sum() / len(y)))


def run_for_benchmark(df_all: pd.DataFrame, flag_col: str, label: str) -> pd.DataFrame:
    df = df_all[df_all[flag_col].astype(int) == 1].reset_index(drop=True)
    y_all = (df["role"].astype(str).str.lower() == "active").to_numpy().astype(int)
    n = len(df)

    score_mat = {pid: pd.to_numeric(df[col], errors="coerce").to_numpy() for pid, col, _ in PROTOCOLS}

    rng = np.random.default_rng(SEED)
    ef_draws = {pid: {f: [] for f in FRACTIONS} for pid, _, _ in PROTOCOLS}
    for _ in range(N_BOOT):
        idx = rng.choice(n, size=n, replace=True)
        y_b = y_all[idx]
        for pid, _, hib in PROTOCOLS:
            s_b = score_mat[pid][idx]
            for frac in FRACTIONS:
                ef_draws[pid][frac].append(ef_at_frac(s_b, y_b, frac, hib))

    rows = []
    ref = ef_draws[REF_PROTOCOL]
    for pid, _, _ in PROTOCOLS:
        row = {"benchmark": label, "protocol": pid}
        for frac in FRACTIONS:
            a = np.array(ef_draws[pid][frac])
            r = np.array(ref[frac])
            ok = np.isfinite(a) & np.isfinite(r)
            diff = a[ok] - r[ok]
            ci_lo, ci_hi = np.percentile(diff, [2.5, 97.5]) if len(diff) else (np.nan, np.nan)
            p_two_sided = 2 * min((diff <= 0).mean(), (diff >= 0).mean()) if len(diff) else np.nan
            row[f"EF{int(frac*100)}pct_mean"] = float(np.nanmean(a))
            row[f"EF{int(frac*100)}pct_diff_vs_P2_mean"] = float(np.mean(diff)) if len(diff) else np.nan
            row[f"EF{int(frac*100)}pct_diff_vs_P2_ci95_low"] = float(ci_lo)
            row[f"EF{int(frac*100)}pct_diff_vs_P2_ci95_high"] = float(ci_hi)
            row[f"EF{int(frac*100)}pct_diff_vs_P2_p_two_sided"] = float(p_two_sided)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_all = pd.read_csv(MOL_SCORES, low_memory=False)

    out_true = run_for_benchmark(df_all, "in_true", "TrueDecoy")
    out_true.to_csv(OUT_DIR / "paired_bootstrap_vs_P2_true_decoy.csv", index=False)
    print("=== TrueDecoy ===")
    print(out_true.to_string(index=False))

    out_random = run_for_benchmark(df_all, "in_random", "RandomDecoy")
    out_random.to_csv(OUT_DIR / "paired_bootstrap_vs_P2_random_decoy.csv", index=False)
    print("\n=== RandomDecoy ===")
    print(out_random.to_string(index=False))

    summary = {
        "n_boot": N_BOOT,
        "reference_protocol": REF_PROTOCOL,
        "method": (
            "Same bootstrap resample (with replacement, molecule rows) applied to all "
            "protocols each draw; EF@frac recomputed per protocol per draw; reported "
            "difference is EF_protocol - EF_P2 on the same draws (paired), with a "
            "two-sided percentile-bootstrap p-value. Interval containing 0 means the "
            "difference from P2 is not statistically distinguishable at alpha=0.05. "
            "Run once on TrueDecoy, once on RandomDecoy."
        ),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
