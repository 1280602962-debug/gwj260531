#!/usr/bin/env python3
"""
Non-docking computational module E — Pareto / multi-objective robustness.

This module quantifies how sensitive the P2 dual-success shortlist is to the ranking rule,
WITHOUT modifying the existing front: it recomputes candidate sets under

  (1) top-k% dual-percentile intersection  (k = 1, 2, 5, 10)
  (2) percentile-threshold gates            (tau = 85, 90, 95 on both S_U and S_N)
  (3) bootstrap Pareto-front membership     (resample pool rows, re-derive front,
      report each compound's front-membership frequency = stability)

Reads the EXISTING dual-docked pool; writes new candidate tables only. The
original pareto_shortlist.csv / pareto_merged_scores.csv are never overwritten.

Input (read-only):
  results/repurposing/pareto_merged_scores.csv

Outputs:
  results/pareto_robustness/topk_intersection.csv
  results/pareto_robustness/threshold_gates.csv
  results/pareto_robustness/bootstrap_front_stability.csv
  results/pareto_robustness/pareto_robustness_summary.json

Usage:
  python3 scripts/13_pareto_robustness.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARETO_DIR = PROJECT_ROOT / "results" / "repurposing"
OUT_DIR = PROJECT_ROOT / "results" / "pareto_robustness"

SU = "s_u_percentile"
SN = "s_n_percentile"
NAME = "name"


def pareto_front_mask(su: np.ndarray, sn: np.ndarray) -> np.ndarray:
    """Non-dominated set maximizing (su, sn)."""
    n = len(su)
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        if dominated[i]:
            continue
        for j in range(n):
            if i == j:
                continue
            if (su[j] >= su[i] and sn[j] >= sn[i]) and (su[j] > su[i] or sn[j] > sn[i]):
                dominated[i] = True
                break
    return ~dominated


def topk_intersection(df: pd.DataFrame, ks=(1, 2, 5, 10)) -> pd.DataFrame:
    rows = []
    n = len(df)
    for k in ks:
        cut = 100 - k
        sel = df[(df[SU] >= cut) & (df[SN] >= cut)]
        rows.append({
            "top_k_pct": k,
            "percentile_cut": cut,
            "n_intersection": int(len(sel)),
            "names": "; ".join(sel[NAME].astype(str).head(30).tolist()),
        })
    return pd.DataFrame(rows)


def threshold_gates(df: pd.DataFrame, taus=(85, 90, 95)) -> pd.DataFrame:
    rows = []
    for tau in taus:
        sel = df[(df[SU] >= tau) & (df[SN] >= tau)]
        rows.append({
            "tau_percentile": tau,
            "n_pass": int(len(sel)),
            "names": "; ".join(sel[NAME].astype(str).head(30).tolist()),
        })
    return pd.DataFrame(rows)


def bootstrap_stability(df: pd.DataFrame, n_boot: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    su = df[SU].to_numpy()
    sn = df[SN].to_numpy()
    names = df[NAME].astype(str).to_numpy()
    n = len(df)
    counts = {nm: 0 for nm in names}
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        mask = pareto_front_mask(su[idx], sn[idx])
        for nm in np.unique(names[idx][mask]):
            counts[nm] += 1
    rows = [{"name": nm, "front_frequency": round(c / n_boot, 3)} for nm, c in counts.items()]
    out = pd.DataFrame(rows).sort_values("front_frequency", ascending=False)
    return out[out["front_frequency"] > 0].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pareto robustness (non-docking module E)")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--n-boot", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pool)
    df = df[df[SU].notna() & df[SN].notna()].reset_index(drop=True)

    topk = topk_intersection(df)
    topk.to_csv(args.output_dir / "topk_intersection.csv", index=False)

    gates = threshold_gates(df)
    gates.to_csv(args.output_dir / "threshold_gates.csv", index=False)

    boot = bootstrap_stability(df, args.n_boot, args.seed)
    boot.to_csv(args.output_dir / "bootstrap_front_stability.csv", index=False)

    # Cross-check: production front size on full pool
    prod_mask = pareto_front_mask(df[SU].to_numpy(), df[SN].to_numpy())
    prod_front = df.loc[prod_mask, NAME].astype(str).tolist()

    summary = {
        "module": "E_pareto_robustness",
        "n_pool": int(len(df)),
        "recomputed_front_size": int(prod_mask.sum()),
        "recomputed_front_names": prod_front,
        "topk_intersection": topk.to_dict(orient="records"),
        "threshold_gates": gates.to_dict(orient="records"),
        "n_bootstrap": args.n_boot,
        "stable_front_members_freq_ge_0.5": boot[boot["front_frequency"] >= 0.5].to_dict(orient="records"),
    }
    with open(args.output_dir / "pareto_robustness_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=== Pareto robustness ===")
    print(f"  pool={len(df)}  recomputed front size={prod_mask.sum()} (production reported 6)")
    print("\nTop-k% dual intersection (expands the thin front):")
    print(topk[["top_k_pct", "percentile_cut", "n_intersection"]].to_string(index=False))
    print("\nThreshold gates:")
    print(gates[["tau_percentile", "n_pass"]].to_string(index=False))
    print("\nBootstrap front stability (freq >= 0.5):")
    print(boot[boot["front_frequency"] >= 0.5].to_string(index=False))


if __name__ == "__main__":
    main()
