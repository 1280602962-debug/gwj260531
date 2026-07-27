#!/usr/bin/env python3
"""Power / sample-size estimation from paired ligand-bootstrap of ΔAUROC.

ΔAUROC = AUROC(rtm_min_z) − AUROC(vina_mean), Dual vs rest.
No re-docking; uses existing per-ligand ablation score tables only.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]  # Dual_Target_Docking/
PAIRS = {
    "EGFR_HER2_panel40": REPO
    / "data/egfr_her2_panel40_v0/tables/ablation_ligand_scores.csv",
    "PIK3CA_mTOR_panel48": REPO
    / "data/pik3ca_mtor_panel48_v0/tables/ablation_ligand_scores.csv",
}
B = 4000
SEED = 20260728


def load(path: Path):
    rows = list(csv.DictReader(open(path)))
    y = np.array([1 if r["class"] == "dual" else 0 for r in rows])
    vina = np.array([float(r["vina_mean"]) for r in rows])
    rtm = np.array([float(r["rtm_min_z"]) for r in rows])
    return y, vina, rtm


def auroc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    c = 0.0
    for p in pos:
        c += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return c / (len(pos) * len(neg))


def main() -> None:
    rng = np.random.default_rng(SEED)
    summary = {}
    for name, path in PAIRS.items():
        y, vina, rtm = load(path)
        n = len(y)
        d_obs = auroc(y, rtm) - auroc(y, vina)
        deltas = []
        for _ in range(B):
            idx = rng.integers(0, n, n)
            yb = y[idx]
            if yb.sum() in (0, len(yb)):
                continue
            deltas.append(auroc(yb, rtm[idx]) - auroc(yb, vina[idx]))
        deltas = np.array(deltas)
        se = deltas.std(ddof=1)
        lo, hi = np.percentile(deltas, [2.5, 97.5])
        summary[name] = dict(n=n, d=d_obs, se=se, lo=lo, hi=hi)
        print(
            f"{name}: N={n} d={d_obs:+.4f} SE={se:.4f} "
            f"CI[{lo:+.4f},{hi:+.4f}] sig={'YES' if lo > 0 else 'no'}"
        )

    print("\nSample size (fixed effect, SE ~ 1/sqrt(N)):")
    for name, v in summary.items():
        for tag, z in [("95% sig", 1.96), ("80% power", 2.80), ("90% power", 3.24)]:
            se_t = abs(v["d"]) / z
            nn = v["n"] * (v["se"] / se_t) ** 2
            print(f"  {name} {tag}: N≈{nn:.0f}")

    print("\nCross-pair random-effects pooled:")
    ds = np.array([summary[k]["d"] for k in summary])
    tau2 = ds.var(ddof=1)
    for k in (2, 3, 5, 8, 10):
        se_pool = math.sqrt(tau2 / k)
        z = ds.mean() / se_pool if se_pool > 0 else float("inf")
        print(f"  K={k}: z≈{z:.2f} sig={'YES' if z >= 1.96 else 'no'}")


if __name__ == "__main__":
    main()
