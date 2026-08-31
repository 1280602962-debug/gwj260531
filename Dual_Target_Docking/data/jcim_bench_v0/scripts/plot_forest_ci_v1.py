#!/usr/bin/env python3
"""Plot DualFourClass-Bench forest figure (summary_min ± 95% CI)."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

_MPL_CACHE = Path(tempfile.gettempdir()) / "dualfourclass-matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_bench_v0" / "tables"
FIG = ROOT / "data" / "jcim_bench_v0" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PAIRS = ["PIK3CA/mTOR", "AChE/BChE", "PIK3CA/PIK3CB", "EGFR/HER2"]
ARMS_ORDER = ["vina_mean", "rtm_min_z", "gnina_cnn_min", "heavy", "mw", "clogp", "tpsa"]
COLORS = {
    "vina_mean": "#1b6ca8",
    "rtm_min_z": "#2a9d8f",
    "gnina_cnn_min": "#264653",
    "heavy": "#b5651d",
    "mw": "#c27c2c",
    "clogp": "#d4a017",
    "tpsa": "#8b4513",
}


def main():
    rows = list(
        csv.DictReader(
            (TAB / "forest_summary_min_ci_v1.csv").open(encoding="utf-8", newline="")
        )
    )
    by = {(r["pair"], r["arm"]): r for r in rows}

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), sharex=True)
    axes = axes.ravel()
    for ax, pair in zip(axes, PAIRS):
        labels = []
        ys = []
        y = 0
        for arm in ARMS_ORDER:
            r = by.get((pair, arm))
            if not r:
                continue
            c = float(r["summary_min"])
            lo, hi = float(r["ci_lo"]), float(r["ci_hi"])
            col = COLORS.get(arm, "#444444")
            fam = "dock" if arm in ("vina_mean", "rtm_min_z", "gnina_cnn_min") else "base"
            ax.errorbar(
                [c],
                [y],
                xerr=[[c - lo], [hi - c]],
                fmt="o",
                color=col,
                ecolor=col,
                elinewidth=1.6,
                capsize=3,
                markersize=5,
            )
            labels.append(f"{arm} [{fam}]")
            ys.append(y)
            y += 1
        ax.axvline(0.5, color="#888888", ls="--", lw=0.9)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(pair, fontsize=11)
        ax.set_xlim(0.05, 0.95)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.25)
    axes[2].set_xlabel("summary_min AUROC (95% ligand bootstrap CI)")
    axes[3].set_xlabel("summary_min AUROC (95% ligand bootstrap CI)")
    fig.suptitle("DualFourClass-Bench: directional summary_min with bootstrap CIs", fontsize=12)
    fig.tight_layout()
    out = FIG / "forest_summary_min_ci_v1.png"
    fig.savefig(out, dpi=180)
    fig.savefig(FIG / "forest_summary_min_ci_v1.pdf")
    print("wrote", out)

    gate = list(csv.DictReader((TAB / "baseline_gate_bootstrap_v1.csv").open()))
    fig2, ax = plt.subplots(figsize=(8.5, 5.2))
    for i, g in enumerate(gate):
        d = float(g["delta_summary_min"])
        lo, hi = float(g["delta_ci_lo"]), float(g["delta_ci_hi"])
        if g["beats_baseline_ci_excl0"] == "True":
            col = "#2a9d8f"
        elif g["loses_baseline_ci_excl0"] == "True":
            col = "#c1121f"
        else:
            col = "#6c757d"
        ax.errorbar(
            [d],
            [i],
            xerr=[[d - lo], [hi - d]],
            fmt="o",
            color=col,
            ecolor=col,
            elinewidth=1.8,
            capsize=3,
            markersize=5,
        )
    ax.axvline(0, color="#222222", lw=1)
    ax.set_yticks(range(len(gate)))
    ax.set_yticklabels(
        [f"{g['pair']} | {g['dock_arm']}−{g['best_baseline_arm']}" for g in gate],
        fontsize=8,
    )
    ax.set_xlabel("Δ summary_min (dock − best trivial baseline), 95% CI")
    ax.set_title("Baseline gate with joint ligand bootstrap")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.25)
    fig2.tight_layout()
    out2 = FIG / "baseline_gate_delta_ci_v1.png"
    fig2.savefig(out2, dpi=180)
    fig2.savefig(FIG / "baseline_gate_delta_ci_v1.pdf")
    print("wrote", out2)


if __name__ == "__main__":
    main()
