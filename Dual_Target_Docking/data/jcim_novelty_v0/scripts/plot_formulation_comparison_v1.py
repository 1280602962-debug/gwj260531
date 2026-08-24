#!/usr/bin/env python3
"""SI figure: conventional Dual-vs-neither vs DualFourClass summary_min."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
OUT = ROOT / "data" / "jcim_novelty_v0" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]


def load(name):
    with (TAB / name).open() as fh:
        return list(csv.DictReader(fh))


def f(v):
    try:
        return float(v) if v not in ("", None) else np.nan
    except ValueError:
        return np.nan


def main():
    rows = load("formulation_conventional_vs_directional_v1.csv")
    by = {}
    for r in rows:
        by.setdefault(r["pair"], {})[r["contrast"]] = r

    x = np.arange(len(ORDER))
    width = 0.24
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=300)

    def series(contrast, lo_key="ci_lo", hi_key="ci_hi"):
        pts, yerr = [], [[], []]
        for p in ORDER:
            r = by[p].get(contrast, {})
            pt = f(r.get("auroc"))
            lo, hi = f(r.get(lo_key)), f(r.get(hi_key))
            pts.append(pt)
            if pt == pt and lo == lo and hi == hi:
                yerr[0].append(pt - lo)
                yerr[1].append(hi - pt)
            else:
                yerr[0].append(0)
                yerr[1].append(0)
        return np.asarray(pts), np.asarray(yerr)

    sm, sm_e = series("summary_min")
    nei, nei_e = series("D_vs_neither_mean")
    alln, alln_e = series("D_vs_A+B+neither_mean")

    ax.bar(x - width, sm, width, yerr=sm_e, capsize=2.5, label="directional summary_min", color="#1f4e79")
    ax.bar(x, nei, width, yerr=nei_e, capsize=2.5, label="conventional Dual vs neither", color="#c45911")
    ax.bar(x + width, alln, width, yerr=alln_e, capsize=2.5, label="Dual vs all non-duals", color="#7f7f7f")
    ax.axhline(0.5, color="0.4", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"], fontsize=8)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.15, 1.02)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(3, 0.22, "neither n=4\nunderpowered", ha="center", va="bottom", fontsize=6.5, color="#c45911")
    fig.tight_layout()
    fig.savefig(OUT / "FigS_formulation_conventional_vs_directional_v1.png")
    fig.savefig(OUT / "FigS_formulation_conventional_vs_directional_v1.pdf")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
