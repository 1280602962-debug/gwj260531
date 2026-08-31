#!/usr/bin/env python3
"""SI heatmap of detectable-effect probabilities and main-text diagnostic workflow."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "data" / "jcim_bench_v0" / "scripts"))
from jcim_figure_style import (  # noqa: E402
    C,
    FONT,
    FS_ANNO,
    FS_AXIS,
    FS_TICK,
    OUT,
    PAIR_ORDER,
    apply_style,
    save_all,
)

TAB = ROOT / "data" / "jcim_novelty_v0" / "tables" / "detectable_effect_simulation_v1.csv"
TRUE = ["0.55", "0.60", "0.65", "0.70", "0.75"]


def load_summary_min() -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {p: {} for p in PAIR_ORDER}
    with TAB.open() as fh:
        for r in csv.DictReader(fh):
            if r["contrast"] != "summary_min":
                continue
            out[r["pair"]][r["true_auroc"]] = float(r["p_ci_excludes_0p5"])
    return out


def fig_s_detectable() -> None:
    apply_style()
    data = load_summary_min()
    mat = np.array([[data[p][a] for a in TRUE] for p in PAIR_ORDER], dtype=float)
    fig, ax = plt.subplots(figsize=(7.00, 3.40))
    im = ax.imshow(mat, cmap="YlGnBu", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(TRUE)))
    ax.set_xticklabels(TRUE)
    ax.set_yticks(range(len(PAIR_ORDER)))
    ax.set_yticklabels(PAIR_ORDER)
    ax.set_xlabel("True AUROC on both directional arms")
    ax.set_ylabel("")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(
                j,
                i,
                f"{v:.2f}",
                ha="center",
                va="center",
                fontsize=FS_ANNO,
                color="white" if v >= 0.55 else C["ink"],
                fontfamily=FONT,
            )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("P(95% CI excludes 0.5)", fontsize=FS_AXIS)
    cbar.ax.tick_params(labelsize=FS_TICK)
    fig.tight_layout()
    save_all(fig, "FigS_detectable_effect")
    plt.close(fig)


def _box(ax, xy, w, h, text, fc="#F4F7FA", ec=None):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.04",
        linewidth=0.9,
        facecolor=fc,
        edgecolor=ec or C["ink"],
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=7.0,
        fontfamily=FONT,
        color=C["ink"],
        wrap=True,
    )
    return x + w / 2, y, y + h


def fig_workflow() -> None:
    apply_style()
    fig, ax = plt.subplots(figsize=(7.00, 5.40))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Dual-target docking diagnostic", fontsize=FS_AXIS, pad=4)

    cx, _, top = _box(ax, (0.22, 0.88), 0.56, 0.09, "Favorable scores in both pockets", fc="#E8F1F8")
    steps = [
        (0.76, "1  Directional hard-negative test\nDual > A-only and Dual > B-only?"),
        (0.58, "2  Ligand-only chemical baseline\nDoes ECFP/property recover the signal?"),
        (0.40, "3  Unused ligand-pool test\nDoes the signal persist off the panel?"),
        (0.22, "4  Receptor-realization test\nStable across crystals?"),
    ]
    prev_bottom = 0.88
    centers = []
    for y, txt in steps:
        cxi, bot, top_i = _box(ax, (0.18, y), 0.44, 0.13, txt)
        centers.append((cxi, bot, top_i))
        ax.annotate(
            "",
            xy=(cxi, top_i + 0.002),
            xytext=(cxi, prev_bottom - 0.002),
            arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=0.9),
        )
        prev_bottom = y
    # Step 2 is inverted relative to the other gates: chemical recovery is the failure.
    fail = [
        (0.76, "No", "Treat as nonselective / suspect"),
        (0.58, "Yes", "Suspect chemical-series confounding"),
        (0.40, "No", "Suspect panel dependence"),
        (0.22, "No", "Receptor-dependent evidence"),
    ]
    pass_lab = ["Yes", "No", "Yes", "Yes"]
    for (cxi, bot, top_i), (y, side_lab, txt), down_lab in zip(centers, fail, pass_lab):
        _box(ax, (0.66, y + 0.015), 0.31, 0.10, txt, fc="#F8EDE8", ec=C["a_only"])
        ax.annotate(
            side_lab,
            xy=(0.655, y + 0.065),
            xytext=(0.62, y + 0.065),
            fontsize=6.5,
            fontfamily=FONT,
            color=C["a_only"],
            arrowprops=dict(arrowstyle="-|>", color=C["a_only"], lw=0.8),
            va="center",
        )
        ax.text(0.40, y - 0.012, down_lab, ha="center", va="top", fontsize=6.5, fontfamily=FONT, color=C["thick"])
    _box(ax, (0.18, 0.04), 0.44, 0.11, "Stronger computational evidence\nunder the present checks", fc="#E7F4EC", ec=C["thick"])
    ax.annotate(
        "",
        xy=(0.40, 0.15),
        xytext=(0.40, 0.22),
        arrowprops=dict(arrowstyle="-|>", color=C["thick"], lw=0.9),
    )
    fig.tight_layout()
    save_all(fig, "Fig8_diagnostic_workflow")
    plt.close(fig)


def main() -> None:
    fig_s_detectable()
    fig_workflow()
    print("wrote", OUT / "FigS_detectable_effect.png")
    print("wrote", OUT / "Fig8_diagnostic_workflow.png")


if __name__ == "__main__":
    main()
