#!/usr/bin/env python3
"""Supporting figures. Values are read from the CSV files written by this run."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from posthoc_common import FIGURE_DIR, OUTPUT_DIR, load_config

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

COLORS = {
    "dual": "#0072B2",
    "B_only": "#D55E00",
    "neither": "#009E73",
    "A_only": "#E69F00",
    "D_vs_A": "#0072B2",
    "D_vs_B": "#D55E00",
    "S1": "#4D4D4D",
    "improve": "#009E73",
    "worsen": "#D55E00",
    "discordant": "#CC79A7",
    "nochange": "#8C8C8C",
}


def _save(fig, stem: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{stem}.pdf")
    fig.savefig(FIGURE_DIR / f"{stem}.svg")
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=600)
    plt.close(fig)


def _stamp(ax, text: str = "Supporting analysis") -> None:
    ax.text(0.0, 1.02, text, transform=ax.transAxes, fontsize=8, fontstyle="italic", color="#333333")


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(OUTPUT_DIR / name)


def figure_a_similarity() -> None:
    cfg = load_config()
    df = _read("A_TEST_NEAREST_TRAIN.csv")
    order = []
    for pair in cfg["pair_order"]:
        for arm in ("D_vs_A", "D_vs_B"):
            order.append((pair, arm))
    rng = np.random.default_rng(int(cfg["random_seed"]))
    fig, ax = plt.subplots(figsize=(14.5, 5.4))
    for i, (pair, arm) in enumerate(order):
        values = df[(df["pair"] == pair) & (df["arm"] == arm)]["max_ecfp4_tanimoto_to_train"].to_numpy(float)
        jitter = rng.uniform(-0.16, 0.16, size=values.size)
        ax.scatter(np.full(values.size, i) + jitter, values, s=12, color=COLORS[arm], alpha=0.75, linewidths=0, zorder=2)
        ax.hlines(np.median(values), i - 0.28, i + 0.28, colors="black", linewidths=1.2, zorder=3)
    ax.axhline(0.70, color="#666666", linestyle="--", linewidth=0.8)
    ax.axhline(0.85, color="#666666", linestyle=":", linewidth=0.8)
    ax.set_ylabel("Maximum ECFP4 Tanimoto similarity to the training fold")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{p}\n{a.replace('_', ' ')}" for p, a in order], rotation=55, ha="right", fontsize=6.5)
    ax.set_ylim(-0.02, 1.05)
    ax.plot([], [], color=COLORS["D_vs_A"], marker="o", linestyle="none", label="D vs A-only")
    ax.plot([], [], color=COLORS["D_vs_B"], marker="o", linestyle="none", label="D vs B-only")
    ax.plot([], [], color="black", linewidth=1.2, label="Median")
    ax.plot([], [], color="#666666", linestyle="--", label="0.70 descriptive cutoff")
    ax.plot([], [], color="#666666", linestyle=":", label="0.85 descriptive cutoff")
    ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.42))
    ax.set_title("Supporting analysis: train–test chemical neighborhood", pad=8)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.36)
    _save(fig, "A_train_test_similarity")


def figure_b_scores() -> None:
    df = _read("B_SCORE_POINTS.csv")
    cfg = load_config()
    rng = np.random.default_rng(int(cfg["random_seed"]))
    labels = ["dual", "B_only", "neither"]
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for i, label in enumerate(labels):
        values = df.loc[df["class"] == label, "M0_score_A"].to_numpy(float)
        jitter = rng.uniform(-0.12, 0.12, size=values.size)
        ax.scatter(np.full(values.size, i) + jitter, values, s=22, color=COLORS[label], alpha=0.85, linewidths=0, zorder=2)
        q1, med, q3 = np.percentile(values, [25, 50, 75], method="linear")
        ax.vlines(i, q1, q3, colors="black", linewidths=3, zorder=3)
        ax.scatter([i], [med], marker="D", s=22, color="black", zorder=4)
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"dual\nn={int((df['class']=='dual').sum())}", f"B-only\nn={int((df['class']=='B_only').sum())}", f"neither\nn={int((df['class']=='neither').sum())}"])
    ax.set_ylabel("S1 score on EGFR (higher is more favorable)")
    ax.set_title("EGFR/HER2 S1 scores", pad=16)
    ax.plot([], [], color="black", linewidth=3, label="IQR")
    ax.plot([], [], color="black", marker="D", linestyle="none", label="Median")
    ax.legend(frameon=False, loc="lower right")
    _stamp(ax)
    fig.tight_layout()
    _save(fig, "B_egfr_s1_score_distributions")


def figure_b_descriptors() -> None:
    points = _read("B_DESCRIPTOR_POINTS.csv")
    assoc = _read("B_EGFR_SCORE_DESCRIPTOR_ASSOCIATIONS.csv")
    names = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "RingCount", "FormalCharge"]
    fig, axes = plt.subplots(2, 4, figsize=(11.2, 6.2))
    for ax, name in zip(axes.ravel(), names):
        for label in ("dual", "B_only"):
            sub = points[points["class"] == label]
            ax.scatter(sub[name], sub["M0_score_A"], s=16, color=COLORS[label], alpha=0.8, linewidths=0)
        rows = assoc[assoc["descriptor"] == name]
        bits = []
        for scope, short in (("pooled_dual_and_B_only", "all"), ("dual", "dual"), ("B_only", "B")):
            row = rows[rows["scope"] == scope].iloc[0]
            if row["status"] == "CONSTANT_VARIABLE":
                bits.append(f"{short}: const")
            else:
                bits.append(f"{short} ρ={float(row['spearman_rho']):.2f}")
        ax.set_title(name + "\n" + ", ".join(bits), fontsize=8)
        ax.set_xlabel(name)
        ax.set_ylabel("S1 EGFR")
    axes[0, 0].scatter([], [], color=COLORS["dual"], label="dual")
    axes[0, 0].scatter([], [], color=COLORS["B_only"], label="B-only")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Supporting analysis: descriptors and EGFR S1 scores", fontsize=11)
    fig.tight_layout()
    _save(fig, "B_descriptor_score_associations")


def figure_b_activity() -> None:
    df = _read("B_EGFR_ACTIVITY_MOLECULES.csv")
    cfg = load_config()
    rng = np.random.default_rng(int(cfg["random_seed"]))
    panels = [
        ("B_only", "pchembl_EGFR", "B-only EGFR"),
        ("dual", "pchembl_EGFR", "dual EGFR"),
        ("dual", "pchembl_HER2", "dual HER2"),
        ("neither", "pchembl_EGFR", "neither EGFR"),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for i, (label, column, title) in enumerate(panels):
        values = df.loc[df["class"] == label, column].to_numpy(float)
        jitter = rng.uniform(-0.12, 0.12, size=values.size)
        ax.scatter(np.full(values.size, i) + jitter, values, s=22, color=COLORS.get(label, "#333333"), alpha=0.85, linewidths=0)
        ax.hlines(np.median(values), i - 0.25, i + 0.25, colors="black", linewidths=1.3)
    ax.axhline(6.0, color="#666666", linestyle="--", linewidth=0.8, label="Formal activity threshold 6.0")
    ax.set_xticks(range(4))
    ax.set_xticklabels([p[2] for p in panels])
    ax.set_ylabel("Verified aggregate pChEMBL")
    ax.set_title("Experimental activity in the EGFR/HER2 diagnostic set", pad=16)
    ax.legend(frameon=False, loc="upper right")
    _stamp(ax)
    fig.tight_layout()
    _save(fig, "B_verified_activity_values")


def figure_c_dual_fraction() -> None:
    df = _read("C_TOPK_METHOD_COMPARISONS.csv")
    cfg = load_config()
    comparisons = ["S2_vs_S1", "S3_vs_S1", "S4_vs_S1", "S5_vs_S1"]
    titles = {
        "S2_vs_S1": "S2 versus S1",
        "S3_vs_S1": "S3 versus S1",
        "S4_vs_S1": "S4 versus S1 (supporting)",
        "S5_vs_S1": "S5 versus S1",
    }
    xmap = {"top_5pct": 5, "top_10pct": 10, "top_20pct": 20}
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.6), sharex=True, sharey=True)
    cmap = plt.get_cmap("tab10")
    for ax, comparison in zip(axes.ravel(), comparisons):
        for i, pair in enumerate(cfg["pair_order"]):
            sub = df[(df["comparison"] == comparison) & (df["pair"] == pair)].copy()
            sub["x"] = sub["cutoff"].map(xmap)
            sub = sub.sort_values("x")
            ax.plot(sub["x"], sub["delta_dual_fraction"], marker="o", color=cmap(i), label=pair, linewidth=1.2)
        ax.axhline(0, color="#666666", linewidth=0.7)
        ax.set_title(titles[comparison])
        ax.set_xticks([5, 10, 20])
    axes[1, 0].set_xlabel("Top-k cutoff (%)")
    axes[1, 1].set_xlabel("Top-k cutoff (%)")
    axes[0, 0].set_ylabel("Change in Top-k dual fraction")
    axes[1, 0].set_ylabel("Change in Top-k dual fraction")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Supporting analysis", fontsize=9, fontstyle="italic", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    _save(fig, "C_topk_dual_fraction_change")


def figure_c_states() -> None:
    df = _read("C_TOPK_CONCORDANCE_COUNTS.csv")
    states = [
        ("CONCORDANT_IMPROVE", "Concordant improvement", COLORS["improve"]),
        ("CONCORDANT_WORSEN", "Concordant worsening", COLORS["worsen"]),
        ("DISCORDANT", "Discordant", COLORS["discordant"]),
        ("NO_CHANGE", "No change", COLORS["nochange"]),
    ]
    cutoffs = ["top_5pct", "top_10pct", "top_20pct"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, scope, title in (
        (axes[0], "all_32", "All 32 comparisons"),
        (axes[1], "exclude_S4", "24 comparisons, S4 excluded"),
    ):
        width = 0.18
        for i, (state, label, color) in enumerate(states):
            heights = [
                int(df[(df["scope"] == scope) & (df["cutoff"] == cutoff) & (df["concordance"] == state)]["n"].iloc[0])
                for cutoff in cutoffs
            ]
            xpos = np.arange(len(cutoffs)) + (i - 1.5) * width
            ax.bar(xpos, heights, width=width, color=color, label=label)
        ax.set_xticks(np.arange(len(cutoffs)))
        ax.set_xticklabels(["Top 5%", "Top 10%", "Top 20%"])
        ax.set_title(title)
    axes[0].set_ylabel("Number of pair–method comparisons")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Supporting analysis", fontsize=9, fontstyle="italic", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    _save(fig, "C_topk_four_state_counts")


def figure_c_cases() -> None:
    df = _read("C_TOPK_METHOD_COMPARISONS.csv")
    cases = [
        ("F2/F10", "S3_vs_S1", "F2/F10, S3 versus S1"),
        ("PIK3CA/mTOR", "S2_vs_S1", "PIK3CA/mTOR, S2 versus S1"),
        ("PPARG/PPARA", "S3_vs_S1", "PPARG/PPARA, S3 versus S1"),
    ]
    cutoffs = ["top_5pct", "top_10pct", "top_20pct"]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 4.0), sharey=True)
    for ax, (pair, comparison, title) in zip(axes, cases):
        sub = df[(df["pair"] == pair) & (df["comparison"] == comparison)].set_index("cutoff")
        x = np.arange(len(cutoffs))
        s1 = [int(sub.loc[c, "S1_n_dual_top"]) for c in cutoffs]
        new = [int(sub.loc[c, "new_n_dual_top"]) for c in cutoffs]
        ax.bar(x - 0.16, s1, width=0.32, color=COLORS["S1"], label="S1")
        scheme = comparison.split("_")[0]
        ax.bar(x + 0.16, new, width=0.32, color="#0072B2", label=scheme)
        ax.set_xticks(x)
        ax.set_xticklabels(["5%", "10%", "20%"])
        ax.set_title(title, fontsize=9)
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    axes[0].set_ylabel("Dual ligands in the Top-k list")
    fig.suptitle("Supporting analysis", fontsize=9, fontstyle="italic", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, "C_representative_cases")


def make_all_figures() -> None:
    figure_a_similarity()
    figure_b_scores()
    figure_b_descriptors()
    figure_b_activity()
    figure_c_dual_fraction()
    figure_c_states()
    figure_c_cases()
