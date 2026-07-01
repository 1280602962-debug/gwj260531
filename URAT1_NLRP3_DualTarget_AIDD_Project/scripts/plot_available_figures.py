#!/usr/bin/env python3
"""Publication-quality figures from existing data (Arial 8 pt, clear axes, no overlap)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve, roc_curve

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "figures"))

from jmm_style import (  # noqa: E402
    FIG_OUT,
    GRID,
    MARGIN_COMPOSITE,
    MARGIN_SINGLE,
    MARGIN_WIDE_X,
    MUTED,
    NLRP3_COLOR,
    NEUTRAL,
    PROJECT_ROOT as ROOT,
    THRESHOLD,
    URAT1_COLOR,
    URAT1_LIGHT,
    apply_margins,
    apply_style,
    figsize_double,
    figsize_single,
    save_figure,
    set_axis_labels,
    style_axes,
    tag_panel,
)

DATA = ROOT / "data"
RESULTS = ROOT / "results"
MANIFEST_PATH = ROOT / "figures" / "generated" / "figure_manifest.json"

NLRP3_CONTROLS = [
    "verinurad",
    "colchicine",
    "lesinurad",
    "benzbromarone",
    "dotinurad",
    "allopurinol",
    "febuxostat",
]


def _name_col(df: pd.DataFrame) -> str:
    return "name" if "name" in df.columns else "pref_name"


def _headroom_ylim(ax, frac: float = 0.18) -> None:
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + (hi - lo) * frac if hi > lo else hi + 1)


def plot_nlrp3_distribution(ml: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(72))
    fig.text(0.5, 0.98, "Target: NLRP3 (assay-conditioned ML)", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    n_ge = int((ml["p_active_nlrp3"] >= 0.5).sum())
    ax.hist(ml["p_active_nlrp3"], bins=50, color=NLRP3_COLOR, alpha=0.8, edgecolor="white", linewidth=0.25)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9, label=f"Threshold 0.5 (n = {n_ge:,})")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Predicted probability P(active)", "Number of compounds")
    ax.set_xlim(0, 1.0)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0))
    apply_margins(fig, MARGIN_SINGLE)
    paths = save_figure(fig, "nlrp3_fig02a_pactive_histogram", "nlrp3")
    return {"id": "nlrp3_fig02a", "target": "NLRP3", "description": "P(active) distribution (n=8319)", **paths}


def plot_nlrp3_controls(ml: pd.DataFrame) -> dict:
    name_col = _name_col(ml)
    rows = []
    for drug in NLRP3_CONTROLS:
        hit = ml[ml[name_col].astype(str).str.upper() == drug.upper()]
        if len(hit):
            rows.append({"drug": drug, "p_active": float(hit["p_active_nlrp3"].iloc[0])})
    cdf = pd.DataFrame(rows).sort_values("p_active", ascending=True)

    fig, ax = plt.subplots(figsize=figsize_single(82))
    fig.text(0.5, 0.98, "Target: NLRP3 (assay-conditioned ML)", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    colors = [NLRP3_COLOR if p >= 0.5 else MUTED for p in cdf["p_active"]]
    bars = ax.barh(cdf["drug"], cdf["p_active"], color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.62)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9, label="Threshold 0.5")
    ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Predicted probability P(active)", "Reference compound")
    ax.set_xlim(0, 1.12)
    ax.legend(loc="lower right")
    apply_margins(fig, {**MARGIN_SINGLE, "left": 0.22})
    paths = save_figure(fig, "nlrp3_fig02b_control_compounds", "nlrp3")
    return {"id": "nlrp3_fig02b", "target": "NLRP3", "description": "Control compound scores", **paths}


def plot_nlrp3_phase(ml: pd.DataFrame) -> dict:
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = ml[ml["max_phase"].isin(phase_keys)]
    data = [sub.loc[sub["max_phase"] == k, "p_active_nlrp3"].values for k in phase_keys]

    fig, ax = plt.subplots(figsize=figsize_single(76))
    fig.text(0.5, 0.98, "Target: NLRP3 (assay-conditioned ML)", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    bp = ax.boxplot(data, tick_labels=phase_labels, patch_artist=True, widths=0.5, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.55)
        patch.set_edgecolor(NEUTRAL)
    ax.axhline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9, label="Threshold 0.5")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Clinical development stage", "Predicted probability P(active)")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.legend(loc="upper right")
    apply_margins(fig, MARGIN_WIDE_X)
    paths = save_figure(fig, "nlrp3_fig02c_phase_boxplot", "nlrp3")
    return {"id": "nlrp3_fig02c", "target": "NLRP3", "description": "P(active) by clinical phase", **paths}


def plot_screening_funnel(summary: dict) -> dict:
    labels = ["Clinical library scored", "P(active) ≥ 0.5 pool", "Dual-target docking"]
    counts = [summary["n_scored"], summary["n_pred_active_ge_threshold"], summary["n_pred_active_ge_threshold"]]
    colors = [NLRP3_COLOR, NLRP3_COLOR, MUTED]

    fig, ax = plt.subplots(figsize=figsize_single(58))
    fig.text(0.5, 0.98, "Target: NLRP3 prescreen funnel", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    y = np.arange(len(labels))
    bars = ax.barh(y, counts, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.55)
    ax.bar_label(bars, labels=[f"{c:,}" if i < 2 else "pending" for i, c in enumerate(counts)], padding=4, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Number of compounds", "Screening stage")
    ax.set_xlim(0, counts[0] * 1.18)
    ax.invert_yaxis()
    apply_margins(fig, {**MARGIN_SINGLE, "left": 0.34})
    paths = save_figure(fig, "nlrp3_fig02d_screening_funnel", "nlrp3")
    return {"id": "nlrp3_fig02d", "target": "NLRP3", "description": "NLRP3 prescreen funnel", **paths}


def plot_urat1_a_vs_d_violin(dock: pd.DataFrame) -> dict:
    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]]
    d_scores = sub.loc[sub["subset"] == "D", "glide_score_xp"].dropna()
    a_scores = sub.loc[sub["subset"] == "A", "glide_score_xp"].dropna()

    fig, ax = plt.subplots(figsize=figsize_single(72))
    fig.text(0.5, 0.98, "Target: URAT1 (SLC22A12, 9DKB XP)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    parts = ax.violinplot([d_scores, a_scores], positions=[0, 1], showmedians=True, widths=0.65)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(URAT1_COLOR if i == 1 else URAT1_LIGHT)
        body.set_alpha(0.75)
        body.set_edgecolor(NEUTRAL)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Decoy subset D (n={len(d_scores):,})", f"Active subset A (n={len(a_scores):,})"])
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "8973 distill subset label", "Glide XP score (kcal/mol)")
    apply_margins(fig, MARGIN_WIDE_X)
    paths = save_figure(fig, "urat1_fig03a_a_vs_d_violin", "urat1")
    return {"id": "urat1_fig03a", "target": "URAT1", "description": "8973 A vs D docking scores", **paths}


def plot_urat1_roc(dock: pd.DataFrame) -> dict:
    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]]
    y = (sub["subset"] == "A").astype(int)
    fpr, tpr, _ = roc_curve(y, sub["s_u_percentile"])
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=figsize_single(72))
    fig.text(0.5, 0.98, "Target: URAT1 (SLC22A12, 9DKB XP)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    ax.plot(fpr, tpr, color=URAT1_COLOR, linewidth=1.0, label=f"Docking score ROC (AUC = {roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.7, label="Random classifier")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "False positive rate (decoy subset D)", "True positive rate (active subset A)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right")
    apply_margins(fig, MARGIN_SINGLE)
    paths = save_figure(fig, "urat1_fig03b_a_vs_d_roc", "urat1")
    return {"id": "urat1_fig03b", "target": "URAT1", "description": "8973 enrichment ROC", **paths}


def plot_urat1_ef(summary: dict) -> dict:
    enr = summary["enrichment"]
    labels = ["EF @ 5%", "EF @ 10%"]
    vals = [enr["ef_5pct_a_vs_d"], enr["ef_10pct_a_vs_d"]]

    fig, ax = plt.subplots(figsize=figsize_single(62))
    fig.text(0.5, 0.98, "Target: URAT1 (SLC22A12, 9DKB XP)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    bars = ax.bar(labels, vals, color=URAT1_COLOR, edgecolor=NEUTRAL, linewidth=0.3, width=0.48)
    ax.axhline(1.0, color=GRID, linestyle="--", linewidth=0.7, label="No enrichment (EF = 1)")
    ax.bar_label(bars, fmt="%.2f", padding=4, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Enrichment metric (8973 A vs D)", "Enrichment factor")
    ax.legend(loc="upper right")
    _headroom_ylim(ax, 0.22)
    apply_margins(fig, MARGIN_SINGLE)
    paths = save_figure(fig, "urat1_fig03c_enrichment_factor", "urat1")
    return {"id": "urat1_fig03c", "target": "URAT1", "description": "EF@5% and EF@10%", **paths}


def plot_urat1_benchmark(bench: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(78))
    fig.text(0.5, 0.98, "Target: URAT1 (SLC22A12, 9DKB XP)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    x = np.arange(len(bench))
    w = 0.36
    b1 = ax.bar(x - w / 2, bench["s_u_percentile"], width=w, label="9DKB XP docking percentile", color=URAT1_COLOR, edgecolor=NEUTRAL, linewidth=0.3)
    b2 = ax.bar(x + w / 2, bench["ml_percentile_vs_8973"], width=w, label="URAT1 ML percentile", color=URAT1_LIGHT, edgecolor=NEUTRAL, linewidth=0.3)
    ax.axhline(90, color=THRESHOLD, linestyle="--", linewidth=0.7, label="90th percentile (top 10%)")
    ax.set_xticks(x)
    ax.set_xticklabels(bench["compound"], rotation=35, ha="right")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Benchmark URAT1 inhibitor", "Percentile rank on 8973 library")
    ax.set_ylim(0, 108)
    ax.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 1.0))
    apply_margins(fig, MARGIN_WIDE_X)
    paths = save_figure(fig, "urat1_fig03d_benchmark_ml_vs_docking", "urat1")
    return {"id": "urat1_fig03d", "target": "URAT1", "description": "Four-drug ML vs docking recovery", **paths}


def plot_nlrp3_roc_pr(oof: pd.DataFrame) -> dict:
    fpr, tpr, _ = roc_curve(oof["y_true"], oof["y_prob"])
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(oof["y_true"], oof["y_prob"])
    pr_auc = auc(rec, prec)

    fig, axes = plt.subplots(1, 2, figsize=figsize_double(72))
    fig.text(0.5, 0.98, "Target: NLRP3 (scaffold-CV out-of-fold)", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    ax = axes[0]
    ax.plot(fpr, tpr, color=NLRP3_COLOR, linewidth=1.0, label=f"ROC (AUC = {roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.7)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "False positive rate", "True positive rate")
    ax.legend(loc="lower right")
    tag_panel(ax, "a")

    ax = axes[1]
    ax.plot(rec, prec, color=NLRP3_COLOR, linewidth=1.0, label=f"PR (AUPRC = {pr_auc:.3f})")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Recall", "Precision")
    ax.legend(loc="upper right")
    tag_panel(ax, "b")

    apply_margins(fig, {**MARGIN_COMPOSITE, "top": 0.84, "bottom": 0.18})
    paths = save_figure(fig, "si_nlrp3_oof_roc_pr", "si")
    return {"id": "si_nlrp3_roc_pr", "target": "NLRP3", "description": "NLRP3 scaffold-CV OOF ROC/PR", **paths}


def plot_urat1_parity(oof: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(72))
    fig.text(0.5, 0.98, "Target: URAT1 (scaffold-CV out-of-fold)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    ax.scatter(oof["y_true"], oof["y_pred"], s=10, alpha=0.4, color=URAT1_COLOR, edgecolors="none", rasterized=True)
    lo = min(oof["y_true"].min(), oof["y_pred"].min()) - 0.2
    hi = max(oof["y_true"].max(), oof["y_pred"].max()) + 0.2
    ax.plot([lo, hi], [lo, hi], color=NEUTRAL, linestyle="--", linewidth=0.7, label="Ideal (y = x)")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Observed pActivity (ChEMBL training set)", "Predicted pActivity (OOF)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right")
    apply_margins(fig, MARGIN_SINGLE)
    paths = save_figure(fig, "si_urat1_oof_parity", "si")
    return {"id": "si_urat1_parity", "target": "URAT1", "description": "URAT1 OOF parity plot", **paths}


def plot_data_asymmetry(summary: dict) -> dict:
    fig, axes = plt.subplots(1, 2, figsize=figsize_double(68))
    fig.text(0.5, 0.98, "Dual-target asymmetric data design", ha="center", va="top", fontsize=8, fontweight="bold", color=NEUTRAL, transform=fig.transFigure)

    ax = axes[0]
    labels = ["URAT1", "NLRP3", "Shared SMILES"]
    vals = [summary["urat1"]["n_compounds"], summary["nlrp3"]["n_compounds"], summary["overlap_smiles"]]
    colors = [URAT1_COLOR, NLRP3_COLOR, MUTED]
    bars = ax.bar(labels, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3, width=0.52)
    ax.bar_label(bars, labels=[str(v) for v in vals], padding=4, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "ChEMBL training set", "Compound count")
    tag_panel(ax, "a")

    ax = axes[1]
    routes = ["URAT1 retrospective\n(8973 @ 9DKB)", "NLRP3 prescreen\n(8319 clinical)"]
    ns = [8973, 8319]
    bars = ax.bar(routes, ns, color=[URAT1_COLOR, NLRP3_COLOR], edgecolor=NEUTRAL, linewidth=0.3, width=0.48)
    ax.bar_label(bars, labels=[f"{v:,}" for v in ns], padding=4, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Primary computational evidence route", "Library size (compounds)")
    tag_panel(ax, "b")

    apply_margins(fig, {**MARGIN_COMPOSITE, "top": 0.84, "bottom": 0.22})
    paths = save_figure(fig, "si_data_asymmetry", "si")
    return {"id": "si_data_asymmetry", "target": "both", "description": "Asymmetric dual-target data design", **paths}


def plot_library_phase(manifest: pd.DataFrame) -> dict:
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = manifest[manifest["max_phase"].isin(phase_keys)]
    counts = sub["max_phase"].value_counts().reindex(phase_keys).fillna(0).astype(int)
    counts.index = phase_labels

    fig, ax = plt.subplots(figsize=figsize_single(68))
    fig.text(0.5, 0.98, "NLRP3 prescreen input library (ChEMBL clinical)", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    bars = ax.bar(counts.index, counts.values, color=NLRP3_COLOR, alpha=0.8, edgecolor=NEUTRAL, linewidth=0.3, width=0.55)
    ax.bar_label(bars, padding=4, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Clinical development stage", "Number of compounds")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    _headroom_ylim(ax, 0.12)
    apply_margins(fig, MARGIN_WIDE_X)
    paths = save_figure(fig, "si_repurposing_library_phase", "si")
    return {"id": "si_library_phase", "target": "NLRP3", "description": "Clinical library phase composition", **paths}


def plot_fig02_composite(ml: pd.DataFrame, summary: dict) -> dict:
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(140))
    fig.text(0.5, 0.99, "Target: NLRP3 — clinical library machine-learning prescreen", ha="center", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR, transform=fig.transFigure)

    # (a) histogram
    ax = axes[0, 0]
    n_ge = int((ml["p_active_nlrp3"] >= 0.5).sum())
    ax.hist(ml["p_active_nlrp3"], bins=45, color=NLRP3_COLOR, alpha=0.8, edgecolor="white", linewidth=0.2)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9, label=f"Threshold 0.5 (n={n_ge:,})")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Predicted probability P(active)", "Number of compounds")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.02, 0.98))
    tag_panel(ax, "a")

    # (b) controls
    ax = axes[0, 1]
    name_col = _name_col(ml)
    rows = []
    for drug in ["verinurad", "colchicine", "lesinurad", "benzbromarone", "allopurinol"]:
        hit = ml[ml[name_col].astype(str).str.upper() == drug.upper()]
        if len(hit):
            rows.append((drug, float(hit["p_active_nlrp3"].iloc[0])))
    rows.sort(key=lambda x: x[1])
    drugs, vals = zip(*rows) if rows else ([], [])
    colors = [NLRP3_COLOR if v >= 0.5 else MUTED for v in vals]
    bars = ax.barh(drugs, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.58)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9)
    ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Predicted probability P(active)", "Reference compound")
    ax.set_xlim(0, 1.15)
    tag_panel(ax, "b")

    # (c) phase boxplot
    ax = axes[1, 0]
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = ml[ml["max_phase"].isin(phase_keys)]
    data = [sub.loc[sub["max_phase"] == k, "p_active_nlrp3"].values for k in phase_keys]
    bp = ax.boxplot(data, tick_labels=phase_labels, patch_artist=True, widths=0.48, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.55)
    ax.axhline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.9)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Clinical development stage", "Predicted probability P(active)")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    tag_panel(ax, "c")

    # (d) funnel bars
    ax = axes[1, 1]
    stages = ["Library scored", "P(active) ≥ 0.5", "Dual docking"]
    counts = [summary["n_scored"], summary["n_pred_active_ge_threshold"], summary["n_pred_active_ge_threshold"]]
    colors = [NLRP3_COLOR, NLRP3_COLOR, MUTED]
    y = np.arange(len(stages))
    bars = ax.barh(y, counts, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.5)
    ax.bar_label(bars, labels=[f"{c:,}" if i < 2 else "pending" for i, c in enumerate(counts)], padding=4, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Number of compounds", "Screening stage")
    ax.set_xlim(0, counts[0] * 1.15)
    ax.invert_yaxis()
    tag_panel(ax, "d")

    apply_margins(fig, {**MARGIN_COMPOSITE, "top": 0.90, "bottom": 0.14, "left": 0.12, "right": 0.98})
    paths = save_figure(fig, "fig02_nlrp3_screening_composite", "main")
    return {"id": "fig02_composite", "target": "NLRP3", "description": "Main Fig 2 composite (NLRP3 screening)", **paths}


def plot_fig03_composite(dock: pd.DataFrame, summary: dict, bench: pd.DataFrame) -> dict:
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(140))
    fig.text(0.5, 0.99, "Target: URAT1 — 8973 retrospective docking validation (9DKB XP)", ha="center", va="top", fontsize=8, fontweight="bold", color=URAT1_COLOR, transform=fig.transFigure)

    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]]
    d_scores = sub.loc[sub["subset"] == "D", "glide_score_xp"]
    a_scores = sub.loc[sub["subset"] == "A", "glide_score_xp"]

    # (a) violin
    ax = axes[0, 0]
    parts = ax.violinplot([d_scores, a_scores], positions=[0, 1], showmedians=True, widths=0.62)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(URAT1_COLOR if i == 1 else URAT1_LIGHT)
        body.set_alpha(0.75)
        body.set_edgecolor(NEUTRAL)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Decoy D", "Active A"])
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "8973 distill subset", "Glide XP score (kcal/mol)")
    tag_panel(ax, "a")

    # (b) ROC
    ax = axes[0, 1]
    y = (sub["subset"] == "A").astype(int)
    fpr, tpr, _ = roc_curve(y, sub["s_u_percentile"])
    roc_auc = auc(fpr, tpr)
    ax.plot(fpr, tpr, color=URAT1_COLOR, linewidth=1.0, label=f"AUC = {roc_auc:.3f}")
    ax.plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.7)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "False positive rate (decoy D)", "True positive rate (active A)")
    ax.legend(loc="lower right", fontsize=8)
    tag_panel(ax, "b")

    # (c) EF
    ax = axes[1, 0]
    enr = summary["enrichment"]
    bars = ax.bar(["EF @ 5%", "EF @ 10%"], [enr["ef_5pct_a_vs_d"], enr["ef_10pct_a_vs_d"]], color=URAT1_COLOR, width=0.48, edgecolor=NEUTRAL, linewidth=0.3)
    ax.axhline(1, color=GRID, linestyle="--", linewidth=0.7)
    ax.bar_label(bars, fmt="%.2f", padding=4, fontsize=8)
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Enrichment metric", "Enrichment factor (A vs D)")
    _headroom_ylim(ax, 0.2)
    tag_panel(ax, "c")

    # (d) benchmark
    ax = axes[1, 1]
    x = np.arange(len(bench))
    w = 0.36
    ax.bar(x - w / 2, bench["s_u_percentile"], width=w, color=URAT1_COLOR, label="Docking", edgecolor=NEUTRAL, linewidth=0.3)
    ax.bar(x + w / 2, bench["ml_percentile_vs_8973"], width=w, color=URAT1_LIGHT, label="ML", edgecolor=NEUTRAL, linewidth=0.3)
    ax.axhline(90, color=THRESHOLD, linestyle="--", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(bench["compound"], rotation=35, ha="right")
    style_axes(ax, grid_y=True)
    set_axis_labels(ax, "Benchmark inhibitor", "Percentile on 8973 library")
    ax.set_ylim(0, 108)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    tag_panel(ax, "d")

    apply_margins(fig, {**MARGIN_COMPOSITE, "top": 0.90, "bottom": 0.16, "left": 0.11, "right": 0.98, "hspace": 0.58, "wspace": 0.45})
    paths = save_figure(fig, "fig03_urat1_retrospective_composite", "main")
    return {"id": "fig03_composite", "target": "URAT1", "description": "Main Fig 3 composite (URAT1 retrospective)", **paths}


def main() -> None:
    apply_style()

    ml = pd.read_csv(DATA / "repurposing" / "screening" / "nlrp3_ml_scores_clinical_all.csv", low_memory=False)
    nlrp3_summary = json.loads((DATA / "repurposing" / "screening" / "nlrp3_screening_summary_clinical_all.json").read_text())
    dock = pd.read_csv(DATA / "docking" / "8973_9DKB_with_manifest.csv", low_memory=False)
    urat1_summary = json.loads((DATA / "docking" / "urat1_docking_vs_ml_summary.json").read_text())
    bench = pd.read_csv(DATA / "docking" / "urat1_benchmark_rankings_docking.csv")
    data_summary = json.loads((DATA / "processed" / "data_summary.json").read_text())
    manifest = pd.read_csv(DATA / "repurposing" / "repurposing_manifest.csv", low_memory=False)

    nlrp3_oof = pd.read_csv(RESULTS / "training" / "nlrp3_oof_predictions.csv")
    urat1_oof = pd.read_csv(RESULTS / "training" / "urat1_oof_predictions.csv")

    entries = [
        plot_nlrp3_distribution(ml),
        plot_nlrp3_controls(ml),
        plot_nlrp3_phase(ml),
        plot_screening_funnel(nlrp3_summary),
        plot_urat1_a_vs_d_violin(dock),
        plot_urat1_roc(dock),
        plot_urat1_ef(urat1_summary),
        plot_urat1_benchmark(bench),
        plot_nlrp3_roc_pr(nlrp3_oof),
        plot_urat1_parity(urat1_oof),
        plot_data_asymmetry(data_summary),
        plot_library_phase(manifest),
        plot_fig02_composite(ml, nlrp3_summary),
        plot_fig03_composite(dock, urat1_summary, bench),
    ]

    manifest_out = {
        "style": {
            "base": "SciencePlots science+no-latex",
            "font": "Arial",
            "size_pt": 8,
            "urat1_color": URAT1_COLOR,
            "nlrp3_color": NLRP3_COLOR,
            "notes": "Target label in figure header; panel tags outside axes; stats in legend",
        },
        "figures": entries,
    }
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest_out, indent=2))
    print(json.dumps(manifest_out, indent=2))
    print(f"\nWrote {len(entries)} figure sets to {FIG_OUT}")


if __name__ == "__main__":
    main()
