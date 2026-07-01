#!/usr/bin/env python3
"""Generate manuscript figures from existing project data (Arial 8 pt, URAT1/NLRP3 coded)."""
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
    NLRP3_COLOR,
    NEUTRAL,
    PROJECT_ROOT as ROOT,
    TARGET_STYLES,
    THRESHOLD,
    URAT1_COLOR,
    add_target_banner,
    apply_style,
    figsize_double,
    figsize_single,
    panel_label,
    save_figure,
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
    "MCC950",
    "GDC-2394",
]
URAT1_BENCHMARKS = ["lesinurad", "benzbromarone", "verinurad", "dotinurad"]


def _name_col(df: pd.DataFrame) -> str:
    return "name" if "name" in df.columns else "pref_name"


def plot_nlrp3_distribution(ml: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(65))
    add_target_banner(ax, "NLRP3")
    ax.hist(ml["p_active_nlrp3"], bins=60, color=NLRP3_COLOR, alpha=0.75, edgecolor="white", linewidth=0.3)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8, label="P(active) = 0.5")
    n_ge = int((ml["p_active_nlrp3"] >= 0.5).sum())
    ax.text(
        0.52,
        ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] else 1,
        f"n = {n_ge}",
        fontsize=8,
        color=THRESHOLD,
    )
    ax.set_xlabel("NLRP3 ML P(active)")
    ax.set_ylabel("Clinical library count")
    ax.set_xlim(-0.02, 1.02)
    paths = save_figure(fig, "nlrp3_fig02a_pactive_histogram", "nlrp3")
    return {"id": "nlrp3_fig02a", "target": "NLRP3", "description": "P(active) distribution (n=8319)", **paths}


def plot_nlrp3_controls(ml: pd.DataFrame) -> dict:
    name_col = _name_col(ml)
    rows = []
    for drug in NLRP3_CONTROLS:
        hit = ml[ml[name_col].astype(str).str.upper() == drug.upper()]
        if len(hit):
            rows.append(
                {
                    "drug": drug,
                    "p_active": float(hit["p_active_nlrp3"].iloc[0]),
                    "percentile": float(hit["nlrp3_percentile"].iloc[0]),
                }
            )
    cdf = pd.DataFrame(rows).sort_values("p_active", ascending=True)

    fig, ax = plt.subplots(figsize=figsize_single(75))
    add_target_banner(ax, "NLRP3")
    colors = [NLRP3_COLOR if p >= 0.5 else "#BBBBBB" for p in cdf["p_active"]]
    ax.barh(cdf["drug"], cdf["p_active"], color=colors, edgecolor=NEUTRAL, linewidth=0.3)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8)
    ax.set_xlabel("NLRP3 ML P(active)")
    ax.set_ylabel("Reference compound")
    ax.set_xlim(0, 1.05)
    paths = save_figure(fig, "nlrp3_fig02b_control_compounds", "nlrp3")
    return {"id": "nlrp3_fig02b", "target": "NLRP3", "description": "Control compound scores", **paths}


def plot_nlrp3_phase(ml: pd.DataFrame) -> dict:
    phase_map = {1.0: "Phase I", 2.0: "Phase II", 3.0: "Phase III", 4.0: "Approved", -1.0: "Unknown"}
    sub = ml[ml["max_phase"].isin(phase_map)].copy()
    sub["phase_label"] = sub["max_phase"].map(phase_map)
    order = ["Phase I", "Phase II", "Phase III", "Approved"]
    data = [sub.loc[sub["phase_label"] == p, "p_active_nlrp3"].values for p in order]

    fig, ax = plt.subplots(figsize=figsize_single(70))
    add_target_banner(ax, "NLRP3")
    bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55)
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.55)
        patch.set_edgecolor(NEUTRAL)
    ax.axhline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8)
    ax.set_ylabel("NLRP3 ML P(active)")
    ax.set_xlabel("Clinical development stage")
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    paths = save_figure(fig, "nlrp3_fig02c_phase_boxplot", "nlrp3")
    return {"id": "nlrp3_fig02c", "target": "NLRP3", "description": "P(active) by clinical phase", **paths}


def plot_screening_funnel(summary: dict) -> dict:
    stages = ["Clinical library\n(NLRP3 ML)", "P(active) ≥ 0.5\n(docking pool)", "Dual docking\n(pending)"]
    counts = [summary["n_scored"], summary["n_pred_active_ge_threshold"], None]
    colors = [NLRP3_COLOR, NLRP3_COLOR, "#BBBBBB"]

    fig, ax = plt.subplots(figsize=figsize_single(55))
    add_target_banner(ax, "NLRP3")
    y = np.arange(len(stages))
    widths = [counts[0], counts[1], counts[1] * 0.35]
    for i, (stage, w, c) in enumerate(zip(stages, widths, colors)):
        ax.barh(i, w, color=c, edgecolor=NEUTRAL, linewidth=0.3, height=0.6)
        label = f"{counts[i]:,}" if counts[i] is not None else "—"
        ax.text(w + counts[0] * 0.02, i, label, va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.set_xlabel("Compound count")
    ax.set_xlim(0, counts[0] * 1.15)
    ax.invert_yaxis()
    paths = save_figure(fig, "nlrp3_fig02d_screening_funnel", "nlrp3")
    return {"id": "nlrp3_fig02d", "target": "NLRP3", "description": "NLRP3 prescreen funnel", **paths}


def plot_urat1_a_vs_d_violin(dock: pd.DataFrame) -> dict:
    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]].copy()
    a = sub.loc[sub["subset"] == "A", "glide_score_xp"].dropna()
    d = sub.loc[sub["subset"] == "D", "glide_score_xp"].dropna()

    fig, ax = plt.subplots(figsize=figsize_single(65))
    add_target_banner(ax, "URAT1")
    parts = ax.violinplot([d, a], positions=[0, 1], showmeans=False, showmedians=True, widths=0.7)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(URAT1_COLOR if i == 1 else "#88BCE4")
        body.set_alpha(0.7)
        body.set_edgecolor(NEUTRAL)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Decoy D\n(n={len(d)})", f"Active A\n(n={len(a)})"])
    ax.set_ylabel("URAT1 Glide XP score (9DKB)")
    ax.set_xlabel("8973 distill subset")
    paths = save_figure(fig, "urat1_fig03a_a_vs_d_violin", "urat1")
    return {"id": "urat1_fig03a", "target": "URAT1", "description": "8973 A vs D docking scores", **paths}


def plot_urat1_roc(dock: pd.DataFrame) -> dict:
    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]].copy()
    y = (sub["subset"] == "A").astype(int)
    scores = sub["s_u_percentile"]
    fpr, tpr, _ = roc_curve(y, scores)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=figsize_single(65))
    add_target_banner(ax, "URAT1")
    ax.plot(fpr, tpr, color=URAT1_COLOR, linewidth=1.0, label=f"AUC = {roc_auc:.3f}")
    ax.plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.6)
    ax.set_xlabel("False positive rate (decoy D)")
    ax.set_ylabel("True positive rate (active A)")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    paths = save_figure(fig, "urat1_fig03b_a_vs_d_roc", "urat1")
    return {"id": "urat1_fig03b", "target": "URAT1", "description": "8973 enrichment ROC", **paths}


def plot_urat1_ef(summary: dict) -> dict:
    enr = summary["enrichment"]
    labels = ["EF @ 5%", "EF @ 10%"]
    vals = [enr["ef_5pct_a_vs_d"], enr["ef_10pct_a_vs_d"]]

    fig, ax = plt.subplots(figsize=figsize_single(50))
    add_target_banner(ax, "URAT1")
    ax.bar(labels, vals, color=URAT1_COLOR, edgecolor=NEUTRAL, linewidth=0.3, width=0.55)
    ax.axhline(1.0, color=GRID, linestyle="--", linewidth=0.6)
    ax.set_ylabel("Enrichment factor")
    ax.set_xlabel("8973 active A vs decoy D")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.08, f"{v:.2f}", ha="center", fontsize=8)
    paths = save_figure(fig, "urat1_fig03c_enrichment_factor", "urat1")
    return {"id": "urat1_fig03c", "target": "URAT1", "description": "EF@5% and EF@10%", **paths}


def plot_urat1_benchmark(bench: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(70))
    add_target_banner(ax, "URAT1")
    x = np.arange(len(bench))
    w = 0.35
    ax.bar(x - w / 2, bench["s_u_percentile"], width=w, label="9DKB XP percentile", color=URAT1_COLOR, edgecolor=NEUTRAL, linewidth=0.3)
    ax.bar(x + w / 2, bench["ml_percentile_vs_8973"], width=w, label="URAT1 ML percentile", color="#88BCE4", edgecolor=NEUTRAL, linewidth=0.3)
    ax.axhline(90, color=THRESHOLD, linestyle="--", linewidth=0.6, label="Top 10%")
    ax.set_xticks(x)
    ax.set_xticklabels(bench["compound"], rotation=20, ha="right")
    ax.set_ylabel("Percentile on 8973 library")
    ax.set_xlabel("URAT1 benchmark inhibitor")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.set_ylim(0, 105)
    paths = save_figure(fig, "urat1_fig03d_benchmark_ml_vs_docking", "urat1")
    return {"id": "urat1_fig03d", "target": "URAT1", "description": "Four-drug ML vs docking recovery", **paths}


def plot_nlrp3_roc_pr(oof: pd.DataFrame) -> dict:
    fpr, tpr, _ = roc_curve(oof["y_true"], oof["y_prob"])
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(oof["y_true"], oof["y_prob"])
    pr_auc = auc(rec, prec)

    fig, axes = plt.subplots(1, 2, figsize=figsize_double(65))
    for ax in axes:
        add_target_banner(ax, "NLRP3")

    axes[0].plot(fpr, tpr, color=NLRP3_COLOR, linewidth=1.0)
    axes[0].plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.6)
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")
    axes[0].text(0.55, 0.15, f"AUC = {roc_auc:.3f}", fontsize=8)
    panel_label(axes[0], "a")

    axes[1].plot(rec, prec, color=NLRP3_COLOR, linewidth=1.0)
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].text(0.45, 0.05, f"AUPRC = {pr_auc:.3f}", fontsize=8)
    panel_label(axes[1], "b")

    fig.subplots_adjust(wspace=0.35)
    paths = save_figure(fig, "si_nlrp3_oof_roc_pr", "si")
    return {"id": "si_nlrp3_roc_pr", "target": "NLRP3", "description": "NLRP3 scaffold-CV OOF ROC/PR", **paths}


def plot_urat1_parity(oof: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(65))
    add_target_banner(ax, "URAT1")
    ax.scatter(oof["y_true"], oof["y_pred"], s=8, alpha=0.45, color=URAT1_COLOR, edgecolors="none")
    lo = min(oof["y_true"].min(), oof["y_pred"].min()) - 0.3
    hi = max(oof["y_true"].max(), oof["y_pred"].max()) + 0.3
    ax.plot([lo, hi], [lo, hi], color=NEUTRAL, linestyle="--", linewidth=0.6)
    ax.set_xlabel("Observed pActivity (ChEMBL)")
    ax.set_ylabel("Predicted pActivity (OOF)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    paths = save_figure(fig, "si_urat1_oof_parity", "si")
    return {"id": "si_urat1_parity", "target": "URAT1", "description": "URAT1 OOF parity plot", **paths}


def plot_data_asymmetry(summary: dict) -> dict:
    fig, axes = plt.subplots(1, 2, figsize=figsize_double(55))

    # Panel a: training set sizes
    ax = axes[0]
    add_target_banner(ax, "URAT1", x=0.02, y=0.98)
    ax.text(0.52, 0.98, TARGET_STYLES["NLRP3"]["label"], transform=ax.transAxes, ha="left", va="top", fontsize=8, fontweight="bold", color=NLRP3_COLOR,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=NLRP3_COLOR, linewidth=0.6))
    labels = ["URAT1\ncompounds", "NLRP3\ncompounds", "Shared\nSMILES"]
    vals = [summary["urat1"]["n_compounds"], summary["nlrp3"]["n_compounds"], summary["overlap_smiles"]]
    colors = [URAT1_COLOR, NLRP3_COLOR, NEUTRAL]
    ax.bar(labels, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3, width=0.55)
    ax.set_ylabel("Count")
    ax.set_xlabel("ChEMBL training sets")
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.03, str(v), ha="center", fontsize=8)
    panel_label(ax, "a")

    # Panel b: evidence route
    ax = axes[1]
    routes = ["URAT1\n(8973 retrospective)", "NLRP3\n(8319 prescreen)"]
    ns = [8973, 8319]
    ax.bar(routes, ns, color=[URAT1_COLOR, NLRP3_COLOR], edgecolor=NEUTRAL, linewidth=0.3, width=0.5)
    ax.set_ylabel("Library size")
    ax.set_xlabel("Primary evidence route")
    for i, v in enumerate(ns):
        ax.text(i, v + 150, f"{v:,}", ha="center", fontsize=8)
    panel_label(ax, "b")

    fig.subplots_adjust(wspace=0.4)
    paths = save_figure(fig, "si_data_asymmetry", "si")
    return {"id": "si_data_asymmetry", "target": "both", "description": "Asymmetric dual-target data design", **paths}


def plot_library_phase(manifest: pd.DataFrame) -> dict:
    phase_map = {1.0: "Phase I", 2.0: "Phase II", 3.0: "Phase III", 4.0: "Approved"}
    sub = manifest[manifest["max_phase"].isin(phase_map)].copy()
    sub["phase_label"] = sub["max_phase"].map(phase_map)
    counts = sub["phase_label"].value_counts().reindex(["Phase I", "Phase II", "Phase III", "Approved"])

    fig, ax = plt.subplots(figsize=figsize_single(55))
    ax.text(0.02, 0.98, "Repurposing library (NLRP3 prescreen input)", transform=ax.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold", color=NLRP3_COLOR,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=NLRP3_COLOR, linewidth=0.6))
    ax.bar(counts.index, counts.values, color=NLRP3_COLOR, alpha=0.75, edgecolor=NEUTRAL, linewidth=0.3)
    ax.set_ylabel("Compound count")
    ax.set_xlabel("Clinical development stage")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    paths = save_figure(fig, "si_repurposing_library_phase", "si")
    return {"id": "si_library_phase", "target": "NLRP3", "description": "Clinical library phase composition", **paths}


def plot_fig02_composite(ml: pd.DataFrame, summary: dict) -> dict:
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(130))
    # a: histogram
    ax = axes[0, 0]
    add_target_banner(ax, "NLRP3")
    ax.hist(ml["p_active_nlrp3"], bins=50, color=NLRP3_COLOR, alpha=0.75, edgecolor="white", linewidth=0.2)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8)
    ax.set_xlabel("NLRP3 ML P(active)")
    ax.set_ylabel("Count")
    panel_label(ax, "a")

    # b: controls
    ax = axes[0, 1]
    add_target_banner(ax, "NLRP3")
    name_col = _name_col(ml)
    rows = []
    for drug in ["verinurad", "colchicine", "lesinurad", "benzbromarone", "allopurinol"]:
        hit = ml[ml[name_col].astype(str).str.upper() == drug.upper()]
        if len(hit):
            rows.append((drug, float(hit["p_active_nlrp3"].iloc[0])))
    rows.sort(key=lambda x: x[1])
    drugs, vals = zip(*rows) if rows else ([], [])
    colors = [NLRP3_COLOR if v >= 0.5 else "#BBBBBB" for v in vals]
    ax.barh(drugs, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3)
    ax.axvline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8)
    ax.set_xlabel("NLRP3 ML P(active)")
    panel_label(ax, "b")

    # c: phase box
    ax = axes[1, 0]
    add_target_banner(ax, "NLRP3")
    phase_map = {1.0: "I", 2.0: "II", 3.0: "III", 4.0: "IV"}
    sub = ml[ml["max_phase"].isin(phase_map)].copy()
    order_keys = [1.0, 2.0, 3.0, 4.0]
    order = ["I", "II", "III", "IV"]
    data = [sub.loc[sub["max_phase"] == k, "p_active_nlrp3"].values for k in order_keys]
    bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55)
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.55)
    ax.axhline(0.5, color=THRESHOLD, linestyle="--", linewidth=0.8)
    ax.set_ylabel("NLRP3 ML P(active)")
    ax.set_xlabel("Clinical phase")
    panel_label(ax, "c")

    # d: funnel text
    ax = axes[1, 1]
    add_target_banner(ax, "NLRP3")
    ax.axis("off")
    lines = [
        f"Clinical library scored: {summary['n_scored']:,}",
        f"P(active) ≥ 0.5 pool: {summary['n_pred_active_ge_threshold']:,}",
        "Dual docking (URAT1 + NLRP3): pending",
        "",
        "URAT1 evidence: separate 8973 retrospective",
    ]
    ax.text(0.05, 0.85, "\n".join(lines), va="top", fontsize=8, family="sans-serif")
    panel_label(ax, "d")

    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    paths = save_figure(fig, "fig02_nlrp3_screening_composite", "main")
    return {"id": "fig02_composite", "target": "NLRP3", "description": "Main Fig 2 composite (NLRP3 screening)", **paths}


def plot_fig03_composite(dock: pd.DataFrame, summary: dict, bench: pd.DataFrame) -> dict:
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(130))

    # a violin
    ax = axes[0, 0]
    add_target_banner(ax, "URAT1")
    sub = dock[(dock["subset"].isin(["A", "D"])) & dock["docked"]]
    a = sub.loc[sub["subset"] == "A", "glide_score_xp"]
    d = sub.loc[sub["subset"] == "D", "glide_score_xp"]
    parts = ax.violinplot([d, a], positions=[0, 1], showmedians=True, widths=0.7)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(URAT1_COLOR if i == 1 else "#88BCE4")
        body.set_alpha(0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Decoy D", "Active A"])
    ax.set_ylabel("Glide XP (9DKB)")
    panel_label(ax, "a")

    # b ROC
    ax = axes[0, 1]
    add_target_banner(ax, "URAT1")
    y = (sub["subset"] == "A").astype(int)
    fpr, tpr, _ = roc_curve(y, sub["s_u_percentile"])
    ax.plot(fpr, tpr, color=URAT1_COLOR, linewidth=1.0)
    ax.plot([0, 1], [0, 1], color=GRID, linestyle="--", linewidth=0.6)
    ax.text(0.5, 0.12, f"AUC = {auc(fpr, tpr):.3f}", fontsize=8)
    ax.set_xlabel("FPR (decoy D)")
    ax.set_ylabel("TPR (active A)")
    panel_label(ax, "b")

    # c EF
    ax = axes[1, 0]
    add_target_banner(ax, "URAT1")
    enr = summary["enrichment"]
    ax.bar(["EF 5%", "EF 10%"], [enr["ef_5pct_a_vs_d"], enr["ef_10pct_a_vs_d"]], color=URAT1_COLOR, width=0.5, edgecolor=NEUTRAL, linewidth=0.3)
    ax.axhline(1, color=GRID, linestyle="--", linewidth=0.6)
    ax.set_ylabel("Enrichment factor")
    panel_label(ax, "c")

    # d benchmark
    ax = axes[1, 1]
    add_target_banner(ax, "URAT1")
    x = np.arange(len(bench))
    w = 0.35
    ax.bar(x - w / 2, bench["s_u_percentile"], width=w, color=URAT1_COLOR, label="Docking")
    ax.bar(x + w / 2, bench["ml_percentile_vs_8973"], width=w, color="#88BCE4", label="ML")
    ax.set_xticks(x)
    ax.set_xticklabels(bench["compound"], rotation=20, ha="right")
    ax.axhline(90, color=THRESHOLD, linestyle="--", linewidth=0.6)
    ax.set_ylabel("Percentile")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    panel_label(ax, "d")

    fig.subplots_adjust(hspace=0.45, wspace=0.35)
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
        "style": {"font": "Arial", "size_pt": 8, "urat1_color": URAT1_COLOR, "nlrp3_color": NLRP3_COLOR},
        "figures": entries,
    }
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest_out, indent=2))
    print(json.dumps(manifest_out, indent=2))
    print(f"\nWrote {len(entries)} figure sets to {FIG_OUT}")


if __name__ == "__main__":
    main()
