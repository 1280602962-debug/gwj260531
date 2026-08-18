#!/usr/bin/env python3
"""Clean publication figures: Arial 8 pt, no grid, labels outside plot elements."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
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
    WARN,
    FONT_SIZE_PT,
    annotate_threshold_hist,
    apply_panel_tags,
    apply_style,
    clean_axes,
    clean_hbar_axes,
    figsize_double,
    figsize_single,
    label_bars_horizontal_outside,
    label_bars_vertical,
    label_hbars_counts,
    legend_below,
    legend_lower_right,
    save_figure,
    set_axis_labels,
    tag_panel,
    target_header,
    ylim_headroom,
)

DATA = ROOT / "data"
RESULTS = ROOT / "results"
MANIFEST_PATH = ROOT / "figures" / "generated" / "figure_manifest.json"

NLRP3_GOUT_COMEDS = ["benzbromarone", "dotinurad", "lesinurad", "allopurinol", "febuxostat"]
# Indirect / off-target — known confounders (discuss in text)
NLRP3_OFFTARGET = ["colchicine", "verinurad"]


def _name_col(df: pd.DataFrame) -> str:
    return "name" if "name" in df.columns else "pref_name"


def _nlrp3_control_rows(ml: pd.DataFrame) -> list[dict]:
    name_col = _name_col(ml)
    rows = []
    for drug in NLRP3_GOUT_COMEDS + NLRP3_OFFTARGET:
        hit = ml[ml[name_col].astype(str).str.upper() == drug.upper()]
        if not len(hit):
            continue
        group = "off_target" if drug in NLRP3_OFFTARGET else "gout_comed"
        rows.append({"drug": drug, "p_active": float(hit["p_active_nlrp3"].iloc[0]), "group": group})
    return sorted(rows, key=lambda r: r["p_active"])


def _funnel_counts(summary: dict) -> tuple[list[int], list[str]]:
    n_scored = summary["n_scored"]
    n_active = summary["n_pred_active_ge_threshold"]
    n_dual = summary.get("n_dual_docked")
    if n_dual is not None:
        return [n_scored, n_active, n_dual], [f"{n_scored:,}", f"{n_active:,}", f"{n_dual:,}"]
    pending_stub = 80
    return [n_scored, n_active, pending_stub], [f"{n_scored:,}", f"{n_active:,}", "pending"]


def _control_color(group: str, p: float) -> str:
    if group == "off_target":
        return WARN
    return MUTED


def plot_nlrp3_distribution(ml: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(76))
    target_header(fig, "Target: NLRP3 (assay-conditioned ML)", NLRP3_COLOR)

    n_ge = int((ml["p_active_nlrp3"] >= 0.5).sum())
    ax.hist(ml["p_active_nlrp3"], bins=45, color=NLRP3_COLOR, alpha=0.85, edgecolor="white", linewidth=0.2)
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    clean_axes(ax)
    set_axis_labels(ax, "Predicted probability P(active)", "Number of compounds")
    ax.set_xlim(0, 1.0)
    ylim_headroom(ax, 0.08)
    annotate_threshold_hist(ax, n_ge)
    fig.subplots_adjust(**MARGIN_SINGLE)
    paths = save_figure(fig, "nlrp3_fig02a_pactive_histogram", "nlrp3")
    return {"id": "nlrp3_fig02a", "target": "NLRP3", "description": "P(active) distribution (n=8319)", **paths}


def plot_nlrp3_controls(ml: pd.DataFrame) -> dict:
    rows = _nlrp3_control_rows(ml)
    drugs = [r["drug"] for r in rows]
    vals = [r["p_active"] for r in rows]
    colors = [_control_color(r["group"], r["p_active"]) for r in rows]

    fig, ax = plt.subplots(figsize=figsize_single(88))
    target_header(fig, "Target: NLRP3 (assay-conditioned ML)", NLRP3_COLOR)

    bars = ax.barh(drugs, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.58)
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    ax.set_xlim(0, 1.05)
    label_bars_horizontal_outside(ax, bars, vals, label_col=1.08)
    clean_axes(ax)
    set_axis_labels(ax, "Predicted probability P(active)", "Reference compound")
    legend_handles = [
        Patch(facecolor=MUTED, edgecolor=NEUTRAL, label="Gout co-medication"),
        Patch(facecolor=WARN, edgecolor=NEUTRAL, label="Indirect / off-target"),
    ]
    legend_below(ax, handles=legend_handles, ncol=2, y=-0.30)
    fig.subplots_adjust(**{**MARGIN_SINGLE, "left": 0.26, "bottom": 0.28})
    paths = save_figure(fig, "nlrp3_fig02b_control_compounds", "nlrp3")
    return {"id": "nlrp3_fig02b", "target": "NLRP3", "description": "Reference compounds by role", **paths}


def plot_nlrp3_phase(ml: pd.DataFrame) -> dict:
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = ml[ml["max_phase"].isin(phase_keys)]
    data = [sub.loc[sub["max_phase"] == k, "p_active_nlrp3"].values for k in phase_keys]

    fig, ax = plt.subplots(figsize=figsize_single(78))
    target_header(fig, "Target: NLRP3 (assay-conditioned ML)", NLRP3_COLOR)

    bp = ax.boxplot(data, tick_labels=phase_labels, patch_artist=True, widths=0.48, showfliers=False, medianprops={"color": NEUTRAL, "linewidth": 0.9})
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.5)
        patch.set_edgecolor(NEUTRAL)
    ax.axhline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    clean_axes(ax)
    set_axis_labels(ax, "Clinical development stage", "Predicted probability P(active)")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    ax.text(0.55, 0.97, "Threshold 0.5", transform=ax.transAxes, ha="center", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)
    fig.subplots_adjust(**MARGIN_WIDE_X)
    paths = save_figure(fig, "nlrp3_fig02c_phase_boxplot", "nlrp3")
    return {"id": "nlrp3_fig02c", "target": "NLRP3", "description": "P(active) by clinical phase", **paths}


def plot_screening_funnel(summary: dict) -> dict:
    stages = ["Library scored", "P(active) ≥ 0.5", "Dual dock 9DKB+7ALV"]
    counts, count_labels = _funnel_counts(summary)
    colors = [NLRP3_COLOR, NLRP3_COLOR, MUTED]

    fig, ax = plt.subplots(figsize=figsize_single(62))
    target_header(fig, "Target: NLRP3 prescreen funnel", NLRP3_COLOR)

    y = np.arange(len(stages))
    bars = ax.barh(y, counts, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.52)
    ax.set_xlim(0, counts[0] * 1.18)
    label_hbars_counts(ax, bars, count_labels)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.invert_yaxis()
    clean_hbar_axes(ax, y_grid=True)
    set_axis_labels(ax, "Number of compounds", "Screening stage")
    fig.subplots_adjust(**{**MARGIN_SINGLE, "left": 0.32})
    paths = save_figure(fig, "nlrp3_fig02d_screening_funnel", "nlrp3")
    return {"id": "nlrp3_fig02d", "target": "NLRP3", "description": "NLRP3 prescreen funnel", **paths}


def plot_nlrp3_roc_pr(oof: pd.DataFrame) -> dict:
    fpr, tpr, _ = roc_curve(oof["y_true"], oof["y_prob"])
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(oof["y_true"], oof["y_prob"])
    pr_auc = auc(rec, prec)

    fig, axes = plt.subplots(1, 2, figsize=figsize_double(76))
    target_header(fig, "Target: NLRP3 (scaffold-CV out-of-fold)", NLRP3_COLOR)

    ax = axes[0]
    ax.plot(fpr, tpr, color=NLRP3_COLOR, linewidth=1.1, label=f"AUC = {roc_auc:.3f}")
    ax.plot([0, 1], [0, 1], color=GRID, linestyle=(0, (4, 3)), linewidth=0.7)
    clean_axes(ax)
    set_axis_labels(ax, "False positive rate", "True positive rate")
    ax.legend(loc="lower right")
    tag_panel(ax, "a")

    ax = axes[1]
    ax.plot(rec, prec, color=NLRP3_COLOR, linewidth=1.1, label=f"AUPRC = {pr_auc:.3f}")
    clean_axes(ax)
    set_axis_labels(ax, "Recall", "Precision")
    ax.legend(loc="upper right")
    tag_panel(ax, "b")

    fig.subplots_adjust(**{**MARGIN_COMPOSITE, "top": 0.82, "bottom": 0.18})
    paths = save_figure(fig, "si_nlrp3_oof_roc_pr", "si")
    return {"id": "si_nlrp3_roc_pr", "target": "NLRP3", "description": "NLRP3 scaffold-CV OOF ROC/PR", **paths}


def plot_urat1_parity(oof: pd.DataFrame) -> dict:
    fig, ax = plt.subplots(figsize=figsize_single(76))
    target_header(fig, "Target: URAT1 (scaffold-CV out-of-fold)", URAT1_COLOR)

    ax.scatter(oof["y_true"], oof["y_pred"], s=9, alpha=0.35, color=URAT1_COLOR, edgecolors="none", rasterized=True)
    lo = min(oof["y_true"].min(), oof["y_pred"].min()) - 0.2
    hi = max(oof["y_true"].max(), oof["y_pred"].max()) + 0.2
    ax.plot([lo, hi], [lo, hi], color=NEUTRAL, linestyle=(0, (4, 3)), linewidth=0.7, label="y = x")
    clean_axes(ax)
    set_axis_labels(ax, "Observed pActivity (ChEMBL)", "Predicted pActivity (OOF)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right")
    fig.subplots_adjust(**MARGIN_SINGLE)
    paths = save_figure(fig, "si_urat1_oof_parity", "si")
    return {"id": "si_urat1_parity", "target": "URAT1", "description": "URAT1 OOF parity plot", **paths}


def plot_data_asymmetry(summary: dict) -> dict:
    fig, axes = plt.subplots(1, 2, figsize=figsize_double(72))
    target_header(fig, "Dual-target asymmetric data design", NEUTRAL)

    ax = axes[0]
    labels = ["URAT1", "NLRP3", "Shared SMILES"]
    vals = [summary["urat1"]["n_compounds"], summary["nlrp3"]["n_compounds"], summary["overlap_smiles"]]
    bars = ax.bar(labels, vals, color=[URAT1_COLOR, NLRP3_COLOR, MUTED], edgecolor=NEUTRAL, linewidth=0.3, width=0.5)
    clean_axes(ax)
    set_axis_labels(ax, "ChEMBL training set", "Compound count")
    ylim_headroom(ax, 0.18)
    label_bars_vertical(ax, bars, fmt="{:.0f}")
    tag_panel(ax, "a")

    ax = axes[1]
    routes = ["URAT1 dock pool\n~1588 @ 9DKB P2", "NLRP3 shrink\n8319 clinical"]
    ns = [1588, 8319]
    bars = ax.bar(routes, ns, color=[URAT1_COLOR, NLRP3_COLOR], edgecolor=NEUTRAL, linewidth=0.3, width=0.46)
    clean_axes(ax)
    set_axis_labels(ax, "Production evidence route", "Library size")
    ylim_headroom(ax, 0.15)
    label_bars_vertical(ax, bars, fmt="{:.0f}")
    tag_panel(ax, "b")

    fig.subplots_adjust(**{**MARGIN_COMPOSITE, "top": 0.82, "bottom": 0.20})
    paths = save_figure(fig, "si_data_asymmetry", "si")
    return {"id": "si_data_asymmetry", "target": "both", "description": "Asymmetric dual-target data design", **paths}


def plot_library_phase(manifest: pd.DataFrame) -> dict:
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = manifest[manifest["max_phase"].isin(phase_keys)]
    counts = sub["max_phase"].value_counts().reindex(phase_keys).fillna(0).astype(int)

    fig, ax = plt.subplots(figsize=figsize_single(72))
    target_header(fig, "NLRP3 prescreen input library (ChEMBL clinical)", NLRP3_COLOR)

    bars = ax.bar(phase_labels, counts.values, color=NLRP3_COLOR, alpha=0.85, edgecolor=NEUTRAL, linewidth=0.3, width=0.52)
    clean_axes(ax)
    set_axis_labels(ax, "Clinical development stage", "Number of compounds")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    ylim_headroom(ax, 0.14)
    label_bars_vertical(ax, bars, fmt="{:.0f}")
    fig.subplots_adjust(**MARGIN_WIDE_X)
    paths = save_figure(fig, "si_repurposing_library_phase", "si")
    return {"id": "si_library_phase", "target": "NLRP3", "description": "Clinical library phase composition", **paths}


def plot_fig02_composite(ml: pd.DataFrame, summary: dict) -> dict:
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(158))
    target_header(fig, "Target: NLRP3 — clinical library machine-learning prescreen", NLRP3_COLOR, y=0.975)

    control_legend = [
        Patch(facecolor=MUTED, edgecolor=NEUTRAL, label="Gout co-medication"),
        Patch(facecolor=WARN, edgecolor=NEUTRAL, label="Indirect / off-target"),
    ]

    # (a) histogram
    ax = axes[0, 0]
    n_ge = int((ml["p_active_nlrp3"] >= 0.5).sum())
    ax.hist(ml["p_active_nlrp3"], bins=40, color=NLRP3_COLOR, alpha=0.85, edgecolor="white", linewidth=0.15)
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    clean_axes(ax)
    set_axis_labels(ax, "Predicted probability P(active)", "Number of compounds", xpad=6, ypad=6)
    ax.set_xlim(0, 1)
    ax.margins(x=0.02)
    ylim_headroom(ax, 0.06)
    annotate_threshold_hist(ax, n_ge)

    # (b) controls
    ax = axes[0, 1]
    rows = _nlrp3_control_rows(ml)
    drugs = [r["drug"] for r in rows]
    vals = [r["p_active"] for r in rows]
    colors = [_control_color(r["group"], r["p_active"]) for r in rows]
    bars = ax.barh(drugs, vals, color=colors, edgecolor=NEUTRAL, linewidth=0.3, height=0.55)
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    ax.set_xlim(0, 1.05)
    label_bars_horizontal_outside(ax, bars, vals, label_col=1.08)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "Predicted probability P(active)", "Reference compound", xpad=6, ypad=6)
    legend_lower_right(ax, handles=control_legend, ncol=1, x=0.98, y=0.03)

    # (c) phase boxplot — wider tick spacing
    ax = axes[1, 0]
    phase_keys = [1.0, 2.0, 3.0, 4.0]
    phase_labels = ["Phase I", "Phase II", "Phase III", "Approved"]
    sub = ml[ml["max_phase"].isin(phase_keys)]
    data = [sub.loc[sub["max_phase"] == k, "p_active_nlrp3"].values for k in phase_keys]
    positions = [1.0, 2.6, 4.2, 5.8]
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": NEUTRAL, "linewidth": 0.9},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor(NLRP3_COLOR)
        patch.set_alpha(0.5)
        patch.set_edgecolor(NEUTRAL)
    ax.axhline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.9)
    ax.set_xticks(positions)
    ax.set_xticklabels(phase_labels, rotation=28, ha="right")
    ax.set_xlim(0.2, 6.6)
    clean_axes(ax)
    set_axis_labels(ax, "Clinical development stage", "Predicted probability P(active)", xpad=10, ypad=6)
    ax.text(0.55, 0.97, "Threshold 0.5", transform=ax.transAxes, ha="center", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)

    # (d) funnel — no vertical grid
    ax = axes[1, 1]
    stages = ["Library scored", "P(active) ≥ 0.5", "Dual dock 9DKB+7ALV"]
    counts, funnel_labels = _funnel_counts(summary)
    y = np.arange(len(stages))
    bars = ax.barh(y, counts, color=[NLRP3_COLOR, NLRP3_COLOR, MUTED], edgecolor=NEUTRAL, linewidth=0.3, height=0.48)
    ax.set_xlim(0, counts[0] * 1.18)
    label_hbars_counts(ax, bars, funnel_labels)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.invert_yaxis()
    clean_hbar_axes(ax, y_grid=True)
    set_axis_labels(ax, "Number of compounds", "Screening stage", xpad=6, ypad=6)

    fig.subplots_adjust(**{**MARGIN_COMPOSITE, "top": 0.86, "bottom": 0.15, "left": 0.15, "right": 0.93, "hspace": 0.80, "wspace": 0.58})
    apply_panel_tags(fig, axes, ("a", "b", "c", "d"))
    paths = save_figure(fig, "fig02_nlrp3_screening_composite", "main", tight=False)
    return {"id": "fig02_composite", "target": "NLRP3", "description": "Main Fig 2 composite (NLRP3 screening)", **paths}


def main() -> None:
    apply_style()

    ml = pd.read_csv(DATA / "repurposing" / "screening" / "nlrp3_ml_scores_clinical_all.csv", low_memory=False)
    nlrp3_summary = json.loads((DATA / "repurposing" / "screening" / "nlrp3_screening_summary_clinical_all.json").read_text())
    # Do not inject Glide-era Pareto n_dual_docked. Funnel shows "pending" until P2 scores are committed.
    data_summary = json.loads((DATA / "processed" / "data_summary.json").read_text())
    manifest = pd.read_csv(DATA / "repurposing" / "repurposing_manifest.csv", low_memory=False)

    oof_nlrp3 = RESULTS / "training" / "nlrp3_oof_predictions.csv"
    oof_urat1 = RESULTS / "training" / "urat1_oof_predictions.csv"
    if not oof_nlrp3.exists():
        oof_nlrp3 = DATA / "models" / "nlrp3_oof_predictions.csv"
    if not oof_urat1.exists():
        oof_urat1 = DATA / "models" / "urat1_oof_predictions.csv"

    entries = [
        plot_nlrp3_distribution(ml),
        plot_nlrp3_controls(ml),
        plot_nlrp3_phase(ml),
        plot_screening_funnel(nlrp3_summary),
    ]
    if oof_nlrp3.exists():
        entries.append(plot_nlrp3_roc_pr(pd.read_csv(oof_nlrp3)))
    if oof_urat1.exists():
        entries.append(plot_urat1_parity(pd.read_csv(oof_urat1)))
    entries += [
        plot_data_asymmetry(data_summary),
        plot_library_phase(manifest),
        plot_fig02_composite(ml, nlrp3_summary),
    ]

    # Keep previously generated OOF SI files when prediction CSVs are absent.
    existing_ids = {e["id"] for e in entries}
    for fid, target, desc, stem in (
        ("si_nlrp3_roc_pr", "NLRP3", "NLRP3 scaffold-CV OOF ROC/PR", "si_nlrp3_oof_roc_pr"),
        ("si_urat1_parity", "URAT1", "URAT1 OOF parity plot", "si_urat1_oof_parity"),
    ):
        if fid in existing_ids:
            continue
        png = FIG_OUT / "si" / f"{stem}.png"
        pdf = FIG_OUT / "si" / f"{stem}.pdf"
        if png.exists() and pdf.exists():
            entries.append({"id": fid, "target": target, "description": desc, "pdf": str(pdf), "png": str(png)})

    manifest_out = {
        "style": {
            "font": "Arial 8 pt",
            "grid": "none",
            "notes": "Current figures only: NLRP3 ML + SI. 8973/Glide Pareto plots were removed; regenerate P2 dual-dock figures after local funnel scores exist.",
        },
        "figures": entries,
    }
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest_out, indent=2))
    print(f"Wrote {len(entries)} figure sets to {FIG_OUT}")


if __name__ == "__main__":
    main()
