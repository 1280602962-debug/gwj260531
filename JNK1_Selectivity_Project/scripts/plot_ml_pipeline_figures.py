#!/usr/bin/env python3
"""
Generate ML pipeline figures from archived results (no model retraining required).

Uses journal style from plot_style.py (Arial, 300 dpi, English axis labels).

Outputs under results/model_comparison/ and results/screening_v2/:
  - model_comparison_r2.png/pdf          (delegates to plot_model_comparison)
  - screening_funnel.png/pdf
  - score_distribution.png/pdf
  - decoy_validation_metrics.png/pdf
  - ml_benchmark_isoform_prediction.png/pdf
  - selectivity_label_scarcity.png/pdf

Usage:
    python3 scripts/plot_ml_pipeline_figures.py
    python3 scripts/plot_ml_pipeline_figures.py --output-dir results/figures/ml
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import FIGSIZE_DOUBLE, FIGSIZE_SINGLE, apply_journal_style, save_figure  # noqa: E402

ISOFORMS = ["JNK1", "JNK2", "JNK3"]
COLORS_ISO = {"JNK1": "#4C72B0", "JNK2": "#DD8452", "JNK3": "#55A868"}


def _save(fig: plt.Figure, stem: Path) -> tuple[Path, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    apply_journal_style()
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


def plot_screening_funnel(report_path: Path, out_stem: Path) -> tuple[Path, Path]:
    with open(report_path) as f:
        stats = json.load(f)["funnel_stats"]
    stages = ["Input", "Drug-like", "F1 p_family", "SA/QED"]
    keys = ["input", "druglike", "f1_pass", "sa_qed_pass"]
    values = [stats[k] for k in keys]

    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    bars = ax.bar(stages, values, color="#55A868", edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Compound count")
    ax.set_title("ML virtual screening funnel (demo library)")
    ax.tick_params(axis="x", rotation=20)
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val:,}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    fig.tight_layout()
    return _save(fig, out_stem)


def plot_score_distribution(hits_path: Path, out_stem: Path) -> tuple[Path, Path] | None:
    df = pd.read_csv(hits_path)
    if df.empty or "final_score" not in df.columns:
        return None

    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    ax.hist(df["final_score"], bins=30, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("final_score (weighted composite)")
    ax.set_ylabel("Hit count")
    ax.set_title("Distribution of final screening scores")
    ax.axvline(df["final_score"].quantile(1 - 500 / len(df)), color="#C44E52", linestyle="--", linewidth=0.8, label="Top-500 cutoff (demo)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return _save(fig, out_stem)


def plot_decoy_metrics(metrics_path: Path, out_stem: Path) -> tuple[Path, Path]:
    with open(metrics_path) as f:
        m = json.load(f)
    cf = m.get("confusion_full", m)

    labels = ["Recall", "Decoy FPR", "Specificity", "ROC-AUC", "EF1%"]
    values = [
        cf.get("sensitivity_recall", cf.get("recall", 0)) * 100,
        cf.get("false_positive_rate", cf.get("decoy_fpr", 0)) * 100,
        cf.get("specificity", 0) * 100,
        m.get("roc_auc", 0),
        m.get("ef1pct", m.get("ef1_percent", 0)),
    ]
    colors = ["#55A868", "#C44E52", "#8172B3", "#4C72B0", "#CCB974"]

    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Value (% or unitless)")
    ax.set_title("External decoy validation (F1 @ p_family ≥ 6.0)")
    ax.tick_params(axis="x", rotation=15)
    for bar, val in zip(bars, values):
        fmt = f"{val:.1f}" if val > 1.5 else f"{val:.2f}"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), fmt, ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    return _save(fig, out_stem)


def plot_benchmark_isoform_prediction(bench_path: Path, out_stem: Path) -> tuple[Path, Path]:
    df = pd.read_csv(bench_path)
    pred_cols = ["pred_pAct_JNK1", "pred_pAct_JNK2", "pred_pAct_JNK3"]
    if not all(c in df.columns for c in pred_cols):
        raise ValueError(f"Missing prediction columns in {bench_path}")

    df = df.copy()
    df["ml_top_isoform"] = df[pred_cols].idxmax(axis=1).str.replace("pred_pAct_", "")

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    x = np.arange(len(df))
    width = 0.25
    for i, iso in enumerate(ISOFORMS):
        ax.bar(x + (i - 1) * width, df[f"pred_pAct_{iso}"], width, label=iso, color=COLORS_ISO[iso], edgecolor="white", linewidth=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels(df["name"], rotation=30, ha="right")
    ax.set_ylabel("Predicted pActivity")
    ax.set_title("ML isoform activity predictions on literature benchmarks")
    ax.legend(ncol=3, fontsize=7, frameon=False)
    fig.tight_layout()
    return _save(fig, out_stem)


def plot_selectivity_labels(counts_path: Path, out_stem: Path) -> tuple[Path, Path]:
    df = pd.read_csv(counts_path)
    label_col = "sel_class" if "sel_class" in df.columns else df.columns[0]
    count_col = "count" if "count" in df.columns else df.columns[1]

    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    colors = ["#C44E52" if "JNK1" in str(l) else "#AAAAAA" for l in df[label_col]]
    bars = ax.bar(df[label_col], df[count_col], color=colors, edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Compound count")
    ax.set_title("Selectivity label scarcity (paired ChEMBL set)")
    ax.tick_params(axis="x", rotation=20)
    for bar, val in zip(bars, df[count_col]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(int(val)), ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    return _save(fig, out_stem)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot ML pipeline figures from archived data")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional extra copy directory")
    args = parser.parse_args()

    written: list[Path] = []

    from plot_model_comparison import generate_comparison_figure

    png, pdf = generate_comparison_figure()
    written.extend([png, pdf])
    print(f"Wrote {png}")

    screening_dir = ROOT / "results" / "screening_v2"
    png, pdf = plot_screening_funnel(
        screening_dir / "screening_report.json",
        screening_dir / "screening_funnel",
    )
    written.extend([png, pdf])
    print(f"Wrote {png}")

    hits = screening_dir / "all_hits.csv"
    if hits.exists():
        out = plot_score_distribution(hits, screening_dir / "score_distribution")
        if out:
            written.extend(out)
            print(f"Wrote {out[0]}")

    decoy_json = ROOT / "results" / "ml_external_validation" / "ml_external_validation_metrics_9bd8.json"
    if decoy_json.exists():
        png, pdf = plot_decoy_metrics(decoy_json, screening_dir / "decoy_validation_metrics")
        written.extend([png, pdf])
        print(f"Wrote {png}")

    bench_csv = screening_dir / "benchmark_validation.csv"
    if bench_csv.exists():
        png, pdf = plot_benchmark_isoform_prediction(
            bench_csv, screening_dir / "ml_benchmark_isoform_prediction"
        )
        written.extend([png, pdf])
        print(f"Wrote {png}")

    sel_csv = ROOT / "results" / "similarity" / "sel_class_counts.csv"
    if sel_csv.exists():
        png, pdf = plot_selectivity_labels(sel_csv, screening_dir / "selectivity_label_scarcity")
        written.extend([png, pdf])
        print(f"Wrote {png}")

    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        import shutil

        for p in written:
            dest = args.output_dir / p.name
            shutil.copy2(p, dest)
            print(f"Copied to {dest}")

    fig_dir = ROOT / "docs" / "popular_science" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    import shutil

    mapping = {
        "model_comparison_r2.png": ROOT / "results" / "model_comparison" / "model_comparison_r2.png",
        "screening_funnel.png": screening_dir / "screening_funnel.png",
        "score_distribution.png": screening_dir / "score_distribution.png",
        "decoy_validation_metrics.png": screening_dir / "decoy_validation_metrics.png",
        "ml_benchmark_isoform_prediction.png": screening_dir / "ml_benchmark_isoform_prediction.png",
        "selectivity_label_scarcity.png": screening_dir / "selectivity_label_scarcity.png",
    }
    for name, src in mapping.items():
        if src.exists():
            shutil.copy2(src, fig_dir / name)
            print(f"Synced {fig_dir / name}")


if __name__ == "__main__":
    main()
