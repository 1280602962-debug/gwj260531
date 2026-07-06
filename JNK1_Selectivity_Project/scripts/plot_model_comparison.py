#!/usr/bin/env python3
"""
Plot XGBoost vs Chemprop model comparison from archived JSON results.

Reads:
  - results/model_comparison/comparison.json  (primary: XGBoost + Chemprop single-target)
  - results/training/training_report.json   (fallback: Chemprop MTL exploratory holdout)

Outputs:
  - results/model_comparison/model_comparison_r2.png
  - results/model_comparison/model_comparison_r2.pdf

Usage:
    python3 scripts/plot_model_comparison.py
    python3 scripts/plot_model_comparison.py --comparison-json path/to/comparison.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys_path = ROOT / "scripts"
import sys

sys.path.insert(0, str(sys_path))
from plot_style import FIGSIZE_DOUBLE, apply_journal_style  # noqa: E402

ISOFORMS = ["JNK1", "JNK2", "JNK3"]
COLORS = {
    "xgboost": "#4C72B0",
    "chemprop": "#DD8452",
    "chemprop_mtl": "#AAAAAA",
}


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    return x


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def extract_xgboost(comparison: dict) -> dict:
    xgb = comparison.get("xgboost", {})
    cv, holdout = {}, {}
    for iso in ISOFORMS:
        cv_row = xgb.get("cv", {}).get(iso, {})
        ho_row = xgb.get("holdout", {}).get(iso, {})
        cv[iso] = {
            "mean_r2": _safe_float(cv_row.get("mean_r2")),
            "std_r2": _safe_float(cv_row.get("std_r2")) or 0.0,
        }
        holdout[iso] = _safe_float(ho_row.get("r2"))
    return {"cv": cv, "holdout": holdout}


def extract_chemprop_single(comparison: dict) -> dict:
    cp = comparison.get("chemprop", {})
    cv, holdout = {}, {}
    for iso in ISOFORMS:
        cv_row = cp.get("cv", {}).get(iso, {})
        ho_row = cp.get("holdout", {}).get(iso, {})
        cv[iso] = {
            "mean_r2": _safe_float(cv_row.get("mean_r2")),
            "std_r2": _safe_float(cv_row.get("std_r2")) or 0.0,
        }
        if isinstance(ho_row, dict):
            holdout[iso] = _safe_float(ho_row.get("r2"))
        else:
            holdout[iso] = _safe_float(ho_row)
    return {"cv": cv, "holdout": holdout}


def extract_chemprop_mtl(training: dict) -> dict[str, float | None]:
    metrics = training.get("stage_b_mtl", {}).get("metrics", {})
    key_map = {"JNK1": "pAct_JNK1", "JNK2": "pAct_JNK2", "JNK3": "pAct_JNK3"}
    out: dict[str, float | None] = {}
    for iso in ISOFORMS:
        row = metrics.get(key_map[iso], {})
        out[iso] = _safe_float(row.get("r2"))
    return out


def has_any(values: dict[str, float | None]) -> bool:
    return any(v is not None for v in values.values())


def build_figure(
    xgb: dict,
    chemprop: dict,
    chemprop_mtl: dict[str, float | None],
) -> plt.Figure:
    apply_journal_style()
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE_DOUBLE, sharey=True)

    cp_holdout_available = has_any(chemprop["holdout"])
    show_mtl = not cp_holdout_available and has_any(chemprop_mtl)

    n_groups = len(ISOFORMS)
    x = np.arange(n_groups)
    width = 0.24

    def _bar_group(ax, metric: str, ylabel: str, title: str):
        bars_xgb, errs_xgb = [], []
        bars_cp, errs_cp = [], []
        for iso in ISOFORMS:
            if metric == "cv":
                bars_xgb.append(xgb["cv"][iso]["mean_r2"])
                errs_xgb.append(xgb["cv"][iso]["std_r2"])
                bars_cp.append(chemprop["cv"][iso]["mean_r2"])
                errs_cp.append(chemprop["cv"][iso]["std_r2"])
            else:
                bars_xgb.append(xgb["holdout"][iso])
                bars_cp.append(chemprop["holdout"][iso])

        ax.bar(
            x - width,
            [v if v is not None else 0 for v in bars_xgb],
            width,
            yerr=errs_xgb if metric == "cv" else None,
            capsize=2,
            color=COLORS["xgboost"],
            label="XGBoost",
            edgecolor="white",
            linewidth=0.5,
        )

        if metric == "cv" and has_any({iso: chemprop["cv"][iso]["mean_r2"] for iso in ISOFORMS}):
            ax.bar(
                x,
                [v if v is not None else 0 for v in bars_cp],
                width,
                yerr=errs_cp,
                capsize=2,
                color=COLORS["chemprop"],
                label="Chemprop 2.0",
                edgecolor="white",
                linewidth=0.5,
            )
        elif metric == "holdout" and cp_holdout_available:
            ax.bar(
                x,
                [v if v is not None else 0 for v in bars_cp],
                width,
                color=COLORS["chemprop"],
                label="Chemprop 2.0",
                edgecolor="white",
                linewidth=0.5,
            )
        elif metric == "holdout" and show_mtl:
            mtl_vals = [chemprop_mtl[iso] if chemprop_mtl[iso] is not None else 0 for iso in ISOFORMS]
            ax.bar(
                x,
                mtl_vals,
                width,
                color=COLORS["chemprop_mtl"],
                label="Chemprop MTL (exploratory)",
                edgecolor="#666666",
                linewidth=0.8,
                hatch="///",
            )
            for i, iso in enumerate(ISOFORMS):
                val = chemprop_mtl[iso]
                if val is not None:
                    ax.text(i, val + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=6, color="#444444")

        ax.set_xticks(x)
        ax.set_xticklabels(ISOFORMS)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(0, 1.0)
        ax.axhline(0.7, color="#999999", linestyle="--", linewidth=0.7, alpha=0.8)

        for i, iso in enumerate(ISOFORMS):
            val = xgb["cv"][iso]["mean_r2"] if metric == "cv" else xgb["holdout"][iso]
            if val is not None:
                ax.text(i - width, val + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=6)

    _bar_group(axes[0], "cv", "R²", "5-fold scaffold CV (mean ± SD)")
    holdout_title = "Holdout test (scaffold split)"
    if show_mtl:
        holdout_title += "\nChemprop: MTL exploratory (single-target N/R)"
    elif not cp_holdout_available:
        holdout_title += "\nChemprop 2.0: not available in archive"
    _bar_group(axes[1], "holdout", "R²", holdout_title)

    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle("JNK1/2/3 activity QSAR: XGBoost vs Chemprop", y=1.08, fontsize=9)
    fig.tight_layout()
    return fig


def generate_comparison_figure(
    comparison_json: Path | None = None,
    training_json: Path | None = None,
    output: Path | None = None,
) -> tuple[Path, Path]:
    comparison_json = comparison_json or ROOT / "results" / "model_comparison" / "comparison.json"
    training_json = training_json or ROOT / "results" / "training" / "training_report.json"
    output = output or ROOT / "results" / "model_comparison" / "model_comparison_r2"

    comparison = load_json(comparison_json)
    training = load_json(training_json) if training_json.exists() else {}

    xgb = extract_xgboost(comparison)
    chemprop = extract_chemprop_single(comparison)
    chemprop_mtl = extract_chemprop_mtl(training)

    fig = build_figure(xgb, chemprop, chemprop_mtl)
    png_path = output.with_suffix(".png")
    pdf_path = output.with_suffix(".pdf")
    apply_journal_style()
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot XGBoost vs Chemprop comparison")
    parser.add_argument(
        "--comparison-json",
        type=Path,
        default=ROOT / "results" / "model_comparison" / "comparison.json",
    )
    parser.add_argument(
        "--training-json",
        type=Path,
        default=ROOT / "results" / "training" / "training_report.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "model_comparison" / "model_comparison_r2",
    )
    args = parser.parse_args()

    png_path, pdf_path = generate_comparison_figure(
        args.comparison_json, args.training_json, args.output
    )
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
