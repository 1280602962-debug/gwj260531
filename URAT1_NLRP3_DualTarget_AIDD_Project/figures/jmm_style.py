"""Unified JMM figure style: Arial 8 pt, URAT1/NLRP3 color coding."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIG_OUT = PROJECT_ROOT / "figures" / "generated"

# Okabe–Ito, colorblind-friendly
URAT1_COLOR = "#0072B2"  # blue
NLRP3_COLOR = "#D55E00"  # vermillion
NEUTRAL = "#333333"
GRID = "#CCCCCC"
THRESHOLD = "#009E73"
CONTROL_HIGHLIGHT = "#000000"

TARGET_STYLES = {
    "URAT1": {"color": URAT1_COLOR, "label": "URAT1 (SLC22A12)"},
    "NLRP3": {"color": NLRP3_COLOR, "label": "NLRP3"},
}

FONT_SIZE_PT = 8
SINGLE_COL_MM = 84
DOUBLE_COL_MM = 174


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.size": FONT_SIZE_PT,
            "axes.labelsize": FONT_SIZE_PT,
            "axes.titlesize": FONT_SIZE_PT,
            "xtick.labelsize": FONT_SIZE_PT,
            "ytick.labelsize": FONT_SIZE_PT,
            "legend.fontsize": FONT_SIZE_PT,
            "axes.linewidth": 0.5,
            "lines.linewidth": 0.8,
            "patch.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.grid": False,
        }
    )


def figsize_single(height_mm: float = 70) -> tuple[float, float]:
    return (mm_to_in(SINGLE_COL_MM), mm_to_in(height_mm))


def figsize_double(height_mm: float = 120) -> tuple[float, float]:
    return (mm_to_in(DOUBLE_COL_MM), mm_to_in(height_mm))


def add_target_banner(ax, target: str, x: float = 0.02, y: float = 0.98) -> None:
    style = TARGET_STYLES[target]
    ax.text(
        x,
        y,
        style["label"],
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_SIZE_PT,
        fontweight="bold",
        color=style["color"],
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=style["color"], linewidth=0.6),
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=FONT_SIZE_PT,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save_figure(fig: plt.Figure, stem: str, subdir: str = "") -> dict[str, str]:
    out_dir = FIG_OUT / subdir if subdir else FIG_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext in ("pdf", "png"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None)
        paths[ext] = str(p)
    plt.close(fig)
    return paths
