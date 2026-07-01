"""Publication figure style: SciencePlots + JMM (Arial 8 pt, URAT1/NLRP3)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    import scienceplots  # noqa: F401

    _STYLE = ["science", "no-latex"]
except ImportError:
    _STYLE = []

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIG_OUT = PROJECT_ROOT / "figures" / "generated"

URAT1_COLOR = "#0072B2"
NLRP3_COLOR = "#D55E00"
URAT1_LIGHT = "#56B4E9"
NEUTRAL = "#333333"
GRID = "#BBBBBB"
THRESHOLD = "#009E73"
MUTED = "#999999"

TARGET_META = {
    "URAT1": {"color": URAT1_COLOR, "short": "URAT1", "full": "URAT1 (SLC22A12, 9DKB XP)"},
    "NLRP3": {"color": NLRP3_COLOR, "short": "NLRP3", "full": "NLRP3 (assay-conditioned ML)"},
}

FONT_SIZE_PT = 8
SINGLE_COL_MM = 84
DOUBLE_COL_MM = 174

# Margins tuned to avoid clipping panel tags / axis labels (SciencePlots-inspired)
MARGIN_SINGLE = dict(left=0.16, right=0.97, top=0.82, bottom=0.20)
MARGIN_WIDE_X = dict(left=0.14, right=0.97, top=0.82, bottom=0.28)
MARGIN_COMPOSITE = dict(left=0.10, right=0.98, top=0.86, bottom=0.12, hspace=0.55, wspace=0.42)


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def apply_style() -> None:
    if _STYLE:
        plt.style.use(_STYLE)
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
            "axes.labelpad": 4.0,
            "xtick.major.pad": 3.0,
            "ytick.major.pad": 3.0,
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "axes.grid": False,
            "legend.frameon": False,
        }
    )


def figsize_single(height_mm: float = 72) -> tuple[float, float]:
    return (mm_to_in(SINGLE_COL_MM), mm_to_in(height_mm))


def figsize_double(height_mm: float = 72) -> tuple[float, float]:
    return (mm_to_in(DOUBLE_COL_MM), mm_to_in(height_mm))


def style_axes(ax, *, grid_y: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_y:
        ax.yaxis.grid(True, color=GRID, linewidth=0.4, alpha=0.8)
        ax.set_axisbelow(True)


def set_axis_labels(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xlabel(xlabel, labelpad=4)
    ax.set_ylabel(ylabel, labelpad=4)


def set_target_title(fig: plt.Figure, target: str, y: float = 0.98) -> None:
    meta = TARGET_META[target]
    fig.text(
        0.5,
        y,
        f"Target: {meta['full']}",
        ha="center",
        va="top",
        fontsize=FONT_SIZE_PT,
        fontweight="bold",
        color=meta["color"],
        transform=fig.transFigure,
    )


def tag_panel(ax, label: str) -> None:
    ax.text(
        -0.16,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=FONT_SIZE_PT,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def apply_margins(fig: plt.Figure, margins: dict) -> None:
    fig.subplots_adjust(**margins)


def save_figure(fig: plt.Figure, stem: str, subdir: str = "") -> dict[str, str]:
    out_dir = FIG_OUT / subdir if subdir else FIG_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for ext in ("pdf", "png"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None, facecolor="white")
        paths[ext] = str(p)
    plt.close(fig)
    return paths
