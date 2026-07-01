"""Publication figure style: clean axes, Arial 8 pt, no grid."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
GRID = "#CCCCCC"
THRESHOLD = "#009E73"
MUTED = "#AAAAAA"
WARN = "#CC6677"

FONT_SIZE_PT = 8
SINGLE_COL_MM = 84
DOUBLE_COL_MM = 174


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
            "axes.labelpad": 5.0,
            "xtick.major.pad": 4.0,
            "ytick.major.pad": 4.0,
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.visible": False,
            "ytick.minor.visible": False,
            "axes.grid": False,
            "grid.alpha": 0.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.pad_inches": 0.06,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def figsize_single(height_mm: float = 76) -> tuple[float, float]:
    return (mm_to_in(SINGLE_COL_MM), mm_to_in(height_mm))


def figsize_double(height_mm: float = 78) -> tuple[float, float]:
    return (mm_to_in(DOUBLE_COL_MM), mm_to_in(height_mm))


def clean_axes(ax) -> None:
    """Minimal axes: left/bottom spines only, no grid, no top/right ticks."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.set_axisbelow(False)
    ax.tick_params(top=False, right=False, which="both")


def set_axis_labels(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xlabel(xlabel, labelpad=5)
    ax.set_ylabel(ylabel, labelpad=5)


def target_header(fig: plt.Figure, text: str, color: str, y: float = 0.97) -> None:
    fig.text(0.5, y, text, ha="center", va="top", fontsize=FONT_SIZE_PT, fontweight="bold", color=color, transform=fig.transFigure)


def tag_panel(ax, label: str) -> None:
    ax.text(-0.13, 1.04, label, transform=ax.transAxes, fontsize=FONT_SIZE_PT, fontweight="bold", va="top", ha="left", clip_on=False)


def ylim_headroom(ax, ratio: float = 0.15) -> None:
    lo, hi = ax.get_ylim()
    if hi > lo:
        ax.set_ylim(lo, hi + (hi - lo) * ratio)


def label_bars_vertical(ax, bars, fmt: str = "{:.2f}") -> None:
    for bar in bars:
        h = bar.get_height()
        if h <= 0:
            continue
        ax.text(bar.get_x() + bar.get_width() / 2, h, fmt.format(h), ha="center", va="bottom", fontsize=FONT_SIZE_PT, clip_on=False)


def label_bars_horizontal_outside(ax, bars, values, pad: float = 0.02, xmin: float = 0.04) -> None:
    """Place numeric labels to the right of horizontal bars; near-zero bars label inside."""
    for bar, val in zip(bars, values):
        y = bar.get_y() + bar.get_height() / 2
        if val < 0.08:
            ax.text(xmin, y, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
        else:
            ax.text(bar.get_width() + pad, y, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)


def label_hbars_counts(ax, bars, labels: list[str], pad_ratio: float = 0.015) -> None:
    xmax = ax.get_xlim()[1]
    pad = xmax * pad_ratio
    for bar, txt in zip(bars, labels):
        ax.text(bar.get_width() + pad, bar.get_y() + bar.get_height() / 2, txt, va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)


def save_figure(fig: plt.Figure, stem: str, subdir: str = "") -> dict[str, str]:
    out_dir = FIG_OUT / subdir if subdir else FIG_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for ext in ("pdf", "png"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight", pad_inches=0.06)
        paths[ext] = str(p)
    plt.close(fig)
    return paths

# Standard margins (fraction of figure)
MARGIN_SINGLE = dict(left=0.17, right=0.96, top=0.76, bottom=0.18)
MARGIN_WIDE_X = dict(left=0.15, right=0.96, top=0.76, bottom=0.26)
MARGIN_COMPOSITE = dict(left=0.11, right=0.97, top=0.88, bottom=0.13, hspace=0.62, wspace=0.48)
