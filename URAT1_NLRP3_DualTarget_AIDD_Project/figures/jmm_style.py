"""Publication figure style: clean axes, Arial 8 pt, no grid."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
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
PANEL_TAG_OFFSET_X = 0.014  # figure coords left of subplot
PANEL_TAG_OFFSET_Y = 0.010  # figure coords above subplot top
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
            "savefig.pad_inches": 0.12,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def figsize_single(height_mm: float = 76) -> tuple[float, float]:
    return (mm_to_in(SINGLE_COL_MM), mm_to_in(height_mm))


def figsize_double(height_mm: float = 78) -> tuple[float, float]:
    return (mm_to_in(DOUBLE_COL_MM), mm_to_in(height_mm))


def set_axis_labels(ax, xlabel: str, ylabel: str, xpad: float = 5.0, ypad: float = 5.0) -> None:
    ax.set_xlabel(xlabel, labelpad=xpad)
    ax.set_ylabel(ylabel, labelpad=ypad)


def clean_axes(ax) -> None:
    """Minimal axes: left/bottom spines only, no grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)
    ax.tick_params(top=False, right=False, which="both")


def clean_hbar_axes(ax, y_grid: bool = False) -> None:
    """Horizontal bar axes: never show vertical (x) grid lines."""
    clean_axes(ax)
    ax.xaxis.grid(False)
    if y_grid:
        ax.yaxis.grid(True, linestyle=(0, (2, 4)), color=GRID, alpha=0.55, linewidth=0.4)
    else:
        ax.yaxis.grid(False)


def target_header(fig: plt.Figure, text: str, color: str, y: float = 0.97) -> None:
    fig.text(0.5, y, text, ha="center", va="top", fontsize=FONT_SIZE_PT, fontweight="bold", color=color, transform=fig.transFigure)


def tag_panel(ax, label: str) -> None:
    """Legacy per-axes tag; prefer apply_panel_tags for composites."""
    ax.text(-0.11, 1.10, label, transform=ax.transAxes, fontsize=FONT_SIZE_PT, fontweight="bold", va="bottom", ha="right", clip_on=False)


def apply_panel_tags(fig: plt.Figure, axes, labels: tuple[str, ...]) -> None:
    """Place a/b/c/d labels at a uniform offset from each subplot corner (figure coords)."""
    fig.canvas.draw()
    for ax, lab in zip(np.ravel(axes), labels):
        bb = ax.get_position()
        fig.text(
            bb.x0 - PANEL_TAG_OFFSET_X,
            bb.y1 + PANEL_TAG_OFFSET_Y,
            lab,
            ha="right",
            va="bottom",
            fontsize=FONT_SIZE_PT,
            fontweight="bold",
            color=NEUTRAL,
        )


def legend_lower_left(ax, handles=None, labels=None, ncol: int = 1, x: float = 0.02, y: float = 0.04) -> None:
    """Unified in-panel legend anchored to lower-left empty corner."""
    kwargs = dict(
        loc="lower left",
        bbox_to_anchor=(x, y),
        ncol=ncol,
        fontsize=FONT_SIZE_PT,
        frameon=False,
        borderaxespad=0.0,
    )
    if handles is not None:
        ax.legend(handles=handles, labels=labels, **kwargs)
    elif labels is not None:
        ax.legend(labels=labels, **kwargs)
    else:
        ax.legend(**kwargs)


def legend_lower_right(ax, handles=None, labels=None, ncol: int = 1, x: float = 0.98, y: float = 0.04) -> None:
    """Unified in-panel legend anchored to lower-right empty corner."""
    kwargs = dict(
        loc="lower right",
        bbox_to_anchor=(x, y),
        ncol=ncol,
        fontsize=FONT_SIZE_PT,
        frameon=False,
        borderaxespad=0.0,
    )
    if handles is not None:
        ax.legend(handles=handles, labels=labels, **kwargs)
    elif labels is not None:
        ax.legend(labels=labels, **kwargs)
    else:
        ax.legend(**kwargs)


def legend_upper_center(ax, handles=None, labels=None, ncol: int = 2, y: float = 0.97) -> None:
    """Legend in the upper margin of the axes, clear of bar tops."""
    kwargs = dict(
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=FONT_SIZE_PT,
        frameon=False,
        borderaxespad=0.0,
    )
    if handles is not None:
        ax.legend(handles=handles, labels=labels, **kwargs)
    elif labels is not None:
        ax.legend(labels=labels, **kwargs)
    else:
        ax.legend(**kwargs)


def subplot_legend_centered(
    fig: plt.Figure, ax, ncol: int = 2, pad: float = 0.088, handles=None, labels=None
) -> None:
    """Place legend centered below a subplot (under the x-axis title)."""
    if ax.get_legend() is not None:
        leg = ax.get_legend()
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()
        leg.remove()
    fig.canvas.draw()
    bb = ax.get_position()
    fig.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=(bb.x0 + bb.width / 2, bb.y0 - pad),
        ncol=ncol,
        fontsize=FONT_SIZE_PT,
        frameon=False,
    )


def subplot_xlabel_centered(
    fig: plt.Figure, ax, text: str, pad: float = 0.028, tick_clearance: float = 0.0
) -> None:
    """Center an x-axis title under one subplot (not under a single tick)."""
    ax.set_xlabel("")
    fig.canvas.draw()
    bb = ax.get_position()
    fig.text(
        bb.x0 + bb.width / 2,
        bb.y0 - pad - tick_clearance,
        text,
        ha="center",
        va="top",
        fontsize=FONT_SIZE_PT,
        color=NEUTRAL,
    )


def legend_below(ax, handles=None, labels=None, ncol: int = 2, y: float = -0.24) -> None:
    """Place legend below the axes, clear of data."""
    kwargs = dict(
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=FONT_SIZE_PT,
        frameon=False,
        borderaxespad=0.0,
    )
    if handles is not None:
        ax.legend(handles=handles, labels=labels, **kwargs)
    elif labels is not None:
        ax.legend(labels=labels, **kwargs)
    else:
        ax.legend(**kwargs)


def annotate_threshold_hist(ax, n_ge: int) -> None:
    """Threshold note in axes coordinates, clear of the P=0.5 line and histogram peaks."""
    ax.text(
        0.54,
        0.94,
        f"Threshold 0.5\nn = {n_ge:,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_SIZE_PT,
        color=THRESHOLD,
        linespacing=1.35,
        clip_on=False,
    )


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


def label_bars_horizontal_outside(
    ax,
    bars,
    values,
    pad: float = 0.03,
    xmin: float = 0.05,
    label_col: float | None = None,
) -> None:
    """Place numeric labels outside horizontal bars; align high values in one column."""
    vals = list(values)
    if label_col is None and vals:
        label_col = min(max(vals) + pad + 0.04, 1.18)
    for bar, val in zip(bars, values):
        y = bar.get_y() + bar.get_height() / 2
        if val < 0.08:
            ax.text(xmin, y, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
        else:
            x = label_col if label_col is not None else bar.get_width() + pad
            ax.text(x, y, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    if vals:
        ax.set_xlim(0, max(ax.get_xlim()[1], (label_col or max(vals) + pad) + 0.06))


def label_hbars_counts(ax, bars, labels: list[str], pad_ratio: float = 0.02) -> None:
    """Place count labels to the right of bars."""
    xmax = ax.get_xlim()[1]
    pad = max(xmax * pad_ratio, 80)
    max_label_x = xmax
    for bar, txt in zip(bars, labels):
        y = bar.get_y() + bar.get_height() / 2
        numeric = txt.replace(",", "").replace(" ", "").isdigit()
        if numeric:
            x = bar.get_width() + pad
            ax.text(x, y, txt, va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
            max_label_x = max(max_label_x, x + xmax * 0.06)
        else:
            x = bar.get_width() + pad * 2.5
            ax.text(
                x,
                y,
                txt,
                va="center",
                ha="left",
                fontsize=FONT_SIZE_PT,
                color=MUTED,
                style="italic",
                clip_on=False,
            )
            max_label_x = max(max_label_x, x + xmax * 0.05)
    ax.set_xlim(0, max_label_x)


def save_figure(fig: plt.Figure, stem: str, subdir: str = "", tight: bool = True) -> dict[str, str]:
    out_dir = FIG_OUT / subdir if subdir else FIG_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for ext in ("pdf", "png"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(
            p,
            dpi=300 if ext == "png" else None,
            facecolor="white",
            bbox_inches="tight" if tight else None,
            pad_inches=0.12,
        )
        paths[ext] = str(p)
    plt.close(fig)
    return paths

# Standard margins (fraction of figure) — extra padding for print (~0.5–1 cm)
MARGIN_SINGLE = dict(left=0.17, right=0.94, top=0.76, bottom=0.18)
MARGIN_WIDE_X = dict(left=0.15, right=0.94, top=0.76, bottom=0.28)
MARGIN_COMPOSITE = dict(left=0.14, right=0.94, top=0.84, bottom=0.16, hspace=0.76, wspace=0.56)
