"""
Journal-style matplotlib settings for pipeline figures.

- Font: Arial (fallback: Helvetica, DejaVu Sans)
- Language: English labels only
- Resolution: 300 dpi
- Sizes aligned with single/double-column journal layouts (~89 mm / ~183 mm)
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# Register Arial when available (e.g. ttf-mscorefonts-installer on Linux)
for _arial_path in (
    Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
    Path("/usr/share/fonts/truetype/msttcorefonts/arial.ttf"),
):
    if _arial_path.exists():
        mpl.font_manager.fontManager.addfont(str(_arial_path))
        break

# Nature / ACS typical dimensions (inches)
FIGSIZE_SINGLE = (3.5, 2.75)   # ~89 mm wide
FIGSIZE_DOUBLE = (7.2, 3.0)   # ~183 mm wide
FIGSIZE_SQUARE = (3.5, 3.5)
DPI = 300

JOURNAL_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.titlesize": 9,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.0,
}


def apply_journal_style() -> None:
    """Apply global rcParams for all pipeline figures."""
    mpl.rcParams.update(JOURNAL_RC)


def save_figure(path, fig=None) -> None:
    """Save figure at 300 dpi with tight bounding box."""
    apply_journal_style()
    if fig is None:
        fig = plt.gcf()
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
