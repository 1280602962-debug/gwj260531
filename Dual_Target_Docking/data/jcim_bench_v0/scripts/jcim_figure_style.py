"""Shared JCIM Article figure style (Liberation Sans ≈ Arial/Helvetica).

ACS/JCIM artwork: double-column figures ≤ 7.00 in wide, ≤ 9.167 in deep;
color 300 dpi; RGB TIFF/PNG; PDF with embedded TrueType (fonttype 42).
TOC graphic: exactly 3.25 × 1.75 in, 300 dpi TIFF, sans-serif ≥ 6 pt.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "figures" / "jcim_article"
OUT.mkdir(parents=True, exist_ok=True)

# ACS/JCIM: sans-serif, embed TrueType in PDF
FONT = "Liberation Sans"
FS_PANEL = 10
FS_AXIS = 8
FS_TICK = 7
FS_LEGEND = 7
FS_ANNO = 7
FS_TOC = 8

# Okabe–Ito (colorblind-safe)
C = {
    "vina": "#0072B2",
    "rtm": "#009E73",
    "gnina": "#56B4E9",
    "desc": "#E69F00",
    "dual": "#0072B2",
    "a_only": "#D55E00",
    "b_only": "#CC79A7",
    "neither": "#999999",
    "chance": "#666666",
    "thick": "#009E73",
    "egfr": "#D55E00",
    "metal": "#B0B0B0",
    "other": "#A6CEE3",
    "holdout": "#E69F00",
    "main": "#0072B2",
    "swap_bad": "#D55E00",
    "swap_ok": "#0072B2",
    "ink": "#222222",
}

PAIR_ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
PAIR_SHORT = {
    "EGFR/HER2": "EGFR/HER2",
    "AChE/BChE": "AChE/BChE",
    "PIK3CA/PIK3CB": "PIK3CA/PIK3CB",
    "PIK3CA/mTOR": "PIK3CA/mTOR",
}
DESC_LABEL = {
    "heavy": "heavy atoms",
    "mw": "MW",
    "clogp": "cLogP",
    "tpsa": "TPSA",
}


def apply_style() -> None:
    found = findfont(FontProperties(family=FONT))
    if "LiberationSans" not in found.replace(" ", "") and "Liberation" not in found:
        raise RuntimeError(f"Liberation Sans not available (got {found})")
    mpl.rcParams.update(
        {
            "font.family": FONT,
            "font.size": FS_TICK,
            "axes.labelsize": FS_AXIS,
            "axes.titlesize": FS_AXIS,
            "xtick.labelsize": FS_TICK,
            "ytick.labelsize": FS_TICK,
            "legend.fontsize": FS_LEGEND,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": C["ink"],
            "text.color": C["ink"],
            "axes.labelcolor": C["ink"],
            "xtick.color": C["ink"],
            "ytick.color": C["ink"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "pdf.use14corefonts": False,
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "standard",
            "savefig.pad_inches": 0.04,
            "legend.frameon": False,
            "legend.handlelength": 1.4,
            "legend.handletextpad": 0.4,
            "legend.borderaxespad": 0.2,
            "legend.columnspacing": 1.0,
        }
    )


def panel_label(ax, letter: str, x: float = -0.12, y: float = 1.08) -> None:
    """Place a panel letter in axes coordinates, outside the data area."""
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=FS_PANEL,
        fontweight="bold",
        fontfamily=FONT,
        va="bottom",
        ha="left",
        clip_on=False,
    )


def _rgb_file(path: Path) -> None:
    """ACS prefers RGB TIFF/PNG, not RGBA."""
    from PIL import Image

    im = Image.open(path)
    dpi = im.info.get("dpi", (300, 300))
    if isinstance(dpi, (int, float)):
        dpi = (float(dpi), float(dpi))
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im.close()
        im = bg
    elif im.mode != "RGB":
        converted = im.convert("RGB")
        im.close()
        im = converted
    suf = path.suffix.lower()
    if suf in {".tif", ".tiff"}:
        im.save(path, format="TIFF", dpi=dpi, compression="tiff_lzw")
    elif suf == ".png":
        im.save(path, format="PNG", dpi=dpi)
    im.close()


def save_all(fig, stem: str, toc: bool = False):
    """Write PDF/PNG/TIF (main figures) or exact-size TOC TIF+PNG."""
    OUT.mkdir(parents=True, exist_ok=True)
    paths = []
    if toc:
        # Exact ACS TOC size; do not use bbox_inches='tight' (would change dimensions)
        kw_toc = dict(dpi=300, bbox_inches=None, pad_inches=0, facecolor="white", edgecolor="none")
        tif = OUT / f"{stem}.tif"
        fig.savefig(
            tif,
            format="tiff",
            pil_kwargs={"compression": "tiff_lzw"},
            **kw_toc,
        )
        paths.append(tif)
        png = OUT / f"{stem}.png"
        fig.savefig(png, **kw_toc)
        paths.append(png)
        _rgb_file(tif)
        _rgb_file(png)
    else:
        for ext in ("pdf", "png", "tif"):
            p = OUT / f"{stem}.{ext}"
            kw = {"dpi": 300, "facecolor": "white", "edgecolor": "none", "bbox_inches": None, "pad_inches": 0.02}
            if ext == "tif":
                kw["format"] = "tiff"
                kw["pil_kwargs"] = {"compression": "tiff_lzw"}
            fig.savefig(p, **kw)
            paths.append(p)
            if ext in ("png", "tif"):
                _rgb_file(p)
    return paths
