"""Shared JCIM Article figure style (ACS Arial).

ACS/JCIM artwork: double-column figures ≤ 7.00 in wide, ≤ 9.167 in deep;
color 300 dpi; RGB TIFF/PNG; PDF with embedded TrueType (fonttype 42).
TOC graphic: exactly 3.25 × 1.75 in, 300 dpi TIFF, sans-serif ≥ 6 pt.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont, fontManager

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "figures" / "jcim_article"
OUT.mkdir(parents=True, exist_ok=True)

# ACS/JCIM prefers Arial or Helvetica. Artwork is generated in Arial.
FONT = "Arial"
ARIAL_PATH: Path | None = None
ARIAL_REGULAR_NAMES = ("Arial.ttf", "arial.ttf", "Arial.TTF")
ARIAL_FAMILY_FILES = (
    ("Arial.ttf", "arial.ttf", "Arial.TTF"),
    ("Arial-Bold.ttf", "arialbd.ttf", "Arialbd.TTF", "Arialbd.ttf"),
    ("Arial-Italic.ttf", "ariali.ttf", "Ariali.TTF", "Ariali.ttf"),
    ("Arial-BoldItalic.ttf", "arialbi.ttf", "Arialbi.TTF", "Arialbi.ttf"),
)
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

# Historical four-pair order (v1/v2 SI). Do not use for the eight-row primary set.
PAIR_ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
# Display order is protein-system grouped (kinase, hydrolase, protease, NR).
# The next two lists only say which frozen CSV holds that pair's Table 2 row.
PRIMARY_PAIRS = [
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
]
UNIFIED_THRESHOLD_PAIRS = [
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
]
COMPARABLE_THETA6_PAIRS = []
HOLDOUT_PAIRS = [p for p in PRIMARY_PAIRS if p != "EGFR/HER2"]
PAIR_SHORT = {
    "EGFR/HER2": "EGFR/HER2",
    "AChE/BChE": "AChE/BChE",
    "PIK3CA/PIK3CB": "PIK3CA/PIK3CB",
    "PIK3CA/mTOR": "PIK3CA/mTOR",
    "F2/F10": "F2/F10",
    "JAK1/TYK2": "JAK1/TYK2",
    "JAK1/JAK2": "JAK1/JAK2",
    "PPARG/PPARA": "PPARG/PPARA",
    "PPARA/PPARD": "PPARA/PPARD",
}
# Okabe–Ito plus two extra colorblind-safe inks for eight pair traces.
PAIR_COLOR = {
    "EGFR/HER2": "#D55E00",
    "AChE/BChE": "#009E73",
    "PIK3CA/mTOR": "#0072B2",
    "F2/F10": "#E69F00",
    "JAK1/TYK2": "#56B4E9",
    "JAK1/JAK2": "#CC79A7",
    "PPARG/PPARA": "#000000",
    "PPARA/PPARD": "#999999",
    "PIK3CA/PIK3CB": "#882255",
}
DESC_LABEL = {
    "heavy": "heavy atoms",
    "mw": "MW",
    "clogp": "cLogP",
    "tpsa": "TPSA",
}


def _arial_search_dirs() -> list[Path]:
    return [
        Path(r"C:\Windows\Fonts"),
        Path("/usr/share/fonts/truetype/msttcorefonts"),
        Path("/usr/share/fonts/truetype/msttcore"),
        Path.home() / ".local" / "share" / "fonts",
        Path("/usr/local/share/fonts"),
        Path(__file__).resolve().parents[1] / "fonts",
    ]


def _register_arial() -> Path:
    """Register Arial TTFs so headless Matplotlib does not fall back to DejaVu."""
    dirs = [d for d in _arial_search_dirs() if d.exists()]
    registered = []
    for names in ARIAL_FAMILY_FILES:
        hit = next((d / name for d in dirs for name in names if (d / name).exists()), None)
        if hit is None:
            continue
        try:
            fontManager.addfont(str(hit))
        except OSError:
            continue
        registered.append(hit)
    regular = next((p for p in registered if p.name.lower() in {n.lower() for n in ARIAL_REGULAR_NAMES}), None)
    if regular is None:
        try:
            found = Path(findfont(FontProperties(family="Arial"), fallback_to_default=False))
        except ValueError as exc:
            raise FileNotFoundError(
                "Arial is required for JCIM artwork. Install Arial "
                "(Windows Fonts, msttcorefonts, or figures/jcim_article/fonts/)."
            ) from exc
        if "arial" not in found.name.lower():
            raise FileNotFoundError(
                f"Arial is required for JCIM artwork; Matplotlib resolved {found}."
            )
        return found
    resolved = Path(findfont(FontProperties(family="Arial"), fallback_to_default=False))
    if "arial" not in resolved.name.lower():
        raise FileNotFoundError(
            f"Arial registration failed; Matplotlib resolved {resolved}."
        )
    return resolved


def apply_style() -> None:
    global FONT, ARIAL_PATH
    FONT = "Arial"
    ARIAL_PATH = _register_arial()
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT,
            "mathtext.it": f"{FONT}:italic",
            "mathtext.bf": f"{FONT}:bold",
            "mathtext.cal": FONT,
            "mathtext.default": "regular",
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
    print(f"figure font: Arial ({ARIAL_PATH})")


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
    from matplotlib.text import Text
    # Presentation labels use the same mathematical typography throughout.
    # Source keys and numeric provenance remain unchanged.
    for item in fig.findobj(match=Text):
        label = item.get_text()
        for raw, display in {
            "summary_min": r"summary$_{\mathrm{min}}$",
            "vina_mean": r"Vina$_{\mathrm{mean}}$",
            "vina_worst": r"Vina$_{\mathrm{worst}}$",
        }.items():
            label = label.replace(raw, display)
        item.set_text(label)
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
