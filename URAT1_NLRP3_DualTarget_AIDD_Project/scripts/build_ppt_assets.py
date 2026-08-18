#!/usr/bin/env python3
"""Build PPT-ready figure assets: PDB images, OA literature pages, schematics, project plots."""

from __future__ import annotations

import json
import shutil
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "figures" / "ppt_assets"
STRUCT = OUT / "structures"
LIT = OUT / "literature"
SCHEM = OUT / "schematics"
PROJ = OUT / "project_results"

RCSB_PDbs = ["9dkb", "7alv", "8etr", "9dka", "9b1k", "9b1l"]
FRONTIERS_LIU = [
    ("g001", "https://www.frontiersin.org/files/Articles/1137822/fimmu-14-1137822-HTML/image_m/fimmu-14-1137822-g001.jpg"),
    ("g002", "https://www.frontiersin.org/files/Articles/1137822/fimmu-14-1137822-HTML/image_m/fimmu-14-1137822-g002.jpg"),
    ("g003", "https://www.frontiersin.org/files/Articles/1137822/fimmu-14-1137822-HTML/image_m/fimmu-14-1137822-g003.jpg"),
]
OA_PDFS = [
    (
        "eurycoma2025",
        "https://www.nature.com/articles/s41467-025-62645-6.pdf",
        "Zhang et al. Nat Commun 2025 (CC BY) — dual-target anti-gout from Eurycoma",
        [2, 3, 4, 5],
    ),
    (
        "fedor2025_urat1",
        "https://www.nature.com/articles/s41467-025-60480-3.pdf",
        "Fedor/Suo et al. Nat Commun 2025 (CC BY) — URAT1 cryo-EM 9DKB",
        [1, 2, 3, 4],
    ),
]
PROJECT_PNGS = [
    ("fig02_nlrp3_screening_composite.png", "NLRP3 ML prescreen (project Fig 2)"),
    ("fig03_urat1_retrospective_composite.png", "8973 URAT1 retrospective (project Fig 3)"),
    ("fig04_pareto_dual_docking_9dkb_7alv.png", "Pareto dual docking (project Fig 4)"),
    ("si_data_asymmetry.png", "ChEMBL data asymmetry (SI)"),
    ("si_nlrp3_oof_roc_pr.png", "NLRP3 ML OOF ROC/PR (SI)"),
]


def download(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, dest)
        return dest.stat().st_size > 0
    except Exception:
        return False


def fetch_rcsb_images() -> list[dict]:
    records = []
    for pdb in RCSB_PDbs:
        dest = STRUCT / f"{pdb}_assembly-1.jpeg"
        url = f"https://cdn.rcsb.org/images/structures/{pdb}_assembly-1.jpeg"
        ok = download(url, dest)
        records.append(
            {
                "file": str(dest.relative_to(PROJECT)),
                "source": f"RCSB PDB ({pdb.upper()})",
                "license": "RCSB PDB free use with citation (CC0 for structure images)",
                "ok": ok,
            }
        )
    return records


def fetch_frontiers() -> list[dict]:
    records = []
    for tag, url in FRONTIERS_LIU:
        dest = LIT / f"frontiers2023_liu_nlrp3_gout_{tag}.jpg"
        ok = download(url, dest)
        records.append(
            {
                "file": str(dest.relative_to(PROJECT)),
                "source": "Liu et al. Front Immunol 2023; doi:10.3389/fimmu.2023.1137822",
                "license": "CC BY 4.0",
                "ppt_use": f"Background: NLRP3 in gout ({tag})",
                "ok": ok,
            }
        )
    return records


def extract_pdf_pages() -> list[dict]:
    try:
        import fitz
    except ImportError:
        return [{"error": "pymupdf not installed; skip PDF extraction"}]

    records = []
    for prefix, url, cite, pages in OA_PDFS:
        pdf_path = LIT / f"{prefix}_article.pdf"
        if not pdf_path.exists():
            download(url, pdf_path)
        if not pdf_path.exists():
            records.append({"prefix": prefix, "ok": False, "error": "PDF download failed"})
            continue
        doc = fitz.open(str(pdf_path))
        for pno in pages:
            if pno - 1 >= len(doc):
                continue
            dest = LIT / f"{prefix}_pdf_page{pno}.png"
            pix = doc[pno - 1].get_pixmap(matrix=fitz.Matrix(2.2, 2.2))
            pix.save(str(dest))
            records.append(
                {
                    "file": str(dest.relative_to(PROJECT)),
                    "source": cite,
                    "license": "CC BY 4.0 (Nature Communications)",
                    "ppt_use": f"Literature figure spread — PDF page {pno}",
                    "ok": True,
                }
            )
        doc.close()
    return records


def copy_project_figures() -> list[dict]:
    records = []
    src_dir = PROJECT / "figures" / "generated"
    for name, desc in PROJECT_PNGS:
        src = None
        for sub in ("main", "si"):
            candidate = src_dir / sub / name
            if candidate.exists():
                src = candidate
                break
        if src is None:
            records.append({"file": name, "ok": False})
            continue
        dest = PROJ / name
        shutil.copy2(src, dest)
        records.append(
            {
                "file": str(dest.relative_to(PROJECT)),
                "source": "This project (plot_available_figures.py)",
                "license": "Project MIT",
                "ppt_use": desc,
                "ok": True,
            }
        )
    return records


def draw_gout_dual_pathway() -> Path:
    SCHEM.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")

    def box(x, y, w, h, text, color):
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02", linewidth=1.2, edgecolor=color, facecolor=color, alpha=0.15
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10, wrap=True)

    box(0.3, 3.2, 2.2, 1.0, "Hyperuricemia\n(renal URAT1 reabsorption)", "#0072B2")
    box(3.0, 3.2, 2.0, 1.0, "MSU crystal\ndeposition", "#666666")
    box(5.7, 3.2, 2.3, 1.0, "NLRP3 inflammasome\nIL-1β release", "#D55E00")
    box(8.3, 3.2, 1.5, 1.0, "Gout flare\n& chronic inflammation", "#CC6677")

    for x1, x2 in [(2.5, 3.0), (5.0, 5.7), (8.0, 8.3)]:
        ax.add_patch(
            FancyArrowPatch((x1, 3.7), (x2, 3.7), arrowstyle="->", mutation_scale=12, linewidth=1.5, color="#333333")
        )

    box(0.5, 1.0, 3.5, 1.2, "Metabolic axis\nULT / URAT1 inhibitors", "#0072B2")
    box(6.0, 1.0, 3.5, 1.2, "Inflammatory axis\nNSAID / colchicine / NLRP3 modulators", "#D55E00")
    ax.text(5, 0.3, "Clinical practice: often combination therapy — not a single validated dual-node drug", ha="center", fontsize=9, style="italic")

    dest = SCHEM / "schematic_gout_dual_pathway.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def draw_workflow() -> Path:
    SCHEM.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 5)
    ax.axis("off")
    steps = [
        (0.2, "ChEMBL clinical library\nn = 8,319", "#444444"),
        (2.0, "NLRP3 ML prescreen\nP≥0.5 → 1,588", "#D55E00"),
        (3.8, "gnina P2 dual dock\n9DKB + 7ALV", "#0072B2"),
        (5.6, "Percentiles +\nPareto audit", "#009E73"),
        (7.4, "Chemistry nomination\n(dual-dock gate)", "#6A3D9A"),
        (9.2, "MD benchmarks\n+ hypotheses", "#666666"),
    ]
    for i, (x, text, c) in enumerate(steps):
        box = FancyBboxPatch((x, 2.0), 1.5, 1.3, boxstyle="round,pad=0.02", linewidth=1.2, edgecolor=c, facecolor=c, alpha=0.12)
        ax.add_patch(box)
        ax.text(x + 0.75, 2.65, text, ha="center", va="center", fontsize=8)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 1.5, 2.65), (steps[i + 1][0], 2.65), arrowstyle="->", mutation_scale=12, color="#333"))
    ax.text(2.0, 0.8, "Parallel track: 8973 distill → URAT1 enrichment only (A vs D)", ha="center", fontsize=9, color="#0072B2")
    dest = SCHEM / "schematic_project_workflow.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def draw_inspiration_comparison() -> Path:
    SCHEM.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.3, 1.0), 4.2, 2.5, boxstyle="round,pad=0.03", linewidth=1.5, edgecolor="#009E73", facecolor="#009E73", alpha=0.08))
    ax.text(2.4, 3.1, "Nat Commun 2025 (Eurycoma)", ha="center", fontsize=11, fontweight="bold")
    ax.text(2.4, 2.3, "Phenotypic screen → 64 derivatives\nURAT1 + NLRP3 wet validation\nCompound 32 (de novo)", ha="center", fontsize=9)

    ax.add_patch(FancyBboxPatch((5.5, 1.0), 4.2, 2.5, boxstyle="round,pad=0.03", linewidth=1.5, edgecolor="#0072B2", facecolor="#0072B2", alpha=0.08))
    ax.text(7.6, 3.1, "This project", ha="center", fontsize=11, fontweight="bold")
    ax.text(7.6, 2.3, "Protocol-first gnina P2\nNLRP3 ML shrink + dual dock\nChemistry nomination (computational)", ha="center", fontsize=9)

    ax.add_patch(FancyArrowPatch((4.6, 2.25), (5.4, 2.25), arrowstyle="<->", mutation_scale=14, linewidth=1.5, color="#333"))
    ax.text(5.0, 2.55, "complementary", ha="center", fontsize=9, style="italic")
    dest = SCHEM / "schematic_inspiration_vs_this_project.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def main() -> None:
    for d in (STRUCT, LIT, SCHEM, PROJ):
        d.mkdir(parents=True, exist_ok=True)

    manifest = {
        "generated_by": "scripts/build_ppt_assets.py",
        "rcsb": fetch_rcsb_images(),
        "frontiers": fetch_frontiers(),
        "literature_pdf_pages": extract_pdf_pages(),
        "project_results": copy_project_figures(),
        "schematics": [
            {"file": str(draw_gout_dual_pathway().relative_to(PROJECT)), "license": "Original (project)", "ppt_use": "Slide: disease background"},
            {"file": str(draw_workflow().relative_to(PROJECT)), "license": "Original (project)", "ppt_use": "Slide: methods workflow"},
            {"file": str(draw_inspiration_comparison().relative_to(PROJECT)), "license": "Original (project)", "ppt_use": "Slide: motivation vs Eurycoma 2025"},
        ],
    }
    manifest_path = OUT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
