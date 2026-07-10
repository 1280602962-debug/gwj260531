#!/usr/bin/env python3
"""Generate figures for docs/report/REPORT.md from project references."""

from __future__ import annotations

import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "docs" / "report" / "figures"

PDB_IMAGES = {
    "8ELC": "https://cdn.rcsb.org/images/structures/8elc_assembly-1.jpeg",
    "3NPC": "https://cdn.rcsb.org/images/structures/3npc_assembly-1.jpeg",
    "4WHZ": "https://cdn.rcsb.org/images/structures/4whz_assembly-1.jpeg",
    "3V6S": "https://cdn.rcsb.org/images/structures/3v6s_assembly-1.jpeg",
    "7N8T": "https://cdn.rcsb.org/images/structures/7n8t_assembly-1.jpeg",
}

COMPOUNDS = [
    ("YL5084", "C[C@H]1CN(C(=O)c2ccc(NC(=O)/C=C/CN(C)C)cc2)C[C@H]1Nc1nccc(-c2c(-c3ccccc3)nn3ccccc23)n1", "JNK2 covalent hit [R1]"),
    ("YL2056", "CN(C)C/C=C/C(=O)Nc1ccc(C(=O)N2CC[C@H](Nc3nccc(-c4c(-c5ccccc5)nn5ccccc45)n3)C2)cc1", "8ELC ligand [R3]"),
    ("JNK-IN-8", "Cc1cc(NC(=O)c2cccc(NC(=O)/C=C/CN(C)C)c2)ccc1Nc1nccc(-c2cccnc2)n1", "Pan-JNK covalent [R2]"),
    ("56d", "C=CC(=O)Nc1cccc(C(=O)Nc2cccc(-n3cc(NC(=O)Nc4cccc5ccccc45)cn3)c2)c1", "Ligand-first covalent [R1]"),
    ("26k", "Clc1ccccc1NC(=O)Nc1cn(nc1)c2cccc(c2)C(=O)Nc2cn(nc2)[C@H]3CCNC3", "4WHZ ligand 26k/3NL [R7b]"),
    ("BIRB796", "CN1CCN(Cc2ccc(NC(=O)Nc3ccc(F)c(C(F)(F)F)c3)cc2C(F)(F)F)CC1=O", "DFG-out Type II [R7]"),
]


def download_pdb_images() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for pdb, url in PDB_IMAGES.items():
        dest = FIG_DIR / f"{pdb}_assembly.jpeg"
        if dest.exists() and dest.stat().st_size > 1000:
            continue
        print(f"Downloading {pdb} ...")
        urllib.request.urlretrieve(url, dest)


def fig01_five_stages() -> None:
    stages = [
        ("Stage 1\nLiterature & PDB", "8ELC / 4WHZ / 3NPC"),
        ("Stage 2\nSeeds & Decoys", "YL5084 + acrylamide decoys"),
        ("Stage 3\nAF3 Gate", "mPAE + EF@1% >= 2 [R15]"),
        ("Stage 4\nCovalent Dock", "DFG-in 8ELC only"),
        ("Stage 5\nSelectivity", "kinact/KI, C116S, cells"),
    ]
    fig, ax = plt.subplots(figsize=(14, 3.2))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 3)
    ax.axis("off")
    colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B"]
    for i, ((title, sub), color) in enumerate(zip(stages, colors)):
        x = 0.4 + i * 2.7
        box = FancyBboxPatch(
            (x, 0.8), 2.3, 1.6,
            boxstyle="round,pad=0.08,rounding_size=0.15",
            facecolor=color, edgecolor="white", linewidth=2, alpha=0.92,
        )
        ax.add_patch(box)
        ax.text(x + 1.15, 1.65, title, ha="center", va="center", fontsize=10,
                fontweight="bold", color="white")
        ax.text(x + 1.15, 1.05, sub, ha="center", va="center", fontsize=7.5,
                color="#f0f0f0", wrap=True)
        if i < len(stages) - 1:
            ax.annotate(
                "", xy=(x + 2.45, 1.6), xytext=(x + 2.25, 1.6),
                arrowprops=dict(arrowstyle="->", color="#555", lw=2),
            )
    ax.set_title("JNK2 Covalent Selectivity — Five-Stage Decision Workflow", fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig01_five_stages_flow.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig02_compounds() -> None:
    mols, legends = [], []
    for name, smi, note in COMPOUNDS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"WARN: skip {name}, invalid SMILES")
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        mols.append(mol)
        legends.append(f"{name}\n{note}")
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(320, 260), legends=legends, useSVG=False)
    img.save(FIG_DIR / "fig02_key_compounds_2d.png")


def fig03_dfg_schematic() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, title, dfg, cys, status, color in [
        (axes[0], "DFG-in (8ELC) — covalent primary", "IN", "Cys116 accessible", "Use for covalent AF3/docking", "#2E7D32"),
        (axes[1], "DFG-out (3NPC) — exclude covalent", "OUT", "Cys116 occluded", "BIRB796 reversible only", "#C62828"),
    ]:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")
        ax.add_patch(FancyBboxPatch((1, 2), 8, 6, boxstyle="round", facecolor="#ECEFF1", edgecolor="#455A64", lw=2))
        ax.add_patch(FancyBboxPatch((2, 5.5), 6, 2, boxstyle="round", facecolor="#BBDEFB", edgecolor="#1565C0"))
        ax.text(5, 6.5, "ATP pocket", ha="center", fontsize=10, fontweight="bold")
        ax.add_patch(mpatches.Circle((3.2, 4.2), 0.35, color="#FF9800"))
        ax.text(3.2, 3.5, "Cys116", ha="center", fontsize=8)
        ax.add_patch(FancyBboxPatch((5.5, 3.2), 2.5, 1.2, boxstyle="round", facecolor=color, alpha=0.85))
        ax.text(6.75, 3.8, f"DFG {dfg}", ha="center", va="center", color="white", fontweight="bold")
        ax.text(5, 1.3, cys, ha="center", fontsize=9)
        ax.text(5, 0.6, status, ha="center", fontsize=9, fontweight="bold", color=color)
        ax.set_title(title, fontsize=11, fontweight="bold")
    fig.suptitle("Why covalent primary screening uses DFG-in (8ELC), not DFG-out (3NPC)", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig03_dfg_in_vs_out_schematic.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig04_leu106() -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    for x, iso, res, pocket, sel in [
        (2.5, "JNK1", "Ile106", "Smaller HR-I pocket", "Weaker covalent fit"),
        (7.5, "JNK2", "Leu106", "Larger HR-I pocket", "148–340× vs JNK1 (reversible)\n~21× kinact/KI (YL5084)"),
    ]:
        ax.add_patch(FancyBboxPatch((x - 1.8, 1), 3.6, 4.2, boxstyle="round", facecolor="#E3F2FD", edgecolor="#1976D2", lw=2))
        ax.text(x, 4.7, iso, ha="center", fontsize=14, fontweight="bold")
        ax.text(x, 3.8, res, ha="center", fontsize=12, color="#D32F2F", fontweight="bold")
        ax.text(x, 2.8, pocket, ha="center", fontsize=9)
        ax.text(x, 1.8, sel, ha="center", fontsize=8.5, style="italic")
    ax.annotate("", xy=(5.5, 3), xytext=(4.5, 3), arrowprops=dict(arrowstyle="<->", color="#555", lw=2))
    ax.text(5, 5.3, "Selectivity axis (Lu 2023; Wydra 2025; Zheng 2014)", ha="center", fontsize=11, fontweight="bold")
    fig.savefig(FIG_DIR / "fig04_leu106_selectivity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig05_pdb_panel() -> None:
    panels = [
        ("8ELC", "YL2056\nDFG-in covalent [R3]"),
        ("3NPC", "BIRB796\nDFG-out [R7]"),
        ("4WHZ", "26k\nDFG-in reversible [R7b]"),
        ("3V6S", "JNK-IN-7\nJNK3 covalent [R2]"),
        ("7N8T", "JNK2–AMP\n1.6 Å [R4]"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()
    for ax, (pdb, caption) in zip(axes, panels):
        path = FIG_DIR / f"{pdb}_assembly.jpeg"
        img = Image.open(path)
        ax.imshow(img)
        ax.set_title(f"{pdb}\n{caption}", fontsize=10, fontweight="bold")
        ax.axis("off")
    axes[-1].axis("off")
    fig.suptitle("Reference PDB structures (RCSB assembly images)", fontsize=13, y=0.98)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig05_pdb_structure_panel.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    download_pdb_images()
    fig01_five_stages()
    fig02_compounds()
    fig03_dfg_schematic()
    fig04_leu106()
    fig05_pdb_panel()
    print(f"Figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
