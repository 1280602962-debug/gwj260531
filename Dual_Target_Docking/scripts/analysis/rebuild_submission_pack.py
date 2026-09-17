#!/usr/bin/env python3
"""Copy current publication files into Dual_Target_Docking/submission_pack/.

Zero-dock. Does not copy remediation, audit reports, or historical packs.
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
ART = ROOT / "figures" / "jcim_article"
CANON = ROOT / "results" / "canonical"
PACK = ROOT / "submission_pack"

STEMS = (
    "Fig1_four_state_and_supply",
    "Fig1_C_chEMBL_supply",
    "Fig2_negative_class_formulation",
    "Fig3_ligand_chemistry",
    "Fig4_mismatched_pocket",
    "Fig5_computational_realization",
    "FigS1_ligand_chemistry_detail",
    "FigS2_protocol_sensitivity",
    "FigS3_cognate_rmsd",
    "FigS4_label_source_robustness",
    "FigS5_external_eligibility",
    "TOC_graphic",
)


def copyfile(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def main() -> int:
    if PACK.exists():
        shutil.rmtree(PACK)
    (PACK / "manuscript").mkdir(parents=True)
    (PACK / "SI").mkdir()
    (PACK / "figures").mkdir()
    (PACK / "tables" / "canonical").mkdir(parents=True)

    copyfile(DOCS / "MANUSCRIPT_JCIM_EN.md", PACK / "manuscript" / "MANUSCRIPT_JCIM_EN.md")
    copyfile(DOCS / "MANUSCRIPT_JCIM_ZH.md", PACK / "manuscript" / "MANUSCRIPT_JCIM_ZH.md")
    copyfile(DOCS / "FIGURE_TABLE_LOCK_POSTFIX_V4.md", PACK / "manuscript" / "FIGURE_TABLE_LOCK_POSTFIX_V4.md")
    copyfile(DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md", PACK / "SI" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md")
    copyfile(
        DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
        PACK / "SI" / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    )
    copyfile(ART / "CAPTIONS.md", PACK / "figures" / "CAPTIONS.md")
    copyfile(ART / "MANUSCRIPT_FIGURE_CAPTIONS.md", PACK / "figures" / "MANUSCRIPT_FIGURE_CAPTIONS.md")
    if (ART / "plotted_values_postfix.json").is_file():
        copyfile(ART / "plotted_values_postfix.json", PACK / "figures" / "plotted_values_postfix.json")

    aliases = {
        "Fig1_four_state_and_supply": "Figure1",
        "Fig2_negative_class_formulation": "Figure2",
        "Fig3_ligand_chemistry": "Figure3",
        "Fig4_mismatched_pocket": "Figure4",
        "Fig5_computational_realization": "Figure5",
        "FigS1_ligand_chemistry_detail": "FigureS1",
        "FigS2_protocol_sensitivity": "FigureS2",
        "FigS3_cognate_rmsd": "FigureS3",
        "FigS4_label_source_robustness": "FigureS4",
        "FigS5_external_eligibility": "FigureS5",
        "TOC_graphic": "TOC",
    }
    for stem in STEMS:
        for ext in (".pdf", ".tif", ".png"):
            src = ART / f"{stem}{ext}"
            if src.is_file():
                copyfile(src, PACK / "figures" / src.name)
                alias = aliases.get(stem)
                if alias:
                    copyfile(src, PACK / "figures" / f"{alias}{ext}")

    for csv in sorted(CANON.glob("*.csv")):
        copyfile(csv, PACK / "tables" / "canonical" / csv.name)
    for box in (
        ROOT / "data/egfr_her2_panel120_v0/boxes/3POZ_box_corrected.json",
        ROOT / "data/egfr_her2_panel120_v0/boxes/3RCD_box_corrected.json",
    ):
        if box.is_file():
            copyfile(box, PACK / "tables" / "source" / box.name)

    readme = PACK / "README.md"
    readme.write_text(
        "# JCIM submission pack\n\n"
        "Current English/Chinese manuscripts, SI, figures, and canonical result tables.\n"
        "Authoritative per-ligand scores: `tables/canonical/current_score_master.csv`.\n"
        "Regenerate from Dual_Target_Docking with "
        "`python3 scripts/analysis/compute_canonical_results.py` then "
        "`python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root Dual_Target_Docking`.\n",
        encoding="utf-8",
    )
    n = sum(1 for p in PACK.rglob("*") if p.is_file())
    print(f"wrote {PACK} ({n} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
