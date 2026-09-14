#!/usr/bin/env python3
"""Copy submission-facing manuscripts, tables, figures, and scripts into submission_pack/."""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACK = ROOT / "submission_pack"

DOCS = [
    "MANUSCRIPT_JCIM_EN.md",
    "MANUSCRIPT_JCIM_ZH.md",
    "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    "STATISTICAL_LOCK_V1.md",
    "FIGURE_PANEL_LOCK_V3.md",
    "SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md",
    "DATA_AND_SOFTWARE_AVAILABILITY_JCIM_EN.md",
    "DATA_AND_SOFTWARE_AVAILABILITY_JCIM_ZH.md",
    "REFERENCES_JCIM.md",
    "PRODUCTION_POSE_ARCHIVE_PREP_V1.md",
    "HIGH_COMPARABILITY_FEASIBILITY_V1.md",
    "PRIMARY_RESULT_INDEX_V1.md",
    "MANUSCRIPT_CONSOLIDATION_REPORT_20260912.md",
    "LANGUAGE_CHANGELOG.md",
    "FACT_CHECK_TABLE.md",
    "UNRESOLVED_ISSUES.md",
]

TABLES = [
    "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv",
    "data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv",
    "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv",
    "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv",
    "data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv",
    "data/jcim_novelty_v0/tables/incremental_information_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv",
    "data/jcim_novelty_v0/tables/ecfp4_docking_scaler_sensitivity_v1.csv",
    "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv",
    "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv",
    "data/jcim_multiseed_v0/tables/multiseed_auroc_aggregate_v2.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/fiveseed_summary_min_aggregate_v1.csv",
    "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv",
    "data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv",
    "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv",
    "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_modes_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/layer3_cognate_rmsd_calcrrms_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/layer3_cognate_rmsd_calcrrms_modes_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/pocket_unidirectional_delta_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/multiseed_fixed_membership_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/multiseed_fixed_membership_ids_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/production_pose_sha256_manifest_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_pose_replay_v1.csv",
    "data/jcim_novelty_v0/tables/DUALFOURCLASS_EVALUATION_CONTRACT_v1.json",
    "data/jcim_novelty_v0/tables/REVISION_CHECKSUM_MANIFEST_v1.csv",
    "data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv",
    "data/jcim_j0j1_v0/tables/j0_strict_label_supply.csv",
    "data/jcim_chembl_universe_v0/tables/track_b_panel_summary_v1.csv",
    "data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv",
    "data/jcim_novelty_v0/tables/primary_directional_intervals_review_v1.csv",
    "data/jcim_novelty_v0/tables/summary_min_stratified_sensitivity_review_v1.csv",
    "data/jcim_novelty_v0/tables/gnina_jak_interval_provenance_review_v1.csv",
    "data/jcim_novelty_v0/tables/operating_point_examples_review_v1.csv",
    "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv",
    "data/jcim_novelty_v0/tables/review_statistics_provenance_v1.json",
]

# Experiment scripts only. Run them from Dual_Target_Docking/, not this flattened folder.
# Do not pack assemble / audit / validate / checksum / freeze / pack helpers.
SCRIPTS = [
    (
        "data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py",
        "Table 2 / S3 θ grid and ligand-chemistry baselines for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR from frozen Vina scores.",
    ),
    (
        "data/jcim_novelty_v0/scripts/benchmark_formulation_v1.py",
        "Table 3, Table S4 fixed-channel Δ, and Table S5 incremental analyses for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR.",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/analyze_five_pair_stack_v1.py",
        "Table 2 / Table 3 / Table S4 / Table S5 for F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, and PPARA/PPARD from frozen scores. Does not redock.",
    ),
    (
        "data/jcim_strengthen_t0t1_v0/scripts/build_p0_missing_tables_v1.py",
        "Matched-versus-mismatched pocket bootstrap for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S6 / Figure 4).",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/analyze_five_pair_local_channels_v1.py",
        "Pocket channels, holdout join, JAK independent GNINA, and five-seed aggregate for F2/F10, JAK, and PPAR pairs (Tables S6–S7 / Figure 5).",
    ),
    (
        "data/jcim_holdout_v0/scripts/analyze_holdout_v1.py",
        "Unused-pool holdout AUROCs in Table S6 / Figure 4. Not external validation.",
    ),
    (
        "data/jcim_independent_dock_v0/scripts/analyze_independent_dock_v1.py",
        "Independent GNINA formulation AUROCs for EGFR/HER2 and PIK3CA/mTOR in Table S7 / Figure 5.",
    ),
    (
        "data/jcim_multiseed_v0/scripts/analyze_multiseed_vina_v2.py",
        "Five-seed AUROC aggregate for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Figure 5C).",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/multiseed_fixed_membership_v1.py",
        "Fixed-membership five-seed summary for F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, and PPARA/PPARD (Figure 5C). Does not redock.",
    ),
    (
        "data/jcim_novelty_v0/scripts/equal_score_cluster_bootstrap_v1.py",
        "Document/scaffold cluster resampling of the two flagship fixed-channel Δ values (Table S4 / Figure 6).",
    ),
    (
        "data/jcim_novelty_v0/scripts/claim_hardening_v1.py",
        "Four-descriptor directional AUROCs for the Table 2 descriptor column and Figure S2.",
    ),
    (
        "data/jcim_novelty_v0/scripts/bindingdb_native_slice_eight_pairs_v1.py",
        "Independence-filtered BindingDB/PubChem eligibility counts for Table S8 / Figure 6. No external docking.",
    ),
    (
        "data/jcim_novelty_v0/scripts/assay_aggregation_max_vs_median_v1.py",
        "API max-versus-median label sensitivity for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S3).",
    ),
    (
        "data/jcim_novelty_v0/scripts/high_confidence_label_rebuild_v1.py",
        "High-confidence field screen for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S3).",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/analyze_five_pair_dump_gated_v1.py",
        "Dump-gated max-versus-median footnote for F2/F10, JAK, and PPAR pairs (Table S3).",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/analyze_eight_pair_dump_gated_v1.py",
        "ChEMBL 37 dump join, max-versus-median, and leftover counts for all eight Table 2 pairs (Table S3).",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py",
        "ChEMBL pair-census summary plotted in Figure 1.",
    ),
    (
        "data/jcim_novelty_v0/scripts/ecfp4_docking_scaler_sensitivity_v1.py",
        "Table S5 StandardScaler sensitivity on the same GroupKFold splits.",
    ),
    (
        "data/jcim_novelty_v0/scripts/reaudit_all14_cognate_rmsd_calcrrms_v1.py",
        "Unified 14-receptor chemically mapped CalcRMS for Table S2 / Figure S4. Does not redock.",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/reaudit_layer3_cognate_rmsd_v1.py",
        "Atom-mapping / CalcRMS helper imported by the all-14 cognate RMSD script.",
    ),
    (
        "data/jcim_chembl_universe_v0/scripts/replay_track_b_vina_mode1_v1.py",
        "Re-reads REMARK VINA RESULT from the committed production pose tree. Does not redock.",
    ),
    (
        "figures/jcim_article/scripts/update_figures_pr32.py",
        "Official generator for Figures 1–6, S1–S4, and the TOC graphic from pinned CSVs.",
    ),
    (
        "figures/jcim_article/scripts/plot_jcim_article_figures_v3.py",
        "Panel functions imported by the official figure generator.",
    ),
    (
        "figures/jcim_article/scripts/jcim_figure_style.py",
        "Shared figure style used by the official figure generator.",
    ),
]


def copy_file(rel: str, dest_root: Path) -> None:
    src = ROOT / rel
    dst = dest_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    if PACK.exists():
        shutil.rmtree(PACK)
    (PACK / "manuscript").mkdir(parents=True)
    (PACK / "tables").mkdir()
    (PACK / "scripts").mkdir()
    (PACK / "figures").mkdir()

    for name in DOCS:
        src = ROOT / "docs" / name
        if src.exists():
            shutil.copy2(src, PACK / "manuscript" / name)
    for rel in TABLES:
        src = ROOT / rel
        dst = PACK / "tables" / Path(rel).name
        if src.exists():
            shutil.copy2(src, dst)
    script_lines = [
        "# Experiment scripts packed for upload",
        "",
        "These scripts reproduce typeset tables and figures from deposited scores or saved poses.",
        "Run them from the `Dual_Target_Docking/` tree. This folder is a flat copy for upload.",
        "Do not upload assemble / audit / validate / checksum / freeze / pack helpers.",
        "",
        "| Script | Repository path | Paper role |",
        "|---|---|---|",
    ]
    for rel, role in SCRIPTS:
        src = ROOT / rel
        if src.exists():
            shutil.copy2(src, PACK / "scripts" / Path(rel).name)
            script_lines.append(f"| `{Path(rel).name}` | `{rel}` | {role} |")
    (PACK / "SCRIPTS.md").write_text("\n".join(script_lines) + "\n", encoding="utf-8")

    fig_src = ROOT / "figures" / "jcim_article"
    allowed_suffix = {".png", ".pdf", ".tif", ".tiff", ".md", ".json"}
    official_stems = {
        "Fig1_four_state_and_supply",
        "Fig1_C_chEMBL_supply",
        "Fig2_negative_class_formulation",
        "Fig3_ligand_chemistry",
        "Fig4_mismatched_pocket",
        "Fig5_computational_realization",
        "Fig6_evidence_boundary",
        "TOC_graphic",
        "FigS1_posthoc_diagnostics",
        "FigS2_pocket_matched_forest",
        "FigS3_protocol_sensitivity",
        "FigS4_cognate_rmsd",
        "MANUSCRIPT_FIGURE_CAPTIONS",
        "ARCHIVED_FIGURES",
        "plotted_values",
        "FIGURE_AUDIT_PR32",
        "FIGURE_SHA256",
    }
    if fig_src.exists():
        for path in sorted(fig_src.iterdir()):
            if path.suffix.lower() not in allowed_suffix:
                continue
            if path.stem in official_stems:
                shutil.copy2(path, PACK / "figures" / path.name)

    readme = """# JCIM submission pack

This directory is the publication-facing slice of DualFourClass-Bench.
It is not a second copy of the docking pose workspaces.

## Use these files for submission

| Path | Role |
|---|---|
| `manuscript/MANUSCRIPT_JCIM_EN.md` | English manuscript (assemble from section drafts; do not hand-edit the assembled file as the source of truth) |
| `manuscript/MANUSCRIPT_JCIM_ZH.md` | Chinese working manuscript |
| `manuscript/SUPPORTING_INFORMATION_JCIM_EN_V1.md` | English SI Tables S1–S9 |
| `manuscript/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` | Chinese SI |
| `manuscript/STATISTICAL_LOCK_V1.md` | Table 2 / Table 3 estimand lock |
| `manuscript/FIGURE_PANEL_LOCK_V3.md` | Figure-to-CSV map |
| `manuscript/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` | Five-round numeric audit |
| `tables/` | Frozen CSVs cited by Tables 1–3 and S1–S9 |
| `figures/` | Regenerated main and SI figures from `figures/jcim_article/scripts/update_figures_pr32.py` |
| `scripts/` | Experiment analysis and official figure generators only. See `SCRIPTS.md` |
| `SCRIPTS.md` | Upload list: each script, repository path, and paper role |

## Do not submit as primary evidence

- Assemble / audit / validate / checksum / freeze / pack helpers (`assemble_manuscript_*.py`, `audit_*.py`, `validate_revision_v1.py`, `build_checksum_manifest_v1.py`, `freeze_submission_v1.py`, `pack_submission_v1.py`, `bootstrap_primary.py`, `build_master_results_table_v1.py`, `review_statistics_sensitivity_v1.py`)
- `pocket_unidirectional_delta_v1.py` (S6b removed from the typeset SI)
- Docking campaign runners (the paper reproduces statistics from deposited score tables)
- `plot_jcim_article_figures_v1.py` / `v2.py` (withdrawn-pair leftovers)
- `plot_jcim_si_composites_v1.py` (historical original-set artwork; may still tick PIK3CA/PIK3CB)
- `data/pik3ca_pik3cb_panel_v0/` (withdrawn pair archive)
- `external_slice_summary_202608_contract_v1.csv` (legacy three-pair BindingDB snapshot)
- Production poses (`local_track_b_v0/poses/`; in git, not copied into this pack slice)
- A minted Zenodo DOI (not created)

## Rebuild / submission freeze

Author-side only. Do not upload the freeze runner.

```bash
python3 scripts/audit/freeze_submission_v1.py
```
"""
    (PACK / "README.md").write_text(readme, encoding="utf-8")
    inventory = ["# Inventory", "", "## root", ""]
    for path in sorted(PACK.iterdir()):
        if path.is_file():
            inventory.append(f"- `{path.name}` ({path.stat().st_size} bytes)")
    inventory.append("")
    for folder in ("manuscript", "tables", "figures", "scripts"):
        inventory.append(f"## {folder}")
        for path in sorted((PACK / folder).iterdir()):
            inventory.append(f"- `{path.name}` ({path.stat().st_size} bytes)")
        inventory.append("")
    (PACK / "INVENTORY.md").write_text("\n".join(inventory), encoding="utf-8")
    print("wrote", PACK)


if __name__ == "__main__":
    main()
