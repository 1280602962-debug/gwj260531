# Figure lock and file mapping (submission)

Branch: `cursor/jcim-pack-consistency-0b1a`  
Artwork generator: `figures/jcim_article/scripts/update_figures_pr32.py`  
Legacy lock script: `figures/jcim_article/scripts/plot_jcim_article_figures_v3.py` (do not use `data/jcim_bench_v0` or `submission_pack/scripts/plot_jcim_article_figures_v3.py` for current numbering)  
Rule: every plotted number is read from the CSV in this table. No hand-typed AUROCs. No AI-drawn figures. The primary set is the eight pairs in Table 1, shown in protein-system order.

Display order (Figures 2–6 and Tables 1–3): EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD.

Article numbering follows the Chinese working manuscript. Disk filenames need not match the final figure number. Figure S5 may cite `FigS6_detectable_effect`.

## Final figure → title → file → plotter → source

| Final no. | Title | Physical stem | Plot function | Data source |
|---|---|---|---|---|
| Figure 1 | Four-state evaluation and ChEMBL supply | `Fig1_four_state_and_supply` | `fig1_design()` | `universe_census_summary_v1.csv` |
| Figure 2 | Docking depends on the experimental-state comparison | `Fig2_negative_class_formulation` | `fig2_negative()` | `formulation_equal_score_negative_v1.csv`; `equal_score_negative_s34_v1.csv`; `unified_threshold_sensitivity_v2.csv`; `table2_comparable_theta6_v1.csv` |
| Figure 3 | Ligand chemistry as a competing explanation | `Fig3_ligand_chemistry` | `fig3_chemistry()` | `ligand_ml_baseline_scaffold_cv_v1.csv`; `ecfp4_incremental_s20s24_v1.csv`; `incremental_information_v1.csv`; `assembled_AChE_BChE.csv` |
| Figure 4 | Matched versus mismatched pocket | `Fig4_mismatched_pocket` | `fig5_pocket()` | `wrong_pocket_paired_delta_bootstrap_v1.csv`; `wrong_pocket_by_channel_v1.csv`; `holdout_pocket_matched_v1.csv` |
| Figure 5 | Computational realization | `Fig5_computational_realization` | `fig4_compute()` | `independent_dock_formulation_v1.csv`; alt-receptor CSVs; `multiseed_auroc_by_seed_v2.csv`; `fiveseed_summary_min_aggregate_v1.csv` |
| Figure 6 | Evidence boundary | `Fig6_evidence_boundary` | `fig6_boundary()` | `unified_threshold_sensitivity_v2.csv`; `equal_score_cluster_bootstrap_v1.csv`; `external_slice_summary_v1.csv` |
| Figure S1 | EGFR/HER2 Top-10 and AND filter | `FigS1_posthoc_diagnostics` | `fig_s7_posthoc()` | `mixed_library_enrichment_v1.csv`; `and_filter_operating_point_v1.csv` |
| Figure S2 | Pocket-matched forest | `FigS2_pocket_matched_forest` | `fig_s4_forest()` | Table-2 sources; `descriptor_all_four_directional_v1.csv` |
| Figure S3 | PIK3CA/mTOR protocol sensitivity | `FigS3_protocol_sensitivity` | `fig_s1_protocol()` | `pm110_vs_pm48_pocket_matched_v1.csv`; `scores_vina_E8_best.csv` |
| Figure S4 | Cognate redocking RMSD | `FigS4_cognate_rmsd` | `fig_s12_cognate()` | Original six: `cognate_rank_rmsd_reaudit_v1.csv`, `pm48_01_rmsd_E16.csv`. Added eight: `layer3_cognate_rmsd_calcrrms_v1.csv`. Historical only: `layer3_cognate_rmsd_v1.csv` |
| Figure S5 | Detectable-effect scenario | `FigS6_detectable_effect` | `fig_s6_sim()` | `detectable_effect_simulation_v1.csv` |
| TOC | Four states → tasks → controls | `TOC_graphic` | `toc_graphic()` | schematic |

## Archived artwork (kept in repo, not packed, not typeset)

These files duplicate a main-text panel. They stay in `figures/jcim_article/` with `ARCHIVED_FIGURES.md`. Packer must not copy them.

| Stem | Why archived |
|---|---|
| `FigS5_unused_pool_holdout` | Same content as Figure 4B |
| `FigS7_bindingdb_native_slice` | Same content as Figure 6C,D |
| `FigS8_cluster_uncertainty` | Same content as Figure 6B |

S2, S3, S9, and S10 remain SI **table** records. Table numbers S1–S13 are unchanged.
