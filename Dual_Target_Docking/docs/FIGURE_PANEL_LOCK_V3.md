# Eight-row panel lock (submission)

Branch: `cursor/jcim-submission-pack-0b1a`  
Artwork generator: `figures/jcim_article/scripts/update_figures_pr32.py`  
Legacy lock script: `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`  
Rule: every plotted number is read from the CSV in this table. No hand-typed AUROCs. No AI-drawn figures. No decorative arrows unrelated to the data. The primary set is the eight pairs in Table 1, shown in protein-system order.

Display order (Figures 2–6 and Tables 1–3): EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD. CSV routing may still use historical file splits; those splits are not drawn on the figures.

| Figure | Panel | Content | Unique source |
|---|---|---|---|
| 1 | A | Four ligand states (schematic) | none |
| 1 | B | Pocket-matched directional tasks (schematic) | none |
| 1 | C | ChEMBL universe census: 2,164,618 / 63,790 / 5,253 / 86 pairs, then 8 primary after panel construction | `universe_census_summary_v1.csv` |
| 2 | A | Fixed-score ΔAUROC (neither − selective), both pockets, eight pairs | `formulation_equal_score_negative_v1.csv`; `equal_score_negative_s34_v1.csv` |
| 2 | B | Directional D/A and D/B AUROC, eight primary pairs | `unified_threshold_sensitivity_v2.csv` θ=6.0; `five_pair_stack_v1/table2_comparable_theta6_v1.csv` |
| 2 | C | Dual-vs-neither (`vina_mean`) vs directional `summary_min` | formulation CSV; same table2 file |
| 3 | A | ECFP4 GroupKFold vs Vina rank AUROC, both arms, eight pairs | `ligand_ml_baseline_scaffold_cv_v1.csv`; `ecfp4_incremental_s20s24_v1.csv` |
| 3 | B | ECFP4 → ECFP4+docking ΔAUROC, 16 contrasts, axis ±0.03 | `incremental_information_v1.csv`; same ECFP4 file |
| 3 | C | Illustrative AChE/BChE TPSA jitter + median/IQR | `assembled_AChE_BChE.csv` |
| 4 | A | Independent GNINA pose generation vs Vina: EGFR/HER2, JAK1/TYK2, PIK3CA/mTOR | `independent_dock_formulation_v1.csv`; `table2_comparable_by_channel_v1.csv` `gnina_independent_jak1_tyk2` |
| 4 | B | PIK3CA 4L23/4JPS/5DXT and mTOR 4JSX on PIK3CA/mTOR only | alt receptor CSVs + unified_threshold |
| 4 | C | Five-seed `summary_min` range, eight pairs (8 × 5) | `multiseed_auroc_by_seed_v2.csv`; `fiveseed_summary_min_aggregate_v1.csv` |
| 5 | A | Main and holdout matched−mismatched Δ on the same rows | `wrong_pocket_paired_delta_bootstrap_v1.csv`; `wrong_pocket_by_channel_v1.csv` |
| 5 | B | Holdout vs main `summary_min` | `holdout_pocket_matched_v1.csv`; `table2_comparable_by_channel_v1.csv` `holdout_vina_20260727` |
| 6 | A | θ-grid `summary_min`, eight pairs | `unified_threshold_sensitivity_v2.csv`; `threshold_grid_v1.csv` |
| 6 | B | EGFR/HER2 and JAK1/TYK2 fixed-score Δ under ligand / scaffold / document resampling | `equal_score_cluster_bootstrap_v1.csv`; ligand-level rows as Figure 2A |
| 6 | C | BindingDB compound counts; color saturates at n=20; 0/8 meet all external criteria | `external_slice_summary_v1.csv` |
| 6 | D | BindingDB source counts; color saturates at n=3 | same BindingDB file |
| S1 | A–B | PIK3CA/mTOR PM48 vs PM110 and exhaustiveness 16 vs 8 | `pm110_vs_pm48_pocket_matched_v1.csv`; `scores_vina_E8_best.csv` |
| S4 | — | Eight-row Vina forest + best single descriptor | same Table-2 sources + `descriptor_all_four_directional_v1.csv` |
| S5 | — | Unused-pool holdout vs main, seven pairs | same holdout sources as Fig 5B |
| S6 | — | Detectable-effect simulation, three available pairs, true AUROC 0.50–0.75 | `detectable_effect_simulation_v1.csv` |
| S7 | — | EGFR/HER2 Top-10 (1/5/4/0) and dual-median AND filter (14/9/24) | `mixed_library_enrichment_v1.csv`; `and_filter_operating_point_v1.csv` |
| S8 | — | BindingDB compound and source counts; eight primary pairs; 0 primary pass | `external_slice_summary_v1.csv` |
| S11 | — | Same cluster intervals as Figure 6B (SI copy) | `equal_score_cluster_bootstrap_v1.csv` |
| S12 | — | Cognate redocking RMSD, 14 primary receptors, top-1 / top-3 / best-of-9 | `cognate_rank_rmsd_reaudit_v1.csv`; `layer3_cognate_rmsd_v1.csv`; `pm48_01_rmsd_E16.csv` |
| TOC | — | Four states → directional and Dual–neither tasks → ligand-only and pocket-correspondence controls | schematic; no AUROCs |

S2, S3, S9, and S10 remain SI **table** records. No artwork for those historical original-set panels is shipped. Table numbers S1–S13 are unchanged.
