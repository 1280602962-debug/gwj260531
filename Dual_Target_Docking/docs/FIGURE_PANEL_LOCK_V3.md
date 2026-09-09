# Eight-row panel lock (submission)

Branch: `cursor/chembl-exhaustive-pair-census-0b1a`  
Script: `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`  
Rule: every plotted number is read from the CSV in this table. No hand-typed AUROCs. No AI-drawn figures. No decorative arrows unrelated to the data. The primary set is the eight pairs in Table 1.

| Figure | Panel | Content | Unique source |
|---|---|---|---|
| 1 | A | Four ligand states (schematic) | none |
| 1 | B | Pocket-matched directional tasks (schematic; no arrows) | none |
| 1 | C | J0 scrape → min selective ≥50 → eight-pair evaluation set; J0 selective-supply bars (HDAC + PIK3CA/mTOR + AChE/BChE + EGFR); later ChEMBL 37 dump selectives for the five census pairs | `j0_strict_label_supply.csv`; `five_pair_crossdb_v1/crossdb_strict_supply_v1.csv` ChEMBL37_dump; `complete_case_usable_pchembl_overlap_v1.csv` (original three maps only) |
| 2 | A | Directional D/A and D/B AUROC, eight primary pairs | original three: `unified_threshold_sensitivity_v2.csv` θ=6.0; five: `five_pair_stack_v1/table2_comparable_theta6_v1.csv` |
| 2 | B | Dual-vs-neither (`vina_mean`) vs directional `summary_min` | original three: formulation CSV; five: same table2 file |
| 2 | C | Fixed pocket-A score, negative-class ΔAUROC | original three: `formulation_equal_score_negative_v1.csv`; five: `equal_score_negative_s34_v1.csv` |
| 3 | A | ECFP4 GroupKFold vs Vina rank AUROC, both arms, eight pairs | original three: `ligand_ml_baseline_scaffold_cv_v1.csv`; five: `ecfp4_incremental_s20s24_v1.csv` |
| 3 | B | ECFP4 → ECFP4+docking ΔAUROC, 16 contrasts | original three: `incremental_information_v1.csv`; five: same ECFP4 file |
| 3 | C | AChE/BChE TPSA jitter + median/IQR | `assembled_AChE_BChE.csv` |
| 4 | A | Independent GNINA pose generation vs Vina: EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2 only | `independent_dock_formulation_v1.csv`; `table2_comparable_by_channel_v1.csv` `gnina_independent_jak1_tyk2` |
| 4 | B | PIK3CA 4L23/4JPS/5DXT and mTOR 4JSX on the receptor-verified PIK3CA/mTOR pair only | alt receptor CSVs + unified_threshold |
| 4 | C | Five-seed `summary_min` range, eight pairs | original three: `multiseed_auroc_by_seed_v2.csv`; five: `fiveseed_summary_min_aggregate_v1.csv` |
| 5 | A | Main matched−mismatched Δ + 95% CI, eight pairs | original three: `wrong_pocket_paired_delta_bootstrap_v1.csv` `set=main_panel`; five: `wrong_pocket_by_channel_v1.csv` `vina_20260727` |
| 5 | B | Holdout matched−mismatched Δ + 95% CI (no EGFR; no withdrawn PIK3CB) | original two: same bootstrap CSV `unused_pool_holdout`; five: `wrong_pocket_by_channel_v1.csv` `holdout_vina_20260727` |
| 5 | C | Holdout vs main `summary_min` | original two: `holdout_pocket_matched_v1.csv`; five: `table2_comparable_by_channel_v1.csv` `holdout_vina_20260727` |
| 6 | A | θ-grid `summary_min`, eight pairs | original three: `unified_threshold_sensitivity_v2.csv`; five: `threshold_grid_v1.csv` |
| 6 | B | PM48 vs PM110 Vina | `pm110_vs_pm48_pocket_matched_v1.csv` |
| 6 | C | PM48 E=16 vs E=8 | E=16 from unified_threshold; E=8 from `scores_vina_E8_best.csv` |
| 6 | D | BindingDB-native gate on all eight primary pairs (0 primary pass). | `external_slice_summary_v1.csv` |
| S4 | — | Eight-row Vina forest + best single descriptor | same Table-2 sources + `descriptor_all_four_directional_v1.csv` |
| S5 | — | Unused-pool holdout vs main, seven pairs | same holdout sources as Fig 5C |
| S7 | — | J0 θ=6.0 candidate-pair census plus current primary n=8 | `theta6_pair_census_v1.csv` |
| S8 | — | BindingDB-native cascade; eight primary pairs; 0 primary pass | `external_slice_summary_v1.csv` |
| TOC | — | Four states and Dual-vs-neither ≠ Dual-vs-selective | schematic; no AUROCs; no arrows |

S1–S3, S9, S10 are SI records for the pairs that have those sensitivity tables. They are not a second primary set.
