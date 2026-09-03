# Six-figure panel lock (submission)

Branch: `cursor/jcim-final-integration-0b1a`  
Script: `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v2.py`  
Rule: every plotted number is read from the CSV in this table. No hand-typed AUROCs. No AI-drawn figures.

| Figure | Panel | Content | Unique source |
|---|---|---|---|
| 1 | A | Four ligand states (schematic) | none |
| 1 | B | Pocket-matched directional tasks (schematic) | none |
| 1 | C | 49-pair strict hard-neg funnel → K=4; complete-case fractions | `j0_strict_label_supply.csv`; `complete_case_usable_pchembl_overlap_v1.csv` |
| 2 | A | Directional D/A and D/B AUROC | `unified_threshold_sensitivity_v2.csv` θ=6.0 |
| 2 | B | Dual-vs-neither (`vina_mean`) vs directional `summary_min` (descriptive) | formulation CSV + unified_threshold CIs |
| 2 | C | Fixed-pocket score, negative-class ΔAUROC | `formulation_equal_score_negative_v1.csv` |
| 3 | A | ECFP4 GroupKFold vs Vina rank AUROC | `ligand_ml_baseline_scaffold_cv_v1.csv` |
| 3 | B | ECFP4 → ECFP4+docking ΔAUROC | `incremental_information_v1.csv` |
| 3 | C | AChE/BChE TPSA jitter + box | `assembled_AChE_BChE.csv` |
| 4 | A | Independent GNINA pose generation vs Vina | `independent_dock_formulation_v1.csv` + unified_threshold + formulation |
| 4 | B | Aligned PIK3CA 4L23/4JPS/5DXT; 4JSX distinct | alt receptor CSVs + unified_threshold |
| 4 | C | Five-seed `summary_min` range | `multiseed_auroc_by_seed_v2.csv` |
| 5 | A | Main matched−mismatched Δ + 95% CI | `wrong_pocket_paired_delta_bootstrap_v1.csv` `set=main_panel` |
| 5 | B | Holdout matched−mismatched Δ + 95% CI | same CSV `set=unused_pool_holdout` |
| 5 | C | Holdout potency/size matching Δ (point) | `holdout_matched_wrong_pocket_summary_v1.csv` |
| 6 | A | θ-grid `summary_min` | `unified_threshold_sensitivity_v2.csv` |
| 6 | B | PM48 vs PM110 Vina | `pm110_vs_pm48_pocket_matched_v1.csv` |
| 6 | C | PM48 E=16 vs E=8 | E=16 from unified_threshold; E=8 from `scores_vina_E8_best.csv` via SI loader |
| 6 | D | BindingDB-native gate | `external_slice_summary_v1.csv` |
| TOC | — | Four states → directional evaluation → Dual-vs-neither ≠ Dual-vs-selective | schematic; no AUROCs |

Demoted (not main figures):

- Former Fig 2B cross-database counts → Figure S2
- Former Fig 6 bar charts and contact-count → not redrawn as main; paired-Δ is Figure 5
- Former Fig 7 B/C/D → Figure S9 (`FigS_ligand_controls`)
- Former Fig 8 workflow → TOC
- `PRIMARY_METRIC_V2.md` is never a figure source
