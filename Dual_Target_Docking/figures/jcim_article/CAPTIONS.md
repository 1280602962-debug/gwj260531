# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V2.md`.
Regenerate: `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v2.py`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimentally labeled ligand states: dual, A-only, B-only, and neither. A-only and B-only are selectivity hard negatives. (B) Primary tasks are pocket-matched directional AUROCs: Dual versus A-only scored in pocket B, Dual versus B-only scored in pocket A. `summary_min` is a descriptive worst-arm summary. (C) J0 ChEMBL audit of 49 candidate pairs (`j0_strict_label_supply.csv`). Four pairs meet a thick hard-negative gate (min ≥50); HDAC1/HDAC6 is excluded as metal-dependent; EGFR/HER2 is retained as a supply-limited K=4 case (strict B-only = 7). Complete-case map coverage is 14.5%–34.0% (`complete_case_usable_pchembl_overlap_v1.csv`). Cross-database counts are Figure S2.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Same frozen AutoDock Vina scores, unified θ = 6.0. (A) Directional Dual versus A-only (pocket B) and Dual versus B-only (pocket A) from `unified_threshold_sensitivity_v2.csv`. (B) Descriptive comparison of directional `summary_min` (Table 2 CIs from the same unified-threshold file) with Dual versus neither using per-ligand `vina_mean` (`formulation_conventional_vs_directional_v1.csv`). These two columns differ in both negative class and score aggregation; the difference is not a paired test of one estimand. PIK3CA/mTOR Dual versus neither is hatched (neither n = 4). (C) Pocket A score held fixed; only the negative class is replaced (B-only versus neither). EGFR/HER2 ΔAUROC = 0.378 [0.205, 0.547] (`formulation_equal_score_negative_v1.csv`). Vertical dashed line in (C), zero.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Scaffold GroupKFold ECFP4 logistic AUROC versus pocket-matched Vina rank AUROC on both directional arms (`ligand_ml_baseline_scaffold_cv_v1.csv`). EGFR/HER2 Dual versus B-only: ECFP4 0.8895 versus Vina 0.4297. (B) Change in GroupKFold AUROC when the pocket-matched Vina score is added to ECFP4 (`incremental_information_v1.csv`). The largest absolute change among the eight contrasts is ≤0.020. (C) AChE/BChE TPSA by class: individual ligands (jittered) with median and IQR (`assembled_AChE_BChE.csv`). n = 27/25/28.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation (not CNN rescoring of Vina poses) on EGFR/HER2 and PIK3CA/mTOR (`independent_dock_formulation_v1.csv`) versus the same-panel Vina values. (B) Replacing PIK3CA 4L23 with 4JPS or 5DXT while holding the second pocket frozen: PIK3CA/mTOR `summary_min` 0.692 → 0.486 / 0.505; PIK3CA/PIK3CB 0.500 → 0.691 / 0.685. 4JSX is an mTOR-pocket swap and is plotted with a distinct marker; it is not applied to PIK3CA/PIK3CB. CIs from the deposited swap tables and Table 2. (C) Directional `summary_min` across five frozen Vina seeds (`multiseed_auroc_by_seed_v2.csv`); diamond, production seed 20260727.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Δ = matched-pocket `summary_min` − mismatched-pocket `summary_min`, ligand bootstrap B = 2000 (`wrong_pocket_paired_delta_bootstrap_v1.csv`). Matched uses Dual versus A-only in pocket B and Dual versus B-only in pocket A; mismatched swaps those score channels. This is a scoring-channel control, not redocking into a physically wrong site. (A) Main panels. EGFR/HER2 and AChE/BChE CIs exclude 0; PIK3CA pairs include 0. (B) Unused-pool holdout. All three CIs include 0; point estimates are negative (mismatched ≥ matched). EGFR/HER2 has no holdout. (C) Holdout potency (|ΔpChEMBL| ≤ 0.5) and size (|Δheavy| ≤ 2) matching: point Δ only (`holdout_matched_wrong_pocket_summary_v1.csv`). Dark, CI excludes 0; gray, CI includes 0.

## Figure 6. Robustness checks and evidence boundary.

(A) Pocket-matched `summary_min` on the unified label-threshold grid (`unified_threshold_sensitivity_v2.csv`). (B) PIK3CA/mTOR PM48 versus PM110 Vina (`pm110_vs_pm48_pocket_matched_v1.csv`). (C) PM48 exhaustiveness 16 versus 8, recomputed from `scores_vina_E8_best.csv` with the same pocket-matched definition. (D) BindingDB-native 202608 slice: zero of four pairs meet the pre-frozen external gate; nothing was docked (`external_slice_summary_v1.csv`).

## Figure S1. Protocol and panel sensitivities.

(A) Unified label-threshold grid. Open markers, underpowered cells. (B) GNINA CNN rescoring of Vina poses (mode-1 versus best-of-9) versus primary Vina; this is not independent GNINA pose generation (Figure 4A / Table S32). (C) PM48 versus PM110 for Vina, RTMScore, and GNINA CNN best-of-9 rescoring. (D) Exhaustiveness 16 versus 8 and single-target enrichment on 4L23/4JT6.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`.

## Figure S3. Additional paired bootstrap differences.

Descriptor and scaffold-versus-random leakage checks. Matched-versus-mismatched main/holdout Δ CIs are now Figure 5.

## Figure S4. Pocket-matched summary_min forest (former main figure).

Vina CIs from `unified_threshold_sensitivity_v2.csv`. GNINA in this figure is CNN rescoring of Vina poses (best-of-9), not independent pose generation.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` only; mismatched-pocket Δ CIs are Figure 5B.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 pair census, AND-like dual filter, and full-map ligand-only ECFP4. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

Filter cascade and remaining four-state counts after literature, structure, and ECFP4 < 0.70 (`external_slice_summary_v1.csv`). Zero of four pairs meet the pre-frozen external gate; nothing was docked.

## Figure S9. Additional ligand-structure controls.

Prespecified descriptors, covariate-adjusted logistic AUROC, and matched-subset weak-arm tests (former main Figure 7B–D). The Vina-only logistic AUROC is not the Table 2 rank AUROC.

## Figure S10. Matched versus mismatched point estimates.

Bar charts of matched versus mismatched `summary_min` on the main panel, unused-pool holdout, potency/size matching, and scoring-free contact counts. Paired Δ CIs are Figure 5. Contact count is exploratory and does not explain holdout reversal.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs.
