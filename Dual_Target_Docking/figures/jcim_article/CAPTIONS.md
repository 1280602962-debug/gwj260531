# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate: `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimentally labeled ligand states: dual, A-only, B-only, and neither. A-only and B-only are selectivity hard negatives. (B) Primary tasks are pocket-matched directional AUROCs: Dual versus A-only scored in pocket B, Dual versus B-only scored in pocket A. `summary_min` is a descriptive worst-arm summary. (C) Left: J0 ChEMBL audit of 49 candidate pairs (`j0_strict_label_supply.csv`). Four pairs meet a thick hard-negative gate (min ≥50). HDAC1/HDAC6 is excluded as metal-dependent. EGFR/HER2 is retained as a supply-limited case (strict B-only = 7), giving four historically docked pairs. PIK3CA/PIK3CB was later withdrawn after a receptor-identity failure. A later ChEMBL 37 dump census added five ordinary pairs, leaving eight primary rows. This is collection completion, not a claim that eight pairs were in the J0 thick set. Right: min strict hard-negatives on the J0 scrape (left group) and on the later ChEMBL 37 dump for the five census pairs (right group; `crossdb_strict_supply_v1.csv`). Complete-case map coverage on the J0-docked pairs is 14.5%–34.0%. Cross-database counts for the original scrape are Figure S2.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Same frozen AutoDock Vina scores, unified θ = 6.0, eight primary rows. PIK3CA/PIK3CB is withdrawn and is not plotted. Horizontal gray rules separate the 2026-07-23 three from the five post-census pairs. (A) Directional Dual versus A-only (pocket B) and Dual versus B-only (pocket A). Original three: `unified_threshold_sensitivity_v2.csv`. Five: `table2_comparable_theta6_v1.csv`. (B) Descriptive comparison of directional `summary_min` with Dual versus neither using per-ligand `vina_mean`. These two columns differ in both negative class and score aggregation. PIK3CA/mTOR Dual versus neither is hatched (neither n = 4). (C) Pocket A score held fixed; only the negative class is replaced (B-only versus neither). EGFR/HER2 ΔAUROC = 0.378 [0.205, 0.547]; JAK1/TYK2 reproduces the gap (0.444 [0.263, 0.620]). Diamond, underpowered neither. Vertical dashed line, zero.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Scaffold GroupKFold ECFP4 logistic AUROC versus pocket-matched Vina rank AUROC on both directional arms for the eight primary pairs. EGFR/HER2 Dual versus B-only: ECFP4 0.8895 versus Vina 0.4297. Five-pair ECFP4 from `ecfp4_incremental_s20s24_v1.csv`. (B) Change in GroupKFold AUROC when the pocket-matched Vina score is added to ECFP4 (16 contrasts). (C) AChE/BChE TPSA by class: individual ligands (jittered) with median and IQR (`assembled_AChE_BChE.csv`). n = 27/25/28.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation (not CNN rescoring of Vina poses) on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2 only. Independent search was not run on F2/F10, JAK1/JAK2, or the PPAR pairs. JAK1/TYK2 independent GNINA: `summary_min` 0.317 [0.183, 0.463], Dual versus neither 0.705. (B) Replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR frozen: PIK3CA/mTOR `summary_min` 0.692 → 0.486 / 0.505. 4JSX is an mTOR-pocket swap. The parallel swap on the withdrawn PIK3CA/PIK3CB pair is not plotted as a peer contrast. (C) Directional `summary_min` across five frozen Vina seeds; diamond, production seed 20260727. No post-census five-seed range crossed 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Δ = matched-pocket `summary_min` − mismatched-pocket `summary_min`, ligand bootstrap B = 2000. Matched uses Dual versus A-only in pocket B and Dual versus B-only in pocket A; mismatched swaps those score channels. This is a scoring-channel control, not redocking into a physically wrong site. Dark, CI excludes 0; gray, CI includes 0. (A) Main panels, eight primary pairs. EGFR/HER2 and AChE/BChE CIs exclude 0; the five post-census production-Vina CIs include 0. (B) Unused-pool holdout Δ. All seven CIs include 0. EGFR/HER2 has no holdout; withdrawn PIK3CA/PIK3CB is omitted. (C) Holdout versus main-panel `summary_min`. PPARG/PPARA holdout is 0.535 [0.350, 0.717]; JAK1/JAK2 stays same-direction (0.619 [0.420, 0.749]; drawn 20/20/18).

## Figure 6. Robustness checks and evidence boundary.

(A) Pocket-matched `summary_min` on the unified label-threshold grid for the eight primary pairs. Solid, 2026-07-23 three; dashed, post-census five. (B) PIK3CA/mTOR PM48 versus PM110 Vina. (C) PM48 exhaustiveness 16 versus 8, recomputed from `scores_vina_E8_best.csv` with the same pocket-matched definition. (D) BindingDB-native 202608 slice under the original four-pair contract: zero pairs meet the pre-frozen external gate; nothing was docked. The five census pairs were not re-opened as a BindingDB external set.

## Figure S1. Protocol and panel sensitivities.

Original-set protocol grid, GNINA CNN rescoring of Vina poses, PM48 versus PM110, and exhaustiveness. This SI record still includes the later-withdrawn PIK3CA/PIK3CB row where that is what the source CSVs contain. Independent GNINA pose generation is Figure 4A.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C (right group).

## Figure S3. Additional paired bootstrap differences.

Descriptor and scaffold-versus-random leakage checks on the original docked set. Matched-versus-mismatched main/holdout Δ CIs are Figure 5.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows. PIK3CA/PIK3CB is omitted.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Mismatched-pocket Δ CIs are Figure 5B.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 J0 pair census (`docked_in_this_paper` = 4 includes the later-withdrawn PIK3CA/PIK3CB row), current primary n = 8, AND-like dual filter on the original three, and ligand-only full-map ECFP4 on the original three. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

Filter cascade and remaining four-state counts after literature, structure, and ECFP4 < 0.70 on the original four-pair contract (`external_slice_summary_v1.csv`). Zero of four pairs meet the pre-frozen external gate; nothing was docked. Five census pairs were not re-opened as BindingDB external.

## Figure S9. Additional ligand-structure controls.

Prespecified descriptors, covariate-adjusted logistic AUROC, and matched-subset weak-arm tests on the original docked set. The Vina-only logistic AUROC is not the Table 2 rank AUROC.

## Figure S10. Matched versus mismatched point estimates.

Bar charts of matched versus mismatched `summary_min` on the original docked set. Paired Δ CIs are Figure 5.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs and no decorative arrows.
