# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate: `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimentally labeled ligand states: dual, A-only, B-only, and neither. A-only and B-only are single-target-selective controls. (B) Primary tasks are pocket-matched directional AUROCs: Dual versus A-only scored in pocket B, Dual versus B-only scored in pocket A. `summary_min` is a descriptive worst-arm summary. (C) Bidirectional selective supply is scarce under a strict 6.5/5.5 rule. HDAC1/HDAC6 is excluded as metal-dependent. EGFR/HER2 is retained as a supply-limited construction case (strict B-only = 7; primary θ = 6.0 n = 28/38/32). The eight primary pairs use one extract protocol with pair-specific candidate pools, quotas, and scaffold caps (Table 1). Right: J0 scrape selectives (HDAC + PIK3CA/mTOR + AChE/BChE + EGFR) and later ChEMBL 37 dump selectives for the five census pairs. Cross-database counts are Figure S2.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Same frozen AutoDock Vina scores, unified θ = 6.0, eight primary rows. (A) Directional Dual versus A-only (pocket B) and Dual versus B-only (pocket A). (B) Descriptive comparison of directional `summary_min` with Dual versus neither using per-ligand `vina_mean`. These two columns differ in both negative class and score aggregation. PIK3CA/mTOR Dual versus neither is hatched (neither n = 4). (C) Pocket A score held fixed; only the negative class is replaced (B-only versus neither). EGFR/HER2 ΔAUROC = 0.378 [0.205, 0.547]; JAK1/TYK2 reproduces the gap (0.444 [0.263, 0.620]). Document-cluster and scaffold-cluster intervals for those two Δ values are in Table S4. Diamond, underpowered neither. Vertical dashed line, zero.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Scaffold GroupKFold ECFP4 logistic AUROC versus pocket-matched Vina rank AUROC on both directional arms for the eight primary pairs. EGFR/HER2 Dual versus B-only: ECFP4 0.8895 versus Vina 0.4297. Five-pair ECFP4 from `ecfp4_incremental_s20s24_v1.csv`. (B) Change in GroupKFold AUROC when the pocket-matched Vina score is added to ECFP4 (16 contrasts). (C) AChE/BChE TPSA by class: individual ligands (jittered) with median and IQR (`assembled_AChE_BChE.csv`). n = 27/25/28.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation (not CNN rescoring of Vina poses) on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2 only. Independent search was not run on F2/F10, JAK1/JAK2, or the PPAR pairs. JAK1/TYK2 independent GNINA: `summary_min` 0.317 [0.183, 0.463], Dual versus neither 0.705. (B) Replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR frozen: PIK3CA/mTOR `summary_min` 0.692 → 0.486 / 0.505. 4JSX is an mTOR-pocket swap. (C) Directional `summary_min` across five frozen Vina seeds; diamond, production seed 20260727. No five-seed range on F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, or PPARA/PPARD crossed 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Δ = matched-pocket `summary_min` − mismatched-pocket `summary_min`, ligand bootstrap B = 2000. Matched uses Dual versus A-only in pocket B and Dual versus B-only in pocket A; mismatched swaps those score channels. This is a scoring-channel control, not redocking into a physically wrong site. Dark, CI excludes 0; gray, CI includes 0. (A) Main panels, eight primary pairs. EGFR/HER2 and AChE/BChE CIs exclude 0. (B) Unused-pool holdout Δ. All seven CIs include 0. EGFR/HER2 has no holdout. (C) Holdout versus main-panel `summary_min`. PPARG/PPARA holdout is 0.535 [0.350, 0.717]; JAK1/JAK2 stays same-direction (0.619 [0.420, 0.749]; drawn 20/20/18).

## Figure 6. Robustness checks and evidence boundary.

(A) Pocket-matched `summary_min` on the unified label-threshold grid for the eight primary pairs. (B) PIK3CA/mTOR PM48 versus PM110 Vina. (C) PM48 exhaustiveness 16 versus 8, recomputed from `scores_vina_E8_best.csv` with the same pocket-matched definition. (D) BindingDB and PubChem counts of all eight pairs plus the external-docking gate: no pair was packaged or docked as an external evaluation set (Table S11).

## Figure S1. Protocol and panel sensitivities.

Protocol grid, GNINA CNN rescoring of Vina poses, PM48 versus PM110, and exhaustiveness. Independent GNINA pose generation is Figure 4A.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C (right group).

## Figure S3. Additional paired bootstrap differences.

Descriptor and scaffold-versus-random leakage checks on the original docked set. Matched-versus-mismatched main/holdout Δ CIs are Figure 5.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Mismatched-pocket Δ CIs are Figure 5B.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 candidate-pair census, current primary n = 8, AND-like dual filter, and ligand-only full-map ECFP4. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

BindingDB / PubChem supply counts for all eight pairs and the remainder after independence filters (`crossdb_strict_supply_v1.csv`; `external_slice_summary_v1.csv`). No pair met the external-docking gate; nothing was docked as an external evaluation set.

## Figure S9. Additional ligand-structure controls.

Prespecified descriptors, covariate-adjusted logistic AUROC, and matched-subset weak-arm tests on the original docked set. The Vina-only logistic AUROC is not the Table 2 rank AUROC.

## Figure S10. Matched versus mismatched point estimates.

Bar charts of matched versus mismatched `summary_min` on the original docked set. Paired Δ CIs are Figure 5.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs and no decorative arrows.
