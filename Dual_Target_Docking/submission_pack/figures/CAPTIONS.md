# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate the current PR32 artwork with `python3 figures/jcim_article/scripts/update_figures_pr32.py`; the source commit, input files, and SHA-256 values are recorded in `figures/jcim_article/plotted_values.json`.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B. A-only and B-only serve as single-target selective controls. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, whereas Dual versus B-only is scored at target A. The lower of the two directional AUROCs is reported as the descriptive summary$_{\mathrm{min}}$. (C) ChEMBL universe census under increasing data requirements: 2,164,618 pairs with at least one ligand measured at both targets, 63,790 with at least ten, 5,253 with at least ten Dual/A-only/B-only compounds at θ=6.0, and 86 meeting the strict 6.5/5.5 supply rule. The eight primary pairs were then selected using pair-specific panel and structural criteria.

For a split-asset submission, the conceptual A–B schematic is supplied as the editable PowerPoint `figures/editable/Figure1_AB_editable_v2.pptx`, and the Python-rendered census panel is supplied as `Fig1_C_chEMBL_supply.*`. The combined `Fig1_four_state_and_supply.*` file is retained as a reference composite.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Results are shown for the eight primary target pairs at θ = 6.0. Horizontal gray rules separate the three initially evaluated pairs from the five pairs added after the census. (A) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (B) Directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). Error bars are ligand-level bootstrap 95% confidence intervals. The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands. (C) Change in AUROC after replacing B-only with neither while holding the target-A score fixed. Error bars are 95% confidence intervals; diamonds indicate comparisons with an underpowered neither class. The dashed vertical line indicates no change. Fixed-score comparator cluster intervals for EGFR/HER2 and JAK1/TYK2 are provided in Figure S11.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Pocket-matched Vina rank AUROC and ECFP4 logistic-regression AUROC under Bemis–Murcko scaffold-grouped cross-validation for both directional contrasts. (B) Change in scaffold-grouped cross-validated AUROC after adding the corresponding Vina score to ECFP4. Positive values favor ECFP4 plus Vina; circles and squares denote Dual versus A-only and Dual versus B-only, respectively. (C) Topological polar surface area (TPSA) of the AChE/BChE ligands by activity class. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range, respectively.

## Figure 4. Computational realization.

(A) Comparison of Vina and independent GNINA 1.3.2 pose generation for the three evaluated target pairs. Circles and squares denote Vina and GNINA, respectively; filled blue symbols show directional summary$_{\mathrm{min}}$, and open orange symbols show Dual versus neither. GNINA points use the available-score subsets (EGFR/HER2: 11 versus 12 for Vina; PIK3CA/mTOR: 13 versus 14 for Vina), so they are not interpreted as paired bootstrap estimates. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary structures (PIK3CA 4L23 and mTOR 4JT6). Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds. Horizontal segments show the range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring in the primary panels. (B) The same paired difference in unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval; blue intervals exclude zero and gray intervals include zero. EGFR/HER2 had no holdout. (C) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. Dashed vertical lines indicate zero in panels A and B and AUROC = 0.5 in panel C.

## Figure 6. Robustness checks and evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and the strict 6.5/5.5 rule. (B) PIK3CA/mTOR PM48/PM110 comparison. (C) PIK3CA/mTOR exhaustiveness 16/8 comparison. Panels B and C show descriptive point estimates. (D) BindingDB class counts after independence filters for all eight pairs; color saturates at the external compound-count criterion of 20. No pair meets all external criteria; source counts are shown in Figure S8.

## Figure S1. Protocol and panel sensitivities.

Historical original-set protocol diagnostics are retained in `archive_before_pr32/` and are excluded from the current top-level figure set.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C (right group).

## Figure S3. Additional paired bootstrap differences.

Paired bootstrap comparisons on the original docked set: matched versus mismatched scoring in the primary and holdout panels, Vina versus the strongest single-descriptor baseline, and ECFP4 estimates under scaffold-grouped and random folds. The later-withdrawn PIK3CA/PIK3CB row is marked with a dagger.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows. PIK3CA/PIK3CB is omitted.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Mismatched-pocket Δ CIs are Figure 5B.

## Figure S6. Detectable-effect simulation.

Probability that a ligand-level 95% confidence interval for `summary_min` excludes 0.5 under the simulated true AUROC values and the observed class sizes. The simulation is available for three pairs only, uses 1,000 Monte Carlo replicates per grid point and 2,000 bootstrap replicates, and is a power-oriented diagnostic rather than an estimate of observed docking performance.

## Figure S7. Post-hoc formulation and screening diagnostics.

EGFR/HER2 operating points: class composition among the top ten compounds ranked by mean Vina score and among compounds passing the median dual worst-target-score threshold. These analyses use different denominators and are exploratory (Table S13).

## Figure S8. BindingDB-native slice.

BindingDB compound and source counts after independence filters for all eight target pairs (`external_slice_summary_v1.csv`). No pair meets the external admission criteria.

## Figure S9. Additional ligand-structure controls.

Additional ligand-structure controls on the original docked set: ECFP4 versus Vina directional AUROCs, Vina and prespecified single-descriptor baselines, covariate-adjusted Dual versus B-only logistic models, and potency- or size-matched subsets. The Vina-only logistic AUROC in panel C is distinct from the rank-based AUROC in the primary analysis.

## Figure S10. Matched versus mismatched point estimates.

Matched and mismatched scoring point estimates for the original primary panels, unused-pool holdouts, potency- and size-matched holdout subsets, and contact-count controls. Paired confidence intervals for the primary and holdout differences are shown in Figure 5.

## Figure S11. Fixed-score comparator differences under cluster resampling.

Target-A score differences between Dual-versus-neither and Dual-versus-B-only for EGFR/HER2 and JAK1/TYK2 under ligand, scaffold-cluster, and document-cluster resampling. Cluster intervals supplement the primary ligand-level intervals. EGFR/HER2 excludes zero under both cluster schemes; the JAK1/TYK2 document-cluster interval includes zero.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs and no decorative arrows.
