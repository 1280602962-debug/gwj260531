# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate the current PR32 artwork with `python3 figures/jcim_article/scripts/update_figures_pr32.py`; the source commit, input files, and SHA-256 values are recorded in `figures/jcim_article/plotted_values.json`.

Target pairs are shown in a protein-system order (kinase, cholinesterase, coagulation protease, nuclear receptor) throughout Figures 2–6 and Tables 1–3. Display order is not a data-collection sequence.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B. A-only and B-only serve as single-target selective controls. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, and Dual versus B-only is scored at target A. The lower of the two directional AUROCs is the descriptive summary$_{\mathrm{min}}$. (C) ChEMBL universe census under increasing data requirements: 2,164,618 pairs with at least one ligand measured at both targets, 63,790 with at least ten, 5,253 with at least ten Dual/A-only/B-only compounds at θ=6.0, and 86 meeting the strict 6.5/5.5 supply rule. The eight primary pairs were then selected using pair-specific panel and structural criteria.

For a split-asset submission, the conceptual A–B schematic is supplied as the editable PowerPoint `figures/editable/Figure1_AB_editable_v2.pptx`, and the Python-rendered census panel is supplied as `Fig1_C_chEMBL_supply.*`. The combined `Fig1_four_state_and_supply.*` file is retained as a reference composite.

## Figure 2. Docking performance depends on the experimental-state comparison.

Results are shown for the eight primary target pairs at θ = 6.0. (A) Change in AUROC after replacing the selective control with neither while holding the score channel fixed. Circles: target-A score, B-only → neither. Squares: target-B score, A-only → neither. Error bars are ligand-level bootstrap 95% confidence intervals. The dashed vertical line indicates no change. Numeric labels mark the EGFR/HER2 (0.378) and JAK1/TYK2 (0.444) target-A differences. Cluster-resampling intervals for those two differences are in Figure 6B. (B) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (C) Directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands (neither n=4).

## Figure 3. Ligand chemistry as a competing explanation.

(A) Pocket-matched Vina rank AUROC and ECFP4 logistic-regression AUROC under Bemis–Murcko scaffold-grouped cross-validation for both directional contrasts. Circles and squares are Vina Dual versus A-only and Dual versus B-only; triangles and diamonds are the corresponding ECFP4 estimates. (B) Change in scaffold-grouped cross-validated AUROC after adding the corresponding Vina score to ECFP4. Positive values favor ECFP4 plus Vina. Circles and squares denote Dual versus A-only and Dual versus B-only. The axis is limited to ±0.03 so that the 16 increments remain visible on their native scale. (C) Illustrative TPSA separation in AChE/BChE. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation versus Vina on the three evaluated target pairs. Filled circles are Vina; open squares are GNINA. Blue markers show directional summary$_{\mathrm{min}}$ and orange markers Dual versus neither. Gray segments connect the two tasks for the same engine. GNINA points use the available-score subsets (EGFR/HER2: 11 versus 12 for Vina; PIK3CA/mTOR: 13 versus 14 for Vina), so they are not paired bootstrap estimates. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary structures (PIK3CA 4L23 and mTOR 4JT6). Blue: primary receptors. Orange: PIK3CA substituted. Gold: mTOR substituted. Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds for all eight pairs (8 × 5). Horizontal segments show the min–max range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring. Blue circles: primary panels. Orange squares: unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval. EGFR/HER2 had no holdout (n/a). Color encodes the ligand set, not whether the interval excludes zero. (B) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. The dashed vertical line indicates AUROC = 0.5.

## Figure 6. Evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and the strict 6.5/5.5 rule. Cells marked † have a class with n < 10. (B) Fixed-score Dual-versus-neither minus Dual-versus-B-only differences on EGFR/HER2 and JAK1/TYK2 under ligand, scaffold-cluster, and document-cluster resampling (target-A score). Circles, squares, and diamonds denote the three resampling levels. (C) BindingDB compound counts after independence filters; color saturates at the external compound-count criterion of 20. (D) BindingDB independent-source counts after the same filters; color saturates at the source criterion of 3. No pair meets both external criteria. Protocol-size and exhaustiveness checks for PIK3CA/mTOR are in Figure S1.

## Figure S1. Protocol sensitivity on PIK3CA/mTOR.

(A) Vina summary$_{\mathrm{min}}$ on the PM48 and PM110 panels. (B) Vina summary$_{\mathrm{min}}$ at exhaustiveness 16 versus 8 on PM48. Both panels are descriptive point estimates.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C.

## Figure S3. Additional paired bootstrap differences.

Paired bootstrap comparisons on the original docked set: matched versus mismatched scoring in the primary and holdout panels, Vina versus the strongest single-descriptor baseline, and ECFP4 estimates under scaffold-grouped and random folds. The later-withdrawn PIK3CA/PIK3CB row is marked with a dagger.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows. PIK3CA/PIK3CB is omitted.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Matched-minus-mismatched Δ CIs are Figure 5A.

## Figure S6. Detectable-effect simulation.

Probability that a ligand-level 95% confidence interval for `summary_min` excludes 0.5 under the simulated true AUROC values and the observed class sizes. The simulation is available for three pairs only, uses 1,000 Monte Carlo replicates per grid point and 2,000 bootstrap replicates, and is a power-oriented diagnostic rather than an estimate of observed docking performance.

## Figure S7. Post-hoc formulation and screening diagnostics.

EGFR/HER2 operating points: class composition among the top ten compounds ranked by mean Vina score and among compounds passing the median dual worst-target-score threshold. These analyses use different denominators and are exploratory (Table S13).

## Figure S8. BindingDB-native slice.

BindingDB compound and source counts after independence filters for all eight target pairs (`external_slice_summary_v1.csv`). The same two matrices appear as Figure 6C and 6D. No pair meets the external admission criteria.

## Figure S9. Additional ligand-structure controls.

Additional ligand-structure controls on the original docked set: ECFP4 versus Vina directional AUROCs, Vina and prespecified single-descriptor baselines, covariate-adjusted Dual versus B-only logistic models, and potency- or size-matched subsets. The Vina-only logistic AUROC in panel C is distinct from the rank-based AUROC in the primary analysis.

## Figure S10. Matched versus mismatched point estimates.

Matched and mismatched scoring point estimates for the original primary panels, unused-pool holdouts, potency- and size-matched holdout subsets, and contact-count controls. Paired confidence intervals for the primary and holdout differences are shown in Figure 5A.

## Figure S11. Fixed-score comparator differences under cluster resampling.

Target-A score differences between Dual-versus-neither and Dual-versus-B-only for EGFR/HER2 and JAK1/TYK2 under ligand, scaffold-cluster, and document-cluster resampling. The same estimates are shown in Figure 6B. Cluster intervals supplement the primary ligand-level intervals. EGFR/HER2 excludes zero under both cluster schemes; the JAK1/TYK2 document-cluster interval includes zero.

## Figure S12. Cognate redocking RMSD for the 14 primary receptors.

Heavy-atom RMSD of the cognate ligand after redocking into each primary receptor. Circles: top-1 pose. Squares: minimum among the top three poses (omitted when the best saved pose is beyond rank 3). Diamonds: minimum among all saved poses (best-of-9). The dashed line marks RMSD = 2 Å. Values are those locked in Table S2.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks plus Dual versus neither, and the requirement that interpretation include ligand-only and pocket-correspondence controls. No numerical AUROCs.
