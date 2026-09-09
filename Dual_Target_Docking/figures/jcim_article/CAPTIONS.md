# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate the current PR32 artwork with `python3 figures/jcim_article/scripts/update_figures_pr32.py`; the source commit, input files, and SHA-256 values are recorded in `figures/jcim_article/plotted_values.json`.

Target pairs are presented in a protein-system order throughout Figures 2–6 and Tables 1–3.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B using a pChEMBL threshold θ. Axes mark ≥ θ versus < θ rather than biologically “active” versus “inactive”. A-only and B-only serve as single-target selective controls. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, and Dual versus B-only is scored at target A. The lower of the two directional AUROCs is the descriptive summary$_{\mathrm{min}}$. (C) ChEMBL universe census under increasing data requirements: 2,164,618 pairs with at least one ligand measured at both targets, 63,790 with at least ten, 5,253 with at least ten Dual/A-only/B-only compounds at θ=6.0, and 86 meeting the strict 6.5/5.5 supply rule. The census summarizes the availability of paired experimental data. The primary evaluation comprised eight target pairs retained using the data-supply, panel-construction, and structural criteria described in Table 1.

For a split-asset submission, the conceptual A–B schematic is supplied as the editable PowerPoint `figures/editable/Figure1_AB_editable_v2.pptx`, and the Python-rendered census panel is supplied as `Fig1_C_chEMBL_supply.*`. The combined `Fig1_four_state_and_supply.*` file is retained as a reference composite.

## Figure 2. Docking performance depends on the experimental-state comparison.

Results are shown for the eight primary target pairs at θ = 6.0. Throughout (A) and (B), blue circles denote Dual versus A-only scored at target B, and orange squares Dual versus B-only scored at target A. (A) Change in AUROC after replacing the selective control with neither while holding the score channel fixed. Circles: target-B score, A-only → neither. Squares: target-A score, B-only → neither. Error bars are ligand-level bootstrap 95% confidence intervals. The dashed vertical line indicates no change. Numeric labels mark the EGFR/HER2 (0.378) and JAK1/TYK2 (0.444) target-A differences. Cluster-resampling intervals for those two differences are in Figure 6B. (B) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (C) Descriptive two-pocket comparison: directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands (neither n=4). Panel C is not a fixed-score contrast.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Primary rank-based Vina AUROC is shown alongside out-of-fold ECFP4 logistic-regression AUROC obtained under Bemis–Murcko scaffold-grouped cross-validation for both directional contrasts. Blue markers are Vina and orange markers ECFP4; circles are Dual versus A-only and squares Dual versus B-only. (B) Change in scaffold-grouped cross-validated AUROC after adding the corresponding Vina score to ECFP4. Positive values favor ECFP4 plus Vina. Circles and squares denote Dual versus A-only and Dual versus B-only. The axis is limited to ±0.03 so that the 16 increments remain visible on their native scale. (C) Illustrative TPSA separation in AChE/BChE. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation versus Vina on the three evaluated target pairs. Filled circles are Vina; open squares are GNINA. Blue markers show directional summary$_{\mathrm{min}}$ and orange markers Dual versus neither. Gray segments connect the two tasks for the same engine. GNINA points use the available-score subsets (EGFR/HER2: 11 versus 12 for Vina; PIK3CA/mTOR: 13 versus 14 for Vina), so they are not paired bootstrap estimates. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary pair (PIK3CA 4L23 and mTOR 4JT6). The first column is the primary receptor pair. Blue: primary receptors. Orange: PIK3CA substituted. Gold: mTOR substituted. Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds for all eight pairs (8 × 5). Horizontal segments show the five-seed min–max range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring. Blue circles: primary panels. Orange squares: unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval. EGFR/HER2† had no unused-pool holdout. Color encodes the ligand set, not whether the interval excludes zero. (B) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. The dashed vertical line indicates AUROC = 0.5. † no unused-pool holdout available.

## Figure 6. Evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and the strict 6.5/5.5 rule. Cells marked † have a class with n < 10. (B) Cluster-resampling sensitivity of the fixed-score Dual-versus-neither minus Dual-versus-B-only differences on EGFR/HER2 and JAK1/TYK2 (target-A score). Circles, squares, and diamonds denote ligand, scaffold-cluster, and document-cluster resampling. (C) BindingDB compound counts after independence filters; color saturates at the compound-count gate of 20 per class. (D) BindingDB independent-source counts after the same filters; color saturates at the source gate of 3 per class. No pair met both compound-count and independent-source criteria (0/8). Protocol-size and exhaustiveness checks for PIK3CA/mTOR are in Figure S1.

## Figure S1. Protocol sensitivity on PIK3CA/mTOR.

(A) Vina summary$_{\mathrm{min}}$ on the PM48 and PM110 panels. (B) Vina summary$_{\mathrm{min}}$ at exhaustiveness 16 versus 8 on PM48. Both panels are descriptive point estimates.

## Figure S2–S3 and S9–S10.

These identifiers belong to historical original-set sensitivity tables. They are **not** part of the current eight-pair submission artwork. The corresponding SI tables remain.

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

## Figure S11. Fixed-score comparator differences under cluster resampling.

Target-A score differences between Dual-versus-neither and Dual-versus-B-only for EGFR/HER2 and JAK1/TYK2 under ligand, scaffold-cluster, and document-cluster resampling. The same estimates are shown in Figure 6B. Cluster intervals supplement the primary ligand-level intervals. EGFR/HER2 excludes zero under both cluster schemes; the JAK1/TYK2 document-cluster interval includes zero.

## Figure S12. Cognate redocking RMSD for the 14 primary receptors.

Heavy-atom RMSD of the cognate ligand after redocking into each primary receptor. Circles: top-1 pose. Squares: minimum among the top three poses (omitted when the best saved pose is beyond rank 3). Diamonds: minimum among all saved poses (best-of-9). The dashed line marks RMSD = 2 Å. The axis is limited to 3.5 Å so that the 0–2 Å majority remains readable; EGFR 3POZ top-1 (9.505 Å) is drawn off-scale. Values are those locked in Table S2.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks plus Dual versus neither on the same docking results, and the requirement that interpretation include ligand-only and pocket-correspondence controls. No numerical AUROCs.
