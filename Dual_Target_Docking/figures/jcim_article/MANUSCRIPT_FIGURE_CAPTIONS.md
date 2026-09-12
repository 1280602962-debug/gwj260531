# Figure captions (manuscript)

Target pairs are presented in a protein-system order throughout Figures 2–6 and Tables 1–3.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B using a pChEMBL threshold θ. Axes mark ≥ θ versus < θ rather than biologically “active” versus “inactive”. A-only and B-only serve as single-target selective controls. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, and Dual versus B-only is scored at target A. The lower of the two directional AUROCs is the descriptive summary$_{\mathrm{min}}$. (C) ChEMBL universe census under increasing data requirements: 2,164,618 pairs with at least one ligand measured at both targets, 63,790 with at least ten, 5,253 with at least ten Dual/A-only/B-only compounds at θ=6.0, and 86 meeting the strict 6.5/5.5 supply rule. The census summarizes the availability of paired experimental data. The primary evaluation comprised eight target pairs retained using the data-supply, panel-construction, and structural criteria described in Table 1.

## Figure 2. Docking performance depends on the experimental-state comparison.

Results are shown for the eight primary target pairs at θ = 6.0. Throughout (A) and (B), blue circles denote Dual versus A-only scored at target B, and orange squares Dual versus B-only scored at target A. (A) Change in AUROC after replacing the selective control with neither while holding the score channel fixed. Circles: target-B score, A-only → neither. Squares: target-A score, B-only → neither. Error bars are ligand-level bootstrap 95% confidence intervals. The dashed vertical line indicates no change. Numeric labels mark the EGFR/HER2 (0.378) and JAK1/TYK2 (0.444) target-A differences. Cluster-resampling intervals for those two differences are in Figure 6B. (B) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (C) Descriptive two-pocket comparison: directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands (neither n=4). Panel C is not a fixed-score contrast.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Primary rank-based Vina AUROC is shown alongside out-of-fold ECFP4 logistic-regression AUROC obtained under Bemis–Murcko scaffold-grouped cross-validation for both directional contrasts. Blue markers are Vina and orange markers ECFP4; circles are Dual versus A-only and squares Dual versus B-only. (B) Change after adding the corresponding Vina score to ECFP4 under the same scaffold-grouped cross-validation. Positive values favor ECFP4 plus Vina. Circles and squares denote Dual versus A-only and Dual versus B-only; grayscale is used so that color remains reserved for method in (A). Connectors join the two methods on the same arm and are not confidence intervals. The axis is limited to ±0.03 so that the 16 changes remain visible on their native scale. (C) Illustrative TPSA separation in AChE/BChE. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range.

## Figure 4. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring. Blue circles: primary panels. Orange squares: unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval. EGFR/HER2† had no unused-pool holdout. Color encodes the ligand set, not whether the interval excludes zero. (B) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. The dashed vertical line indicates AUROC = 0.5. † no unused-pool holdout available.

## Figure 5. Computational realization.

(A) Independent GNINA 1.3.2 pose generation versus Vina on the three evaluated target pairs. Filled circles are Vina; open squares are GNINA. Blue markers show directional summary$_{\mathrm{min}}$ and orange markers Dual versus neither. Gray segments connect the two tasks for the same engine and are not confidence intervals. The legend is placed outside the plotting area. GNINA points use the available-score subsets (EGFR/HER2: 11 versus 12 for Vina; PIK3CA/mTOR: 13 versus 14 for Vina), so they are not paired bootstrap estimates. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary pair (PIK3CA 4L23 and mTOR 4JT6). The first column is the primary receptor pair. Blue: primary receptors. Orange: PIK3CA substituted. Gold: mTOR substituted. Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds for all eight pairs (8 × 5). Horizontal segments show the five-seed min–max range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5. Panel C does not show the dual-versus-neither minus summary$_{\mathrm{min}}$ task difference; those values are in Table S9.

## Figure 6. Evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and the strict 6.5/5.5 rule. Cells marked † have a class with n < 10. (B) Cluster-resampling sensitivity of the fixed-score Dual-versus-neither minus Dual-versus-B-only differences on EGFR/HER2 and JAK1/TYK2 (target-A score). Circles, squares, and diamonds denote ligand, scaffold-cluster, and document-cluster resampling. (C) BindingDB compound counts after independence filters; color saturates at n = 20 per class. (D) BindingDB independent-source counts after the same filters; color saturates at 3 sources per class. After the independence filters, 0/8 pairs met the full external-evaluation admission criteria. Protocol-size and exhaustiveness checks for PIK3CA/mTOR are in Figure S3.

## Figure S1. Post-hoc formulation and screening diagnostics.

EGFR/HER2 operating points: class composition among the top ten compounds ranked by mean Vina score (denominator 110) and among compounds passing the median dual worst-target-score threshold on Dual+A-only+B-only (denominator 98, neither excluded). These analyses use different denominators and are exploratory (Table S13).

## Figure S2. Pocket-matched summary_min forest.

Vina ligand-level bootstrap 95% CIs and the best single-descriptor point estimates on the eight primary rows. Difference CIs of Vina minus the best descriptor are in Table S5, not in this figure.

## Figure S3. Protocol sensitivity on PIK3CA/mTOR.

(A) Vina summary$_{\mathrm{min}}$ on the PM48 and PM110 panels. PM48 is the primary PIK3CA/mTOR panel (quota 18/14/12/4; n_scored dual/A/B = 18/14/12; exhaustiveness = 16). PM110 is a larger protocol-sensitivity panel on the same pair. (B) Vina summary$_{\mathrm{min}}$ at exhaustiveness 16 versus 8 on PM48. Both panels are descriptive point estimates.

## Figure S4. Cognate redocking RMSD for the 14 primary receptors.

Heavy-atom RMSD of the cognate ligand after redocking into each primary receptor. Circles: top-1 pose. Diamonds: lowest heavy-atom RMSD among all saved poses. AChE 4EY7 and TYK2 3LXP deposited 8 poses; the other primary receptors deposited 9. The dashed line marks RMSD = 2 Å. Values above 3.5 Å are indicated off-scale; the EGFR 3POZ top-1 RMSD was 9.505 Å. All 14 primary receptors use chemically mapped CalcRMS from Table S2b. Historical coordinate assignment is in Table S2c.

## Figure S5. Detectable-effect simulation.

Probability that a simulated `summary_min` percentile 95% interval has lower bound > 0.5 or upper bound < 0.5, under the simulated true AUROC values and the observed class sizes. The simulation is available for three pairs only, uses 1,000 Monte Carlo replicates per grid point and 2,000 bootstrap replicates, and is an independent scenario diagnostic. At true AUROC = 0.50 the probabilities are 0.109, 0.119, and 0.128, indicating a coverage artifact. The figure is not an estimate of observed docking performance and does not by itself explain Table 2 uncertainty. Artwork file: `FigS6_detectable_effect`. Holdout, BindingDB-matrix, and cluster-copy figures remain in the repository and are not repeated in the submission SI.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks plus Dual versus neither on the same docking results, and the requirement that interpretation include ligand-only and pocket-correspondence controls. No numerical AUROCs.
