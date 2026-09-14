# Figure captions (manuscript)

Target pairs are presented in a protein-system order throughout Figures 2–6 and Tables 1–3.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B using a pChEMBL threshold θ. Axes mark ≥ θ versus < θ rather than biologically “active” versus “inactive”. A-only and B-only are pair-relative experimental states at the chosen threshold, not proteome-wide selective ligands. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, and Dual versus B-only is scored at target A. The descriptive weaker-arm summary is $\mathrm{summary}_{\mathrm{min}}=\min[\mathrm{AUROC}_{D/A}(B),\,\mathrm{AUROC}_{D/B}(A)]$. (C) Compact supply check: paired experimental coverage decreases under four-state requirements (2,164,618 → 63,790 → 5,253 → 86). The primary evaluation comprised eight target pairs selected using the panel-construction and structural criteria in Table 1. That eight-pair set is not a direct continuation of the census counts.

## Figure 2. Docking performance depends on the experimental-state comparison.

Results are shown for the eight primary target pairs at θ = 6.0. Throughout (A) and (B), blue circles denote Dual versus A-only scored at target B, and orange squares Dual versus B-only scored at target A. (A) Change in AUROC after replacing the single-target-active control with neither while holding the score channel fixed. Circles: target-B score, A-only → neither. Squares: target-A score, B-only → neither. Error bars are ligand-level bootstrap 95% confidence intervals. The dashed vertical line indicates no change. Numeric labels mark the EGFR/HER2 (0.378) and JAK1/TYK2 (0.444) target-A differences. Cluster-resampling intervals for those two differences are in Figure 6B. (B) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (C) Descriptive two-pocket comparison: directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands (neither n=4). Panel C is not a fixed-score contrast. (D) JAK1/TYK2 operating point: class composition among the top ten compounds ranked by the two-pocket mean Vina score (denominator 109). The higher dual-versus-neither AUROC does not correspond to exclusion of single-target-active ligands from the top ranks.

## Figure 3. Ligand chemistry as a competing explanation, not a head-to-head predictive benchmark.

(A) Vina points are raw primary-score rankings, whereas ECFP4 points are scaffold-grouped out-of-fold logistic-regression predictions. Blue markers are Vina and orange markers ECFP4; circles are Dual versus A-only and squares Dual versus B-only. The panel is not a head-to-head algorithm comparison. (B) Change after adding the corresponding Vina score to ECFP4 under the same scaffold-grouped cross-validation. Positive values favor ECFP4 plus Vina. Circles and squares denote Dual versus A-only and Dual versus B-only; grayscale is used so that color remains reserved for method in (A). Connectors join the two methods on the same arm and are not confidence intervals. The axis is limited to ±0.03 so that the 16 changes remain visible on their native scale. (C) Illustrative TPSA separation in AChE/BChE. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range.

## Figure 4. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring. Blue circles: primary panels. Orange squares: unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval. EGFR/HER2† had no unused-pool holdout. Color encodes the ligand set, not whether the interval excludes zero. (B) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. The dashed vertical line indicates AUROC = 0.5. † no unused-pool holdout available.

## Figure 5. Computational realization.

(A) Independent GNINA 1.3.2 pose generation versus Vina on the three evaluated target pairs. Filled circles are Vina; open squares are GNINA. Blue markers show directional summary$_{\mathrm{min}}$ and orange markers Dual versus neither. Gray segments connect the two tasks for the same engine and are not confidence intervals. GNINA points use the available-score subsets (EGFR/HER2: 11 versus 12 for Vina; PIK3CA/mTOR: 13 versus 14 for Vina), so they are not paired bootstrap estimates. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary pair (PIK3CA 4L23 and mTOR 4JT6). The first column is the primary receptor pair. Blue: primary receptors. Orange: PIK3CA substituted. Gold: mTOR substituted. Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds for all eight pairs (8 × 5). Horizontal segments show the five-seed min–max range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5. Panel C does not show the dual-versus-neither minus summary$_{\mathrm{min}}$ task difference; those seed-level values are archived with the five-seed score tables.

## Figure 6. Evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and the strict 6.5/5.5 rule. Cells marked † have a class with n < 10. (B) Cluster-resampling sensitivity of the fixed-score Dual-versus-neither minus Dual-versus-B-only differences on EGFR/HER2 and JAK1/TYK2 (target-A score). Circles, squares, and diamonds denote ligand, scaffold-cluster, and document-cluster resampling. (C, D) BindingDB eligibility counts after independence filters; color saturates at n = 20 per class in (C) and at 3 sources per class in (D). PubChem was an additional paired-data availability check and is not merged into these panels (Table S8). After the independence filters, 0/8 pairs met the full external-evaluation admission criteria. Protocol-size and exhaustiveness checks for PIK3CA/mTOR are in Figure S3.

## Figure S1. JAK1/TYK2 AND-type two-pocket filter.

Class composition among JAK1/TYK2 compounds passing the median dual worst-target-score threshold on Dual+A-only+B-only (denominator 95, neither excluded). The corresponding Top-10 ranking is in Figure 2D. This analysis is exploratory; the counts are reported in Results 3.2.

## Figure S2. Pocket-matched summary_min forest.

Vina ligand-level bootstrap 95% CIs and the best single-descriptor point estimates on the eight primary rows. Difference CIs of Vina minus the best descriptor are in Table S5, not in this figure.

## Figure S3. Protocol sensitivity on PIK3CA/mTOR.

(A) Vina summary$_{\mathrm{min}}$ on the PM48 and PM110 panels. PM48 is the primary PIK3CA/mTOR panel (quota 18/14/12/4; n_scored dual/A/B = 18/14/12; exhaustiveness = 16). PM110 is a larger protocol-sensitivity panel on the same pair. (B) Vina summary$_{\mathrm{min}}$ at exhaustiveness 16 versus 8 on PM48. Both panels are descriptive point estimates.

## Figure S4. Cognate redocking RMSD for the 14 primary receptors.

Heavy-atom RMSD of the cognate ligand after redocking into each primary receptor. Circles: top-1 pose. Diamonds: lowest heavy-atom RMSD among all saved poses. AChE 4EY7 and TYK2 3LXP deposited 8 poses; the other primary receptors deposited 9. The dashed line marks RMSD = 2 Å. Values above 3.5 Å are indicated off-scale; the EGFR 3POZ top-1 RMSD was 9.505 Å. All 14 slots use the unified chemically mapped CalcRMS table (Table S2). The detectable-effect simulation (`FigS6_detectable_effect`) and holdout, BindingDB-matrix, and cluster-copy figures remain in the repository and are not repeated in the submission SI.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks plus Dual versus neither on the same docking results, and the requirement that interpretation include ligand-only and pocket-correspondence controls. No numerical AUROCs.
