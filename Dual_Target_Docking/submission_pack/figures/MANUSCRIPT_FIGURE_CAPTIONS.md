# SUPERSEDED as a numbering lock.

Publication-facing figure/table numbering and source provenance are controlled only by `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`. Captions here match V4.

Target pairs are presented in a protein-system order throughout Figures 2–5 and Tables 1–3. There is no main-text Figure 6.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B using a pChEMBL threshold θ. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, and Dual versus B-only is scored at target A. (C) Compact supply check: 2,164,618 → 63,790 → 5,253 → 86. The primary evaluation comprised eight target pairs selected using the panel-construction and structural criteria in Table 1. That eight-pair set is not a direct continuation of the census counts.

## Figure 2. Core dual-target candidate-ranking result.

(A) Change in AUROC after replacing the single-target-active control with neither while holding the score channel fixed. Error bars are ligand-level bootstrap 95% confidence intervals. The EGFR/HER2 pocket-A difference is the post-fix value 0.462. (B) Directional AUROCs. (C) Directional summary_min versus Dual versus neither using the mean Vina score. Panel C is not a fixed-score contrast. The n=4 annotation marks PIK3CA/mTOR neither. (D) Eight-pair top-10% class composition under two-pocket mean Vina ranking, k=ceil(0.10 n). Bars are 100% stacked; k is annotated at right. EFdual,10% is in Table 3.

## Figure 3. Ligand chemistry and incremental docking information.

(A) Raw Vina rank AUROC versus scaffold-grouped out-of-fold ECFP4 AUROC, 8 pairs / 16 directions. Not a head-to-head algorithm comparison. (B) Change after adding the corresponding Vina score to ECFP4 under the same splits. AChE/BChE TPSA is Figure S1.

## Figure 4. Pocket correspondence.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Main-panel matched-minus-mismatched Δsummary_min, eight pairs. (B) Unused-pool internal holdout, seven pairs. Color encodes the ligand set, not whether the interval excludes zero. The AChE/BChE and EGFR/HER2 main-panel intervals exclude zero; the other six main-panel intervals include zero; all seven holdout intervals include zero.

## Figure 5. Computational robustness.

(A) Independent GNINA 1.3.2 pose generation versus Vina on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. Uniform EGFR GNINA: dual-versus-B-only ≈ 0.227; dual-versus-neither ≈ 0.705. (B) PIK3CA/mTOR receptor substitution of primary 4L23 / 4JT6 by 4JPS, 5DXT, or 4JSX. (C) Five-seed summary_min range for all eight pairs, same protocol as Table 2 (EGFR/HER2 uniform RDKit/Meeko, corrected boxes). AChE/BChE seed realizations do not use the primary membership (`same_protocol_as_primary=1`, `same_membership_as_primary=0`; primary n=95, available-case n=94, intersection n=88). Fixed-membership sensitivity (`results/qa/five_seed_fixed_membership_sensitivity.csv`) does not change the qualitative conclusion. RTMScore / CNN rescoring is Table S7.

## Figure S1. Ligand-chemistry detail.

(A) Vina CIs and best single-descriptor points. (B) AChE/BChE TPSA.

## Figure S2. Protocol sensitivity on PIK3CA/mTOR.

(A) PM48 versus PM110. (B) Exhaustiveness 16 versus 8.

## Figure S3. Cognate redocking RMSD for the 14 primary receptors.

Top-1 and lowest saved-pose heavy-atom RMSD. EGFR 3POZ top-1 = 1.019 Å; HER2 3RCD top-1 = 1.947 Å. The historical 9.505 Å value is not current.

## Figure S4. Label and source robustness.

(A) Activity-threshold grid. (B) Cluster resampling for EGFR/HER2 and JAK1/TYK2. EGFR/HER2 cluster bootstrap was recomputed with corrected-box scores and the frozen scaffold/document groupings.

## Figure S5. External-data eligibility.

BindingDB remainder after independence filters. 0/8 pairs met the full external docking gate. Eligibility screen, not external validation.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks, and ligand-only / pocket-correspondence controls. No numerical AUROCs.
