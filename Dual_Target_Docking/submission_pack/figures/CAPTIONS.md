# SUPERSEDED as a numbering lock.

Publication-facing figure/table numbering and source provenance are controlled only by:

`docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`

The captions below match V4. Historical pre-fix tokens (0.378, 9.505 Å, JAK1-only Figure 2D, main-text Figure 6) are not current results.

Regenerate: `python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root Dual_Target_Docking`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimental states defined by threshold θ: dual / A-only / B-only / neither. (B) Directional evaluation: dual versus A-only uses the target B score; dual versus B-only uses the target A score. (C) Paired experimental coverage 2,164,618 / 63,790 / 5,253 / 86. Primary evaluation: 8 target pairs selected using panel-construction and structural criteria. The eight-pair set is not a direct continuation of the census funnel. No performance results.

## Figure 2. Core dual-target candidate-ranking result.

(A) Fixed-score ΔAUROC (dual–neither minus dual–matched single-target-active) for eight pairs and two score directions; 95% CI. EGFR pocket-A difference is 0.446. (B) Directional AUROCs: dual versus A-only using pocket B; dual versus B-only using pocket A. (C) Directional summary_min versus two-pocket mean dual-versus-neither AUROC. Panel C is not a fixed-score comparison. n = 4 marks the PIK3CA/mTOR neither sample. (D) Eight-pair top-10% class composition under two-pocket mean Vina ranking, k=ceil(0.10 n). Bars are 100% stacked dual / A-only / B-only / neither; k is at right. EFdual,10% is in Table 3, not on the figure.

## Figure 3. Ligand chemistry and incremental docking information.

(A) Vina directional rank AUROC versus ECFP4 scaffold-GroupKFold OOF AUROC, 8 pairs / 16 directions. Not a formal algorithm leaderboard. (B) ΔAUROC = ECFP4+docking − ECFP4 on the same folds and ligands. AChE/BChE TPSA is Figure S1.

## Figure 4. Pocket correspondence.

(A) Main-panel matched-minus-mismatched Δsummary_min, eight pairs, 95% CI. (B) Unused-pool internal holdout matched-minus-mismatched Δsummary_min, seven pairs, 95% CI. EGFR/HER2 has no holdout. Color distinguishes main versus holdout, not significance. Only the AChE/BChE main-panel interval excludes zero; EGFR/HER2 includes zero; all seven holdout intervals include zero.

## Figure 5. Computational robustness.

(A) Independent GNINA versus Vina on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2: directional summary_min and dual-versus-neither. Uniform EGFR GNINA: D-versus-B-only ≈ 0.227; D-versus-neither ≈ 0.705. (B) PIK3CA/mTOR receptor substitution of primary 4L23 / 4JT6 by 4JPS, 5DXT, or 4JSX. (C) Five-seed summary_min range, eight pairs, same boxes as Table 2 (EGFR/HER2 uniform RDKit/Meeko, corrected boxes). AChE/BChE seed realizations have changing available-case membership; fixed-membership sensitivity gives the same qualitative conclusion. The figure asks whether results depend on a single docking realization. RTMScore / CNN rescoring is Table S7.

## Figure S1. Ligand-chemistry detail.

(A) Best single descriptor versus Vina summary_min forest. (B) AChE/BChE TPSA distribution.

## Figure S2. PIK3CA/mTOR protocol sensitivity.

(A) PM48 versus PM110. (B) E=16 versus E=8.

## Figure S3. Cognate redocking RMSD.

Fourteen primary receptor slots; top-1 and lowest saved pose; 2 Å reference. Post-fix EGFR 3POZ top-1 = 1.019 Å; HER2 3RCD top-1 = 1.947 Å. The historical 9.505 Å value is not current.

## Figure S4. Label and source robustness.

(A) Theta threshold sensitivity. (B) Cluster-resampling intervals for EGFR/HER2 and JAK1/TYK2 after substituting current scores into the frozen scaffold and document groupings. Cluster intervals are a sensitivity analysis and do not replace Table 2. Does not repeat Figure 4 holdout.

## Figure S5. External-data eligibility.

(A) BindingDB independent-filter compound counts. (B) Independent-source counts. 0/8 passed the full external docking gate. Eligibility screen, not external validation.

## TOC graphic (For Table of Contents Only).

Four experimental states, two directional tasks, and ligand-only / pocket-correspondence controls. No numerical AUROCs.
