# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match the files in `figures/jcim_article/`.
All numerical values are those plotted from the frozen CSVs (unrounded). Table 2 in the text may round to three decimals.

Regenerate (from `Dual_Target_Docking/`):
`python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`
The script re-reads the CSVs, writes PDF/PNG/TIF, and fails if any plotted value disagrees with its source (`plotted_values.json`).

## Figure 1. DualFourClass task and pocket-matched readout.

(A) Dual-target evaluation requires four ligand classes measured on both targets: dual-actives and experimental single-end hard negatives (A_only, B_only), plus neither. (B) The primary metric scores dual versus A_only in pocket B and dual versus B_only in pocket A; summary_min is the smaller of the two AUROCs so that a weak arm cannot be hidden by pooling.

## Figure 2. Public-data supply of strict hard negatives.

(A) Minimum of the two strict hard-negative counts (A_only, B_only) for every target pair in the J0 ChEMBL audit (`j0_strict_label_supply.csv`). Dashed line, thick-panel gate (≥50); dotted line, thin-panel gate (≥20). Highlighted: the three thick pairs used as K=4 main panels, EGFR/HER2 (7 B_only; supply-limited case), and HDAC1/HDAC6 (metal enzyme; excluded). (B) Count-level comparison of the same four pairs in ChEMBL pChEMBL versus BindingDB equal-relation measurements (Table S12). No docking.

## Figure 3. Pocket-matched summary_min on the frozen K=4 set.

Vina (primary), RTMScore, GNINA CNN best-of-9, and the strongest trivial descriptor (heavy-atom count, MW, cLogP, or TPSA) with 95% ligand-bootstrap CIs. Vina CIs are the θ = 6.0 values in `unified_threshold_sensitivity_v2.csv` (Table 2). Best descriptor (right column, from `forest_summary_min_ci_v1.csv`): EGFR/HER2 cLogP; AChE/BChE TPSA; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count. Vertical dashed line, chance (0.5). GNINA is a single CNN channel, not a three-engine competition.

## Figure 4. Weak-arm asymmetry and physicochemical confounding.

(A) Directional Vina AUROCs at θ = 6.0: dual versus A_only (pocket B) and dual versus B_only (pocket A). (B) Vina pocket-matched summary_min versus the strongest trivial descriptor, with 95% CIs. (C) TPSA on the AChE/BChE panel by class (individual ligands from `assembled_AChE_BChE.csv`; horizontal line, median). Dual ligands are more polar than either hard-negative class, matching the TPSA baseline that exceeds Vina on this pair.

## Figure 5. Ligand-set holdout versus receptor swap.

(A) Pocket-matched summary_min on the main panel versus the unused-pool holdout (20/20/20; seed 20260731) for the three pairs with unused-pool supply. EGFR/HER2 has no holdout. (B) PM48 crystal swap: replacing PIK3CA 4L23 with 4JPS or 5DXT (mTOR held at 4JT6) or replacing mTOR 4JT6 with 4JSX (PIK3CA held at 4L23). Error bars are 95% ligand-bootstrap CIs.

## Figure 6. Pocket-matched versus wrong-pocket control.

(A) Main K=4 panels: pocket-matched Vina summary_min versus the wrong-pocket control (`pocket_matched_directional_v1.csv`). Matched exceeds wrong-pocket on all four pairs. (B) Unused-pool holdout: the inequality reverses (wrong-pocket ≥ matched) on all three pairs with holdout supply (`holdout_pocket_matched_v1.csv`). EGFR/HER2 has no holdout. (C) The reversal remains after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) (`holdout_matched_wrong_pocket_summary_v1.csv`). Wrong-pocket remains ≥ matched on all nine cells; matching does not restore the main-panel inequality. (D) Scoring-free contact-count AUROC on pocket A (D vs A_only) and pocket B (D vs B_only) versus Vina wrong-pocket summary_min (`wrong_pocket_contact_v1_output.txt`; not a PLIF). B-arm contact is above chance; the magnitude does not reproduce Vina on PIK3CA/mTOR.

## Figure 7. Ligand-structure association and matched-subset tests.

(A) ECFP4 logistic regression under scaffold GroupKFold versus pocket-matched Vina on both directional arms (`ligand_ml_baseline_scaffold_cv_v1.csv`). Fingerprint AUROCs are chemotype–label association, not evidence of pocket physics. (B) Pocket-matched Vina versus all four trivial descriptors (heavy-atom count, MW, cLogP, TPSA) with 95% CIs. Descriptor CIs are from `forest_summary_min_ci_v1.csv`; Vina from θ = 6.0. Figure 4 reports only the strongest descriptor per pair; this panel shows all four. (C) Weak-arm (D vs B_only) logistic AUROC of Vina alone versus Vina plus heavy-atom count and TPSA, with the Vina odds ratio (`covariate_adjusted_v1.csv`). EGFR/HER2 score-only in that table is 0.5703 (the table’s logistic AUROC of feature `vina_A`), which is not the rank AUROC 0.4297 in Table 2. (D) D vs B_only after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) versus the unmatched full-panel contrast (`matched_subset_directional_v1.csv`). Error bars are the table’s single-contrast 95% CIs.

## Figure S1. Protocol knobs that do not change the ranking.

(A) Pocket-matched summary_min across the unified label-threshold grid (`unified_threshold_sensitivity_v2.csv`). Open markers are underpowered cells (EGFR/HER2 strict, n_B_only=7; PIK3CA/mTOR at θ=5.5, n_B_only=5, and strict, n_B_only=4). At the primary θ=6.0, PIK3CA/mTOR is the highest point estimate; AChE/BChE is flat at 0.6058 across the grid. The underpowered θ=5.5 PIK3CA/mTOR cell (0.5017) is not a ranking contradiction. (B) GNINA CNN mode01 versus best-of-9 versus the same-panel Vina reference (`gnina_pocket_matched_mode01_vs_best9_k4_v1.csv`). Best-of-9 versus mode01 moves summary_min by −0.04 to +0.08. EGFR/HER2 and AChE/BChE remain below chance on both GNINA channels. PIK3CA/PIK3CB GNINA best-of-9 is 0.533 versus Vina 0.500 (near chance). (C) PIK3CA/mTOR PM48 versus the PM110 expansion for Vina, RTMScore, and GNINA best-of-9 (`pm110_vs_pm48_pocket_matched_v1.csv`). (D) PM48 Vina at exhaustiveness 16 versus 8, computed from `scores_vina_E8_best.csv` (empty affinities skipped; ligands labeled neither were excluded) with the same pocket-matched definition, beside single-target enrichment AUROC and EF1% on 4L23 and 4JT6 (`single_target_enrichment_v1.csv`).

## Figure S2. Equal-relation supply and holdout sampling shift.

(A) Minimum strict hard-negative counts for the K=4 pairs in ChEMBL pChEMBL, BindingDB/PubChem `equal_only`, and BindingDB/PubChem `as_is` (`crossdb_strict_supply_v1.csv`). Count-level only; no docking. `as_is` lets EGFR/HER2 pass ≥50 because censored `>` values are treated as point estimates; `equal_only` does not. (B) Holdout minus main-panel mean pChEMBL for dual pA, A_only pA, and B_only pB (`holdout_vs_main_potency_size_v1.csv`). Sampling shift is real, especially on PIK3CA/mTOR, but does not reverse Figure 6C.

## TOC graphic (For Table of Contents Only).

DualFourClass-Bench evaluates whether docking ranks dual-target ligands above experimental single-end hard negatives in both pockets. The graphic does not report numerical AUROCs and is not a reuse of Figure 1.
