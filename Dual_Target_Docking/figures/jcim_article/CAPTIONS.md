# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match the files in `figures/jcim_article/`.
All numerical values are those plotted from the frozen CSVs (unrounded). Table 2 in the text may round to three decimals.

Regenerate (from `Dual_Target_Docking/`):
`python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`
The script re-reads the CSVs, writes PDF/PNG/TIF, and fails if any plotted value disagrees with its source (`plotted_values.json`).

## Figure 1. Dual-target docking as dual-versus-selective discrimination.

(A) A strict dual-target benchmark distinguishes four experimentally labeled ligand states: dual-active (D), A-selective (A_only), B-selective (B_only), and neither. A_only and B_only are selectivity hard negatives: they are active on one target and can produce plausible docking scores, yet they lack activity on the other. Neither is curated as part of the four-state panel but is not used in the primary AUROCs. (B) The prespecified primary readout is two directional pairwise discriminations, not a four-class classifier and not a pooled score. Dual versus A_only is scored in pocket B; dual versus B_only is scored in pocket A. The pair-level summary is the weaker arm (summary_min), so a favorable score on one target cannot hide directional failure on the other. Pooled docking scores are retained only as a control.

## Figure 2. Public-data supply of strict hard negatives.

(A) Minimum of the two strict hard-negative counts (A_only, B_only) for every target pair in the J0 ChEMBL audit (`j0_strict_label_supply.csv`). Dashed line, thick-panel gate (≥50); dotted line, thin-panel gate (≥20). Highlighted: the three thick pairs used as K=4 main panels, EGFR/HER2 (7 B_only; supply-limited case), and HDAC1/HDAC6 (metal enzyme; excluded). (B) Count-level comparison of the same four pairs in ChEMBL pChEMBL versus BindingDB equal-relation measurements (Table S12). No docking.

## Figure 3. Benchmark formulation changes the apparent evidence for dual-target recognition.

Same frozen AutoDock Vina scores under two task formulations (unified θ = 6.0). Dark bars: directional pocket-matched `summary_min` with 95% ligand-bootstrap CIs from `unified_threshold_sensitivity_v2.csv` (Table 2). Orange bars: Dual-versus-neither comparator using pooled `vina_mean`, with CIs from `formulation_conventional_vs_directional_v1.csv` (Table 3). Dual versus neither is a nonselectivity-controlled comparator, not “the conventional dual-target benchmark.” EGFR/HER2 is the proof-of-principle gap (0.756 versus 0.430). AChE/BChE and PIK3CA/PIK3CB increments are small and overlapping. PIK3CA/mTOR Dual versus neither is hatched as underpowered (neither n = 4) and is not a reverse-overestimation result. Dual versus all non-duals is reported in Table 3, not here. Vertical dashed line, chance (0.5). The former pocket-matched forest (Vina / RTM / GNINA / descriptor) is Figure S4.

## Figure 4. Weak-arm asymmetry and physicochemical confounding.

(A) Directional Vina AUROCs at θ = 6.0: dual versus A_only (pocket B) and dual versus B_only (pocket A). (B) Vina pocket-matched summary_min versus the best single-descriptor reference, with 95% CIs. (C) TPSA on the AChE/BChE panel by class (individual ligands from `assembled_AChE_BChE.csv`; horizontal line, median). Dual ligands are more polar than either hard-negative class, matching the TPSA reference that exceeds Vina on this pair.

## Figure 5. Receptor realization can raise or lower apparent dual-target discrimination.

(A) PIK3CA/mTOR (PM48): pocket-matched summary_min after replacing one receptor at a time. 4L23 is the original PIK3CA structure (B = 4JT6 frozen); 4JPS and 5DXT replace PIK3CA; 4JSX replaces mTOR (A = 4L23 frozen). Point estimates: 0.692 (4L23), 0.486 (4JPS), 0.505 (5DXT), 0.639 (4JSX). (B) PIK3CA/PIK3CB: the same PIK3CA crystals with 2WXF held frozen. Point estimates: 0.500 (4L23), 0.691 (4JPS), 0.685 (5DXT). Error bars are 95% ligand-bootstrap CIs from the deposited swap tables and Table 2. 4JSX is an mTOR swap and is not applied to PIK3CA/PIK3CB. Receptor replacement is a realization effect, not a unidirectional collapse and not a robustness certificate. PIK3CA/PIK3CB uses the same 99-ligand set as Table 2 (PAB_034 timeout on original 4L23 and on both 4JPS and 5DXT). Unused-pool holdout is Figure S5.

## Figure 6. Wrong-pocket controls reveal an unresolved out-of-panel failure mode.

(A) Main K=4 panels: pocket-matched Vina summary_min versus the wrong-pocket control (`pocket_matched_directional_v1.csv`). Matched exceeds wrong-pocket on all four pairs. (B) Unused-pool holdout: the inequality reverses (wrong-pocket ≥ matched) on all three pairs with holdout supply (`holdout_pocket_matched_v1.csv`). EGFR/HER2 has no holdout. (C) The reversal remains after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) (`holdout_matched_wrong_pocket_summary_v1.csv`). Wrong-pocket remains ≥ matched on all nine cells; matching does not restore the main-panel inequality. (D) Scoring-free contact-count AUROC on pocket A (D vs A_only) and pocket B (D vs B_only) versus Vina wrong-pocket summary_min (`wrong_pocket_contact_v1_output.txt`; not a PLIF). B-arm contact is above chance; the magnitude does not reproduce Vina on PIK3CA/mTOR.

## Figure 7. Ligand-structure association and matched-subset tests.

(A) ECFP4 logistic regression under scaffold GroupKFold versus pocket-matched Vina on both directional arms (`ligand_ml_baseline_scaffold_cv_v1.csv`). Fingerprint AUROCs are chemotype–label association, not evidence of pocket physics. (B) Pocket-matched Vina versus all four prespecified physicochemical descriptors (heavy-atom count, MW, cLogP, TPSA) with 95% CIs. Descriptor CIs are from `forest_summary_min_ci_v1.csv`; Vina from θ = 6.0. Figure 4 reports only the best single-descriptor reference per pair; this panel shows all four. (C) Weak-arm (D vs B_only) logistic AUROC of Vina alone versus Vina plus heavy-atom count and TPSA, with the Vina odds ratio (`covariate_adjusted_v1.csv`). EGFR/HER2 score-only in that table is 0.5703 (the table’s logistic AUROC of feature `vina_A`), which is not the rank AUROC 0.4297 in Table 2. (D) D vs B_only after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) versus the unmatched full-panel contrast (`matched_subset_directional_v1.csv`). Error bars are the table’s single-contrast 95% CIs.

## Figure S1. Protocol knobs that do not change the ranking.

(A) Pocket-matched summary_min across the unified label-threshold grid (`unified_threshold_sensitivity_v2.csv`). Open markers are underpowered cells (EGFR/HER2 strict, n_B_only=7; PIK3CA/mTOR at θ=5.5, n_B_only=5, and strict, n_B_only=4). At the primary θ=6.0, PIK3CA/mTOR is the highest point estimate; AChE/BChE is flat at 0.6058 across the grid. The underpowered θ=5.5 PIK3CA/mTOR cell (0.5017) is not a ranking contradiction. (B) GNINA CNN mode01 versus best-of-9 versus the same-panel Vina reference (`gnina_pocket_matched_mode01_vs_best9_k4_v1.csv`). Best-of-9 versus mode01 moves summary_min by −0.04 to +0.08. EGFR/HER2 and AChE/BChE remain below chance on both GNINA channels. PIK3CA/PIK3CB GNINA best-of-9 is 0.533 versus Vina 0.500 (near chance). (C) PIK3CA/mTOR PM48 versus the PM110 expansion for Vina, RTMScore, and GNINA best-of-9 (`pm110_vs_pm48_pocket_matched_v1.csv`). (D) PM48 Vina at exhaustiveness 16 versus 8, computed from `scores_vina_E8_best.csv` (empty affinities skipped; ligands labeled neither were excluded) with the same pocket-matched definition, beside single-target enrichment AUROC and EF1% on 4L23 and 4JT6 (`single_target_enrichment_v1.csv`).

## Figure S2. Equal-relation supply and holdout sampling shift.

(A) Minimum strict hard-negative counts for the K=4 pairs in ChEMBL pChEMBL, BindingDB/PubChem `equal_only`, and BindingDB/PubChem `as_is` (`crossdb_strict_supply_v1.csv`). Count-level only; no docking. `as_is` lets EGFR/HER2 pass ≥50 because censored `>` values are treated as point estimates; `equal_only` does not. (B) Holdout minus main-panel mean pChEMBL for dual pA, A_only pA, and B_only pB (`holdout_vs_main_potency_size_v1.csv`). Sampling shift is real, especially on PIK3CA/mTOR, but does not reverse Figure 6C.

## Figure S3. Paired bootstrap differences that Figure 6 does not show.

All values are from `wrong_pocket_paired_delta_bootstrap_v1.csv` and `pocket_matched_vs_best_descriptor_delta_v1.csv` (B = 2000 ligand resamples, seed 20260729). Point Δ equals the rounded Table 2 / Figure 6 AUROCs subtracted at four decimals, not a separately rounded difference. Blue, 95% CI excludes 0; gray, CI includes 0. (A) Main K=4 panels: Δ = pocket-matched − wrong-pocket summary_min. Point Δ is positive on all four pairs (EGFR/HER2 0.1697, AChE/BChE 0.1614, PIK3CA/PIK3CB 0.1511, PIK3CA/mTOR 0.0902). Only EGFR/HER2 and AChE/BChE have CIs that exclude 0; PIK3CA/PIK3CB and PIK3CA/mTOR CIs include 0. (B) Unused-pool holdout: point Δ is negative on all three eligible pairs (wrong-pocket ≥ matched), and every CI includes 0. EGFR/HER2 has no holdout. This panel is the interval on the Figure 6B reversal, not a new docking experiment. (C) Pocket-matched Vina minus the best single-descriptor reference (EGFR/HER2 cLogP 0.4821; AChE/BChE TPSA 0.7333; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count). All four CIs include 0, including PIK3CA/mTOR +0.2291 [−0.0105, 0.4352]. This is not the pooled `vina_mean` gate (EGFR/HER2 0.2824). (D) ECFP4 logistic AUROC under scaffold GroupKFold versus random StratifiedKFold (`ligand_ml_scaffold_vs_random_v1.csv`). Mean (random − scaffold) across eight directional contrasts is 0.0112. Scaffold split remains the primary ML readout; this is a leakage check, not a search for a leakier split.

## Figure S4. Pocket-matched summary_min on the frozen K=4 set (former main Figure 3).

Vina (primary), RTMScore, GNINA CNN best-of-9, and the best single-descriptor reference (heavy-atom count, MW, cLogP, or TPSA) with 95% ligand-bootstrap CIs. Vina CIs are the θ = 6.0 values in `unified_threshold_sensitivity_v2.csv` (Table 2). Best descriptor (right column, from `forest_summary_min_ci_v1.csv`): EGFR/HER2 cLogP; AChE/BChE TPSA; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count. Vertical dashed line, chance (0.5). GNINA is a single CNN channel, not a three-engine competition. The main-text formulation comparison is Figure 3.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched summary_min on the main panel versus the unused-pool holdout (20/20/20; seed 20260731) for the three pairs with unused-pool supply. EGFR/HER2 has no holdout. PM110 is a same-family stability check shown in Figure S1C, not a third independent validation trajectory. This panel was formerly Figure 5A; the main-text Figure 5 is receptor realization only.

## TOC graphic (For Table of Contents Only).

DualFourClass-Bench asks whether docking can distinguish experimentally labeled dual-active ligands from single-target selective hard negatives in both pockets, rather than whether both docking scores are merely favorable. The graphic does not report numerical AUROCs and is not a reuse of Figure 1.
