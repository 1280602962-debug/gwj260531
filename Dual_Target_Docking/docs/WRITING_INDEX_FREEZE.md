# DualFourClass writing index (freeze)

Source of numbers: `results/canonical/` after promotion from a verified `/tmp` freeze rebuild (repository `results/freeze_rebuild/` is historical / unavailable and must be absent).  
Config: `scripts/analysis/analysis_config.py`.  
Do not quote `results/canonical` siblings that predate the freeze ENV.

## What may be claimed

1. **Eight pairs, four experimental states, θ=6.0.** Primary n is activity-eligible complete-case. Table 1 construction quotas can still show EGFR 28/38/32/12; Table 2 n_scored is 28/37/31.
2. **Score channel.** Vina mode-1 stored as higher-better (`score_S = −energy`). Dual vs A-only uses pocket B; dual vs B-only uses pocket A. `summary_min` is the weaker of those two arms. Two-pocket mean is only for dual-vs-neither and ranking/EF.
3. **Scheme B.** Class-stratified percentile bootstrap, B=2000, seed=20260729. Point estimates are on the original sample.
4. **EGFR/HER2 directional.** D vs A (pocket B) 0.663 [0.521, 0.792]; D vs B (pocket A) = `summary_min` 0.334 [0.197, 0.471]. Fixed-score Δ on the EGFR pocket: 0.446 [0.259, 0.632].
5. **Control-class dependence is pair-specific.** JAK1/TYK2 Δ 0.444 [0.261, 0.630]. Other pairs smaller or intervals include 0.
6. **Matched−mismatched.** AChE/BChE and EGFR/HER2 main-panel 95% CIs exclude 0. Seven holdouts include 0. Not multiplicity-adjusted.
7. **ECFP incremental.** Freeze max |Δ| = 0.0234 at PPARA/PPARD D vs B (three-decimal claim 0.023). Figure 3A/B is redrawn from `results/canonical/ecfp4_incremental_information.csv`. Sole plotted-value authority: `figures/jcim_article/plotted_values_postfix.json`. The old 0.0112 value is a superseded historical-environment result; background is the superseded snapshot `docs/archive/PR39_SCIENTIFIC_DATA_AUDIT.md` (not current truth). Current authority: `results/canonical/ecfp4_incremental_information.csv` and `docs/audit/FINAL_FULL_PROJECT_AUDIT.md`.
8. **Max vs median** is a sensitivity on the qualified-record intersection, not Table 2. AChE 1/94 (CHEMBL659); EGFR 5/109. AB_056 stays in Table 2.
9. **GNINA** is another pose-generation realization on deposited poses, not a Vina replacement. EGFR `summary_min` 0.227 [0.104, 0.373], n=89.
10. **Cognate RMSD** from `all14_cognate_rmsd_calcrrms_v1.csv`. EGFR 3POZ CalcRMS top-1 = best = 1.019 Å. Historical 9.505 Å is not current.

## What must not be claimed

- Pre-fix EGFR AUROCs 0.430 / 0.808 / 0.378 as current results.
- PR #35 reconstructed EGFR 0.760 / 9.505 Å as current SI Table S2.
- Five-seed complete-case D vs A 0.6607 (historical 28/38/32) as Table 2. Current five-seed production D vs A is 0.6631 on n = 28/37/31.
- JAK1/TYK2 document-cluster CI (status `unresolved_mapping_unavailable`).
- Detectable-effect probabilities as observed power.
- PIK3CA/PIK3CB as a current primary pair.
- That deposited EGFR/HER2 production or corrected-box five-seed scores were prepared with RDKit/Meeko or LigPrep. Ligand PDBQT for those tables are not recoverable (`docs/EGFR_LIGAND_PREP_PROVENANCE_AUDIT.md`).
- That all 14 receptors were prepared with Meeko or Protein Preparation Wizard.
- That the repository can fully rerun all original docking from SMILES.
- That all eight pairs are assay-level adjudicated.
- That AChE five-seed available-case membership equals Table 2 (it does not; use fixed-membership n=88).
- That full-panel best-descriptor AUROC is a selection-adjusted predictive estimate.
- That checksum manifests are a scientific PASS/FAIL gate.
- Track-B `do_not_now` as a current instruction to skip analyses that were later completed as sensitivity.

## Files to write from

| Claim | File |
|---|---|
| Table 2 | `results/canonical/primary_summary_min.csv`, `results/canonical/primary_directional_auroc.csv`, `results/canonical/class_counts.csv` |
| Table 3 / ranking | `results/canonical/two_pocket_mean_ranking.csv`, `results/canonical/top10_operating_points.csv` |
| Fixed-score Δ | `results/canonical/fixed_score_negative_class_delta.csv` |
| ECFP incremental | `results/canonical/ecfp4_incremental_information.csv` |
| Matched−mismatched | `results/canonical/matched_minus_mismatched.csv`, `results/canonical/holdout_metrics.csv` |
| Descriptors | `results/canonical/descriptor_baselines.csv` (full-panel best = descriptive; nested OOF = predictive) |
| Nested descriptor CV | `results/canonical/descriptor_nested_scaffold_cv.csv` |
| Five-seed available-case | `results/canonical/five_seed_summary_min.csv` |
| Five-seed fixed membership | `results/canonical/five_seed_fixed_membership_sensitivity.csv` |
| Receptor registry | `data/provenance/receptor_input_registry.csv`; `docs/RECEPTOR_PREP_PROVENANCE.md` |
| Ligand input | `docs/LIGAND_INPUT_REPRODUCIBILITY.md` |
| Protocol levels | `docs/PROTOCOL_LEVELS.md`; `docs/TRACK_B_FINAL_EXECUTION_STATUS.md` |
| Labels | `docs/LABEL_PROVENANCE.md` |
| Path audit | `docs/CURRENT_PATH_AUDIT.md` |
| Current forensic audit | `docs/audit/FINAL_POST_REMEDIATION_AUDIT.md` (prior scan: `docs/audit/FINAL_FULL_PROJECT_AUDIT.md`) |
| Historical freeze memos (superseded; not current truth) | `docs/archive/PR39_SCIENTIFIC_DATA_AUDIT.md`, `docs/archive/PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md`, `docs/archive/SUBMISSION_ONLY_CLEANUP_REPORT.md` |
| GNINA / receptor swap | `results/canonical/computational_robustness.csv`, `results/canonical/receptor_substitution.csv` |
| RMSD | `results/canonical/cognate_rmsd.csv` |
| Detectable-effect | `results/canonical/detectable_effect_simulation.csv` |
| Max/median | `results/canonical/max_vs_median_sensitivity.csv` |
| Plotted values | `figures/jcim_article/plotted_values_postfix.json` |
