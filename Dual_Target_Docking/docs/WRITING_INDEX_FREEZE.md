# DualFourClass writing index (freeze)

Source of numbers: `results/canonical/` after promotion from `results/freeze_rebuild/`.  
Config: `scripts/analysis/analysis_config.py`.  
Do not quote `results/canonical` siblings that predate the freeze ENV.

## What may be claimed

1. **Eight pairs, four experimental states, θ=6.0.** Primary n is activity-eligible complete-case. Table 1 construction quotas can still show EGFR 28/38/32/12; Table 2 n_scored is 28/37/32.
2. **Score channel.** Vina mode-1 stored as higher-better (`score_S = −energy`). Dual vs A-only uses pocket B; dual vs B-only uses pocket A. `summary_min` is the weaker of those two arms. Two-pocket mean is only for dual-vs-neither and ranking/EF.
3. **Scheme B.** Class-stratified percentile bootstrap, B=2000, seed=20260729. Point estimates are on the original sample.
4. **EGFR/HER2 directional.** D vs A (pocket B) 0.656 [0.520, 0.791]; D vs B (pocket A) = `summary_min` 0.324 [0.195, 0.474]. Fixed-score Δ on the EGFR pocket: 0.462 [0.260, 0.641].
5. **Control-class dependence is pair-specific.** JAK1/TYK2 Δ 0.444 [0.261, 0.630]. Other pairs smaller or intervals include 0.
6. **Matched−mismatched.** Only AChE/BChE main-panel 95% CI excludes 0. EGFR includes 0. Seven holdouts include 0. Not multiplicity-adjusted.
7. **ECFP incremental.** Freeze max |Δ| = 0.0234 at PPARA/PPARD D vs B (three-decimal claim 0.023). Figure 3A/B and `plotted_values.json` `fig3B_max_abs` are redrawn from `results/canonical/ecfp4_incremental_information.csv`. Files that still store 0.0112 are historical leftovers (`docs/HISTORICAL_LEFTOVER_FILES.md`).
8. **Max vs median** is a sensitivity on the qualified-record intersection, not Table 2. AChE 1/94 (CHEMBL659); EGFR 5/109. AB_056 stays in Table 2.
9. **GNINA** is another pose-generation realization on deposited poses, not a Vina replacement. EGFR `summary_min` 0.265 [0.142, 0.402], n=108.
10. **Cognate RMSD** from `all14_cognate_rmsd_calcrrms_v1.csv`. EGFR 3POZ CalcRMS top-1 = best = 1.019 Å. Historical 9.505 Å is not current.

## What must not be claimed

- Pre-fix EGFR AUROCs 0.430 / 0.808 / 0.378 as current results.
- PR #35 reconstructed EGFR 0.760 / 9.505 Å as current SI Table S2.
- Five-seed archive D vs A 0.6607 (complete-case 28/38/32) as Table 2.
- JAK1/TYK2 document-cluster CI (status `unresolved_mapping_unavailable`).
- Detectable-effect probabilities as observed power.
- PIK3CA/PIK3CB as a current primary pair.
- That dual-vs-neither discrimination implies discrimination against single-target selectives.

## Files to write from

| Claim | File |
|---|---|
| Table 2 | `primary_summary_min.csv`, `primary_directional_auroc.csv`, `class_counts.csv` |
| Table 3 / ranking | `two_pocket_mean_ranking.csv`, `top10_operating_points.csv` |
| Fixed-score Δ | `fixed_score_negative_class_delta.csv` |
| ECFP incremental | `ecfp4_incremental_information.csv` |
| Matched−mismatched | `matched_minus_mismatched.csv`, `holdout_metrics.csv` |
| Descriptors | `descriptor_baselines.csv` |
| GNINA / receptor swap | `computational_robustness.csv`, `receptor_substitution.csv` |
| RMSD | `cognate_rmsd.csv` |
| Detectable-effect | `detectable_effect_simulation.csv` |
| Max/median | `max_vs_median_sensitivity.csv` |
