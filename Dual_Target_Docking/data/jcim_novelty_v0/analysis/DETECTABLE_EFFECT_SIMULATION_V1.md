# Detectable-effect simulation v1

Zero docking. Binormal scores; ligand-level class-preserving bootstrap with fixed class sizes as part of the simulation design.
N_MC = 1000; N_BOOT = 2000; seed = 20260729.

This is **not** observed power on the empirical AUROCs.

## Probability that the 95% CI excludes 0.5 (`summary_min`)

| Pair | n_scored | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |
|---|---:|---:|---:|---:|---:|---:|
| EGFR/HER2 | 28/38/32 | 0.025 | 0.065 | 0.268 | 0.621 | 0.907 |
| JAK1/JAK2 | 32/32/32 | 0.018 | 0.066 | 0.285 | 0.639 | 0.911 |
| JAK1/TYK2 | 31/32/32 | 0.033 | 0.053 | 0.280 | 0.623 | 0.904 |
| PIK3CA/mTOR | 18/14/12 | 0.037 | 0.027 | 0.081 | 0.233 | 0.453 |
| AChE/BChE | 27/26/28 | 0.024 | 0.056 | 0.202 | 0.540 | 0.817 |
| F2/F10 | 31/32/32 | 0.029 | 0.062 | 0.318 | 0.649 | 0.903 |
| PPARG/PPARA | 32/31/32 | 0.025 | 0.068 | 0.261 | 0.654 | 0.911 |
| PPARA/PPARD | 32/32/32 | 0.023 | 0.062 | 0.285 | 0.674 | 0.901 |

## Interpretation freeze

- Current class sizes resolve **large** directional effects more readily than moderate ones.
- Failure of an observed CI to exclude 0.5 does **not** establish equivalence to chance.
- Dual versus neither uses a smaller negative set than the directional B-only/A-only arms on some pairs;
  detectable-effect probabilities are therefore not interchangeable across formulations.
- Current eight-pair class sizes from `results/canonical/class_counts.csv`.
- This is a simulation-based detectable-effect analysis, not observed/post-hoc power.

