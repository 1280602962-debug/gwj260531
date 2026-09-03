# Detectable-effect simulation v1

Zero docking. Binormal scores; ligand-level class-preserving bootstrap with fixed class sizes as part of the simulation design.
N_MC = 1000; N_BOOT = 2000; seed = 20260729.

This is **not** observed power on the empirical AUROCs.

## Probability that the 95% CI excludes 0.5 (`summary_min`)

| Pair | n_scored | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |
|---|---:|---:|---:|---:|---:|---:|
| EGFR/HER2 | 28/38/32 | 0.025 | 0.065 | 0.268 | 0.621 | 0.907 |
| AChE/BChE | 27/25/28 | 0.020 | 0.049 | 0.225 | 0.504 | 0.828 |
| PIK3CA/PIK3CB | 28/27/28 | 0.032 | 0.041 | 0.226 | 0.564 | 0.849 |
| PIK3CA/mTOR | 18/14/12 | 0.037 | 0.025 | 0.072 | 0.219 | 0.452 |

## Interpretation freeze

- Current class sizes resolve **large** directional effects more readily than moderate ones.
- Failure of an observed CI to exclude 0.5 does **not** establish equivalence to chance.
- Dual versus neither uses a smaller negative set than the directional B-only/A-only arms on some pairs;
  detectable-effect probabilities are therefore not interchangeable across formulations.

