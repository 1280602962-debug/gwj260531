# M3 — Trivial baselines vs docking arms

## Rule (frozen)

- Mandatory baselines: `heavy_atoms`, `MW`, `cLogP`, `TPSA` (plus optional `morgan_dual_medsim` = LOO median Tanimoto to other duals).
- Arm comparison uses **summary = min(D/A, D/B)**.
- Docking arm **fail_baseline** if `summary_min ≤ best_volume_summary_min` on that pair.
- Report pairs **separately**; do not average EGFR with PIK3CA/mTOR.

## Results

| pair        | arm                | family   |   auroc_D_vs_A |   auroc_D_vs_B |   summary_min | best_volume_arm   |   best_volume_summary_min | fail_baseline   |
|:------------|:-------------------|:---------|---------------:|---------------:|--------------:|:------------------|--------------------------:|:----------------|
| EGFR_HER2   | vina_mean          | docking  |         0.6889 |         0.3114 |        0.3114 | cLogP             |                    0.4821 | True            |
| EGFR_HER2   | vina_min           | docking  |         0.6711 |         0.2712 |        0.2712 | cLogP             |                    0.4821 | True            |
| EGFR_HER2   | rtm_mean           | docking  |         0.5686 |         0.3203 |        0.3203 | cLogP             |                    0.4821 | True            |
| EGFR_HER2   | rtm_min            | docking  |         0.5912 |         0.3237 |        0.3237 | cLogP             |                    0.4821 | True            |
| EGFR_HER2   | rtm_min_z          | docking  |         0.5874 |         0.317  |        0.317  | cLogP             |                    0.4821 | True            |
| EGFR_HER2   | heavy_atoms        | baseline |         0.6997 |         0.3694 |        0.3694 | cLogP             |                    0.4821 |                 |
| EGFR_HER2   | MW                 | baseline |         0.6476 |         0.4163 |        0.4163 | cLogP             |                    0.4821 |                 |
| EGFR_HER2   | cLogP              | baseline |         0.5263 |         0.4821 |        0.4821 | cLogP             |                    0.4821 |                 |
| EGFR_HER2   | TPSA               | baseline |         0.703  |         0.4275 |        0.4275 | cLogP             |                    0.4821 |                 |
| EGFR_HER2   | morgan_dual_medsim | baseline |         0.555  |         0.5848 |        0.555  | cLogP             |                    0.4821 |                 |
| PIK3CA_mTOR | vina_mean          | docking  |         0.6984 |         0.5972 |        0.5972 | heavy_atoms       |                    0.463  | False           |
| PIK3CA_mTOR | vina_min           | docking  |         0.6825 |         0.588  |        0.588  | heavy_atoms       |                    0.463  | False           |
| PIK3CA_mTOR | rtm_mean           | docking  |         0.619  |         0.7639 |        0.619  | heavy_atoms       |                    0.463  | False           |
| PIK3CA_mTOR | rtm_min            | docking  |         0.623  |         0.7546 |        0.623  | heavy_atoms       |                    0.463  | False           |
| PIK3CA_mTOR | rtm_min_z          | docking  |         0.6111 |         0.7917 |        0.6111 | heavy_atoms       |                    0.463  | False           |
| PIK3CA_mTOR | heavy_atoms        | baseline |         0.4643 |         0.463  |        0.463  | heavy_atoms       |                    0.463  |                 |
| PIK3CA_mTOR | MW                 | baseline |         0.4484 |         0.4722 |        0.4484 | heavy_atoms       |                    0.463  |                 |
| PIK3CA_mTOR | cLogP              | baseline |         0.5198 |         0.3102 |        0.3102 | heavy_atoms       |                    0.463  |                 |
| PIK3CA_mTOR | TPSA               | baseline |         0.2599 |         0.4514 |        0.2599 | heavy_atoms       |                    0.463  |                 |
| PIK3CA_mTOR | morgan_dual_medsim | baseline |         0.6389 |         0.6944 |        0.6389 | heavy_atoms       |                    0.463  |                 |

## Gate (per pair)

- **EGFR_HER2: M3 = No-Go**
- **PIK3CA_mTOR: M3 = Go**

Overall note: EGFR expected No-Go (volume ≥ docking on weak end); PIK3CA/mTOR expected Go if docking exceeds volume on both directions.
