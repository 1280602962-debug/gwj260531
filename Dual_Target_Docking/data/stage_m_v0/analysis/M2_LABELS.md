# M2 — Label validity

## Rules (frozen)

- `dual_strict`: both pChEMBL ≥ 6.5
- `A_only_strict`: A≥6.5 and B≤5.5 (measured)
- `B_only_strict`: B≥6.5 and A≤5.5
- `neither_strict`: both ≤5.5
- `gray`: both measured, else → excluded from strict directional analysis
- Unmeasured end ≠ negative (incomplete excluded from re-binarization)

## Margin counts

| pair        | label          |   n |   n_both_measured |   frac_of_panel |
|:------------|:---------------|----:|------------------:|----------------:|
| EGFR_HER2   | dual_strict    |  26 |               110 |           0.236 |
| EGFR_HER2   | A_only_strict  |  17 |               110 |           0.155 |
| EGFR_HER2   | B_only_strict  |   7 |               110 |           0.064 |
| EGFR_HER2   | neither_strict |   9 |               110 |           0.082 |
| EGFR_HER2   | gray           |  51 |               110 |           0.464 |
| EGFR_HER2   | incomplete     |   0 |               110 |           0     |
| PIK3CA_mTOR | dual_strict    |  17 |                48 |           0.354 |
| PIK3CA_mTOR | A_only_strict  |   7 |                48 |           0.146 |
| PIK3CA_mTOR | B_only_strict  |   4 |                48 |           0.083 |
| PIK3CA_mTOR | neither_strict |   2 |                48 |           0.042 |
| PIK3CA_mTOR | gray           |  18 |                48 |           0.375 |
| PIK3CA_mTOR | incomplete     |   0 |                48 |           0     |

## Directional on strict margin

| pair        | arm         |   n_dual |   n_A_only |   n_B_only |   auroc_D_vs_A |   auroc_D_vs_B | underpowered   |
|:------------|:------------|---------:|-----------:|-----------:|---------------:|---------------:|:---------------|
| EGFR_HER2   | vina_mean   |       26 |         17 |          7 |         0.8122 |         0.1264 | True           |
| EGFR_HER2   | rtm_min_z   |       26 |         17 |          7 |         0.5633 |         0.2363 | True           |
| EGFR_HER2   | heavy_atoms |       26 |         17 |          7 |         0.7285 |         0.228  | True           |
| PIK3CA_mTOR | vina_mean   |       17 |          7 |          4 |         0.605  |         0.6324 | True           |
| PIK3CA_mTOR | rtm_min_z   |       17 |          7 |          4 |         0.5966 |         0.7941 | True           |
| PIK3CA_mTOR | heavy_atoms |       17 |          7 |          4 |         0.4538 |         0.4265 | True           |

## Threshold sensitivity (vina_mean)

| pair        |   cutoff |   n_dual |   n_A_only |   n_B_only |   auroc_D_vs_A |   auroc_D_vs_B |
|:------------|---------:|---------:|-----------:|-----------:|---------------:|---------------:|
| EGFR_HER2   |      5.5 |       69 |         22 |         10 |         0.7885 |         0.3609 |
| EGFR_HER2   |      6   |       28 |         38 |         32 |         0.6889 |         0.3114 |
| EGFR_HER2   |      6.5 |       26 |         29 |         29 |         0.7533 |         0.3249 |
| PIK3CA_mTOR |      5.5 |       33 |          9 |          5 |         0.4949 |         0.5636 |
| PIK3CA_mTOR |      6   |       18 |         14 |         12 |         0.6984 |         0.5972 |
| PIK3CA_mTOR |      6.5 |       17 |         15 |         12 |         0.6745 |         0.5931 |

## Noise / oracle ceiling (σ=0.5 highlight)

| pair        | score                 |   median_D_vs_A |   median_D_vs_B |   median_min_DA_DB |
|:------------|:----------------------|----------------:|----------------:|-------------------:|
| EGFR_HER2   | vina_mean_true_scores |          0.7218 |          0.3451 |             0.3451 |
| EGFR_HER2   | oracle_min_pchembl    |          0.9415 |          0.9055 |             0.9053 |
| PIK3CA_mTOR | vina_mean_true_scores |          0.6519 |          0.5909 |             0.5838 |
| PIK3CA_mTOR | oracle_min_pchembl    |          0.9568 |          0.9335 |             0.9308 |

## Gate

**M2 = Weak**

- Margin same-sign vs main panel (vina): True
- Oracle distinguishable @σ=0.5 (≥1 direction median≥0.65): True
- Details: EGFR_HER2: oracle@σ0.5 D/A=0.942 D/B=0.905; PIK3CA_mTOR: oracle@σ0.5 D/A=0.957 D/B=0.933
- Cutoff flips D/B sign: False
- Margin underpowered flag: True
