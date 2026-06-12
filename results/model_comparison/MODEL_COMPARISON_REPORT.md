# JNK1/2/3 Model Comparison Report

## Data Summary

| Metric | Value |
|--------|-------|
| JNK1_compounds | 1576 |
| JNK2_compounds | 803 |
| JNK3_compounds | 1492 |
| paired_total | 2562 |
| paired_ge2_isoforms | 950 |
| jnk1_selective | 0 |
| train_size | 2063 |
| val_size | 276 |
| test_size | 223 |

## Test Set Performance (Scaffold Split)

### XGBoost MTL

| Target | R² | RMSE | MAE | Spearman | n |
|--------|-----|------|-----|----------|---|
| pAct_JNK1 | 0.566 | 0.676 | 0.517 | 0.741 | 129 |
| pAct_JNK2 | 0.442 | 0.784 | 0.612 | 0.665 | 78 |
| pAct_JNK3 | 0.376 | 0.910 | 0.743 | 0.648 | 130 |
| **Mean** | **0.461** | **0.790** | — | **0.684** | 337 |

### Chemprop 2.0 MTL

| Target | R² | RMSE | MAE | Spearman | n |
|--------|-----|------|-----|----------|---|
| pAct_JNK1 | 0.500 | 0.726 | 0.581 | 0.691 | 129 |
| pAct_JNK2 | 0.360 | 0.840 | 0.653 | 0.624 | 78 |
| pAct_JNK3 | 0.408 | 0.886 | 0.693 | 0.668 | 130 |
| **Mean** | **0.423** | **0.817** | — | **0.661** | 337 |

## Model Selection

**Winner: XGBoost MTL**

Reason: Higher mean R² and/or Spearman on scaffold-test set

Recommendation: Use **XGBoost MTL** as primary activity predictor for virtual screening.

## Notes

- Both models trained on identical scaffold-based train/val/test splits.
- Missing JNK isoform labels handled natively (Chemprop mask; XGBoost per-task training).
- JNK1 has fewer data points; compare JNK1 task performance carefully.
- For selectivity modeling + SHAP, continue using XGBoost selective models (script 04/05).