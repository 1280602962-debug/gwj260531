# Preliminary Data Analysis Summary

Generated from `docs/JNK1.csv`, `JNK2.csv`, `JNK3.csv`.

## Dataset Size (after curation)

| Isoform | Unique compounds | Median pActivity | Active fraction (pAct ≥ 6.5) |
|---------|------------------|------------------|------------------------------|
| JNK1 | 1576 | 6.46 | 49.4% |
| JNK2 | 803 | 6.50 | 50.6% |
| JNK3 | 1492 | 6.42 | 47.2% |

- **Paired multitask molecules**: 2562 total SMILES in merged table
- **Molecules with ≥2 isoform labels**: 950
- **Scaffold split**: train 2063 / val 276 / test 223

## Chemical Space Similarity

- Mean cross-dataset Tanimoto (Morgan FP): **~0.14–0.16** (moderate overlap)
- Shared Murcko scaffolds (all three isoforms): **221**
- JNK1–JNK2 scaffold Jaccard: **0.39**
- JNK1–JNK3 scaffold Jaccard: **0.29**

## Activity Distribution

- JNK1 vs JNK2: KS test p = 0.077 (not significantly different)
- JNK1 vs JNK3: KS test p < 1e-9 (**significantly different**)

## Model Comparison (Scaffold Test Set)

| Model | Mean R² | Mean RMSE | Mean Spearman |
|-------|---------|-----------|---------------|
| **XGBoost MTL** | **0.461** | **0.790** | **0.684** |
| Chemprop 2.0 MTL | 0.423 | 0.817 | 0.661 |

### Per-target (R² / Spearman)

| Target | XGBoost R² | Chemprop R² | XGBoost ρ | Chemprop ρ |
|--------|------------|-------------|-----------|------------|
| JNK1 | **0.566** | 0.500 | **0.741** | 0.691 |
| JNK2 | **0.442** | 0.360 | **0.665** | 0.624 |
| JNK3 | 0.376 | **0.408** | 0.648 | **0.668** |

## Recommendation

**Use XGBoost MTL as the primary activity predictor** for the next screening stage.

Reasons:
1. Higher mean R² and Spearman on held-out scaffold test set
2. Best on JNK1 (the target of interest for selectivity)
3. Better integration with downstream SHAP interpretability

Keep Chemprop as optional benchmark; consider ensemble if both predictions are needed.

See full report: [MODEL_COMPARISON_REPORT.md](model_comparison/MODEL_COMPARISON_REPORT.md)
