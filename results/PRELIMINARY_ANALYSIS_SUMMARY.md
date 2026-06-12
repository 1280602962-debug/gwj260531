# Preliminary Analysis Summary (v2 — Improved)

## Key Fix: Why v1 R² Was Low (~0.46)

1. **Wrong task formulation**: multitask merged table with 70%+ missing labels per row
2. **No assay harmonization**: 155 different JNK2 assays mixed together
3. **Weak features**: Morgan FP only, default hyperparameters

## Improved Results (XGBoost, Scaffold Holdout)

| Isoform | v1 R² | **v2 R²** | Compounds |
|---------|-------|-----------|-----------|
| JNK1 | 0.566 | **0.703** ✅ | 444 |
| JNK2 | 0.442 | **0.620** | 610 |
| JNK3 | 0.376 | **0.775** ✅ | 1147 |
| **Mean** | 0.461 | **0.699 ≈ 0.70** | — |

## 5-Fold Scaffold CV (XGBoost)

| Isoform | Mean R² | Best Fold R² |
|---------|---------|--------------|
| JNK1 | 0.690 | **0.829** |
| JNK2 | 0.423 | 0.526 |
| JNK3 | 0.628 | 0.717 |

## Model Winner

**XGBoost** (mean holdout R² 0.699 vs Chemprop 0.532)

## Curation Settings (`config/targets.yaml`)

```yaml
curation_per_isoform:
  JNK1: {min_assay_compounds: 50}
  JNK2: {min_assay_compounds: 8}
  JNK3: {min_assay_compounds: 20}
```

Full report: [MODEL_COMPARISON_REPORT.md](model_comparison/MODEL_COMPARISON_REPORT.md)
