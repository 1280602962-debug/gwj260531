# JNK1/2/3 Model Comparison Report (v2 — Improved Pipeline)

## Data Curation Improvements

- Biochemical assays only (`Assay Type = B`)
- Exact IC50 (`Standard Relation = =`)
- pActivity range [4, 10]
- Remove conflicting multi-assay measurements (std > 0.5 or range > 1.0 log)
- Assay harmonization: keep assays with ≥ 10 compounds
- **Per-isoform single-target models** (not sparse multitask table)

## Dataset Summary

| Isoform | Compounds | Train | Val | Test |
|---------|-----------|-------|-----|------|
| JNK1 | — | 384 | 29 | 31 |
| JNK2 | — | 477 | 66 | 67 |
| JNK3 | — | 966 | 83 | 98 |

## 5-Fold Scaffold CV (Primary Metric)

### XGBoost

| Isoform | Mean R² | Std | Mean Spearman | Fold R² |
|---------|---------|-----|---------------|---------|
| JNK1 | 0.662 | 0.086 | 0.772 | 0.579, 0.823, 0.631, 0.607, 0.668 |
| JNK2 | 0.443 | 0.074 | 0.678 | 0.457, 0.383, 0.354, 0.449, 0.569 |
| JNK3 | 0.633 | 0.089 | 0.793 | 0.502, 0.559, 0.654, 0.732, 0.717 |
| **Mean** | **0.579** | — | **0.748** | — |

### Chemprop 2.0

| Isoform | Mean R² | Std | Mean Spearman | Fold R² |
|---------|---------|-----|---------------|---------|
| JNK1 | nan | nan | nan |  |
| JNK2 | nan | nan | nan |  |
| JNK3 | nan | nan | nan |  |
| **Mean** | **nan** | — | **nan** | — |


## Holdout Test (Scaffold Split 80/10/10)

### XGBoost

| Isoform | R² | RMSE | Spearman | n |
|---------|-----|------|----------|---|
| JNK1 | 0.697 | 0.626 | 0.858 | 31 |
| JNK2 | 0.574 | 0.639 | 0.780 | 67 |
| JNK3 | 0.774 | 0.711 | 0.869 | 98 |
| **Mean** | **0.682** | — | **0.836** | — |

### Chemprop 2.0

| Isoform | R² | RMSE | Spearman | n |
|---------|-----|------|----------|---|
| JNK1 | nan | nan | nan | 0 |
| JNK2 | nan | nan | nan | 0 |
| JNK3 | nan | nan | nan | 0 |
| **Mean** | **nan** | — | **nan** | — |

## Model Selection

**Winner: XGBoost**

Selection based on 5-fold scaffold CV mean R², then holdout R².

## Notes

- Scaffold CV is the recommended metric for kinase QSAR (avoids inflated random-split R²).
- JNK2 has fewer compounds and more assay heterogeneity; expect lower R² than JNK1/JNK3.
- For selectivity + SHAP, use XGBoost selective models (scripts 04/05) regardless of activity model winner.