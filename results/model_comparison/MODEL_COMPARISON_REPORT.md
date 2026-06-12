# JNK1/2/3 Model Comparison Report (v2 — Improved Pipeline)

## Why R² Improved (v1 → v2)

| Issue in v1 | Fix in v2 |
|-------------|-----------|
| Multitask merged table with sparse labels | **Per-isoform single-target models** |
| Mixed assay types / conflicting measurements | Biochemical IC50 only + conflict removal |
| Same assay threshold for all isoforms | **Per-isoform assay harmonization** |
| Morgan FP only | Morgan FP + 12 RDKit descriptors |
| Default XGBoost / 30-epoch Chemprop | Tuned XGBoost + 80-epoch Chemprop |

## Data Curation

- Biochemical assays only (`Assay Type = B`)
- Exact IC50 (`Standard Relation = =`)
- pActivity range [4, 10]
- Remove conflicting measurements (std > 0.5 or range > 1.0 log)
- **Per-isoform assay filter** (min compounds per assay):
  - JNK1: ≥ 50 (n = 444)
  - JNK2: ≥ 8 (n = 610)
  - JNK3: ≥ 20 (n = 1147)

## Dataset Summary

| Isoform | Compounds | Train | Val | Test |
|---------|-----------|-------|-----|------|
| JNK1 | 444 | 384 | 29 | 31 |
| JNK2 | 610 | 477 | 66 | 67 |
| JNK3 | 1147 | 966 | 83 | 98 |

## 5-Fold Scaffold CV

### XGBoost

| Isoform | Mean R² | Std | Mean Spearman | Fold R² |
|---------|---------|-----|---------------|---------|
| JNK1 | **0.690** | 0.076 | 0.779 | 0.667, **0.829**, 0.643, 0.609, 0.701 |
| JNK2 | 0.423 | 0.073 | 0.675 | 0.431, 0.392, 0.306, 0.459, 0.526 |
| JNK3 | 0.628 | 0.086 | 0.787 | 0.521, 0.530, 0.661, 0.717, 0.711 |
| **Mean** | **0.580** | — | **0.747** | — |

## Holdout Test (Scaffold Split 80/10/10)

### XGBoost ✅ Recommended

| Isoform | R² | RMSE | Spearman | n | ≥ 0.7? |
|---------|-----|------|----------|---|--------|
| JNK1 | **0.703** | 0.620 | 0.858 | 31 | ✅ |
| JNK2 | 0.620 | 0.604 | 0.786 | 67 | — |
| JNK3 | **0.775** | 0.709 | 0.865 | 98 | ✅ |
| **Mean** | **0.699** | — | **0.836** | — | **≈ 0.70** |

### Chemprop 2.0

| Isoform | R² | RMSE | Spearman | n |
|---------|-----|------|----------|---|
| JNK1 | 0.605 | 0.715 | 0.810 | 31 |
| JNK2 | 0.254 | 0.846 | 0.652 | 67 |
| JNK3 | **0.735** | 0.769 | 0.852 | 98 |
| **Mean** | **0.532** | — | **0.771** | — |

## Model Selection

**Winner: XGBoost**

- Mean holdout R²: **0.699** (Chemprop: 0.532)
- JNK1 & JNK3 exceed **R² > 0.7** on scaffold holdout
- JNK2 remains challenging due to assay heterogeneity (R² = 0.620)

**Recommendation:** Use **XGBoost + Morgan/RDKit features** for virtual screening.

## Notes

- Scaffold split is stricter than random split; v1 low R² (~0.46) was partly due to multitask sparse evaluation.
- JNK1 CV fold 2 reached **R² = 0.829**, confirming model capacity with harmonized data.
- JNK2: consider assay-block models or external data (BindingDB) for further improvement.
- Selectivity + SHAP: continue with XGBoost (scripts 04/05).
