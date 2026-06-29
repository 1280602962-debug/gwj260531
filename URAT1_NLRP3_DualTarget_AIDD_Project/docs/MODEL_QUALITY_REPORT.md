# Model Quality & Benchmark Backtest Report

**Overall verdict: CONDITIONAL_GO**

Internal CV is strong; benchmark recovery is partial. Use ML as coarse filter; rely more on docking for URAT1 novel scaffolds.

## 1. Cross-validation (scaffold GroupKFold, 5 folds)

### URAT1 regression + conformal UQ
- RMSE: 0.661
- R²: 0.442
- Spearman: 0.650
- EF@10%: 1.82
- CV screening suitable: True

### NLRP3 assay-conditioned classifier
- AUROC: 0.893
- AUPRC: 0.914
- EF@10%: 1.57
- CV screening suitable: True

## 2. Benchmark backtest

### URAT1 (predicted pActivity ≥ 6 = pass)

| Compound | In train | Max Tc | Pred | Lit. | Pass |
|----------|----------|--------|------|------|------|
| lesinurad | False | 0.689 | 5.642 | 7.0 | False |
| benzbromarone | False | 0.75 | 6.559 | 7.5 | True |
| verinurad | True | 1.0 | 7.083 | 8.0 | True |
| dotinurad | False | 0.837 | 5.071 | 8.2 | False |
| allopurinol | False | 0.22 | 5.416 | — | True |

URAT1 must-recover binary pass: 2/4
- In training set: 1/1
- Scaffold-novel (excluded from curation): 1/3

### NLRP3 (P(active) ≥ 0.5 = pass)

| Compound | In train | Max Tc | P(active) | Pass |
|----------|----------|--------|-----------|------|
| MCC950 | True | 1.0 | 1.0 | True |
| GDC-2394 | True | 1.0 | 0.917 | True |

NLRP3 must-recover binary pass: 2/2

## 3. Interpretation notes

- **lesinurad / benzbromarone** were dropped during ChEMBL curation due to >1 log assay conflict; ChEMBL median pActivity (~5.1–6.5) is lower than literature references used in benchmarks.
- **verinurad** is in the training set; model prediction is consistent with held-in data.
- **MCC950** is in NLRP3 training data; high P(active) confirms correct class assignment.
- For scaffold-novel benchmarks, prioritize **conformational ensemble docking** ($S_{trap}$) over ML rank.
