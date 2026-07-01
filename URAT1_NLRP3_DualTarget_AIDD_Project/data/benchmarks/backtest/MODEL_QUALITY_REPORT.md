# Model Quality & Benchmark Backtest Report

**Overall verdict: URAT1_NO_GO**

NLRP3 model is screening-ready; URAT1 model fails strict CV and/or benchmark recovery. URAT1 library filtering must NOT rely on ML alone — use $S_trap$ conformational ensemble docking as primary evidence.

## 1. Cross-validation (scaffold GroupKFold, 5 folds)

### URAT1 regression + conformal UQ
- RMSE (OOF): 0.663
- R² (OOF): 0.508
- Spearman (OOF): 0.726
- ROC-AUC (p≥7): 0.852
- EF@5% (p≥7, strong actives): 3.29
- EF@10% (p≥6): 1.71 — **misleading** (theoretical max ≈1.75 at 57% base rate)
- Strict CV pass: True

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

## 3. Why the previous URAT1 table was misleading

1. **EF@10% at p≥6 is capped near 1.75** when 57% of training compounds are already actives — even a perfect ranker cannot exceed ~1.75.
2. **Thresholds were too lenient** (R²≥0.25, EF≥1.5), allowing a mediocre model to show all green checks.
3. **Fold-averaged R² (0.44) understates OOF R² (0.51)** but both are only moderate for prospective screening.
4. **Benchmark backtest contradicts** the pass table: lesinurad/dotinurad fail despite CV pass.

## 4. Interpretation notes

- **lesinurad / benzbromarone** were dropped during ChEMBL curation due to >1 log assay conflict; ChEMBL median pActivity (~5.1–6.5) is lower than literature references used in benchmarks.
- **verinurad** is in the training set; model prediction is consistent with held-in data.
- **MCC950** is in NLRP3 training data; high P(active) confirms correct class assignment.
- For scaffold-novel benchmarks, prioritize **conformational ensemble docking** ($S_{trap}$) over ML rank.
