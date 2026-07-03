# ML External Validation Report

## Design (anti-circular)

| Component | Source | Role |
|-----------|--------|------|
| Decoys | Taosu library (n=10,000) | Assumed inactive (label=0) |
| Benchmarks | literature_benchmarks.csv (n=9) | Known active (label=1) |
| ChEMBL actives | pActivity ≥ 6.0 (n=1210) | Known active (label=1) |

**Excluded from Taosu sample:** docked top-5000 + ChEMBL demo/training library (1835 paired compounds).

**NOT used:** Demo ChEMBL library circular validation (84% F1 pass rate).

## Metrics (full evaluation set)

| Metric | Value |
|--------|-------|
| N total | 11,219 |
| N actives | 1,219 |
| N decoys | 10,000 |
| ROC-AUC (p_family) | 0.8759 |
| EF1% | 9.20 |

## Threshold: p_family ≥ 6.0

| | Predicted active | Predicted inactive |
|--|------------------|-------------------|
| True active | TP=1211 | FN=8 |
| True decoy | FP=9528 | TN=472 |

- **Sensitivity (recall):** 99.3% — fraction of actives passing
- **Specificity:** 4.7% — fraction of decoys correctly rejected
- **Precision:** 11.3%
- **Decoy pass rate (FPR):** 95.3% (9528/10000 decoys pass)

## Benchmark-only (recall-only baseline)

- Benchmarks passing p_family ≥ 6.0: **9/9** (100%)
- This confirms recall only; specificity requires external decoys.

## Interpretation

With 10,000 external Taosu decoys, specificity at p_family ≥ 6.0 is **4.7%** (not measurable with benchmarks alone).
EF1% = 9.20 indicates enrichment of actives in the top 1% by p_family.
