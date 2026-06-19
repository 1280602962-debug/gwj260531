# JNK Isoform Benchmark Panel

Curated reference inhibitors for **F1 threshold calibration** and **F3 docking direction validation**.

## Files

| File | Purpose |
|------|---------|
| `literature_benchmarks.csv` | Primary benchmark panel (9 compounds) |

## Usage

```bash
# After training/loading xgboost_jnk{1,2,3}.joblib
python3 scripts/calibrate_threshold.py \
  --benchmarks data/benchmarks/literature_benchmarks.csv \
  --output results/calibration
```

## Field guide

- **expected_profile**: Intended isoform behavior class (not a training label from ChEMBL).
- **must_pass_F1**: If `yes`, compound must survive ML pre-filter at the calibrated threshold.
- **use_f3_validation**: If `yes`, use for tri-structure docking direction checks.

## Data sources

IC50 values follow the user reference figure / primary literature; SMILES from ChEMBL export in `docs/JNK*.csv`.
