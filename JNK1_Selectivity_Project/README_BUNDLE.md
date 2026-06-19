# JNK1 Selectivity CADD Project (Full Bundle)

This folder is a **complete snapshot** of the JNK1/2/3 selectivity virtual screening project, including:

- `data/benchmarks/literature_benchmarks.csv` — 8 reference inhibitors for threshold calibration
- `scripts/calibrate_threshold.py` — Benchmark-calibrated F1 threshold scanner
- `scripts/` — Full pipeline (00–07, plot_style, run_selectivity_pipeline)
- `config/`, `docs/`, `data/processed/`, `results/`

## Quick start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Train per-isoform XGBoost models (if models/ missing)
python3 scripts/07_compare_models.py

# 3. Calibrate F1 threshold on benchmark panel
python3 scripts/calibrate_threshold.py

# 4. Run selectivity pipeline (when docking module is ready)
python3 scripts/run_selectivity_pipeline.py
```

## Benchmark panel

See `data/benchmarks/README.md` for the 8-compound isoform reference set (SP600125, CC-90001, CC-930, JNK-IN-8, TCS JNK 6O, Q63, etc.).

## Note

The git repository root is the parent directory (`gwj260531`). This subfolder mirrors the full project for local organization and sharing.
