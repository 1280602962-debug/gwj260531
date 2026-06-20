# JNK1 Selectivity CADD Project (Full Bundle)

This folder is a **standalone snapshot** of the JNK1/2/3 selectivity virtual screening project.

**Run all commands from this directory** (`JNK1_Selectivity_Project/`), not the parent repo root.

## Contents

- `scripts/06_virtual_screening.py` — F1 `p_family ≥ 6.0` + drug-like + SA/QED (v2, supports million-scale CSV)
- `scripts/07_compare_models.py` — Train per-isoform XGBoost models
- `scripts/calibrate_threshold.py` — Benchmark F1 threshold calibration
- `data/benchmarks/literature_benchmarks.csv` — 9 reference inhibitors
- `config/targets.yaml` — Includes `screening.p_family_threshold: 6.0`

## Quick start

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt

# Train models (first time)
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop

# Calibrate F1 on benchmarks (optional)
python3 scripts/calibrate_threshold.py

# Demo screening
python3 scripts/build_demo_library.py
python3 scripts/06_virtual_screening.py \
  --library data/libraries/screening_demo.smi \
  --output results/screening_v2

# Million-compound CSV (e.g. Taosu library)
python3 scripts/06_virtual_screening.py \
  --library data/libraries/taosu_100w.csv \
  --output results/screening_taosu_1M \
  --batch-size 50000 \
  --top-n 5000
```

See `data/libraries/README.md` for copying large libraries from Windows/WSL paths.

## Note

The git repository root is the parent directory (`gwj260531`). This subfolder is kept in sync with root `scripts/` and `config/` for offline/local use.
