# Compound libraries for virtual screening (Step 06)

Supported formats: `.smi`, `.csv`, `.smi.gz`, `.csv.gz`

## Taosu 1M library (local Windows path)

If your file is on your PC, e.g.:

```
D:\CADD paper exercise\Document_PLK1 and NLRP3\ML\taosu_20210823_100w_asteroid_murcko_protonized.csv
```

**Do not commit this file to GitHub** (already in `.gitignore`). Copy it into this folder:

```
<project>/data/libraries/taosu_100w.csv
```

### Option A — Cursor / VS Code

1. Open the project folder locally.
2. Drag `taosu_20210823_100w_asteroid_murcko_protonized.csv` into `data/libraries/`.
3. Rename to `taosu_100w.csv` (optional, shorter path).

### Option B — Command line (Windows PowerShell)

```powershell
Copy-Item "D:\CADD paper exercise\Document_PLK1 and NLRP3\ML\taosu_20210823_100w_asteroid_murcko_protonized.csv" `
  ".\data\libraries\taosu_100w.csv"
```

### Run screening

Ensure models exist (`python3 scripts/07_compare_models.py` once).

```bash
# Smoke test (first 10k rows, ~1 min)
python3 scripts/06_virtual_screening.py \
  --library data/libraries/taosu_100w.csv \
  --output results/screening_taosu_test \
  --max-rows 10000 \
  --batch-size 50000

# Full ~1M run (expect ~1–3 hours depending on CPU)
python3 scripts/06_virtual_screening.py \
  --library data/libraries/taosu_100w.csv \
  --output results/screening_taosu_1M \
  --batch-size 50000 \
  --top-n 5000 \
  --diverse-n 500
```

If SMILES column is not auto-detected:

```bash
python3 scripts/06_virtual_screening.py \
  --library data/libraries/taosu_100w.csv \
  --smiles-column YOUR_COLUMN_NAME \
  --output results/screening_taosu_1M
```

### Outputs

- `results/screening_taosu_1M/all_hits.csv` — all F1+SA/QED passes
- `top5000.csv` / `top500_diverse.csv` — ranked subsets for docking (F3)

Funnel: preprocess → Lipinski → **p_family ≥ 6.0** → SA/QED → rank.
