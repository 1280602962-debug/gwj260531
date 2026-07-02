#!/usr/bin/env bash
# GNINA CPU batch: P05 pool @ 9DKB + 7ALV → Pareto
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

POOL="${POOL:-data/repurposing/screening/docking_pool_p05.csv}"
LIG_DIR="${LIG_DIR:-results/repurposing/ligands_p05}"
JOBS="${JOBS:-2}"  # WSL CPU: 2 parallel jobs recommended (no GPU)

export PATH="$ROOT/tools:$PATH"

python3 scripts/prepare_ligands_vina.py --input "$POOL" --output-dir "$LIG_DIR"

python3 scripts/run_gnina_batch.py --target urat1_9dkb \
  --manifest "$LIG_DIR/ligand_manifest.csv" \
  --output-dir results/repurposing/docking_gnina/9dkb --jobs "$JOBS"

python3 scripts/run_gnina_batch.py --target nlrp3_7alv \
  --manifest "$LIG_DIR/ligand_manifest.csv" \
  --output-dir results/repurposing/docking_gnina/7alv --jobs "$JOBS"

python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_gnina/9dkb/docking_9dkb_gnina.csv \
  --pdb 9DKB --engine gnina \
  --output results/repurposing/docking_raw/urat1_9dkb_p05.csv

python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_gnina/7alv/docking_7alv_gnina.csv \
  --pdb 7ALV --engine gnina \
  --output results/repurposing/docking_raw/nlrp3_7alv_p05.csv

python3 scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_7alv_p05.csv \
  --pool "$POOL" --sn-mode both

python3 scripts/analyze_pareto_benchmarks.py
python3 scripts/plot_available_figures.py
echo "GNINA pipeline done."
