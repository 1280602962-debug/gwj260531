#!/usr/bin/env bash
# Open-source dual-target docking pipeline (Vina 1.2.5)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

POOL="${POOL:-data/repurposing/screening/docking_pool_p05.csv}"
LIG_DIR="${LIG_DIR:-results/repurposing/ligands_p05}"
JOBS="${JOBS:-8}"

echo "== Receptor prep =="
python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
python3 scripts/prepare_receptor_vina.py --target nlrp3_7alv

echo "== Ligand prep =="
python3 scripts/prepare_ligands_vina.py --input "$POOL" --output-dir "$LIG_DIR"

echo "== URAT1 @ 9DKB =="
python3 scripts/run_vina_batch.py --target urat1_9dkb \
  --manifest "$LIG_DIR/ligand_manifest.csv" \
  --output-dir results/repurposing/docking_vina/9dkb --jobs "$JOBS"

echo "== NLRP3 @ 7ALV =="
python3 scripts/run_vina_batch.py --target nlrp3_7alv \
  --manifest "$LIG_DIR/ligand_manifest.csv" \
  --output-dir results/repurposing/docking_vina/7alv --jobs "$JOBS"

echo "== Normalize =="
python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_vina/9dkb/docking_9dkb_vina.csv \
  --pdb 9DKB --engine vina \
  --output results/repurposing/docking_raw/urat1_9dkb_p05.csv
python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_vina/7alv/docking_7alv_vina.csv \
  --pdb 7ALV --engine vina \
  --output results/repurposing/docking_raw/nlrp3_7alv_p05.csv

echo "== Pareto =="
python3 scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_7alv_p05.csv \
  --pool "$POOL" --sn-mode both

python3 scripts/analyze_pareto_benchmarks.py
python3 scripts/plot_available_figures.py
echo "Done."
