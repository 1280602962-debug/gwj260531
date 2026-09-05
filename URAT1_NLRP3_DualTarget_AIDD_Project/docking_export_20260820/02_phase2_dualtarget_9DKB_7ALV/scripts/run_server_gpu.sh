#!/usr/bin/env bash
# P2 dual-target docking on GPU server (9DKB + 7ALV).
# Resume-safe: skips ligands that already have pose SDF+log under results/.../poses/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

# ---- tune these ----
JOBS="${JOBS:-1}"          # parallel ligands; on 1 GPU usually 1–2 is best
CFG="${CFG:-config/docking_production_p2_gpu.yaml}"
MANIFEST="${MANIFEST:-results/repurposing/ligands_p05/ligand_manifest.csv}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
export PATH="/opt/anaconda3/bin:${PATH}"
PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
GNINA_BIN="${GNINA_BIN:-/home/hww/gwj/NLRP3_URAT1/software/gnina_5090.sh}"
# --------------------

if [[ ! -x "$GNINA_BIN" ]] && ! command -v gnina >/dev/null 2>&1; then
  echo "ERROR: gnina not found ($GNINA_BIN)." >&2
  exit 1
fi
# Ensure `gnina` name resolves for version print / sanity
mkdir -p "$ROOT/tools"
ln -sfn "$GNINA_BIN" "$ROOT/tools/gnina"
export PATH="$ROOT/tools:$PATH"

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: python not found: $PYTHON_BIN" >&2
  exit 1
fi

echo "=== Server P2 GPU funnel ==="
echo "root=$ROOT"
echo "gnina=$("$GNINA_BIN" --version 2>&1 | head -1)"
echo "python=$PYTHON_BIN"
echo "jobs=$JOBS  config=$CFG  CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L 2>/dev/null || echo "(nvidia-smi not found — still trying GPU mode from config)"

mkdir -p results/repurposing/logs results/repurposing/docking_p2/9dkb/poses results/repurposing/docking_p2/7alv/poses

echo ""
echo ">>> [1/3] Dock URAT1 9DKB (resume skips finished poses)"
"$PYTHON_BIN" scripts/run_gnina_batch.py \
  --config "$CFG" \
  --target urat1_9dkb \
  --manifest "$MANIFEST" \
  --output-dir results/repurposing/docking_p2/9dkb \
  --jobs "$JOBS"

echo ""
echo ">>> [2/3] Dock NLRP3 7ALV"
"$PYTHON_BIN" scripts/run_gnina_batch.py \
  --config "$CFG" \
  --target nlrp3_7alv \
  --manifest "$MANIFEST" \
  --output-dir results/repurposing/docking_p2/7alv \
  --jobs "$JOBS"

echo ""
echo ">>> [3/3] Merge + Pareto shortlist"
PARETO_OUT="${PARETO_OUT:-data/repurposing/pareto}"
"$PYTHON_BIN" scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_p2/9dkb/docking_9dkb_gnina.csv \
  --nlrp3-dock results/repurposing/docking_p2/7alv/docking_7alv_gnina.csv \
  --nlrp3-pdb 7ALV \
  --pool data/repurposing/screening/docking_pool_p05.csv \
  --sn-mode both \
  --output-dir "$PARETO_OUT"

echo ""
echo "Done. Key outputs:"
echo "  results/repurposing/docking_p2/9dkb/docking_9dkb_gnina.csv"
echo "  results/repurposing/docking_p2/7alv/docking_7alv_gnina.csv"
echo "  $PARETO_OUT/pareto_merged_scores.csv"
echo "  $PARETO_OUT/pareto_shortlist.csv"
