#!/usr/bin/env bash
# Start gnina alongside an already-running vina job (does not touch vina).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
# Force Blackwell build (ignore stale GNINA_BIN from parent env).
export GNINA_BIN=/home/hww/gwj/NLRP3_URAT1/software/gnina_5090.sh
export GNINA_NO_GPU=0

export GNINA_NO_GPU="${GNINA_NO_GPU:-0}"
export GNINA_MAX_JOBS="${GNINA_MAX_JOBS:-2}"
export GNINA_DEVICES="${GNINA_DEVICES:-0}"
export CPU_PER_TASK="${CPU_PER_TASK:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPU_PER_TASK}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${CPU_PER_TASK}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${CPU_PER_TASK}}"
export TORCH_NUM_THREADS="${TORCH_NUM_THREADS:-${CPU_PER_TASK}}"

mkdir -p "${LOG_DIR}" "${WORK_DIR}/logs/gnina"
find "${WORK_DIR}/gnina" -name '*_out.sdf' -size 0 -delete 2>/dev/null || true

ts() { date '+%F %T'; }
gd=$(find "${WORK_DIR}/gnina" -name '*_out.sdf' -size +0 2>/dev/null | wc -l)
vd=$(find "${WORK_DIR}/vina" -name '*_out.pdbqt' -size +0 2>/dev/null | wc -l)
echo "[$(ts)] start gnina PARALLEL with vina: GNINA_MAX_JOBS=${GNINA_MAX_JOBS} devices=${GNINA_DEVICES} no_gpu=${GNINA_NO_GPU}"
echo "[$(ts)] current vina=${vd}/9838 gnina=${gd}/9839 bin=${GNINA_BIN}"

bash "${ROOT}/scripts/03_run_gnina_local.sh" > "${LOG_DIR}/gnina_all.log" 2>&1
echo "[$(ts)] gnina finished: $(find "${WORK_DIR}/gnina" -name '*_out.sdf' -size +0 2>/dev/null | wc -l)/9839"
wc -l "${WORK_DIR}/logs/gnina/skipped.txt" "${WORK_DIR}/logs/gnina/timeouts.txt" "${WORK_DIR}/logs/gnina/failed.txt" 2>/dev/null || true
