#!/usr/bin/env bash
# Sequential docking + overlapping RTMScore:
#   1) Vina (CPU) — RTM scores vina poses as they finish
#   2) After Vina: give CPUs to gnina (GPU CNN) — RTM scores gnina poses too
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"

export VINA_MAX_JOBS="${VINA_MAX_JOBS:-8}"
export CPU_PER_TASK="${CPU_PER_TASK:-4}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPU_PER_TASK}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${CPU_PER_TASK}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${CPU_PER_TASK}}"
export TORCH_NUM_THREADS="${TORCH_NUM_THREADS:-${CPU_PER_TASK}}"
# Stage-1: do not run gnina yet (CPU reserved for vina + light RTM).
# Stage-2 overrides below after vina finishes.
export GNINA_BIN="${GNINA_BIN:-/home/hww/gwj/NLRP3_URAT1/software/gnina_5090.sh}"
export GNINA_NO_GPU="${GNINA_NO_GPU:-0}"
export GNINA_DEVICES="${GNINA_DEVICES:-0,1}"

mkdir -p "${LOG_DIR}" "${WORK_DIR}/logs/vina" "${WORK_DIR}/logs/gnina"
find "${WORK_DIR}/vina" -name '*_out.pdbqt' -size 0 -delete 2>/dev/null || true
find "${WORK_DIR}/gnina" -name '*_out.sdf' -size 0 -delete 2>/dev/null || true

ts() { date '+%F %T'; }
count_done() {
  local kind="$1"
  if [[ "${kind}" == vina ]]; then
    find "${WORK_DIR}/vina" -name '*_out.pdbqt' -size +0 2>/dev/null | wc -l
  else
    find "${WORK_DIR}/gnina" -name '*_out.sdf' -size +0 2>/dev/null | wc -l
  fi
}
count_rtm() {
  local kind="$1"
  local f="${WORK_DIR}/rtmscore_${kind}/scored_mols.txt"
  if [[ -f "${f}" ]]; then sort -u "${f}" | wc -l; else echo 0; fi
}

echo "[$(ts)] launch SEQUENTIAL: vina(timeout=${VINA_TIMEOUT_SEC}s,jobs=${VINA_MAX_JOBS}) -> gnina(timeout=${GNINA_TIMEOUT_SEC}s) + incremental RTM"
echo "[$(ts)] vina done=$(count_done vina)/9838  gnina done=$(count_done gnina)/9839"

(
  while true; do
    vd=$(count_done vina)
    gd=$(count_done gnina)
    rv=$(count_rtm vina)
    rg=$(count_rtm gnina)
    load=$(cut -d' ' -f1-3 /proc/loadavg)
    echo "[$(ts)] progress vina=${vd}/9838 gnina=${gd}/9839 rtm_v=${rv} rtm_g=${rg} load=${load}"
    sleep 120
  done
) >> "${LOG_DIR}/progress.log" 2>&1 &
MON_PID=$!
echo "${MON_PID}" > "${LOG_DIR}/progress_monitor.pid"

# RTM on vina poses while vina runs (light CPU).
rm -f "${WORK_DIR}/rtmscore_vina/STOP"
nice -n 10 bash "${ROOT}/scripts/06_rtmscore_incremental.sh" vina \
  > "${LOG_DIR}/rtm_vina_incremental.log" 2>&1 &
RTM_V_PID=$!
echo "${RTM_V_PID}" > "${LOG_DIR}/rtm_vina.pid"
echo "[$(ts)] started incremental RTM(vina) pid=${RTM_V_PID}"

echo "[$(ts)] === STAGE 1/2: Vina (CPU) ==="
bash "${ROOT}/scripts/02_run_vina_local.sh" > "${LOG_DIR}/vina_all.log" 2>&1
echo "[$(ts)] vina finished: $(count_done vina)/9838"

# Flush remaining vina RTM, then stop the long-running watcher.
touch "${WORK_DIR}/rtmscore_vina/STOP"
wait "${RTM_V_PID}" 2>/dev/null || true
nice -n 5 bash "${ROOT}/scripts/06_rtmscore_incremental.sh" vina --until-idle \
  > "${LOG_DIR}/rtm_vina_flush.log" 2>&1 || true
echo "[$(ts)] rtm vina scored=$(count_rtm vina)"

echo "[$(ts)] === STAGE 2/2: gnina (CPUs freed from vina; GPU CNN) ==="
export GNINA_MAX_JOBS="${GNINA_MAX_JOBS:-6}"
export GNINA_NO_GPU=0
export GNINA_DEVICES="${GNINA_DEVICES:-0,1}"
export CPU_PER_TASK="${CPU_PER_TASK:-4}"
export OMP_NUM_THREADS="${CPU_PER_TASK}"
export MKL_NUM_THREADS="${CPU_PER_TASK}"
export OPENBLAS_NUM_THREADS="${CPU_PER_TASK}"
export TORCH_NUM_THREADS="${CPU_PER_TASK}"
echo "[$(ts)] gnina jobs=${GNINA_MAX_JOBS} devices=${GNINA_DEVICES} bin=${GNINA_BIN}"

rm -f "${WORK_DIR}/rtmscore_gnina/STOP"
nice -n 10 bash "${ROOT}/scripts/06_rtmscore_incremental.sh" gnina \
  > "${LOG_DIR}/rtm_gnina_incremental.log" 2>&1 &
RTM_G_PID=$!
echo "${RTM_G_PID}" > "${LOG_DIR}/rtm_gnina.pid"
echo "[$(ts)] started incremental RTM(gnina) pid=${RTM_G_PID}"

bash "${ROOT}/scripts/03_run_gnina_local.sh" > "${LOG_DIR}/gnina_all.log" 2>&1
echo "[$(ts)] gnina finished: $(count_done gnina)/9839"

touch "${WORK_DIR}/rtmscore_gnina/STOP"
wait "${RTM_G_PID}" 2>/dev/null || true
nice -n 5 bash "${ROOT}/scripts/06_rtmscore_incremental.sh" gnina --until-idle \
  > "${LOG_DIR}/rtm_gnina_flush.log" 2>&1 || true

kill "${MON_PID}" 2>/dev/null || true
echo "[$(ts)] all docking + RTM finished"
echo "[$(ts)] final vina=$(count_done vina)/9838 gnina=$(count_done gnina)/9839 rtm_v=$(count_rtm vina) rtm_g=$(count_rtm gnina)"
wc -l "${WORK_DIR}/logs/vina/skipped.txt" "${WORK_DIR}/logs/vina/timeouts.txt" "${WORK_DIR}/logs/vina/failed.txt" 2>/dev/null || true
wc -l "${WORK_DIR}/logs/gnina/skipped.txt" "${WORK_DIR}/logs/gnina/timeouts.txt" "${WORK_DIR}/logs/gnina/failed.txt" 2>/dev/null || true
