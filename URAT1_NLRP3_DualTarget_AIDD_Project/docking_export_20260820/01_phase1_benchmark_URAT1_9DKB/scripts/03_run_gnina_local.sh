#!/usr/bin/env bash
# Local parallel gnina (cap concurrency; per-mol timeout inside shard script).
# Uses a PID pool so concurrency works in non-interactive shells (no job control).
# With GPU: round-robin GNINA_DEVICES (e.g. "0" or "0,1") via CUDA_VISIBLE_DEVICES.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
NPROC="${NPROC:-$(nproc)}"
MAX_JOBS=$(( NPROC / CPU_PER_TASK ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi
if [[ -n "${GNINA_MAX_JOBS:-}" ]]; then
  MAX_JOBS="${GNINA_MAX_JOBS}"
elif (( MAX_JOBS > 2 )); then
  MAX_JOBS=2
fi
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

IFS=',' read -r -a GPU_DEVS <<< "${GNINA_DEVICES:-${GNINA_DEVICE:-0}}"
if ((${#GPU_DEVS[@]} == 0)); then GPU_DEVS=(0); fi

mapfile -t SHARDS < <(ls "${SHARD_DIR}"/shard_*.txt | sort)
echo "Running ${#SHARDS[@]} gnina shards, max concurrent=${MAX_JOBS}, timeout/mol=${GNINA_TIMEOUT_SEC}s no_gpu=${GNINA_NO_GPU} devices=${GNINA_DEVICES:-${GNINA_DEVICE:-0}} bin=${GNINA_BIN}"

running=0
idx=0
declare -a PIDS=()
for shard in "${SHARDS[@]}"; do
  while (( running >= MAX_JOBS )); do
    new_pids=()
    running=0
    for pid in "${PIDS[@]}"; do
      if kill -0 "${pid}" 2>/dev/null; then
        new_pids+=("${pid}")
        running=$((running + 1))
      fi
    done
    PIDS=("${new_pids[@]}")
    if (( running >= MAX_JOBS )); then
      sleep 2
    fi
  done
  dev="${GPU_DEVS[$((idx % ${#GPU_DEVS[@]}))]}"
  idx=$((idx + 1))
  if [[ "${GNINA_NO_GPU}" == "1" ]]; then
    bash "${ROOT}/scripts/run_gnina_shard.sh" "${shard}" &
  else
    GNINA_DEVICE="${dev}" CUDA_VISIBLE_DEVICES="${dev}" \
      bash "${ROOT}/scripts/run_gnina_shard.sh" "${shard}" &
  fi
  PIDS+=("$!")
  running=$((running + 1))
done

for pid in "${PIDS[@]}"; do
  wait "${pid}" || true
done
echo "All gnina shards finished"
