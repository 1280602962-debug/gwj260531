#!/usr/bin/env bash
# Local parallel Vina over all shards (per-mol timeout inside shard script).
# Uses a PID-file pool so concurrency works in non-interactive shells (no job control).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
NPROC="${NPROC:-$(nproc)}"
MAX_JOBS=$(( NPROC / CPU_PER_TASK ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi
if [[ -n "${VINA_MAX_JOBS:-}" ]]; then
  MAX_JOBS="${VINA_MAX_JOBS}"
fi
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

mapfile -t SHARDS < <(ls "${SHARD_DIR}"/shard_*.txt | sort)
echo "Running ${#SHARDS[@]} vina shards, max concurrent=${MAX_JOBS}, timeout/mol=${VINA_TIMEOUT_SEC}s"

running=0
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
  bash "${ROOT}/scripts/run_vina_shard.sh" "${shard}" &
  PIDS+=("$!")
  running=$((running + 1))
done

# Wait for all remaining
for pid in "${PIDS[@]}"; do
  wait "${pid}" || true
done
echo "All vina shards finished"
