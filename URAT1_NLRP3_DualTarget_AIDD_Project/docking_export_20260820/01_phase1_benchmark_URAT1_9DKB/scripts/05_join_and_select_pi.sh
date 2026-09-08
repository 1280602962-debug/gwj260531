#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
mkdir -p "${RESULT_DIR}" "${WORK_DIR}/scores"
if [[ ! -f "${WORK_DIR}/mol_index_map.csv" ]]; then
  "${PYTHON_BIN}" "${ROOT}/scripts/build_index_map.py" \
    --sdf "${INPUT_DIR}/ligands_all.sdf" \
    --pool "${POOL_CSV}" \
    --out "${WORK_DIR}/mol_index_map.csv"
fi
"${PYTHON_BIN}" "${ROOT}/scripts/collect_dock_scores.py" \
  --work "${WORK_DIR}" --pool "${POOL_CSV}" \
  --index-map "${WORK_DIR}/mol_index_map.csv"
"${PYTHON_BIN}" "${ROOT}/scripts/05_join_and_select_pi.py" \
  --work "${WORK_DIR}" --pool "${POOL_CSV}" \
  --true-bench "${TRUE_BENCH_CSV}" --random-bench "${RANDOM_BENCH_CSV}" \
  --outdir "${RESULT_DIR}"
echo "See ${RESULT_DIR}/protocol_metrics.csv and selected_pi.json"
# also list timeouts if any
echo "Timeouts vina: $(wc -l < "${WORK_DIR}/logs/vina/timeouts.txt" 2>/dev/null || echo 0)"
echo "Timeouts gnina: $(wc -l < "${WORK_DIR}/logs/gnina/timeouts.txt" 2>/dev/null || echo 0)"
