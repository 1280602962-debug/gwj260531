#!/usr/bin/env bash
# Run Vina on one shard with per-molecule timeout skip.
# Usage: bash run_vina_shard.sh shards/shard_0000.txt
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
# shellcheck disable=SC1091
source "${ROOT}/scripts/run_with_timeout.sh"

SHARD_FILE="${1:?shard file}"
read -r START END < "${SHARD_FILE}"
OUT_DIR="${WORK_DIR}/vina"
LOG_DIR="${WORK_DIR}/logs/vina"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"
TOUCH="${LOG_DIR}/timeouts.txt"
FAILF="${LOG_DIR}/failed.txt"
SKIPF="${LOG_DIR}/skipped.txt"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPU_PER_TASK}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${CPU_PER_TASK}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${CPU_PER_TASK}}"

echo "[vina] shard $(basename "${SHARD_FILE}") rows [${START},${END}) timeout=${VINA_TIMEOUT_SEC}s"
for ((i=START; i<END; i++)); do
  mol_id=$(printf "mol_%05d" "$i")
  lig="${WORK_DIR}/ligands_pdbqt/${mol_id}.pdbqt"
  out="${OUT_DIR}/${mol_id}_out.pdbqt"
  log="${LOG_DIR}/${mol_id}.log"
  if [[ ! -f "${lig}" ]]; then
    echo "${mol_id} missing_ligand" >> "${SKIPF}"
    continue
  fi
  # Resume only if a non-empty result exists (crashed runs may leave 0-byte files).
  if [[ -s "${out}" ]]; then
    continue
  fi
  # Already timed out previously → do not retry (keeps workers on remaining ligands).
  if [[ -f "${TOUCH}" ]] && grep -q "^${mol_id} timeout_" "${TOUCH}" 2>/dev/null; then
    continue
  fi
  rm -f "${out}"
  set +e
  run_with_timeout "${VINA_TIMEOUT_SEC}" "${TOUCH}" \
    "${VINA_BIN}" \
      --receptor "${RECEPTOR_PDBQT}" \
      --ligand "${lig}" \
      --center_x "${CENTER_X}" --center_y "${CENTER_Y}" --center_z "${CENTER_Z}" \
      --size_x "${SIZE_X}" --size_y "${SIZE_Y}" --size_z "${SIZE_Z}" \
      --exhaustiveness "${EXHAUSTIVENESS}" \
      --num_modes "${NUM_MODES}" \
      --cpu "${CPU_PER_TASK}" \
      --seed "${SEED}" \
      --out "${out}" \
      --verbosity 1 > "${log}" 2>&1
  rc=$?
  set -e
  if (( rc == 124 || rc == 137 )); then
    echo "${mol_id} timeout_${VINA_TIMEOUT_SEC}s" >> "${TOUCH}"
    rm -f "${out}"
    continue
  fi
  if (( rc != 0 )) || [[ ! -f "${out}" ]]; then
    echo "${mol_id} rc=${rc}" >> "${FAILF}"
  fi
done
echo "[vina] done $(basename "${SHARD_FILE}")"
