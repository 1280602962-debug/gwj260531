#!/usr/bin/env bash
# Run gnina on one shard with per-molecule timeout skip.
# Uses the SAME explicit box as Vina (not autobox) for fair protocol comparison.
# Usage: bash run_gnina_shard.sh shards/shard_0000.txt
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
# shellcheck disable=SC1091
source "${ROOT}/scripts/run_with_timeout.sh"

SHARD_FILE="${1:?shard file}"
read -r START END < "${SHARD_FILE}"
OUT_DIR="${WORK_DIR}/gnina"
LOG_DIR="${WORK_DIR}/logs/gnina"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"
TOUCH="${LOG_DIR}/timeouts.txt"
FAILF="${LOG_DIR}/failed.txt"
SKIPF="${LOG_DIR}/skipped.txt"

EXTRA=()
if [[ "${GNINA_NO_GPU}" == "1" ]]; then
  EXTRA+=(--no_gpu)
  export CUDA_VISIBLE_DEVICES=""
else
  # Torch backend ignores --device; pin via CUDA_VISIBLE_DEVICES (single visible GPU → --device 0).
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-${GNINA_DEVICE:-0}}"
  EXTRA+=(--device 0)
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPU_PER_TASK}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${CPU_PER_TASK}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${CPU_PER_TASK}}"
export TORCH_NUM_THREADS="${TORCH_NUM_THREADS:-${CPU_PER_TASK}}"

echo "[gnina] shard $(basename "${SHARD_FILE}") rows [${START},${END}) timeout=${GNINA_TIMEOUT_SEC}s box=(${CENTER_X},${CENTER_Y},${CENTER_Z}) size=(${SIZE_X},${SIZE_Y},${SIZE_Z}) gpu_no=${GNINA_NO_GPU} cuda_visible=${CUDA_VISIBLE_DEVICES:-}"
for ((i=START; i<END; i++)); do
  mol_id=$(printf "mol_%05d" "$i")
  lig="${WORK_DIR}/ligands_sdf/${mol_id}.sdf"
  out="${OUT_DIR}/${mol_id}_out.sdf"
  log="${LOG_DIR}/${mol_id}.log"
  if [[ ! -f "${lig}" ]]; then
    echo "${mol_id} missing_ligand" >> "${SKIPF}"
    continue
  fi
  # Resume only if a non-empty result exists (crashed runs may leave 0-byte files).
  if [[ -s "${out}" ]]; then
    continue
  fi
  # Already timed out previously → do not retry.
  if [[ -f "${TOUCH}" ]] && grep -q "^${mol_id} timeout_" "${TOUCH}" 2>/dev/null; then
    continue
  fi
  rm -f "${out}"
  set +e
  run_with_timeout "${GNINA_TIMEOUT_SEC}" "${TOUCH}" \
    "${GNINA_BIN}" \
      -r "${RECEPTOR_PDB}" \
      -l "${lig}" \
      --center_x "${CENTER_X}" --center_y "${CENTER_Y}" --center_z "${CENTER_Z}" \
      --size_x "${SIZE_X}" --size_y "${SIZE_Y}" --size_z "${SIZE_Z}" \
      -o "${out}" \
      --exhaustiveness "${EXHAUSTIVENESS}" \
      --num_modes "${NUM_MODES}" \
      --cnn_scoring rescore \
      --seed "${SEED}" \
      --cpu "${CPU_PER_TASK}" \
      "${EXTRA[@]}" > "${log}" 2>&1
  rc=$?
  set -e
  if (( rc == 124 || rc == 137 )); then
    echo "${mol_id} timeout_${GNINA_TIMEOUT_SEC}s" >> "${TOUCH}"
    rm -f "${out}"
    continue
  fi
  if (( rc != 0 )) || [[ ! -f "${out}" ]]; then
    echo "${mol_id} rc=${rc}" >> "${FAILF}"
  fi
done
echo "[gnina] done $(basename "${SHARD_FILE}")"
