#!/usr/bin/env bash
# Incremental RTMScore: score newly finished vina|gnina poses while docking runs.
# Usage:
#   bash 06_rtmscore_incremental.sh vina|gnina [--until-idle]
#   --until-idle: exit when no pending poses remain (after docking mostly done)
# Stop early by: touch "${WORK_DIR}/rtmscore_${ENGINE}/STOP"
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
# shellcheck disable=SC1091
source "${ROOT}/scripts/run_with_timeout.sh"

ENGINE="${1:?vina or gnina}"
UNTIL_IDLE=0
if [[ "${2:-}" == "--until-idle" ]]; then UNTIL_IDLE=1; fi

CHUNK_SIZE="${RTM_CHUNK_SIZE:-50}"
POLL_SEC="${RTM_POLL_SEC:-60}"
PY="${RTM_PYTHON_BIN:-${PYTHON_BIN}}"

if [[ -n "${RTMSCORE_ENV_ACTIVATE}" ]]; then
  # shellcheck disable=SC2086
  eval "${RTMSCORE_ENV_ACTIVATE}"
fi
[[ -x "${PY}" || -f "${PY}" ]] || { echo "RTM python missing: ${PY}"; exit 1; }
[[ -f "${RTMSCORE_PY}" ]] || { echo "RTMSCORE_PY missing"; exit 1; }
[[ -f "${RTMSCORE_MODEL}" ]] || { echo "RTMSCORE_MODEL missing"; exit 1; }
[[ -f "${RTM_POCKET_PDB}" ]] || { echo "RTM_POCKET_PDB missing"; exit 1; }

POSE_DIR="${WORK_DIR}/${ENGINE}"
OUT_DIR="${WORK_DIR}/rtmscore_${ENGINE}"
mkdir -p "${OUT_DIR}/chunks" "${OUT_DIR}/sdf_chunks" "${OUT_DIR}/per_mol_sdf"
DONE_FILE="${OUT_DIR}/scored_mols.txt"
FAILF="${OUT_DIR}/failed_chunks.txt"
TOUCH="${OUT_DIR}/timeouts.txt"
STOPF="${OUT_DIR}/STOP"
LOCKF="${OUT_DIR}/.rtm.lock"
: >> "${DONE_FILE}"
: >> "${FAILF}"
: >> "${TOUCH}"

export OMP_NUM_THREADS="${RTM_OMP_THREADS:-2}"
export MKL_NUM_THREADS="${RTM_OMP_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${RTM_OMP_THREADS:-2}"
export TORCH_NUM_THREADS="${RTM_OMP_THREADS:-2}"

ts() { date '+%F %T'; }

list_pose_mols() {
  if [[ "${ENGINE}" == "vina" ]]; then
    find "${POSE_DIR}" -maxdepth 1 -name '*_out.pdbqt' -size +0 -printf '%f\n' 2>/dev/null \
      | sed 's/_out\.pdbqt$//' | sort
  else
    find "${POSE_DIR}" -maxdepth 1 -name '*_out.sdf' -size +0 -printf '%f\n' 2>/dev/null \
      | sed 's/_out\.sdf$//' | sort
  fi
}

rebuild_scores_csv() {
  local first=1
  : > "${OUT_DIR}/scores_all.csv"
  # shellcheck disable=SC2012
  mapfile -t csvs < <(ls -1 "${OUT_DIR}/chunks"/chunk_*.csv 2>/dev/null | sort || true)
  for csv in "${csvs[@]:-}"; do
    [[ -s "${csv}" ]] || continue
    if (( first == 1 )); then
      cat "${csv}" > "${OUT_DIR}/scores_all.csv"
      first=0
    else
      tail -n +2 "${csv}" >> "${OUT_DIR}/scores_all.csv"
    fi
  done
  if [[ -s "${OUT_DIR}/scores_all.csv" ]]; then
    cp -f "${OUT_DIR}/scores_all.csv" "${OUT_DIR}/scores.csv"
  fi
}

score_chunk() {
  local -a mols=("$@")
  (( ${#mols[@]} > 0 )) || return 0
  local chunk_id chunk_tag merged score_prefix mol f sdf
  chunk_id=$(printf "%04d" "$(find "${OUT_DIR}/chunks" -name 'chunk_*.log' 2>/dev/null | wc -l)")
  chunk_tag="chunk_${chunk_id}"
  merged="${OUT_DIR}/sdf_chunks/${chunk_tag}.sdf"
  score_prefix="${OUT_DIR}/chunks/${chunk_tag}"
  : > "${merged}"

  for mol in "${mols[@]}"; do
    if [[ "${ENGINE}" == "vina" ]]; then
      f="${POSE_DIR}/${mol}_out.pdbqt"
      sdf="${OUT_DIR}/per_mol_sdf/${mol}_poses.sdf"
      if [[ ! -s "${sdf}" ]]; then
        if ! "${OBABEL_BIN}" "${f}" -O "${sdf}" >/dev/null 2>&1; then
          echo "${mol} obabel_fail" >> "${FAILF}"
          continue
        fi
      fi
      [[ -s "${sdf}" ]] && cat "${sdf}" >> "${merged}"
    else
      f="${POSE_DIR}/${mol}_out.sdf"
      [[ -s "${f}" ]] && cat "${f}" >> "${merged}"
    fi
  done

  if [[ ! -s "${merged}" ]]; then
    echo "${chunk_tag} empty" >> "${FAILF}"
    return 0
  fi

  echo "[$(ts)] [rtm/${ENGINE}] scoring ${#mols[@]} mols -> ${chunk_tag}"
  set +e
  set -o pipefail
  run_with_timeout "${RTM_TIMEOUT_SEC}" "${TOUCH}" \
    "${PY}" "${RTMSCORE_PY}" \
      -p "${RTM_POCKET_PDB}" \
      -l "${merged}" \
      -m "${RTMSCORE_MODEL}" \
      -o "${score_prefix}" \
      2>&1 | tee "${OUT_DIR}/chunks/${chunk_tag}.log"
  local rc=$?
  set +o pipefail
  set -e

  if (( rc == 124 || rc == 137 )); then
    echo "${chunk_tag} timeout_${RTM_TIMEOUT_SEC}s" >> "${TOUCH}"
    return 0
  fi
  if (( rc != 0 )); then
    echo "${chunk_tag} rc=${rc}" >> "${FAILF}"
    return 0
  fi

  local csv=""
  for cand in "${score_prefix}.csv" "${score_prefix}_score.csv"; do
    [[ -f "${cand}" ]] && csv="${cand}" && break
  done
  if [[ -z "${csv}" ]]; then
    csv=$(ls -t "${OUT_DIR}/chunks/${chunk_tag}"*.csv 2>/dev/null | head -1 || true)
  fi
  if [[ -z "${csv}" || ! -f "${csv}" ]]; then
    echo "${chunk_tag} missing_csv" >> "${FAILF}"
    return 0
  fi

  printf '%s\n' "${mols[@]}" >> "${DONE_FILE}"
  rebuild_scores_csv
  local n_done n_pose
  n_done=$(sort -u "${DONE_FILE}" | wc -l)
  n_pose=$(list_pose_mols | wc -l)
  echo "[$(ts)] [rtm/${ENGINE}] ok ${chunk_tag}; scored_unique=${n_done} poses_ready=${n_pose}"
}

# Single-instance guard
exec 9>"${LOCKF}"
if ! flock -n 9; then
  echo "[$(ts)] [rtm/${ENGINE}] another incremental scorer holds ${LOCKF}; exit"
  exit 0
fi

echo "[$(ts)] [rtm/${ENGINE}] start incremental chunk=${CHUNK_SIZE} poll=${POLL_SEC}s py=${PY}"
idle_rounds=0
while true; do
  if [[ -f "${STOPF}" ]]; then
    echo "[$(ts)] [rtm/${ENGINE}] STOP file present; exiting"
    break
  fi

  mapfile -t pending < <(
    comm -23 <(list_pose_mols) <(sort -u "${DONE_FILE}")
  )
  if (( ${#pending[@]} >= CHUNK_SIZE )); then
    score_chunk "${pending[@]:0:${CHUNK_SIZE}}"
    idle_rounds=0
    continue
  fi

  if (( ${#pending[@]} > 0 )); then
    # Flush smaller leftovers if docking is no longer producing quickly.
    if (( UNTIL_IDLE == 1 || idle_rounds >= 2 )); then
      score_chunk "${pending[@]}"
      idle_rounds=0
      continue
    fi
  fi

  if (( UNTIL_IDLE == 1 && ${#pending[@]} == 0 )); then
    echo "[$(ts)] [rtm/${ENGINE}] caught up; exiting"
    break
  fi

  idle_rounds=$((idle_rounds + 1))
  sleep "${POLL_SEC}"
done

rebuild_scores_csv
echo "[$(ts)] [rtm/${ENGINE}] finished scored_unique=$(sort -u "${DONE_FILE}" | wc -l)"
