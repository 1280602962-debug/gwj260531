#!/usr/bin/env bash
# RTMScore rescoring for vina|gnina poses, in chunks (avoids one giant SDF).
# Usage: bash 04_rtmscore_rescore.sh vina|gnina
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/config.sh"
# shellcheck disable=SC1091
source "${ROOT}/scripts/run_with_timeout.sh"

ENGINE="${1:?vina or gnina}"
CHUNK_SIZE="${RTM_CHUNK_SIZE:-500}"   # molecules per RTMScore call
PY="${RTM_PYTHON_BIN:-${PYTHON_BIN}}"

if [[ -n "${RTMSCORE_ENV_ACTIVATE}" ]]; then
  # shellcheck disable=SC2086
  eval "${RTMSCORE_ENV_ACTIVATE}"
fi
[[ -x "${PY}" || -f "${PY}" ]] || { echo "RTM python missing: ${PY}"; exit 1; }
[[ -f "${RTMSCORE_PY}" ]] || { echo "RTMSCORE_PY missing — set in config.sh"; exit 1; }
[[ -f "${RTMSCORE_MODEL}" ]] || { echo "RTMSCORE_MODEL missing — set in config.sh"; exit 1; }
[[ -f "${RTM_POCKET_PDB}" ]] || { echo "RTM_POCKET_PDB missing"; exit 1; }

POSE_DIR="${WORK_DIR}/${ENGINE}"
OUT_DIR="${WORK_DIR}/rtmscore_${ENGINE}"
mkdir -p "${OUT_DIR}/chunks" "${OUT_DIR}/sdf_chunks"
TOUCH="${OUT_DIR}/timeouts.txt"
FAILF="${OUT_DIR}/failed_chunks.txt"
: > "${OUT_DIR}/scores_all.csv"
HEADER_WRITTEN=0

mapfile -t POSES < <(
  if [[ "${ENGINE}" == "vina" ]]; then
    ls "${POSE_DIR}"/*_out.pdbqt 2>/dev/null | sort
  elif [[ "${ENGINE}" == "gnina" ]]; then
    ls "${POSE_DIR}"/*_out.sdf 2>/dev/null | sort
  else
    echo "engine must be vina or gnina" >&2
    exit 1
  fi
)
echo "[rtmscore] engine=${ENGINE} n_pose_files=${#POSES[@]} chunk_size=${CHUNK_SIZE}"
if (( ${#POSES[@]} == 0 )); then
  echo "No pose files under ${POSE_DIR}"; exit 1
fi

chunk_id=0
for ((offset=0; offset<${#POSES[@]}; offset+=CHUNK_SIZE)); do
  chunk_tag=$(printf "chunk_%04d" "${chunk_id}")
  merged="${OUT_DIR}/sdf_chunks/${chunk_tag}.sdf"
  score_prefix="${OUT_DIR}/chunks/${chunk_tag}"
  : > "${merged}"

  end=$(( offset + CHUNK_SIZE ))
  if (( end > ${#POSES[@]} )); then end=${#POSES[@]}; fi
  for ((j=offset; j<end; j++)); do
    f="${POSES[$j]}"
    if [[ "${ENGINE}" == "vina" ]]; then
      base=$(basename "$f" _out.pdbqt)
      sdf="${OUT_DIR}/sdf_chunks/${base}_poses.sdf"
      if [[ ! -f "${sdf}" ]]; then
        if ! command -v "${OBABEL_BIN}" >/dev/null 2>&1; then
          echo "obabel required to convert vina pdbqt->sdf" >&2
          exit 1
        fi
        "${OBABEL_BIN}" "$f" -O "${sdf}" >/dev/null 2>&1 || {
          echo "${base} obabel_fail" >> "${FAILF}"
          continue
        }
      fi
      [[ -f "${sdf}" ]] && cat "${sdf}" >> "${merged}"
    else
      cat "$f" >> "${merged}"
    fi
  done

  if [[ ! -s "${merged}" ]]; then
    echo "${chunk_tag} empty" >> "${FAILF}"
    chunk_id=$((chunk_id + 1))
    continue
  fi

  set +e
  set -o pipefail
  run_with_timeout "${RTM_TIMEOUT_SEC}" "${TOUCH}" \
    "${PY}" "${RTMSCORE_PY}" \
      -p "${RTM_POCKET_PDB}" \
      -l "${merged}" \
      -m "${RTMSCORE_MODEL}" \
      -o "${score_prefix}" \
      2>&1 | tee "${OUT_DIR}/chunks/${chunk_tag}.log"
  rc=$?
  set +o pipefail
  set -e

  if (( rc == 124 || rc == 137 )); then
    echo "${chunk_tag} timeout_${RTM_TIMEOUT_SEC}s" >> "${TOUCH}"
    chunk_id=$((chunk_id + 1))
    continue
  fi
  if (( rc != 0 )); then
    echo "${chunk_tag} rc=${rc}" >> "${FAILF}"
    chunk_id=$((chunk_id + 1))
    continue
  fi

  # RTMScore writes ${prefix}.csv or ${prefix}_*.csv — pick newest csv in chunks
  csv=""
  for cand in "${score_prefix}.csv" "${score_prefix}_score.csv"; do
    [[ -f "${cand}" ]] && csv="${cand}" && break
  done
  if [[ -z "${csv}" ]]; then
    csv=$(ls -t "${OUT_DIR}/chunks/${chunk_tag}"*.csv 2>/dev/null | head -1 || true)
  fi
  if [[ -z "${csv}" || ! -f "${csv}" ]]; then
    echo "${chunk_tag} missing_csv" >> "${FAILF}"
  else
    if (( HEADER_WRITTEN == 0 )); then
      cat "${csv}" > "${OUT_DIR}/scores_all.csv"
      HEADER_WRITTEN=1
    else
      tail -n +2 "${csv}" >> "${OUT_DIR}/scores_all.csv"
    fi
  fi
  chunk_id=$((chunk_id + 1))
done

# Convenience copy expected by collect_dock_scores.py
cp -f "${OUT_DIR}/scores_all.csv" "${OUT_DIR}/scores.csv"
echo "RTMScore done for ${ENGINE}: ${OUT_DIR}/scores.csv (chunks=${chunk_id})"
