#!/usr/bin/env bash
# Sourced helper: run_with_timeout <seconds> <timeout_log> <cmd...>
# TERM at T, then KILL 15s later. Normalize timeout exit to 124.
run_with_timeout() {
  local SEC="$1"; shift
  local TOLOG="$1"; shift
  local TB="${TIMEOUT_BIN:-timeout}"
  if ! command -v "${TB}" >/dev/null 2>&1 && [[ ! -x "${TB}" ]]; then
    "$@"
    return $?
  fi
  "${TB}" -k 15 "${SEC}" "$@"
  local rc=$?
  # 124 = GNU timeout; 137 = 128+SIGKILL (seen with --signal=KILL)
  if (( rc == 124 || rc == 137 )); then
    echo "TIMEOUT after ${SEC}s (rc=${rc}): $*" >> "${TOLOG}"
    return 124
  fi
  return "${rc}"
}
