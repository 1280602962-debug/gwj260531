#!/usr/bin/env bash
# Wait for JCIM four-pair multi-seed Vina to finish, then resume paused C1 GNINA docks.
set -u
MS_ROOT="/tmp/gwj_check/Dual_Target_Docking/data/jcim_multiseed_v0"
LOG="$MS_ROOT/logs/resume_c1_after_multiseed.log"
PIDFILE="$MS_ROOT/logs/multiseed_vina.pid"
DONE_MARKER="$MS_ROOT/tables/multiseed_scores_long_v1.csv"
STDOUT_LOG="$MS_ROOT/logs/multiseed_vina.stdout"

# Explicit PIDs that were SIGSTOP'd for this JCIM run (still verified alive+T at resume time).
RESUME_PIDS=(12634 12635 18951 18971)

mkdir -p "$MS_ROOT/logs"
exec >>"$LOG" 2>&1
echo "[$(date -Is)] watcher start; waiting for multi-seed finish"

multiseed_running() {
  local pid
  if [[ -f "$PIDFILE" ]]; then
    pid=$(cat "$PIDFILE" 2>/dev/null || true)
    if [[ -n "${pid:-}" ]] && kill -0 "$pid" 2>/dev/null; then
      return 0
    fi
  fi
  pgrep -f 'run_four_pair_multiseed_vina_v1.py' >/dev/null 2>&1
}

multiseed_complete() {
  # Prefer final table + no runner process
  if [[ -f "$DONE_MARKER" ]] && ! multiseed_running; then
    return 0
  fi
  if grep -q 'DONE rows=' "$STDOUT_LOG" 2>/dev/null && ! multiseed_running; then
    return 0
  fi
  return 1
}

while true; do
  if multiseed_complete; then
    echo "[$(date -Is)] multi-seed complete"
    break
  fi
  if ! multiseed_running && [[ ! -f "$DONE_MARKER" ]]; then
    # Runner died without final table — wait a bit more in case of restart, then still resume C1
    echo "[$(date -Is)] multi-seed runner not found and no final table yet; recheck in 120s"
    sleep 120
    if ! multiseed_running; then
      echo "[$(date -Is)] still no runner; proceeding to resume C1 anyway"
      break
    fi
    continue
  fi
  sleep 120
done

echo "[$(date -Is)] resuming paused C1 / GNINA processes"
resumed=0
for pid in "${RESUME_PIDS[@]}"; do
  if kill -0 "$pid" 2>/dev/null; then
    st=$(ps -o stat= -p "$pid" 2>/dev/null | tr -d ' ')
    echo "  CONT $pid (stat=$st)"
    kill -CONT "$pid" 2>/dev/null && resumed=$((resumed + 1)) || echo "  CONT failed $pid"
  else
    echo "  skip $pid (not alive)"
  fi
done

# Also CONT any remaining stopped gnina / c1 dockers (pattern-based safety net)
while read -r pid; do
  [[ -z "$pid" ]] && continue
  st=$(ps -o stat= -p "$pid" 2>/dev/null | tr -d ' ')
  if [[ "$st" == *T* ]]; then
    echo "  CONT pattern-match $pid (stat=$st)"
    kill -CONT "$pid" 2>/dev/null && resumed=$((resumed + 1)) || true
  fi
done < <(pgrep -f 'run_c1_acid_dual_dock|gnina.*URAT1|gnina.*c1-campaign|gnina.*7ALV' 2>/dev/null || true)

sleep 3
echo "[$(date -Is)] post-resume snapshot:"
ps -eo pid,stat,etime,cmd | awk '$2 ~ /T/ || /run_c1_acid|gnina.*7ALV|gnina.*URAT1/ {print}' | head -40
echo "[$(date -Is)] resumed_count≈$resumed ; watcher exit"
