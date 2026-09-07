#!/usr/bin/env bash
# Submit all-pairs stack. NEVER stop the chain because a step returned non-zero
# (timeout→skip and ligand parse fails are expected to be recorded, not fatal).
set +e
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
WORKERS="${WORKERS:-6}"
TIMEOUT="${TIMEOUT:-600}"
START_FROM="${START_FROM:-A}"  # A|B|C|D|E
export LD_LIBRARY_PATH="/mnt/d/CADD paper exercise/gnina/conda_env/lib:${LD_LIBRARY_PATH:-}"

LOGDIR="$ROOT/local_track_b_v0/allpairs_stack/logs"
mkdir -p "$LOGDIR"
TS="$(date +%Y%m%d_%H%M%S)"
LOG="$LOGDIR/submit_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== ALLPAIRS START $(date -Is) workers=$WORKERS timeout=${TIMEOUT}s start_from=$START_FROM (continuous; no abort) ==="

run_step() {
  local name="$1"; shift
  echo "=== ${name} BEGIN $(date -Is) ==="
  "$@"
  local rc=$?
  echo "=== ${name} END rc=${rc} $(date -Is) (continuing) ==="
  return 0
}

should_run() {
  # START_FROM gate: A runs all; B skips A; C skips A+B; etc.
  local step="$1"
  case "$START_FROM" in
    A) return 0 ;;
    B) [[ "$step" != A ]] ;;
    C) [[ "$step" == C || "$step" == D || "$step" == E ]] ;;
    D) [[ "$step" == D || "$step" == E ]] ;;
    E) [[ "$step" == E ]] ;;
    *) return 0 ;;
  esac
}

if should_run A; then
  run_step "A GNINA independent JAK1/TYK2" \
    python3 scripts/dock_track_b_gnina_independent_v1.py --workers "$WORKERS" --timeout "$TIMEOUT"
fi
if should_run B; then
  run_step "B five-seed Vina" \
    python3 scripts/dock_track_b_fiveseed_v1.py --workers "$WORKERS" --timeout "$TIMEOUT"
fi
if should_run C; then
  run_step "C RTM rescore" \
    python3 scripts/rescore_track_b_rtm_v1.py
fi
if should_run D; then
  run_step "D GNINA CNN rescore" \
    python3 scripts/rescore_track_b_gnina_cnn_v1.py --workers "$WORKERS" --timeout 180
fi
if should_run E; then
  run_step "E Holdout Vina" \
    python3 scripts/dock_track_b_holdout_v1.py --workers "$WORKERS" --timeout "$TIMEOUT"
fi

echo "=== ALLPAIRS FINISHED $(date -Is) log=$LOG ==="
exit 0
