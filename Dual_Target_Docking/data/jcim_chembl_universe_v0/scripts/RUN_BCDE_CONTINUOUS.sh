#!/usr/bin/env bash
# Run B→C→D→E without stopping on non-zero exit (timeout/skip/fail recorded in CSVs).
# Idempotent enough: fiveseed skips exists; rescores may re-run safely.
set +e
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
WORKERS="${WORKERS:-6}"
TIMEOUT="${TIMEOUT:-600}"
export LD_LIBRARY_PATH="/mnt/d/CADD paper exercise/gnina/conda_env/lib:${LD_LIBRARY_PATH:-}"

LOGDIR="$ROOT/local_track_b_v0/allpairs_stack/logs"
mkdir -p "$LOGDIR"
TS="$(date +%Y%m%d_%H%M%S)"
LOG="$LOGDIR/run_bcde_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== BCDE CONTINUOUS START $(date -Is) workers=$WORKERS timeout=${TIMEOUT}s (never abort on step rc) log=$LOG ==="

run_step() {
  local name="$1"; shift
  echo "=== ${name} BEGIN $(date -Is) ==="
  "$@"
  local rc=$?
  echo "=== ${name} END rc=${rc} $(date -Is) (continuing) ==="
  return 0
}

run_step "B five-seed Vina" \
  python3 scripts/dock_track_b_fiveseed_v1.py --workers "$WORKERS" --timeout "$TIMEOUT"

run_step "C RTM rescore" \
  python3 scripts/rescore_track_b_rtm_v1.py

run_step "D GNINA CNN rescore" \
  python3 scripts/rescore_track_b_gnina_cnn_v1.py --workers "$WORKERS" --timeout 180

run_step "E Holdout Vina" \
  python3 scripts/dock_track_b_holdout_v1.py --workers "$WORKERS" --timeout "$TIMEOUT"

echo "=== BCDE CONTINUOUS FINISHED $(date -Is) ==="
exit 0
