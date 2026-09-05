#!/usr/bin/env bash
# C5 docking handoff after Track B Vina finishes.
set -uo pipefail

PROJECT="/home/gwj/work/c1-campaign/URAT1_NLRP3_DualTarget_AIDD_Project"
TRACK_B_PID="${TRACK_B_PID:-34308}"
TRACK_B_PGREP='dock_track_b_production_v1.py'
LOG_DIR="$PROJECT/data/campaigns/c5/00_verification"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/c5_handoff_$(date +%Y%m%d_%H%M%S).log"
CPU="${C5_CPU:-4}"

exec > >(tee -a "$MASTER_LOG") 2>&1

echo "==== C5 handoff start $(date -Is) ===="
echo "project=$PROJECT"
echo "track_b_pid=$TRACK_B_PID"
echo "gnina=$("$PROJECT/tools/gnina" --version 2>&1 | head -1)"
echo "NOTE: worklist locks GNINA 1.3.1; local binary version logged above."

wait_track_b() {
  echo "Waiting for Track B docking to finish..."
  while true; do
    alive=0
    if kill -0 "$TRACK_B_PID" 2>/dev/null; then alive=1; fi
    if pgrep -f "$TRACK_B_PGREP" >/dev/null 2>&1; then alive=1; fi
    if pgrep -af 'jcim_chembl_universe_v0/local_track_b_v0' | grep -q '[v]ina'; then alive=1; fi
    if [[ "$alive" -eq 0 ]]; then
      sleep 15
      if ! pgrep -af 'jcim_chembl_universe_v0/local_track_b_v0' | grep -q '[v]ina' \
         && ! pgrep -f "$TRACK_B_PGREP" >/dev/null 2>&1; then
        echo "Track B clear at $(date -Is)"
        return 0
      fi
    fi
    prog=$(rg -o '\[[0-9]+/1100\]' \
      "/mnt/d/CADD paper exercise/dual target docking/repo_Dual_Target_Docking/Dual_Target_Docking/data/jcim_chembl_universe_v0/local_track_b_v0/logs/production_vina_run.log" \
      2>/dev/null | tail -1 || true)
    echo "$(date +%H:%M:%S) still waiting... last_progress=$prog vina=$(pgrep -c vina || echo 0)"
    sleep 120
  done
}

cd "$PROJECT"
export PATH="/home/gwj/miniconda3/bin:$PATH"
PY=/home/gwj/miniconda3/bin/python

wait_track_b

echo "==== Task1 W1 gate ===="
set +e
$PY scripts/run_c5_w1_crossdock.py --phase gate --cpu "$CPU"
gate_rc=$?
set -e
if [[ $gate_rc -ne 0 ]]; then
  echo "GATE failed/aborted (rc=$gate_rc). STOP. No Task2/Task3."
  exit 2
fi
$PY - <<'PY2'
import json, sys
from pathlib import Path
g = json.loads(Path('data/campaigns/c5/01_crossdock/gate_benzbromarone_9dka.json').read_text())
print('gate.pass=', g.get('pass'))
sys.exit(0 if g.get('pass') else 2)
PY2

echo "==== Task1 W1 remaining cells ===="
$PY scripts/run_c5_w1_crossdock.py --phase rest --cpu "$CPU"

echo "==== Task2 W4 decoy lock ===="
if [[ ! -f data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv ]]; then
  $PY scripts/sample_c5_w4_decoys.py --n 40 --seed 0xC5DEC0
else
  echo "decoys already locked"
fi

echo "==== Task2 W4 NLRP3 panel docking ===="
$PY scripts/run_c5_w4_nlrp3_panel.py --cpu "$CPU"

echo "==== Task3 W2 IFP rescoring (0 new dock) ===="
if [[ -f scripts/run_acid_gate_benchmark.py ]]; then
  set +e
  $PY scripts/run_acid_gate_benchmark.py 2>&1 | tee "$LOG_DIR/w2_acid_gate_rescore.log"
  set -e
fi

echo "==== C5 handoff complete $(date -Is) ===="
echo "master_log=$MASTER_LOG"
