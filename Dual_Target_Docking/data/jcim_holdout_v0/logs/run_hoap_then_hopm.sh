#!/usr/bin/env bash
set -u
REPO=/home/gwj/repos/gwj260531/Dual_Target_Docking
LOGDIR=$REPO/data/jcim_holdout_v0/logs
cd "$REPO"
echo "=== HOAP resume $(date -Is) workers=8 ===" | tee -a "$LOGDIR/orchestrator.log"
python3 data/jcim_holdout_v0/scripts/dock_holdout_v1.py --prefix HOAP --workers 8 \
  2>&1 | tee -a "$LOGDIR/dock_HOAP_local.log" | tee -a "$LOGDIR/orchestrator.log"
ec_hoap=${PIPESTATUS[0]}
echo "=== HOAP exit=$ec_hoap $(date -Is) ===" | tee -a "$LOGDIR/orchestrator.log"
echo "=== HOPM start $(date -Is) workers=8 ===" | tee -a "$LOGDIR/orchestrator.log"
python3 data/jcim_holdout_v0/scripts/dock_holdout_v1.py --prefix HOPM --workers 8 \
  2>&1 | tee -a "$LOGDIR/dock_HOPM_local.log" | tee -a "$LOGDIR/orchestrator.log"
ec_hopm=${PIPESTATUS[0]}
echo "=== HOPM exit=$ec_hopm $(date -Is) ===" | tee -a "$LOGDIR/orchestrator.log"
echo "ALL_DONE hoap=$ec_hoap hopm=$ec_hopm $(date -Is)" | tee -a "$LOGDIR/orchestrator.log"
