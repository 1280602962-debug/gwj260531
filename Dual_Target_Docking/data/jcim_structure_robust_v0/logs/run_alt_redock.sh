#!/usr/bin/env bash
set -u
REPO=/home/gwj/repos/gwj260531/Dual_Target_Docking
LOG=$REPO/data/jcim_structure_robust_v0/logs
cd "$REPO"
echo "=== 4JPS redock $(date -Is) ===" | tee -a "$LOG/alt_redock.log"
python3 data/jcim_structure_robust_v0/scripts/redock_pm48_alt_pik3ca_v1.py --alt 4JPS --workers 8 \
  2>&1 | tee -a "$LOG/redock_4JPS.log" | tee -a "$LOG/alt_redock.log"
echo "=== 5DXT redock $(date -Is) ===" | tee -a "$LOG/alt_redock.log"
python3 data/jcim_structure_robust_v0/scripts/redock_pm48_alt_pik3ca_v1.py --alt 5DXT --workers 8 \
  2>&1 | tee -a "$LOG/redock_5DXT.log" | tee -a "$LOG/alt_redock.log"
echo "=== ALT_REDOCK_ALL_DONE $(date -Is) ===" | tee -a "$LOG/alt_redock.log"
