#!/usr/bin/env bash
# GNINA fair best-of-9 rescore (all Vina modes) — fixes RTM vs GNINA pose asymmetry
set -euo pipefail
RES="/mnt/d/CADD paper exercise/dual target docking/results"
REPO="/home/gwj/repos/gwj260531/Dual_Target_Docking/data"
PY="/home/gwj/miniconda3/bin/python"
GNINA_PY="$REPO/jcim_bench_v0/scripts/gnina_rescore_panel.py"
LOG="$RES/jcim_gnina_best9.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== START best9 $(date -Iseconds) ==="
source "/mnt/d/CADD paper exercise/gnina/activate.sh"

run_pack() {
  local name="$1"; shift
  echo "=== GNINA best9 $name $(date -Iseconds) ==="
  "$PY" -u "$GNINA_PY" "$@"
  echo "=== GNINA best9 $name DONE $(date -Iseconds) ==="
}

WORKERS=7
TIMEOUT=180

run_pack ache_bche \
  --root "$RES/ache_bche_panel_v0" \
  --targets ACHE BCHE \
  --receptor-map \
    ACHE="$RES/ache_bche_panel_v0/receptors/ACHE_protein.pdb" \
    BCHE="$RES/ache_bche_panel_v0/receptors/BCHE_protein.pdb" \
  --workers "$WORKERS" --modes all --timeout "$TIMEOUT"

run_pack pik3ca_pik3cb \
  --root "$RES/pik3ca_pik3cb_panel_v0" \
  --targets PIK3CA PIK3CB \
  --receptor-map \
    PIK3CA="$RES/pik3ca_pik3cb_panel_v0/receptors/PIK3CA_protein.pdb" \
    PIK3CB="$RES/pik3ca_pik3cb_panel_v0/receptors/PIK3CB_protein.pdb" \
  --workers "$WORKERS" --modes all --timeout "$TIMEOUT"

run_pack pm48 \
  --root "$RES/pik3ca_mtor_panel48_rdkit_v0" \
  --targets 4L23 4JT6 \
  --receptor-map \
    4L23="$RES/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_protein.pdb" \
    4JT6="$RES/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_protein.pdb" \
  --workers "$WORKERS" --modes all --timeout "$TIMEOUT"

EH="$RES/egfr_her2_panel120_v0"
run_pack egfr \
  --root "$EH" \
  --targets 3POZ 3RCD \
  --receptor-map \
    3POZ="$EH/receptors/3POZ_protein.pdb" \
    3RCD="$EH/receptors/3RCD_protein.pdb" \
  --workers "$WORKERS" --modes all --timeout "$TIMEOUT"

# Optional PM110 stability pack
if [[ -d "$RES/pik3ca_mtor_panel110_rdkit_v0/poses" ]]; then
  run_pack pm110 \
    --root "$RES/pik3ca_mtor_panel110_rdkit_v0" \
    --targets 4L23 4JT6 \
    --receptor-map \
      4L23="$RES/pik3ca_mtor_panel110_rdkit_v0/receptors/4L23_protein.pdb" \
      4JT6="$RES/pik3ca_mtor_panel110_rdkit_v0/receptors/4JT6_protein.pdb" \
    --workers "$WORKERS" --modes all --timeout "$TIMEOUT"
fi

# Sync tables into repo packs
for pack in ache_bche_panel_v0 pik3ca_pik3cb_panel_v0 pik3ca_mtor_panel48_rdkit_v0 egfr_her2_panel120_v0 pik3ca_mtor_panel110_rdkit_v0; do
  mkdir -p "$REPO/$pack/tables"
  for f in scores_gnina_long.csv scores_gnina_best.csv scores_gnina_long_mode01_backup.csv scores_gnina_best_mode01_backup.csv; do
    [[ -f "$RES/$pack/tables/$f" ]] && cp -f "$RES/$pack/tables/$f" "$REPO/$pack/tables/" && echo "synced $pack/$f"
  done
done

echo "=== JCIM GNINA BEST9 ALL DONE $(date -Iseconds) ==="
echo JCIM_GNINA_BEST9_DONE > "$RES/JCIM_GNINA_BEST9_DONE.flag"
