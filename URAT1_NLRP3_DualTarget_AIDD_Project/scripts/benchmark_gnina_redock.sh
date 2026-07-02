#!/usr/bin/env bash
# GNINA benchmark redock @ 9DKB + 7ALV (uses local WSL gnina install)
# Usage:
#   export GNINA_ROOT="/mnt/d/CADD paper exercise/gnina"
#   export PROJECT_ROOT="/mnt/d/CADD paper exercise/gwj260531/URAT1_NLRP3_DualTarget_AIDD_Project"
#   bash scripts/benchmark_gnina_redock.sh
set -euo pipefail

GNINA_ROOT="${GNINA_ROOT:-/mnt/d/CADD paper exercise/gnina}"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
OUT_BASE="${OUT_BASE:-$GNINA_ROOT/output/benchmark}"
SUMMARY="${SUMMARY:-$PROJECT_ROOT/results/gnina_benchmark/benchmark_redock_summary.csv}"
# Glide-XP–like GNINA: thorough search + single best pose (see docs/GNINA_BENCHMARK_REDOCK_WSL.md)
EXHAUST="${EXHAUST:-32}"
NUM_MODES="${NUM_MODES:-1}"
CNN_SCORING="${CNN_SCORING:-rescore}"

source "$GNINA_ROOT/activate.sh"

mkdir -p "$(dirname "$SUMMARY")" "$OUT_BASE"

# Convert 9DKB CIF -> PDB if needed
CIF="$PROJECT_ROOT/data/structures/pdb/9DKB.cif"
PDB_9DKB="$OUT_BASE/9DKB.pdb"
PDB_7ALV="$PROJECT_ROOT/data/structures/pdb/7ALV.pdb"
if [[ -f "$CIF" && ! -f "$PDB_9DKB" ]]; then
  obabel "$CIF" -O "$PDB_9DKB" 2>/dev/null || python3 -c "
import gemmi
st = gemmi.read_structure('$CIF')
st.write_pdb('$PDB_9DKB')
"
fi

echo "target,pdb_id,compound,affinity_kcal_mol,cnn_pose_score,cnn_affinity,pose_file,log_file,status" > "$SUMMARY"

dock_one() {
  local target="$1" pdb="$2" compound="$3" resname="$4" smiles="$5" receptor_pdb="$6"
  local prep="$OUT_BASE/${pdb,,}/prepare_${compound}"
  local dockdir="$OUT_BASE/${pdb,,}/dock_${compound}"
  mkdir -p "$prep" "$dockdir"
  local log="$dockdir/dock.log"
  local out_sdf="$dockdir/docked.sdf"

  echo "=== $pdb / $compound ==="
  python "$GNINA_ROOT/scripts/prepare_docking.py" \
    --receptor "$receptor_pdb" \
    --ligand-resname "$resname" \
    --ligand-smiles "$smiles" \
    --out-dir "$prep" || { echo "$target,$pdb,$compound,,,,,$out_sdf,$log,prepare_fail" >> "$SUMMARY"; return; }

  local r="$prep/receptor.pdb"
  local l="$prep/query_ligand.sdf"
  local ref="$prep/ref_ligand.sdf"
  for f in "$r" "$l" "$ref"; do
    [[ -f "$f" ]] || { echo "$target,$pdb,$compound,,,,,$out_sdf,$log,missing_$f" >> "$SUMMARY"; return; }
  done

  if [[ -x "$GNINA_ROOT/bin/gnina" ]]; then
    GNINA_BIN="$GNINA_ROOT/bin/gnina"
  else
    GNINA_BIN="gnina"
  fi

  "$GNINA_BIN" -r "$r" -l "$l" \
    --autobox_ligand "$ref" \
    -o "$out_sdf" \
    --exhaustiveness "$EXHAUST" --num_modes "$NUM_MODES" \
    --cnn_scoring "$CNN_SCORING" --no_gpu \
    --log "$log" 2>&1 | tee "$dockdir/stdout.txt" || true

  # Parse mode 1 from log (best effort)
  aff="$(awk '/^[[:space:]]*1[[:space:]]+-/{print $2; exit}' "$log" 2>/dev/null || true)"
  cnn="$(awk '/^[[:space:]]*1[[:space:]]+-/{print $3; exit}' "$log" 2>/dev/null || true)"
  cnnaff="$(awk '/^[[:space:]]*1[[:space:]]+-/{print $4; exit}' "$log" 2>/dev/null || true)"
  local st="docked"
  [[ -n "$aff" ]] || st="parse_fail"
  echo "$target,$pdb,$compound,$aff,$cnn,$cnnaff,$out_sdf,$log,$st" >> "$SUMMARY"
}

# URAT1 @ 9DKB (autobox A1AIL)
SMILES_LES='O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12'
SMILES_BENZ='CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1'
SMILES_DOT='O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21'
SMILES_EGCG='O=C(O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1'
SMILES_MCC='CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1'

dock_one URAT1 9DKB lesinurad A1AIL "$SMILES_LES" "$PDB_9DKB"
dock_one URAT1 9DKB benzbromarone A1AIL "$SMILES_BENZ" "$PDB_9DKB"
dock_one URAT1 9DKB dotinurad A1AIL "$SMILES_DOT" "$PDB_9DKB"
dock_one URAT1 9DKB EGCG A1AIL "$SMILES_EGCG" "$PDB_9DKB"

# NLRP3 @ 7ALV (autobox RM5)
dock_one NLRP3 7ALV MCC950 RM5 "$SMILES_MCC" "$PDB_7ALV"
dock_one NLRP3 7ALV EGCG RM5 "$SMILES_EGCG" "$PDB_7ALV"

echo "Done. Summary: $SUMMARY"
column -t -s, "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
