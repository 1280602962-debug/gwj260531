#!/usr/bin/env bash
# Production dual-target funnel with locked protocol Pi* = P2 (gnina CNNaffinity).
# Docks the NLRP3-ML-reduced clinical pool (docking_pool_p05.csv, ~1588) at
# URAT1 9DKB and NLRP3 7ALV, then merges + builds the Pareto shortlist.
#
# Prereqs: gnina in PATH (or tools/gnina); Python deps installed (rdkit, meeko,
#          gemmi/openbabel for receptor prep, pandas, numpy, pyyaml).
#
# Usage:
#   bash scripts/run_funnel_p2.sh
#   JOBS=8 bash scripts/run_funnel_p2.sh          # more parallel workers
#   POOL=data/repurposing/screening/docking_pool_p05_phase_ge3.csv bash scripts/run_funnel_p2.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

CFG="${CFG:-config/docking_production_p2.yaml}"
POOL="${POOL:-data/repurposing/screening/docking_pool_p05.csv}"
JOBS="${JOBS:-4}"
OUT="${OUT:-results/repurposing}"
LIGDIR="$OUT/ligands_p05"
MLSCORES="${MLSCORES:-data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv}"

echo "=== Funnel P2 (gnina CNNaffinity) ==="
echo "config=$CFG  pool=$POOL  jobs=$JOBS"

# 0) Receptors (idempotent)
[ -f data/structures/prepared/9DKB_receptor.pdbqt ] || python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
[ -f data/structures/prepared/7ALV_receptor.pdbqt ] || python3 scripts/prepare_receptor_vina.py --target nlrp3_7alv

# 1) Ligands (once; reused for both targets)
if [ ! -f "$LIGDIR/ligand_manifest.csv" ]; then
  python3 scripts/prepare_ligands_vina.py --input "$POOL" --output-dir "$LIGDIR"
fi

# 2) Dock at URAT1 9DKB (P2)
python3 scripts/run_gnina_batch.py \
  --config "$CFG" --target urat1_9dkb \
  --manifest "$LIGDIR/ligand_manifest.csv" \
  --output-dir "$OUT/docking_p2/9dkb" --jobs "$JOBS"

# 3) Dock at NLRP3 7ALV (P2)
python3 scripts/run_gnina_batch.py \
  --config "$CFG" --target nlrp3_7alv \
  --manifest "$LIGDIR/ligand_manifest.csv" \
  --output-dir "$OUT/docking_p2/7alv" --jobs "$JOBS"

# 4) Merge -> Pareto shortlist (S_U vs S_N = max(ML, dock))
#    Write into data/repurposing/pareto/ so downstream audit scripts
#    (10_admet, 11_chemical_space, 13_pareto_robustness, 14_candidate_nomination)
#    pick it up via their default --pool/--shortlist paths.
PARETO_OUT="${PARETO_OUT:-data/repurposing/pareto}"
python3 scripts/merge_docking_pareto.py \
  --ml-scores "$MLSCORES" \
  --urat1-dock "$OUT/docking_p2/9dkb/docking_9dkb_gnina.csv" \
  --nlrp3-dock "$OUT/docking_p2/7alv/docking_7alv_gnina.csv" \
  --nlrp3-pdb 7ALV \
  --pool "$POOL" \
  --sn-mode both \
  --output-dir "$PARETO_OUT"

echo ""
echo "Done. Key outputs:"
echo "  $OUT/docking_p2/9dkb/docking_9dkb_gnina.csv"
echo "  $OUT/docking_p2/7alv/docking_7alv_gnina.csv"
echo "  $PARETO_OUT/pareto_merged_scores.csv"
echo "  $PARETO_OUT/pareto_shortlist.csv"
echo "Next: python3 scripts/10_admet_druglikeness.py && python3 scripts/11_chemical_space_novelty.py"
echo "      python3 scripts/13_pareto_robustness.py && python3 scripts/14_candidate_nomination.py"
