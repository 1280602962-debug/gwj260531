#!/usr/bin/env bash
set -euo pipefail
source "/mnt/d/CADD paper exercise/gnina/activate.sh"
ROOT="/mnt/d/CADD paper exercise/NLRP3_URAT1/URAT1 docking methond test"
RD="$ROOT/redock_9DKB_lesinurad"
PREP="$RD/prep"
"${GNINA_BIN}" --no_gpu \
  -r "$PREP/9DKB_receptor.pdb" \
  -l "$ROOT/lesinurad_prepared.sdf" \
  --autobox_ligand "$PREP/lesinurad_crystal_ref.sdf" \
  -o "$RD/gnina/lesinurad_exh8_out.sdf" \
  --exhaustiveness 8 \
  --num_modes 9 \
  --cnn_scoring rescore \
  --seed 42 \
  --cpu 8
