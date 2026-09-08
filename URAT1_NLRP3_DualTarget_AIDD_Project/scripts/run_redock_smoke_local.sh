#!/usr/bin/env bash
# Local redock smoke test: lesinurad (+ optional controls) @ 9DKB with Vina and gnina.
# Usage:
#   bash scripts/run_redock_smoke_local.sh
#   EXHAUST=32 COMPOUNDS=lesinurad bash scripts/run_redock_smoke_local.sh
#   SKIP_GNINA=1 bash scripts/run_redock_smoke_local.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

EXHAUST="${EXHAUST:-8}"
NUM_MODES="${NUM_MODES:-9}"
JOBS="${JOBS:-2}"
COMPOUNDS="${COMPOUNDS:-lesinurad}"   # comma-separated names from pool, or "all"
SKIP_GNINA="${SKIP_GNINA:-0}"
OUT="${OUT:-$ROOT/results/redock_smoke}"
POOL="$ROOT/data/redock_smoke/redock_pool.csv"
LIGDIR="$OUT/ligands"
CFG="$ROOT/config/docking_open_source.yaml"

mkdir -p "$OUT" "$LIGDIR"

echo "=== Redock smoke @ 9DKB ==="
echo "exhaustiveness=$EXHAUST  num_modes=$NUM_MODES  compounds=$COMPOUNDS"
echo "output: $OUT"

# Filter pool to URAT1 compounds (and optionally subset by name)
FILTERED="$OUT/pool_run.csv"
python3 - <<PY
import csv
from pathlib import Path
names = [x.strip().lower() for x in "$COMPOUNDS".split(",") if x.strip()]
rows = []
with open("$POOL", newline="") as f:
    for r in csv.DictReader(f):
        if r.get("receptor_hint") != "urat1_9dkb":
            continue
        if names != ["all"] and r.get("name", "").lower() not in names:
            continue
        rows.append(r)
if not rows:
    raise SystemExit("No compounds selected; check COMPOUNDS= or pool CSV")
with open("$FILTERED", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r["repurposing_id"], r["name"], r["role"])
PY

# Receptor
if [[ ! -f data/structures/prepared/9DKB_receptor.pdbqt ]]; then
  python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
fi

# Ligands
python3 scripts/prepare_ligands_vina.py \
  --input "$FILTERED" \
  --output-dir "$LIGDIR"

# --- Vina (P1 poses; also input for P4) ---
if ! command -v vina >/dev/null 2>&1; then
  echo "WARNING: 'vina' not in PATH. Install AutoDock Vina 1.2.x or set PATH."
  echo "Skipping Vina. See docs/REDOCK_SMOKE_TEST_SHEET.md"
else
  echo "=== Vina batch ==="
  # Temporarily override exhaustiveness / num_modes via env-friendly edit of a run-local cfg
  RUN_CFG="$OUT/docking_smoke.yaml"
  python3 - <<PY
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path("$CFG").read_text())
cfg["vina"]["exhaustiveness"] = int("$EXHAUST")
cfg["vina"]["num_modes"] = int("$NUM_MODES")
cfg["gnina"]["exhaustiveness"] = int("$EXHAUST")
cfg["gnina"]["num_modes"] = int("$NUM_MODES")
cfg["gnina"]["score_mode"] = "affinity"
Path("$RUN_CFG").write_text(yaml.dump(cfg, sort_keys=False))
print("wrote", "$RUN_CFG")
PY
  python3 scripts/run_vina_batch.py \
    --config "$RUN_CFG" \
    --target urat1_9dkb \
    --manifest "$LIGDIR/ligand_manifest.csv" \
    --output-dir "$OUT/vina" \
    --jobs "$JOBS"
  echo "Vina CSV: $OUT/vina/docking_9dkb_vina.csv"
fi

# --- gnina (P0/P2/P3 poses; also input for P5) ---
if [[ "$SKIP_GNINA" == "1" ]]; then
  echo "SKIP_GNINA=1 — skipping gnina"
elif [[ -x "$ROOT/tools/gnina" ]] || command -v gnina >/dev/null 2>&1; then
  echo "=== gnina batch (CPU, CNN rescore) ==="
  RUN_CFG="${RUN_CFG:-$OUT/docking_smoke.yaml}"
  if [[ ! -f "$RUN_CFG" ]]; then
    python3 - <<PY
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path("$CFG").read_text())
cfg["vina"]["exhaustiveness"] = int("$EXHAUST")
cfg["vina"]["num_modes"] = int("$NUM_MODES")
cfg["gnina"]["exhaustiveness"] = int("$EXHAUST")
cfg["gnina"]["num_modes"] = int("$NUM_MODES")
Path("$RUN_CFG").write_text(yaml.dump(cfg, sort_keys=False))
PY
  fi
  python3 scripts/run_gnina_batch.py \
    --config "$RUN_CFG" \
    --target urat1_9dkb \
    --manifest "$LIGDIR/ligand_manifest.csv" \
    --output-dir "$OUT/gnina" \
    --jobs "$JOBS"
  echo "gnina CSV: $OUT/gnina/docking_9dkb_gnina.csv"
  echo "Tip: open the .log under $OUT/gnina/poses/ for CNNscore / CNNaffinity columns (P0/P2/P3)."
else
  echo "WARNING: gnina not found (tools/gnina or PATH). Skipping."
  echo "Install: bash scripts/setup_gnina_wsl_cpu.sh"
fi

# Merge score summary for the sheet
python3 - <<'PY'
import csv
from pathlib import Path

out = Path("results/redock_smoke")
rows = []

def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

for r in read_csv(out / "vina" / "docking_9dkb_vina.csv"):
    rows.append({
        "protocol_hint": "P1",
        "engine": "vina",
        "compound_id": r.get("repurposing_id", ""),
        "score_name": "vina_affinity",
        "score": r.get("dock_score") or r.get("minimizedAffinity", ""),
        "status": r.get("status", ""),
        "pose_file": r.get("pose_pdbqt") or r.get("out_pdbqt", ""),
    })
for r in read_csv(out / "gnina" / "docking_9dkb_gnina.csv"):
    rows.append({
        "protocol_hint": "P3_default_csv",
        "engine": "gnina",
        "compound_id": r.get("repurposing_id", ""),
        "score_name": "gnina_primary_csv",
        "score": r.get("dock_score", ""),
        "status": r.get("status", ""),
        "pose_file": r.get("pose_sdf") or r.get("out_sdf", ""),
    })
summary = out / "scores_summary.csv"
fields = ["protocol_hint", "engine", "compound_id", "score_name", "score", "status", "pose_file"]
with summary.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {summary} ({len(rows)} rows)")
print("Fill RMSD: cp data/redock_smoke/redock_results_template.csv results/redock_smoke/redock_results_filled.csv")
print("See docs/REDOCK_SMOKE_TEST_SHEET.md")
PY

echo ""
echo "Done. Next:"
echo "  1) cp data/redock_smoke/redock_results_template.csv results/redock_smoke/redock_results_filled.csv"
echo "  2) Measure Top-1 / Best RMSD vs 9DKB A1AIL in PyMOL (heavy atoms)"
echo "  3) For formal gate, re-run: EXHAUST=32 bash scripts/run_redock_smoke_local.sh"
