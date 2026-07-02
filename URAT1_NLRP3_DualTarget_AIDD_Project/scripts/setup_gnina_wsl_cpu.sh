#!/usr/bin/env bash
# Setup GNINA (CPU-only) on WSL2 Ubuntu for URAT1_NLRP3_DualTarget_AIDD_Project
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "=== [1/6] Check WSL / Ubuntu ==="
uname -a
if ! grep -qi microsoft /proc/version 2>/dev/null; then
  echo "Note: not detected as WSL; script still works on native Linux."
fi

echo "=== [2/6] System packages ==="
sudo apt-get update -qq
sudo apt-get install -y -qq curl wget ca-certificates python3 python3-pip python3-venv \
  build-essential git

echo "=== [3/6] Python deps ==="
pip3 install --user -r requirements.txt

echo "=== [4/6] Download GNINA binary (Ubuntu 22.04 x86_64) ==="
mkdir -p tools
GNINA_URL="https://github.com/gnina/gnina/releases/download/v1.3.1/gnina-1.3.1-x86_64-ubuntu22.04"
if [ ! -x tools/gnina ]; then
  curl -fsSL -o tools/gnina "$GNINA_URL" || {
    echo "Direct download failed. Try manually from https://github.com/gnina/gnina/releases"
    exit 1
  }
  chmod +x tools/gnina
fi
tools/gnina --help | head -5 || true

echo "=== [5/6] Prepare receptors ==="
python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
python3 scripts/prepare_receptor_vina.py --target nlrp3_7alv

echo "=== [6/6] Smoke test (lesinurad, CPU) ==="
mkdir -p results/repurposing/smoke_gnina/ligands results/repurposing/smoke_gnina/dock
cat > results/repurposing/smoke_gnina/pool.csv <<'CSV'
repurposing_id,canonical_smiles,name
BENCH_LES,CC(C)(C)OC(=O)N[C@@H](CS)c1nc(-c2cccc3ccccc23)no1,lesinurad
CSV
python3 scripts/prepare_ligands_vina.py \
  --input results/repurposing/smoke_gnina/pool.csv \
  --output-dir results/repurposing/smoke_gnina/ligands
python3 scripts/run_gnina_batch.py \
  --target urat1_9dkb \
  --manifest results/repurposing/smoke_gnina/ligands/ligand_manifest.csv \
  --output-dir results/repurposing/smoke_gnina/dock \
  --jobs 1 --limit 1
cat results/repurposing/smoke_gnina/dock/docking_9dkb_gnina.csv

echo ""
echo "Setup complete. GNINA binary: $ROOT/tools/gnina"
echo "Full pipeline: bash scripts/run_gnina_docking_pipeline.sh"
