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
GNINA_URL="https://github.com/gnina/gnina/releases/download/v1.3.1/gnina1.3.1"
if [ ! -x tools/gnina ]; then
  if curl -fsSL -o tools/gnina "$GNINA_URL"; then
    chmod +x tools/gnina
  else
    echo ""
    echo ">>> curl 下载 GNINA 失败（可能无法访问 GitHub）。"
    echo ">>> 请在 Windows 浏览器打开："
    echo ">>>   $GNINA_URL"
    echo ">>> 下载后复制到 WSL："
    echo ">>>   mkdir -p tools"
    echo ">>>   cp /mnt/c/Users/你的用户名/Downloads/gnina1.3.1 tools/gnina"
    echo ">>>   chmod +x tools/gnina"
    echo ">>> 然后重新运行: bash scripts/setup_gnina_wsl_cpu.sh"
    exit 1
  fi
fi
tools/gnina --help | head -5 || true

echo "=== [5/6] Prepare receptors ==="
python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
python3 scripts/prepare_receptor_vina.py --target nlrp3_7alv

echo "=== [6/6] Smoke test (lesinurad, CPU) ==="
mkdir -p results/repurposing/smoke_gnina/ligands results/repurposing/smoke_gnina/dock
cat > results/repurposing/smoke_gnina/pool.csv <<'CSV'
repurposing_id,canonical_smiles,name
BENCH_LES,O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12,lesinurad
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
