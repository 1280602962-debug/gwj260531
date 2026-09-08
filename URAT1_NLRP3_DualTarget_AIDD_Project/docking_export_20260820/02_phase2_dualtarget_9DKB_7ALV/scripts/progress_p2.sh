#!/usr/bin/env bash
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
n9=$(find "$ROOT/results/repurposing/docking_p2/9dkb/poses" -name '*_out.sdf' -size +0 2>/dev/null | wc -l)
n7=$(find "$ROOT/results/repurposing/docking_p2/7alv/poses" -name '*_out.sdf' -size +0 2>/dev/null | wc -l)
echo "p2_9DKB=$n9/1583  p2_7ALV=$n7/1583"
pgrep -af 'run_gnina_batch|run_server_gpu|gnina_5090|gnina_install/bin/gnina' | grep -v pgrep | head -8 || true
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader
uptime | sed 's/^/load: /'
