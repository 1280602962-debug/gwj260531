#!/usr/bin/env bash
# Local B5 runner: PIK3CA/PIK3CB panel into 4JPS then 5DXT (keep 2WXF scores).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
PYTHON="${PYTHON:-python3}"
cd "$ROOT"
echo "repo=$ROOT vina=$(command -v vina || true)"
"$PYTHON" data/jcim_structure_robust_v0/scripts/redock_pik3cb_alt_pik3ca_v1.py --alt 4JPS --workers "${WORKERS:-6}"
"$PYTHON" data/jcim_structure_robust_v0/scripts/redock_pik3cb_alt_pik3ca_v1.py --alt 5DXT --workers "${WORKERS:-6}"
echo "B5 PIK3CA/PIK3CB receptor swap finished"
