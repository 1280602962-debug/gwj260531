#!/usr/bin/env bash
# Minimal local reproducibility entry point for DualFourClass-Bench.
#
# Rebuilds the zero-dock analysis artifacts described in the top-level README.
# It does not attempt docking, GNINA rescoring, or RTM rescoring.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-python3}"

echo "== DualFourClass-Bench local reproducibility =="
echo "repo: ${ROOT}"
echo "python: ${PYTHON_BIN}"
echo

cd "${ROOT}"

echo "[1/4] Checking local Python environment"
"${PYTHON_BIN}" scripts/check_local_env.py

echo
echo "[2/4] Rebuilding pocket-matched diagnostics"
"${PYTHON_BIN}" data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py

echo
echo "[3/4] Rebuilding T0 strengthen analysis"
"${PYTHON_BIN}" data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py

echo
echo "[4/4] Rebuilding forest / gate figures"
"${PYTHON_BIN}" data/jcim_bench_v0/scripts/plot_forest_ci_v1.py

echo
echo "Optional phase-1 revision analyses: bash scripts/run_phase1_revision.sh"

echo
echo "Done."
echo "Key outputs:"
echo "  - data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv  (Table 2 CIs)"
echo "  - data/jcim_strengthen_t0t1_v0/analysis/PRIMARY_METRIC_V2.md  (deprecated CI source; do not cite)"
echo "  - data/jcim_bench_v0/tables/forest_summary_min_ci_v1.csv  (vina_mean forest, not Table 2)"
