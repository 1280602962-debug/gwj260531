#!/usr/bin/env bash
# Environment checker for Maestro-prepped docking package.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -f "${ROOT}/config.sh" ]]; then
  # shellcheck disable=SC1091
  source "${ROOT}/config.sh"
else
  # shellcheck disable=SC1091
  source "${ROOT}/config.example.sh"
fi

ok=0; warn=0; fail=0
pass() { echo "[PASS] $*"; ok=$((ok+1)); }
warn_() { echo "[WARN] $*"; warn=$((warn+1)); }
fail_() { echo "[FAIL] $*"; fail=$((fail+1)); }

echo "=== Host ==="; hostname; date; echo "CPU(s): $(nproc 2>/dev/null || echo ?)"
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -3 || echo "nvidia-smi: none"
echo
echo "=== Timeout helper ==="
if command -v "${TIMEOUT_BIN}" >/dev/null 2>&1; then
  pass "timeout binary: ${TIMEOUT_BIN} (VINA=${VINA_TIMEOUT_SEC}s GNINA=${GNINA_TIMEOUT_SEC}s)"
else
  fail_ "GNU timeout not found (TIMEOUT_BIN=${TIMEOUT_BIN}) — required for skip-on-timeout"
fi

echo; echo "=== Binaries ==="
if command -v "${VINA_BIN}" >/dev/null 2>&1 || [[ -x "${VINA_BIN}" ]]; then
  pass "vina: $("${VINA_BIN}" --version 2>&1 | head -1)"
else fail_ "vina missing"; fi
if command -v "${GNINA_BIN}" >/dev/null 2>&1 || [[ -x "${GNINA_BIN}" ]]; then
  pass "gnina present"
else fail_ "gnina missing"; fi
command -v "${OBABEL_BIN}" >/dev/null 2>&1 && pass "obabel" || warn_ "obabel missing (needed for vina→RTM SDF)"

echo; echo "=== Inputs / prepared ligands ==="
for f in RECEPTOR_PDB RECEPTOR_PDBQT AUTOBOX_LIGAND; do
  [[ -f "${!f}" ]] && pass "$f" || fail_ "$f missing: ${!f}"
done
nsdf=$(ls "${WORK_DIR}/ligands_sdf"/mol_*.sdf 2>/dev/null | wc -l)
npdbqt=$(ls "${WORK_DIR}/ligands_pdbqt"/mol_*.pdbqt 2>/dev/null | wc -l)
echo "ligands_sdf=${nsdf} ligands_pdbqt=${npdbqt}"
(( nsdf > 0 )) && pass "per-mol SDF ready" || fail_ "no ligands_sdf — run scripts/split_ligands_meeko.py"
(( npdbqt > 0 )) && pass "per-mol PDBQT ready" || fail_ "no ligands_pdbqt"

echo; echo "=== Summary PASS=${ok} WARN=${warn} FAIL=${fail} ==="
(( fail > 0 )) && exit 1 || exit 0
