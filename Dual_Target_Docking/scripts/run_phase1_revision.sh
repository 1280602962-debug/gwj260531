#!/usr/bin/env bash
# Phase-1 zero-dock revision analyses (no Vina / GNINA).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-python3}"

cd "${ROOT}"
echo "== Phase-1 revision analyses =="
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/document_blocked_cv_v1.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/assay_context_audit_v1.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/time_split_validation_v1.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/cognate_artifact_inventory_v1.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/build_master_results_table_v1.py
"${PYTHON_BIN}" docs/assemble_manuscript_en.py
"${PYTHON_BIN}" docs/assemble_manuscript_zh.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py
"${PYTHON_BIN}" data/jcim_novelty_v0/scripts/validate_revision_v1.py
echo "Done."
