#!/usr/bin/env python3
"""Rerun ECFP4+docking incremental tables for EGFR/HER2 and AChE/BChE.

Ligand-only ECFP4/descriptor models are recomputed only for AChE/BChE because
panel membership changed. EGFR ligand-only chemistry is unchanged; docking
features are updated.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SCRIPT = ROOT / "data/jcim_novelty_v0/scripts/benchmark_formulation_v1.py"
TAB = ROOT / "remediation_outputs/canonical_tables/jcim_novelty_v0"


def main() -> int:
    TAB.mkdir(parents=True, exist_ok=True)
    code = SCRIPT.read_text(encoding="utf-8")
    old = 'TAB = OUT / "tables"'
    new = f'TAB = Path("{TAB}"); TAB.mkdir(parents=True, exist_ok=True)'
    if old not in code:
        raise SystemExit(f"TAB assignment not found in {SCRIPT}")
    ns = {"__name__": "__main__", "__file__": str(SCRIPT)}
    exec(compile(code.replace(old, new, 1), str(SCRIPT), "exec"), ns)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
