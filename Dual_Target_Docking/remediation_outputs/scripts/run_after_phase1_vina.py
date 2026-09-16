#!/usr/bin/env python3
"""Run Phase 1 metrics + AChE holdout dock after corrected-box Vina finishes."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SCRIPTS = ROOT / "remediation_outputs/scripts"
SCORES = ROOT / "remediation_outputs/phase1_vina/scores_vina_mode1_corrected_box.csv"


def run(script: str) -> int:
    cmd = [sys.executable, str(SCRIPTS / script)]
    print("RUN", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def main() -> int:
    if not SCORES.exists():
        print("waiting for", SCORES)
        return 2
    rc = run("phase1_egfr_her2_metrics.py")
    if rc != 0:
        return rc
    rc = run("phase1_cognate_redock.py")
    # holdout dock is independent
    rc2 = run("phase2_dock_new_holdout_ligands.py")
    if rc2 == 0:
        rc2 = run("phase2_ache_holdout_metrics.py")
    return rc2 or rc


if __name__ == "__main__":
    raise SystemExit(main())
