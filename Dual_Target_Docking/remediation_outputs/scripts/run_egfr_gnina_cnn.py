#!/usr/bin/env python3
"""GNINA CNN rescore of corrected-box Vina poses (best-of-K)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SCRIPT = ROOT / "data/jcim_bench_v0/scripts/gnina_rescore_panel.py"
POSE_ROOT = ROOT / "remediation_outputs/phase_rtm"
REC = ROOT / "data/egfr_her2_panel40_v0/receptors"


def main() -> int:
    (POSE_ROOT / "tables").mkdir(parents=True, exist_ok=True)
    rec_a = REC / "3POZ_protein.pdb"
    rec_b = REC / "3RCD_protein.pdb"
    if not rec_a.exists():
        rec_a = REC / "3POZ_receptor.pdbqt"
    if not rec_b.exists():
        rec_b = REC / "3RCD_receptor.pdbqt"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--root",
        str(POSE_ROOT),
        "--targets",
        "3POZ",
        "3RCD",
        "--receptor-map",
        f"3POZ={rec_a}",
        f"3RCD={rec_b}",
        "--workers",
        "2",
        "--modes",
        "all",
        "--timeout",
        "180",
    ]
    print(" ".join(cmd), flush=True)
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
