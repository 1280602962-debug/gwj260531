#!/usr/bin/env python3
"""LEGACY / PRE-REMEDIATION pointer.

Official analysis environment: scripts/check_analysis_env.py
Official docking environment: scripts/check_docking_env.py

This file no longer treats /home/gwj or /mnt/d/... as official lookup logic.
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    print("LEGACY check_local_env.py — delegating to scripts/check_analysis_env.py")
    print("Docking binaries: scripts/check_docking_env.py (PATH or VINA_BIN/GNINA_BIN/RTMSCORE_PYTHON)")
    sys.argv = [str(ROOT / "scripts" / "check_analysis_env.py"), *sys.argv[1:]]
    try:
        runpy.run_path(str(ROOT / "scripts" / "check_analysis_env.py"), run_name="__main__")
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
