#!/usr/bin/env python3
"""Score existing D1 SDFs without launching GNINA."""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_urat1_capability_dock import main as run_main  # noqa: E402


if __name__ == "__main__":
    sys.argv[1:1] = ["--score-only"]
    run_main()
