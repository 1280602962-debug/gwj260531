#!/usr/bin/env python3
"""Deprecated wrapper.

The old C5 W1 runner stopped after benzbromarone Top-1 RMSD > 2 Å and included
9DK9 apo. Local capability docking now uses scripts/run_urat1_capability_dock.py
(D1: 9DKA/B/C only, 20 new jobs, no gate).
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_urat1_capability_dock import main as capability_main  # noqa: E402


def main() -> None:
    print(
        "WARNING: run_c5_w1_crossdock.py now forwards to "
        "run_urat1_capability_dock.py (no Top-1 gate, no 9DK9).",
        flush=True,
    )
    # Drop the old --phase flag if present; D1 has no gate/rest split.
    cleaned: list[str] = []
    args = sys.argv[1:]
    skip_next = False
    for i, item in enumerate(args):
        if skip_next:
            skip_next = False
            continue
        if item == "--phase":
            skip_next = True
            continue
        if item.startswith("--phase="):
            continue
        cleaned.append(item)
    sys.argv = [sys.argv[0], *cleaned]
    capability_main()


if __name__ == "__main__":
    main()
