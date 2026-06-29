#!/usr/bin/env python3
"""STAD-AIDD v1.0 compatibility wrapper — delegates to TAPE-GATE v2.0 pipeline."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

TAPE_GATE = Path(__file__).resolve().parent / "run_tape_gate_pipeline.py"


def main() -> None:
    print("Note: STAD-AIDD v1.0 is superseded by TAPE-GATE v2.0.")
    print("Delegating to run_tape_gate_pipeline.py ...\n")
    cmd = [sys.executable, str(TAPE_GATE)] + sys.argv[1:]
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
