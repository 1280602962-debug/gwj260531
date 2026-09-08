#!/usr/bin/env python3
"""Build URAT1/NLRP3 models, assess quality, and run benchmark backtest."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
PROJECT = SCRIPTS.parent


def run(script: str) -> None:
    print(f"\n>>> {script}")
    subprocess.run([sys.executable, str(SCRIPTS / script)], check=True, cwd=SCRIPTS)


def main() -> None:
    run("00_prepare_data.py")
    run("02_train_asymmetric_models.py")
    run("07_benchmark_backtest.py")
    print("\nDone. See results/training/ and results/benchmark_backtest/")
    print("Summary docs: docs/MODEL_QUALITY_REPORT.md")


if __name__ == "__main__":
    main()
