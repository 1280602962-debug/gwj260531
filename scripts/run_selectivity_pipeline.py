#!/usr/bin/env python3
"""
Run the JNK1 selectivity pipeline: data prep → train → SHAP → screening.

Steps:
  00  Prepare curated datasets + paired labels
  04  Train MTL + selectivity models (v2 features)
  05  SHAP interpretation
  06  Virtual screening on demo library

Usage:
    python scripts/run_selectivity_pipeline.py
    python scripts/run_selectivity_pipeline.py --skip-data-prep
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def run(cmd: list[str], desc: str) -> None:
    logger.info("=== %s ===", desc)
    logger.info("Command: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        raise SystemExit(f"Step failed ({desc}): exit code {result.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run JNK selectivity pipeline (04→05→06)")
    parser.add_argument("--skip-data-prep", action="store_true", help="Skip 00_prepare_user_data.py")
    parser.add_argument("--skip-shap", action="store_true", help="Skip SHAP step")
    parser.add_argument("--library", type=Path, default=ROOT / "data" / "libraries" / "screening_demo.smi")
    args = parser.parse_args()

    py = sys.executable

    if not args.skip_data_prep:
        run([py, str(SCRIPTS / "00_prepare_user_data.py")], "Data preparation (00)")

    run([py, str(SCRIPTS / "build_demo_library.py")], "Build demo screening library")

    run(
        [
            py,
            str(SCRIPTS / "04_train_selectivity_model.py"),
            "--input",
            str(ROOT / "data" / "processed"),
            "--output",
            str(ROOT / "models"),
            "--plots",
            str(ROOT / "results" / "training"),
        ],
        "Selectivity model training (04)",
    )

    if not args.skip_shap:
        run(
            [
                py,
                str(SCRIPTS / "05_model_interpretation.py"),
                "--model",
                str(ROOT / "models" / "best_model.joblib"),
                "--data",
                str(ROOT / "data" / "processed" / "paired_set.csv"),
                "--output",
                str(ROOT / "results" / "shap"),
            ],
            "SHAP interpretation (05)",
        )

    run(
        [
            py,
            str(SCRIPTS / "06_virtual_screening.py"),
            "--models-dir",
            str(ROOT / "models" / "xgboost"),
            "--library",
            str(args.library),
            "--output",
            str(ROOT / "results" / "screening_v2"),
        ],
        "Virtual screening (06)",
    )

    logger.info("Selectivity pipeline complete.")


if __name__ == "__main__":
    main()
