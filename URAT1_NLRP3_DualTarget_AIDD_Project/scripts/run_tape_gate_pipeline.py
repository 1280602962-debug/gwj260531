#!/usr/bin/env python3
"""TAPE-GATE end-to-end pipeline: dual-path library + generative screening."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"


def run_step(name: str, script: str, extra_args: list[str] | None = None) -> None:
    cmd = [sys.executable, str(SCRIPTS / script)] + (extra_args or [])
    print(f"\n{'='*60}\n>>> {name}\n{'='*60}")
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run TAPE-GATE dual-path pipeline")
    parser.add_argument("--skip-generative", action="store_true", help="Path A only (library screening)")
    parser.add_argument("--library", type=Path, help="Screening library for Path A")
    parser.add_argument("--run-plk1-baseline", action="store_true", help="Also run PLK1-style ablation baseline")
    args = parser.parse_args()

    run_step("Stage 0: Data preparation", "00_prepare_data.py")
    run_step("Stage 1: Dataset analysis", "01_dataset_analysis.py")
    run_step("Stage 2: Asymmetric dual-evidence models", "02_train_asymmetric_models.py")

    screen_args: list[str] = []
    if args.library:
        screen_args = ["--library", str(args.library)]
    run_step("Stage 3: Path A — Library screening", "03_library_screening.py", screen_args)

    if not args.skip_generative:
        run_step("Stage 4: Path B — Generative optimization", "04_generative_optimization.py")

    run_step("Stage 5: Reliability fusion + Pareto ranking", "05_fusion_and_ranking.py")

    val_args = []
    if args.run_plk1_baseline:
        val_args = ["--run-plk1-baseline"]
    run_step("Stage 6: Retrospective validation", "06_retrospective_validation.py", val_args)

    print("\nTAPE-GATE pipeline complete.")
    print("See results/ and docs/MANUSCRIPT_OUTLINE.md")
    print("Differentiation vs PLK1/NLRP3: docs/DIFFERENTIATION_VS_PLK1_NLRP3.md")


if __name__ == "__main__":
    main()
