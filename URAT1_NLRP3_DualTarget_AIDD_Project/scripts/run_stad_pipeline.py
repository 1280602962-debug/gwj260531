#!/usr/bin/env python3
"""STAD-AIDD end-to-end pipeline orchestrator."""
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
    parser = argparse.ArgumentParser(description="Run STAD-AIDD pipeline")
    parser.add_argument("--skip-generative", action="store_true")
    parser.add_argument("--library", type=Path, help="Screening library for stage 3")
    args = parser.parse_args()

    run_step("Stage 0: Data preparation", "00_prepare_data.py")
    run_step("Stage 1: Dataset analysis", "01_dataset_analysis.py")
    run_step("Stage 2: MTL training", "02_train_mtl_models.py")

    screen_args = []
    if args.library:
        screen_args = ["--library", str(args.library)]
    run_step("Stage 3: Structure screening", "03_structure_screening.py", screen_args)

    if not args.skip_generative:
        run_step("Stage 4: Generative optimization", "04_generative_optimization.py")

    run_step("Stage 5: Retrospective validation", "05_retrospective_validation.py")
    print("\nPipeline complete. See results/ and docs/MANUSCRIPT_OUTLINE.md")


if __name__ == "__main__":
    main()
