#!/usr/bin/env python3
"""
STAD-AIDD Stage 4: RL-guided dual-target molecular generation.

Reward: MTL predictions + ensemble docking + QED + SA + novelty
Reference: POLYGON (Nat Commun 2024), CLM dual-target (Nat Commun 2024)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=5000)
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "generation")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    report = {
        "status": "skeleton",
        "method": "CLM + REINFORCE",
        "reward_components": [
            "predicted_pActivity_URAT1",
            "predicted_pActivity_NLRP3",
            "S_trap_URAT1",
            "S_struct_NLRP3",
            "QED",
            "SA_score",
            "novelty",
        ],
        "compute_tip": "Compute docking reward every 500 RL steps to save CPU",
    }
    out = args.output / "generation_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Generation config snapshot: {out}")


if __name__ == "__main__":
    main()
