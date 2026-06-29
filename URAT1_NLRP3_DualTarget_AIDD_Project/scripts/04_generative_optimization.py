#!/usr/bin/env python3
"""
TAPE-GATE Stage 4 / Path B: RL-guided dual-target molecular generation.

Reward: URAT1 regression + NLRP3 assay-conditioned prob + S_trap + NLRP3 struct
        + QED + SA + novelty + conformal penalty

Reference: POLYGON (Nat Commun 2024), CLM dual-target (Nat Commun 2024)
See config/dual_path.yaml path_b_generative
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
        "path": "B_generative",
        "method": "CLM cross-fine-tune + REINFORCE",
        "reward_components": [
            "urat1_predicted_pactivity_conformal",
            "nlrp3_assay_conditioned_P_active",
            "S_trap_URAT1",
            "S_struct_NLRP3",
            "QED",
            "SA_score",
            "novelty",
            "conformal_width_penalty",
        ],
        "compute_tip": "Compute docking reward every 100 RL steps (see dual_path.yaml)",
        "differentiation": "PLK1/NLRP3 paper has no generative path — this is TAPE-GATE innovation",
    }
    out = args.output / "generation_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Generation config snapshot: {out}")


if __name__ == "__main__":
    main()
