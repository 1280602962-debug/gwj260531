#!/usr/bin/env python3
"""
TAPE-GATE Stage 5: Merge Path A/B candidates + reliability-weighted Pareto fusion.

NOT fixed 0.5/0.5 linear fusion (see DIFFERENTIATION_VS_PLK1_NLRP3.md).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-a", type=Path, default=PROJECT_ROOT / "results" / "path_a_library")
    parser.add_argument("--path-b", type=Path, default=PROJECT_ROOT / "results" / "generation")
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "model_hierarchy.yaml")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "fusion")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        model_cfg = yaml.safe_load(f)

    fusion = model_cfg.get("fusion", {})
    report = {
        "status": "skeleton",
        "framework": "TAPE-GATE",
        "inputs": {"path_a": str(args.path_a), "path_b": str(args.path_b)},
        "fusion_method": fusion.get("method", "reliability_weighted_pareto"),
        "ml_weights": fusion.get("ml_weights"),
        "pareto_objectives": fusion.get("pareto_objectives"),
        "formula": {
            "omega_U": "proportional to 1/conformal_width",
            "omega_N": "proportional to assay_confidence",
            "S_dual": "omega_U * y_U + omega_N * P_active + gamma * sqrt(S_U * S_N)",
        },
        "ablation_fixed_fusion": fusion.get("ablation_fixed_fusion", "equal_weight_0.5"),
    }
    out = args.output / "fusion_ranking_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Fusion ranking config snapshot: {out}")


if __name__ == "__main__":
    main()
