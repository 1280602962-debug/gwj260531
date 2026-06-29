#!/usr/bin/env python3
"""
TAPE-GATE Stage 2: Asymmetric dual-evidence modeling.

URAT1 arm: regression + split conformal UQ + SLC22 transfer
NLRP3 arm: assay-conditioned classification (NOT anchor similarity)

See docs/ALGORITHM_FRAMEWORK.md §3 and config/model_hierarchy.yaml
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "model_hierarchy.yaml")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "training")
    parser.add_argument("--run-mtl-ablation", action="store_true", help="Also train MTL for Abl-6")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        model_cfg = yaml.safe_load(f)

    report = {
        "status": "skeleton",
        "framework": "TAPE-GATE",
        "strategy": model_cfg.get("modeling_strategy", "independent_dual_evidence"),
        "urat1_arm": model_cfg.get("urat1_arm"),
        "nlrp3_arm": model_cfg.get("nlrp3_arm"),
        "fusion": model_cfg.get("fusion"),
        "cv_protocol": model_cfg.get("cv_protocol"),
        "baselines": model_cfg.get("baselines"),
        "implementation_notes": [
            "URAT1: MiniMol/Chemprop regression + split conformal (alpha=0.1)",
            "NLRP3: assay-conditioned classification; preserve assay_id metadata",
            "Train independent models by default (0 overlapping SMILES)",
            "Run PLK1-style baseline (SVR + anchor similarity) for ablation only",
            "Murcko GroupKFold (n=5); Wilcoxon vs XGBoost and PLK1-style",
        ],
        "avoid_as_primary": model_cfg.get("nlrp3_arm", {}).get("deprecated_primary", []),
    }
    if args.run_mtl_ablation:
        report["mtl_ablation"] = model_cfg.get("mtl_ablation")

    out = args.output / "asymmetric_training_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Asymmetric model config snapshot: {out}")


if __name__ == "__main__":
    main()
