#!/usr/bin/env python3
"""
TAPE-GATE ablation only: Multi-task learning (Abl-6).

Primary training path is 02_train_asymmetric_models.py (independent dual evidence).
MTL is optional because ChEMBL data has 0 overlapping URAT1/NLRP3 SMILES.

Architecture (see docs/ALGORITHM_FRAMEWORK.md §3.3):
  MiniMol fingerprint (frozen) + MLP heads (urat1, nlrp3, dual)
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
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        model_cfg = yaml.safe_load(f)

    report = {
        "status": "skeleton",
        "mtl_architecture": model_cfg.get("mtl_architecture"),
        "cv_protocol": model_cfg.get("cv_protocol"),
        "baselines": model_cfg.get("baselines"),
        "implementation_notes": [
            "Use Murcko GroupKFold (n=5)",
            "Report RMSE, MAE, R2, Spearman, EF@1%",
            "Wilcoxon test vs XGBoost baseline",
        ],
    }
    out = args.output / "training_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Training config snapshot: {out}")


if __name__ == "__main__":
    main()
