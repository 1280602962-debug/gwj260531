#!/usr/bin/env python3
"""
STAD-AIDD Stage 3: Conformational ensemble docking + dual-target funnel.

URAT1: S_trap conformation-trapping score (transporter-aware)
NLRP3: NACHT domain docking + optional MM-GBSA

See config/docking_ensemble.yaml and docs/URAT1_TRANSPORTER_VALIDATION.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=False, help="SMILES library CSV")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "screening")
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    dock_cfg_path = PROJECT_ROOT / "config" / "docking_ensemble.yaml"
    with open(dock_cfg_path) as f:
        dock_cfg = yaml.safe_load(f)

    report = {
        "status": "skeleton",
        "urat1_ensemble": [s["pdb_id"] for s in dock_cfg["urat1_ensemble"]["pdb_structures"]],
        "nlrp3_ensemble": [s["pdb_id"] for s in dock_cfg["nlrp3_ensemble"]["pdb_structures"]],
        "funnel": dock_cfg["dual_target_funnel"],
        "scoring_formula": "S_dual = sqrt(S_U * S_N) + 0.2 * min(S_U, S_N)",
        "transporter_note": "URAT1 requires S_trap; do NOT use single-structure enzyme-style docking",
    }
    out = args.output / "screening_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Docking config snapshot: {out}")


if __name__ == "__main__":
    main()
