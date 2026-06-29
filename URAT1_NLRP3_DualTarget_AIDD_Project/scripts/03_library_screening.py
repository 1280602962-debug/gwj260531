#!/usr/bin/env python3
"""
TAPE-GATE Stage 3 / Path A: Library screening funnel.

Flow: ADMET → URAT1 conformal filter → NLRP3 assay-conditioned filter
      → ensemble docking (S_trap + NLRP3 struct) → diversity selection

See config/dual_path.yaml path_a_library_screening
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=False, help="SMILES library CSV/SMI")
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "dual_path.yaml")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "path_a_library")
    parser.add_argument("--top-n", type=int, default=500)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        path_cfg = yaml.safe_load(f)
    with open(PROJECT_ROOT / "config" / "docking_ensemble.yaml") as f:
        dock_cfg = yaml.safe_load(f)

    path_a = path_cfg.get("path_a_library_screening", {})
    report = {
        "status": "skeleton",
        "path": "A_library_screening",
        "library": str(args.library) if args.library else path_a.get("library", {}).get("primary"),
        "ml_filters": path_a.get("ml_filters"),
        "urat1_scoring": "S_trap conformation-trapping (transporter-aware)",
        "nlrp3_scoring": "assay-conditioned P_active + NACHT docking",
        "urat1_ensemble": [s["pdb_id"] for s in dock_cfg["urat1_ensemble"]["pdb_structures"]],
        "nlrp3_ensemble": [s["pdb_id"] for s in dock_cfg["nlrp3_ensemble"]["pdb_structures"]],
        "output_max": path_a.get("output_max", args.top_n),
        "not_used": "anchor ECFP max-pooling (reserved for PLK1-style baseline)",
    }
    out = args.output / "path_a_screening_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Path A screening config snapshot: {out}")


if __name__ == "__main__":
    main()
