#!/usr/bin/env python3
"""
STAD-AIDD Stage 0: Download and curate URAT1 / NLRP3 bioactivity data.

Sources:
  - ChEMBL API (CHEMBL6120, CHEMBL1741208)
  - Patent supplements (manual CSV in data/raw/)
  - SLC22 auxiliary targets for transfer learning

Output:
  data/processed/urat1_curated.csv
  data/processed/nlrp3_curated.csv
  data/processed/splits/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "targets.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare URAT1/NLRP3 datasets")
    parser.add_argument("--skip-chembl", action="store_true", help="Skip API download")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "data" / "processed")
    args = parser.parse_args()

    config = load_config()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "status": "skeleton",
        "message": (
            "Implement ChEMBL download using chembl_webresource_client. "
            "See docs/PREPARATION_CHECKLIST.md and config/targets.yaml for curation rules."
        ),
        "targets": {
            "urat1": config["targets"]["urat1"]["chembl_id"],
            "nlrp3": config["targets"]["nlrp3"]["chembl_id"],
        },
        "curation": config["data_curation"],
    }

    summary_path = args.output_dir / "data_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Config loaded. Summary written to {summary_path}")
    print("Next: implement ChEMBL fetch + curate_isoform_raw (see JNK1 scripts/utils_ml.py)")


if __name__ == "__main__":
    main()
