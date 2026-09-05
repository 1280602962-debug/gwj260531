#!/usr/bin/env python3
"""Curate URAT1 / NLRP3 bioactivity data from ChEMBL exports.

Input (default): data/raw/URAT1_CHEMBL.csv, data/raw/NLRP3_CHEMBL.csv
Output:
  data/processed/urat1_curated.csv
  data/processed/nlrp3_records.csv
  data/processed/data_summary.json
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from utils_ml import canonicalize, curate_nlrp3_records, curate_urat1_raw

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "targets.yaml"
DEFAULT_UPLOADS = Path("/home/ubuntu/.cursor/projects/workspace/uploads")


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def resolve_input(path: Path | None, fallback_name: str) -> Path:
    if path and path.exists():
        return path
    local = PROJECT_ROOT / "data" / "raw" / fallback_name
    if local.exists():
        return local
    upload = DEFAULT_UPLOADS / fallback_name
    if upload.exists():
        return upload
    raise FileNotFoundError(f"Could not find {fallback_name}; pass --urat1-csv / --nlrp3-csv")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare URAT1/NLRP3 datasets")
    parser.add_argument("--urat1-csv", type=Path, default=None)
    parser.add_argument("--nlrp3-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "data" / "processed")
    parser.add_argument("--copy-raw", action="store_true", help="Copy source CSVs into data/raw/")
    args = parser.parse_args()

    config = load_config()
    curation = config["data_curation"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = PROJECT_ROOT / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    urat1_src = resolve_input(args.urat1_csv, "URAT1_CHEMBL_cf12.csv")
    nlrp3_src = resolve_input(args.nlrp3_csv, "NLRP3_CHEMBL_4807.csv")

    if args.copy_raw:
        shutil.copy2(urat1_src, raw_dir / urat1_src.name)
        shutil.copy2(nlrp3_src, raw_dir / nlrp3_src.name)

    urat1 = curate_urat1_raw(
        urat1_src,
        pactivity_min=curation["pactivity_range"][0],
        pactivity_max=curation["pactivity_range"][1],
        max_std=curation["conflict_std_threshold"],
        max_range=curation["conflict_range_threshold"],
    )
    nlrp3 = curate_nlrp3_records(
        nlrp3_src,
        pactivity_min=curation["pactivity_range"][0],
        pactivity_max=curation["pactivity_range"][1],
        min_assay_compounds=5,
        active_threshold=6.0,
    )

    urat1_path = args.output_dir / "urat1_curated.csv"
    nlrp3_path = args.output_dir / "nlrp3_records.csv"
    urat1.to_csv(urat1_path, index=False)
    nlrp3.to_csv(nlrp3_path, index=False)

    overlap = set(urat1["canonical_smiles"]) & set(nlrp3["canonical_smiles"].unique())
    nlrp3_by_smiles = nlrp3.groupby("canonical_smiles")["pActivity"].agg(["count", "min", "max"])
    nlrp3_multi = nlrp3_by_smiles[nlrp3_by_smiles["count"] > 1]
    nlrp3_conflict = nlrp3_multi[(nlrp3_multi["max"] - nlrp3_multi["min"]) > 1.0]
    n_compounds = int(nlrp3["canonical_smiles"].nunique())
    thp_mask = nlrp3["assay_cell_type"].astype(str).str.contains("THP", case=False, na=False)
    thp_df = nlrp3[thp_mask]
    summary = {
        "urat1": {
            "n_compounds": int(len(urat1)),
            "pactivity_mean": float(urat1["pActivity"].mean()),
            "pactivity_std": float(urat1["pActivity"].std()),
            "n_scaffolds": int(urat1["scaffold"].nunique()),
            "source": str(urat1_src),
        },
        "nlrp3": {
            "n_records": int(len(nlrp3)),
            "n_compounds": n_compounds,
            "n_assays": int(nlrp3["Assay ChEMBL ID"].nunique()),
            "n_multi_assay_compounds": int(len(nlrp3_multi)),
            "n_conflict_gt_1log_compounds": int(len(nlrp3_conflict)),
            "pct_conflict_gt_1log": round(100.0 * len(nlrp3_conflict) / n_compounds, 1),
            "active_rate": float(nlrp3["active"].mean()),
            "thp1_n_records": int(len(thp_df)),
            "thp1_n_compounds": int(thp_df["canonical_smiles"].nunique()),
            "source": str(nlrp3_src),
        },
        "overlap_smiles": int(len(overlap)),
        "curation": {
            "nlrp3_filter": "IL-1beta endpoint + Assay Type B + assay>=5 compounds",
            "urat1_aggregation": "median per SMILES with conflict discard",
        },
    }
    summary_path = args.output_dir / "data_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"URAT1 curated: {len(urat1)} compounds -> {urat1_path}")
    print(f"NLRP3 records: {len(nlrp3)} ({nlrp3['canonical_smiles'].nunique()} compounds) -> {nlrp3_path}")
    print(f"Overlap SMILES: {len(overlap)}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
