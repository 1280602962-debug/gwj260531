#!/usr/bin/env python3
"""
Prepare user-uploaded ChEMBL CSV exports with improved curation for QSAR.

Key improvements vs v1:
  - Biochemical assays only (Assay Type B)
  - IC50 with exact relation (=)
  - pActivity range filter [4, 10]
  - Remove conflicting multi-assay measurements
  - Keep assays with >= N compounds (assay harmonization)
  - Per-isoform single-target datasets + splits (better R²)
  - Multitask merged table retained for selectivity analysis

Usage:
    python3 scripts/00_prepare_user_data.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from utils_ml import curate_isoform_raw, murcko_scaffold, scaffold_holdout_split  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CONFIG_PATH = ROOT / "config" / "targets.yaml"
ISO_MAP = {"JNK1": "JNK1.csv", "JNK2": "JNK2.csv", "JNK3": "JNK3.csv"}


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def build_paired(datasets: dict[str, pd.DataFrame], config: dict) -> pd.DataFrame:
    threshold_active = config["selectivity"]["jnk1_active_threshold_pchembl"]
    threshold_inactive = config["selectivity"]["jnk23_inactive_threshold_pchembl"]
    delta_threshold = config["selectivity"]["delta_log_threshold"]

    merged = None
    for iso, df in datasets.items():
        sub = df[["molecule_chembl_id", "canonical_smiles", "pActivity"]].rename(
            columns={"pActivity": f"pAct_{iso}"}
        )
        merged = sub if merged is None else merged.merge(
            sub, on=["molecule_chembl_id", "canonical_smiles"], how="outer"
        )

    merged["n_isoforms"] = merged[["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]].notna().sum(axis=1)
    merged["delta_12"] = merged["pAct_JNK1"] - merged["pAct_JNK2"]
    merged["delta_13"] = merged["pAct_JNK1"] - merged["pAct_JNK3"]
    merged["delta_min"] = merged["pAct_JNK1"] - merged[["pAct_JNK2", "pAct_JNK3"]].max(axis=1)

    def classify(row):
        if pd.isna(row["pAct_JNK1"]):
            return "unknown"
        jnk1_active = row["pAct_JNK1"] >= threshold_active
        jnk2_inactive = pd.isna(row["pAct_JNK2"]) or row["pAct_JNK2"] < threshold_inactive
        jnk3_inactive = pd.isna(row["pAct_JNK3"]) or row["pAct_JNK3"] < threshold_inactive
        delta_ok = (
            (not pd.isna(row["delta_min"]) and row["delta_min"] >= delta_threshold)
            or (not pd.isna(row["delta_12"]) and row["delta_12"] >= delta_threshold)
            or (not pd.isna(row["delta_13"]) and row["delta_13"] >= delta_threshold)
        )
        if jnk1_active and jnk2_inactive and jnk3_inactive and delta_ok:
            return "JNK1-selective"
        if jnk1_active and not jnk2_inactive and not jnk3_inactive:
            return "pan-JNK"
        return "non-selective"

    merged["sel_class"] = merged.apply(classify, axis=1)
    return merged


def save_isoform_splits(df: pd.DataFrame, isoform: str, splits_dir: Path, target_col: str = "pActivity"):
    iso_dir = splits_dir / isoform.lower()
    iso_dir.mkdir(parents=True, exist_ok=True)
    chem = df.rename(columns={"canonical_smiles": "smiles", target_col: f"pAct_{isoform}"})
    chem = chem[["smiles", f"pAct_{isoform}"]]

    train_idx, val_idx, test_idx = scaffold_holdout_split(chem["smiles"].tolist())
    chem.iloc[train_idx].to_csv(iso_dir / "train.csv", index=False)
    chem.iloc[val_idx].to_csv(iso_dir / "val.csv", index=False)
    chem.iloc[test_idx].to_csv(iso_dir / "test.csv", index=False)
    chem.to_csv(iso_dir / "full.csv", index=False)
    return len(train_idx), len(val_idx), len(test_idx)


def main():
    parser = argparse.ArgumentParser(description="Prepare uploaded JNK CSV data (v2 curation)")
    parser.add_argument("--docs-dir", type=Path, default=ROOT / "docs")
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "processed")
    parser.add_argument("--min-assay-compounds", type=int, default=10)
    args = parser.parse_args()
    config = load_config()
    args.output.mkdir(parents=True, exist_ok=True)
    splits_dir = args.output / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)

    curation = config.get("curation", {})
    per_iso = config.get("curation_per_isoform", {})

    datasets = {}
    split_info = {}
    for iso, fname in ISO_MAP.items():
        path = args.docs_dir / fname
        iso_cur = {**curation, **per_iso.get(iso, {})}
        if "min_assay_compounds" not in iso_cur:
            iso_cur["min_assay_compounds"] = args.min_assay_compounds
        df = curate_isoform_raw(path, **iso_cur)
        datasets[iso] = df
        out = args.output / f"{iso.lower()}_curated.csv"
        df.to_csv(out, index=False)
        tr, va, te = save_isoform_splits(df, iso, splits_dir)
        split_info[iso] = {"train": tr, "val": va, "test": te, "total": len(df)}
        logger.info("%s curated: %d compounds (train/val/test = %d/%d/%d)", iso, len(df), tr, va, te)

    paired = build_paired(datasets, config)
    paired.to_csv(args.output / "paired_set.csv", index=False)

    # Legacy multitask merged file
    mtl = paired[["canonical_smiles", "pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]].rename(columns={"canonical_smiles": "smiles"})
    mtl = mtl.drop_duplicates(subset=["smiles"])
    mtl.to_csv(args.output / "chemprop_mtl.csv", index=False)

    summary = {
        "curation": curation,
        "curation_per_isoform": per_iso,
        "JNK1_compounds": len(datasets["JNK1"]),
        "JNK2_compounds": len(datasets["JNK2"]),
        "JNK3_compounds": len(datasets["JNK3"]),
        "paired_total": len(paired),
        "paired_ge2_isoforms": int((paired["n_isoforms"] >= 2).sum()),
        "jnk1_selective": int((paired["sel_class"] == "JNK1-selective").sum()),
        "splits": split_info,
    }
    with open(args.output / "data_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Data preparation complete.")


if __name__ == "__main__":
    main()
