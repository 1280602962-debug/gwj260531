#!/usr/bin/env python3
"""
Phase 0: Download and curate JNK1/2/3 bioactivity data from ChEMBL.

Usage:
    python scripts/01_download_chembl_data.py --output data/raw
    python scripts/01_download_chembl_data.py --output data/raw --build-paired
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "targets.yaml"

VALID_TYPES = {"IC50", "Ki", "Kd", "EC50"}


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def fetch_activities(target_chembl_id: str) -> pd.DataFrame:
    """Fetch activities from ChEMBL webresource client."""
    try:
        from chembl_webresource_client.new_client import new_client
    except ImportError as exc:
        raise ImportError(
            "Install chembl-webresource-client: pip install chembl-webresource-client"
        ) from exc

    activity = new_client.activity
    logger.info("Fetching activities for %s ...", target_chembl_id)
    records = activity.filter(target_chembl_id=target_chembl_id).only(
        [
            "molecule_chembl_id",
            "canonical_smiles",
            "standard_type",
            "standard_value",
            "standard_units",
            "pchembl_value",
            "assay_type",
            "assay_description",
            "assay_confidence_score",
            "bao_label",
            "target_chembl_id",
        ]
    )
    df = pd.DataFrame(list(records))
    logger.info("  Raw records: %d", len(df))
    return df


def curate(df: pd.DataFrame, min_confidence: int = 6) -> pd.DataFrame:
    """Apply standard curation filters."""
    out = df.copy()

    # Activity type filter
    out = out[out["standard_type"].isin(VALID_TYPES)]

    # Assay confidence
    out = out[out["assay_confidence_score"].notna()]
    out = out[out["assay_confidence_score"] >= min_confidence]

    # Valid pchembl / standard value
    out = out[out["pchembl_value"].notna() | out["standard_value"].notna()]

    # Compute pActivity
    if "pchembl_value" in out.columns:
        out["pActivity"] = out["pchembl_value"]
    else:
        import numpy as np

        out["pActivity"] = out.apply(
            lambda r: -1 * np.log10(r["standard_value"] * 1e-9)
            if r["standard_units"] == "nM"
            else None,
            axis=1,
        )

    # Valid SMILES
    try:
        from rdkit import Chem

        out = out[out["canonical_smiles"].notna()]
        out = out[
            out["canonical_smiles"].apply(lambda s: Chem.MolFromSmiles(s) is not None)
        ]
    except ImportError:
        logger.warning("RDKit not installed; skipping SMILES validation")

    # Aggregate duplicates: geometric mean of standard_value per molecule
    agg = (
        out.groupby(["molecule_chembl_id", "canonical_smiles"], as_index=False)
        .agg(
            pActivity=("pActivity", "mean"),
            n_measurements=("pActivity", "count"),
            assay_types=("assay_type", lambda x: "|".join(sorted(set(x.dropna())))),
        )
    )
    logger.info("  Curated unique compounds: %d", len(agg))
    return agg


def build_paired(datasets: dict[str, pd.DataFrame], config: dict) -> pd.DataFrame:
    """Build paired dataset for molecules measured on multiple isoforms."""
    threshold_active = config["selectivity"]["jnk1_active_threshold_pchembl"]
    threshold_inactive = config["selectivity"]["jnk23_inactive_threshold_pchembl"]
    delta_threshold = config["selectivity"]["delta_log_threshold"]

    merged = None
    for isoform, df in datasets.items():
        sub = df[["molecule_chembl_id", "canonical_smiles", "pActivity"]].rename(
            columns={"pActivity": f"pAct_{isoform}"}
        )
        merged = sub if merged is None else merged.merge(sub, on=["molecule_chembl_id", "canonical_smiles"], how="outer")

    merged["n_isoforms"] = merged[["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]].notna().sum(axis=1)
    paired = merged[merged["n_isoforms"] >= 2].copy()

    paired["delta_12"] = paired["pAct_JNK1"] - paired["pAct_JNK2"]
    paired["delta_13"] = paired["pAct_JNK1"] - paired["pAct_JNK3"]
    paired["delta_min"] = paired["pAct_JNK1"] - paired[["pAct_JNK2", "pAct_JNK3"]].max(axis=1)

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
        if not pd.isna(row["delta_12"]) and row["delta_12"] < -delta_threshold:
            return "JNK2-biased"
        return "non-selective"

    paired["sel_class"] = paired.apply(classify, axis=1)
    logger.info("Paired set: %d molecules (≥2 isoforms)", len(paired))
    logger.info("  JNK1-selective: %d", (paired["sel_class"] == "JNK1-selective").sum())
    return paired


def main():
    parser = argparse.ArgumentParser(description="Download JNK1/2/3 ChEMBL data")
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--build-paired", action="store_true")
    parser.add_argument("--skip-download", action="store_true", help="Use existing raw CSVs")
    args = parser.parse_args()

    config = load_config()
    args.output.mkdir(parents=True, exist_ok=True)
    processed_dir = ROOT / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    curated = {}
    for isoform, meta in config["targets"].items():
        chembl_id = meta["chembl_id"]
        out_raw = args.output / f"{isoform.lower()}_raw.csv"
        out_cur = processed_dir / f"{isoform.lower()}_curated.csv"

        if args.skip_download and out_raw.exists():
            df_raw = pd.read_csv(out_raw)
        else:
            df_raw = fetch_activities(chembl_id)
            df_raw.to_csv(out_raw, index=False)

        df_cur = curate(df_raw, min_confidence=config["min_assay_confidence"])
        df_cur.to_csv(out_cur, index=False)
        curated[isoform] = df_cur
        logger.info("Saved %s", out_cur)

    if args.build_paired:
        paired = build_paired(curated, config)
        paired_path = processed_dir / "paired_set.csv"
        paired.to_csv(paired_path, index=False)
        logger.info("Saved paired set: %s", paired_path)


if __name__ == "__main__":
    main()
