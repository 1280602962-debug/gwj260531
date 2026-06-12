#!/usr/bin/env python3
"""
Prepare user-uploaded ChEMBL CSV exports (docs/JNK1/2/3.csv) for the pipeline.

Outputs:
  data/processed/jnk{1,2,3}_curated.csv
  data/processed/paired_set.csv
  data/processed/chemprop_mtl.csv          # merged multitask (missing = blank)
  data/processed/splits/{train,val,test}.csv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "targets.yaml"

VALID_TYPES = {"IC50", "Ki", "Kd", "EC50"}
ISO_MAP = {"JNK1": "docs/JNK1.csv", "JNK2": "docs/JNK2.csv", "JNK3": "docs/JNK3.csv"}


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def canonicalize(smiles: str) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def curate_chembl_export(path: Path, isoform: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df[df["Standard Type"].isin(VALID_TYPES)]
    df = df[df["Standard Relation"].astype(str).str.strip("'\"") == "="]
    df = df[df["Smiles"].notna()]

    # Prefer pChEMBL Value
    df["pActivity"] = pd.to_numeric(df["pChEMBL Value"], errors="coerce")
    missing = df["pActivity"].isna()
    if missing.any():
        sv = pd.to_numeric(df.loc[missing, "Standard Value"], errors="coerce")
        units = df.loc[missing, "Standard Units"].astype(str).str.lower()
        nm_mask = units == "nm"
        df.loc[missing & nm_mask, "pActivity"] = -np.log10(sv[nm_mask] * 1e-9)
        um_mask = units == "um"
        df.loc[missing & um_mask, "pActivity"] = -np.log10(sv[um_mask] * 1e-6)

    df = df[df["pActivity"].notna()]
    df["canonical_smiles"] = df["Smiles"].apply(canonicalize)
    df = df[df["canonical_smiles"].notna()]

    agg = (
        df.groupby(["Molecule ChEMBL ID", "canonical_smiles"], as_index=False)
        .agg(
            pActivity=("pActivity", "mean"),
            n_measurements=("pActivity", "count"),
        )
        .rename(columns={"Molecule ChEMBL ID": "molecule_chembl_id"})
    )
    logger.info("%s: %d raw rows → %d unique compounds", isoform, len(df), len(agg))
    return agg


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


def build_chemprop_mtl(paired: pd.DataFrame) -> pd.DataFrame:
    out = paired[["canonical_smiles", "pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]].copy()
    out = out.rename(columns={"canonical_smiles": "smiles"})
    return out.drop_duplicates(subset=["smiles"])


def scaffold_split_df(df: pd.DataFrame, test_frac=0.1, val_frac=0.1, seed=42):
    from sklearn.model_selection import GroupShuffleSplit

    smiles = df["smiles"].tolist()
    scaffolds = [murcko_scaffold(s) for s in smiles]
    groups = np.array(scaffolds)
    idx = np.arange(len(smiles))

    gss_test = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    train_val_idx, test_idx = next(gss_test.split(idx, groups=groups))

    sub_groups = groups[train_val_idx]
    val_size = val_frac / (1 - test_frac)
    gss_val = GroupShuffleSplit(n_splits=1, test_size=val_size, random_state=seed)
    tr_rel, va_rel = next(gss_val.split(train_val_idx, groups=sub_groups))

    train_idx = train_val_idx[tr_rel]
    val_idx = train_val_idx[va_rel]
    return (
        df.iloc[train_idx].reset_index(drop=True),
        df.iloc[val_idx].reset_index(drop=True),
        df.iloc[test_idx].reset_index(drop=True),
    )


def main():
    parser = argparse.ArgumentParser(description="Prepare uploaded JNK CSV data")
    parser.add_argument("--docs-dir", type=Path, default=ROOT / "docs")
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "processed")
    args = parser.parse_args()
    config = load_config()
    args.output.mkdir(parents=True, exist_ok=True)
    splits_dir = args.output / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)

    datasets = {}
    for iso, rel in ISO_MAP.items():
        path = args.docs_dir / Path(rel).name
        if not path.exists():
            path = ROOT / rel
        datasets[iso] = curate_chembl_export(path, iso)
        out = args.output / f"{iso.lower()}_curated.csv"
        datasets[iso].to_csv(out, index=False)

    paired = build_paired(datasets, config)
    paired.to_csv(args.output / "paired_set.csv", index=False)
    logger.info(
        "Paired set: %d molecules, JNK1-selective: %d",
        len(paired),
        (paired["sel_class"] == "JNK1-selective").sum(),
    )

    chemprop_df = build_chemprop_mtl(paired)
    chemprop_df.to_csv(args.output / "chemprop_mtl.csv", index=False)

    train, val, test = scaffold_split_df(chemprop_df)
    train.to_csv(splits_dir / "train.csv", index=False)
    val.to_csv(splits_dir / "val.csv", index=False)
    test.to_csv(splits_dir / "test.csv", index=False)
    logger.info("Splits: train=%d, val=%d, test=%d", len(train), len(val), len(test))

    summary = {
        "JNK1_compounds": len(datasets["JNK1"]),
        "JNK2_compounds": len(datasets["JNK2"]),
        "JNK3_compounds": len(datasets["JNK3"]),
        "paired_total": len(paired),
        "paired_ge2_isoforms": int((paired["n_isoforms"] >= 2).sum()),
        "jnk1_selective": int((paired["sel_class"] == "JNK1-selective").sum()),
        "train_size": len(train),
        "val_size": len(val),
        "test_size": len(test),
    }
    import json

    with open(args.output / "data_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Summary: %s", summary)


if __name__ == "__main__":
    main()
