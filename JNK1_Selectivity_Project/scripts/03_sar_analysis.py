#!/usr/bin/env python3
"""
Phase 2 (Method 1): Classical SAR / MMP / activity cliff analysis.

Usage:
    python scripts/03_sar_analysis.py --input data/processed --output results/sar
"""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]


def morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def tanimoto_smiles(s1: str, s2: str) -> float:
    from rdkit import DataStructs

    fp1, fp2 = morgan_fp(s1), morgan_fp(s2)
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def find_activity_cliffs(df: pd.DataFrame, tanimoto_min: float = 0.85, delta_min: float = 1.0) -> pd.DataFrame:
    """Find activity cliffs within a dataset."""
    smiles = df["canonical_smiles"].dropna().unique()
    pactivity = df.set_index("canonical_smiles")["pActivity"].to_dict()
    cliffs = []
    n = len(smiles)
    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = smiles[i], smiles[j]
            sim = tanimoto_smiles(s1, s2)
            if sim >= tanimoto_min:
                delta = abs(pactivity.get(s1, 0) - pactivity.get(s2, 0))
                if delta >= delta_min:
                    cliffs.append(
                        {
                            "smiles_1": s1,
                            "smiles_2": s2,
                            "tanimoto": sim,
                            "delta_pActivity": delta,
                            "pAct_1": pactivity.get(s1),
                            "pAct_2": pactivity.get(s2),
                        }
                    )
    return pd.DataFrame(cliffs)


def find_selectivity_cliffs(paired: pd.DataFrame, tanimoto_min: float = 0.85) -> pd.DataFrame:
    """Find pairs with high similarity but large selectivity difference."""
    if "delta_min" not in paired.columns:
        return pd.DataFrame()
    smiles = paired["canonical_smiles"].dropna().tolist()
    delta_map = paired.set_index("canonical_smiles")["delta_min"].to_dict()
    cliffs = []
    for i in range(len(smiles)):
        for j in range(i + 1, len(smiles)):
            s1, s2 = smiles[i], smiles[j]
            sim = tanimoto_smiles(s1, s2)
            if sim >= tanimoto_min:
                d1, d2 = delta_map.get(s1), delta_map.get(s2)
                if d1 is not None and d2 is not None and not (np.isnan(d1) or np.isnan(d2)):
                    if abs(d1 - d2) >= 1.0:
                        cliffs.append(
                            {
                                "smiles_1": s1,
                                "smiles_2": s2,
                                "tanimoto": sim,
                                "delta_min_1": d1,
                                "delta_min_2": d2,
                                "sel_cliff": abs(d1 - d2),
                            }
                        )
    return pd.DataFrame(cliffs)


def simple_mmp_analysis(paired: pd.DataFrame) -> pd.DataFrame:
    """
    Simplified MMP: find pairs differing by <= 1 heavy atom change
    (exact atom-level MMP requires mmpdb; this is a lightweight proxy).
    """
    from rdkit import Chem

    def num_heavy(s):
        m = Chem.MolFromSmiles(s)
        return m.GetNumHeavyAtoms() if m else 0

    records = []
    smiles_list = paired["canonical_smiles"].dropna().tolist()
    for i in range(len(smiles_list)):
        for j in range(i + 1, len(smiles_list)):
            s1, s2 = smiles_list[i], smiles_list[j]
            if abs(num_heavy(s1) - num_heavy(s2)) > 2:
                continue
            sim = tanimoto_smiles(s1, s2)
            if 0.7 <= sim < 1.0:
                row1 = paired[paired["canonical_smiles"] == s1].iloc[0]
                row2 = paired[paired["canonical_smiles"] == s2].iloc[0]
                if "delta_12" in paired.columns:
                    d12 = row1.get("delta_12", np.nan) - row2.get("delta_12", np.nan)
                    if not np.isnan(d12) and abs(d12) >= 0.5:
                        records.append(
                            {
                                "smiles_1": s1,
                                "smiles_2": s2,
                                "tanimoto": sim,
                                "delta_delta_12": d12,
                                "delta_12_1": row1.get("delta_12"),
                                "delta_12_2": row2.get("delta_12"),
                            }
                        )
    return pd.DataFrame(records).sort_values("delta_delta_12", key=abs, ascending=False)


def murcko_series_analysis(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    def scaffold(s):
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)

    rows = []
    all_scaffolds = set()
    for iso, df in datasets.items():
        for s in df["canonical_smiles"].dropna():
            sc = scaffold(s)
            if sc:
                all_scaffolds.add(sc)

    for sc in all_scaffolds:
        row = {"scaffold": sc}
        for iso, df in datasets.items():
            mask = df["canonical_smiles"].apply(lambda x: scaffold(x) == sc if scaffold(x) else False)
            sub = df[mask]
            row[f"n_{iso}"] = len(sub)
            row[f"median_pAct_{iso}"] = sub["pActivity"].median() if len(sub) else np.nan
        if row.get("n_JNK1", 0) >= 3:
            d12 = row.get("median_pAct_JNK1", np.nan) - row.get("median_pAct_JNK2", np.nan)
            row["series_delta_12"] = d12
            rows.append(row)

    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("series_delta_12", ascending=False, na_position="last")
    return out


def main():
    parser = argparse.ArgumentParser(description="SAR / MMP / activity cliff analysis")
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "processed")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "sar")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    datasets = {}
    for iso in ["jnk1", "jnk2", "jnk3"]:
        path = args.input / f"{iso}_curated.csv"
        if path.exists():
            datasets[iso.upper()] = pd.read_csv(path)

    # Activity cliffs per isoform
    for iso, df in datasets.items():
        cliffs = find_activity_cliffs(df)
        cliffs.to_csv(args.output / f"activity_cliffs_{iso}.csv", index=False)
        logger.info("%s activity cliffs: %d", iso, len(cliffs))

    # Paired selectivity analysis
    paired_path = args.input / "paired_set.csv"
    if paired_path.exists():
        paired = pd.read_csv(paired_path)
        sel_cliffs = find_selectivity_cliffs(paired)
        sel_cliffs.to_csv(args.output / "selectivity_cliffs.csv", index=False)
        logger.info("Selectivity cliffs: %d", len(sel_cliffs))

        mmp = simple_mmp_analysis(paired)
        mmp.to_csv(args.output / "mmp_selectivity_pairs.csv", index=False)
        logger.info("MMP-like pairs: %d", len(mmp))

    # Scaffold series
    if datasets:
        series = murcko_series_analysis(datasets)
        series.to_csv(args.output / "scaffold_series_analysis.csv", index=False)
        jnk1_sel_series = series[series["series_delta_12"] >= 1.0] if "series_delta_12" in series.columns else series
        jnk1_sel_series.to_csv(args.output / "jnk1_selective_series.csv", index=False)
        logger.info("JNK1-selective series: %d", len(jnk1_sel_series))

    logger.info("SAR analysis complete → %s", args.output)


if __name__ == "__main__":
    main()
