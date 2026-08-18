#!/usr/bin/env python3
"""Sample distill subset D: unlabeled diversity negatives from a large library.

Subset D is for URAT1 retrospective only. It is not the TrueDecoy protocol-selection pool.

Output:
  data/distill/distill_subset_d.csv
  data/distill/distill_manifest.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from utils_ml import canonicalize, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DISTILL_DIR = PROJECT_ROOT / "data" / "distill"

EXCLUDE_SOURCES = [
    PROJECT_ROOT / "data" / "processed" / "urat1_curated.csv",
    PROJECT_ROOT / "data" / "processed" / "nlrp3_records.csv",
    PROJECT_ROOT / "data" / "auxiliary" / "oat_combined_transfer.csv",
    PROJECT_ROOT / "data" / "auxiliary" / "oat1_chembl_curated.csv",
    PROJECT_ROOT / "data" / "auxiliary" / "oat3_chembl_curated.csv",
    PROJECT_ROOT / "data" / "auxiliary" / "oct1_chembl_curated.csv",
    PROJECT_ROOT / "data" / "auxiliary" / "oct2_chembl_curated.csv",
    PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks.csv",
]


def load_exclude_smiles() -> set[str]:
    excluded: set[str] = set()
    for path in EXCLUDE_SOURCES:
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        col = "canonical_smiles" if "canonical_smiles" in df.columns else "Smiles"
        if col not in df.columns:
            continue
        for smi in df[col].dropna().astype(str):
            c = canonicalize(smi) if col == "Smiles" else smi
            if c:
                excluded.add(c)
    return excluded


def read_library_smiles(path: Path, smiles_col: str | None) -> pd.Series:
    path = Path(path)
    if path.suffix.lower() == ".smi":
        lines = [ln.split()[0] for ln in path.read_text().splitlines() if ln.strip()]
        return pd.Series(lines, name="smiles_raw")

    df = pd.read_csv(path, low_memory=False)
    if smiles_col:
        if smiles_col not in df.columns:
            raise ValueError(f"Column {smiles_col!r} not in {path}")
        return df[smiles_col].astype(str)
    for c in ("smiles", "SMILES", "Smiles", "canonical_smiles", "mol"):
        if c in df.columns:
            return df[c].astype(str)
    raise ValueError(f"No SMILES column found in {path}; pass --smiles-col")


def lipinski_pass(mol) -> bool:
    from rdkit.Chem import Descriptors

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    return mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10


def pains_pass(mol) -> bool:
    try:
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        return catalog.GetFirstMatch(mol) is None
    except Exception:
        return True


def prefilt_smiles(raw: str) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    if not lipinski_pass(mol):
        return None
    if not pains_pass(mol):
        return None
    return canonicalize(raw)


def diversity_pick(df: pd.DataFrame, n: int, tanimoto: float = 0.4) -> pd.DataFrame:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.ML.Cluster import Butina

    mols = [Chem.MolFromSmiles(s) for s in df["canonical_smiles"]]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols if m]
    if len(fps) < n:
        return df.head(n)
    dists = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[: i + 1])
        dists.append(1 - sims[-1])
    clusters = Butina.ClusterData(dists, len(fps), tanimoto, isDistData=True)
    picked_idx = [cluster[0] for cluster in clusters]
    if len(picked_idx) >= n:
        return df.iloc[picked_idx[:n]].reset_index(drop=True)
    return df.iloc[picked_idx].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample MASFL distill subset D from a large library")
    parser.add_argument("--library", type=Path, required=True, help="Enamine/ChEMBL CSV or .smi")
    parser.add_argument("--smiles-col", type=str, default=None, help="SMILES column name if CSV")
    parser.add_argument("--n", type=int, default=8000, help="Target sample size (3000–10000)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chunk-size", type=int, default=200_000, help="Rows per read chunk for huge CSV")
    parser.add_argument("--pool-multiplier", type=int, default=5, help="Random pool = n * multiplier before prefilt")
    parser.add_argument("--diversity", action="store_true", help="Butina diversity pick after random pool")
    parser.add_argument("--tanimoto", type=float, default=0.4)
    parser.add_argument("--output-dir", type=Path, default=DISTILL_DIR)
    args = parser.parse_args()

    if not 3000 <= args.n <= 10000:
        print(f"Warning: n={args.n} outside MASFL recommended 3000–10000")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    excluded = load_exclude_smiles()
    print(f"Exclusion set: {len(excluded)} known project SMILES")

    raw_series = read_library_smiles(args.library, args.smiles_col)
    print(f"Library rows read: {len(raw_series)}")

    # Reservoir-style random subsample before expensive RDKit
    pool_n = min(len(raw_series), args.n * args.pool_multiplier)
    pool = raw_series.sample(n=pool_n, random_state=args.seed).tolist()

    rows = []
    for raw in pool:
        if len(rows) >= args.n * 2:
            break
        canon = prefilt_smiles(raw)
        if canon is None or canon in excluded:
            continue
        rows.append({"smiles_raw": raw, "canonical_smiles": canon, "scaffold": murcko_scaffold(canon)})

    if len(rows) < args.n:
        print(f"Warning: only {len(rows)} candidates after prefilt/exclusion; increase --pool-multiplier")

    df = pd.DataFrame(rows).drop_duplicates("canonical_smiles")
    if args.diversity and len(df) > args.n:
        df = diversity_pick(df, args.n, args.tanimoto)
    else:
        df = df.sample(n=min(args.n, len(df)), random_state=args.seed).reset_index(drop=True)

    df["subset"] = "D"
    df["label_type"] = "unlabeled_negative"
    df["source_library"] = args.library.name
    df["has_bioactivity_label"] = False

    out_d = args.output_dir / "distill_subset_d.csv"
    df.to_csv(out_d, index=False)
    print(f"Wrote {len(df)} compounds -> {out_d}")

    manifest_path = args.output_dir / "distill_manifest.csv"
    manifest_row = df[
        ["canonical_smiles", "scaffold", "subset", "label_type", "source_library", "has_bioactivity_label"]
    ].copy()
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path)
        manifest = pd.concat([existing, manifest_row], ignore_index=True).drop_duplicates("canonical_smiles")
    else:
        manifest = manifest_row
    manifest.to_csv(manifest_path, index=False)
    print(f"Manifest -> {manifest_path} ({len(manifest)} total rows)")

    summary = {
        "subset": "D",
        "n_sampled": int(len(df)),
        "n_scaffolds": int(df["scaffold"].nunique()),
        "source_library": str(args.library),
        "excluded_known_smiles": len(excluded),
        "seed": args.seed,
        "diversity": args.diversity,
    }
    with open(args.output_dir / "distill_subset_d_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
