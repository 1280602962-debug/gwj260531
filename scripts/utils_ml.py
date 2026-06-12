"""Shared ML utilities for JNK QSAR pipeline."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

# RDKit descriptors used alongside Morgan fingerprints
RDKIT_DESC_FUNCS = None


def _get_desc_funcs():
    global RDKIT_DESC_FUNCS
    if RDKIT_DESC_FUNCS is None:
        from rdkit.Chem import Descriptors

        RDKIT_DESC_FUNCS = [
            Descriptors.MolWt,
            Descriptors.MolLogP,
            Descriptors.TPSA,
            Descriptors.NumHDonors,
            Descriptors.NumHAcceptors,
            Descriptors.NumRotatableBonds,
            Descriptors.RingCount,
            Descriptors.NumAromaticRings,
            Descriptors.FractionCSP3,
            Descriptors.BertzCT,
            Descriptors.NumHeteroatoms,
            Descriptors.qed,
        ]
    return RDKIT_DESC_FUNCS


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


def murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def curate_isoform_raw(
    path,
    *,
    min_assay_compounds: int = 10,
    pactivity_min: float = 4.0,
    pactivity_max: float = 10.0,
    max_std: float = 0.5,
    max_range: float = 1.0,
) -> pd.DataFrame:
    """Curate one ChEMBL export CSV for single-target QSAR."""
    df = pd.read_csv(path, low_memory=False)
    df = df[df["Standard Type"] == "IC50"]
    df = df[df["Standard Relation"].astype(str).str.strip("'\"") == "="]
    df = df[df["Assay Type"] == "B"]  # biochemical only
    df = df[df["Smiles"].notna()]

    df["pActivity"] = pd.to_numeric(df["pChEMBL Value"], errors="coerce")
    df = df[df["pActivity"].notna()]
    df = df[(df["pActivity"] >= pactivity_min) & (df["pActivity"] <= pactivity_max)]

    if min_assay_compounds > 0:
        assay_counts = df.groupby("Assay ChEMBL ID").size()
        keep_assays = set(assay_counts[assay_counts >= min_assay_compounds].index)
        df = df[df["Assay ChEMBL ID"].isin(keep_assays)]

    rows = []
    for smi, g in df.groupby("Smiles"):
        vals = g["pActivity"].values.astype(float)
        if len(vals) > 1 and (np.std(vals) > max_std or (vals.max() - vals.min()) > max_range):
            continue
        canon = canonicalize(smi)
        if canon is None:
            continue
        rows.append(
            {
                "molecule_chembl_id": g["Molecule ChEMBL ID"].iloc[0],
                "canonical_smiles": canon,
                "pActivity": float(np.median(vals)),
                "n_measurements": len(vals),
            }
        )
    return pd.DataFrame(rows)


def featurize_smiles(smiles_list: list[str], morgan_bits: int = 2048) -> np.ndarray:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    desc_funcs = _get_desc_funcs()
    n_desc = len(desc_funcs)
    X = np.zeros((len(smiles_list), morgan_bits + n_desc), dtype=np.float32)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=morgan_bits)
        DataStructs.ConvertToNumpyArray(fp, X[i, :morgan_bits])
        X[i, morgan_bits:] = [f(mol) for f in desc_funcs]
    return X


def regression_metrics(y_true, y_pred) -> dict:
    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
    if mask.sum() == 0:
        return {"rmse": np.nan, "mae": np.nan, "r2": np.nan, "spearman": np.nan, "n": 0}
    yt, yp = y_true[mask], y_pred[mask]
    rmse = float(np.sqrt(np.mean((yt - yp) ** 2)))
    mae = float(np.mean(np.abs(yt - yp)))
    ss_res = np.sum((yt - yp) ** 2)
    ss_tot = np.sum((yt - yt.mean()) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot else np.nan
    rho, _ = spearmanr(yt, yp)
    return {"rmse": rmse, "mae": mae, "r2": r2, "spearman": float(rho), "n": int(mask.sum())}


def scaffold_holdout_split(smiles: list[str], test_frac=0.1, val_frac=0.1, seed=42):
    scaffolds = [murcko_scaffold(s) for s in smiles]
    groups = np.array(scaffolds)
    idx = np.arange(len(smiles))

    gss_test = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    train_val_idx, test_idx = next(gss_test.split(idx, groups=groups))

    sub_groups = groups[train_val_idx]
    val_size = val_frac / (1 - test_frac)
    gss_val = GroupShuffleSplit(n_splits=1, test_size=val_size, random_state=seed)
    tr_rel, va_rel = next(gss_val.split(train_val_idx, groups=sub_groups))

    return train_val_idx[tr_rel], train_val_idx[va_rel], test_idx


def scaffold_cv_indices(smiles: list[str], n_splits: int = 5, seed: int = 42):
    groups = [murcko_scaffold(s) for s in smiles]
    gkf = GroupKFold(n_splits=n_splits)
    return list(gkf.split(np.arange(len(smiles)), groups=groups))
