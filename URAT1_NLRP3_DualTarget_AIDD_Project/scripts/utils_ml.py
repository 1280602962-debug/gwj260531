"""Shared ML utilities for URAT1/NLRP3 (asymmetric dual-evidence models)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

RDKIT_DESC_NAMES = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "NumAromaticRings",
    "FractionCSP3",
    "BertzCT",
    "NumHeteroatoms",
    "QED",
]

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


def _clean_relation(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip("'\"")


def fill_pactivity_column(df: pd.DataFrame) -> pd.Series:
    """Use pChEMBL Value; if missing, derive from Standard Value + Units (nM/µM)."""
    pact = pd.to_numeric(df["pChEMBL Value"], errors="coerce")
    miss = pact.isna()
    if not miss.any():
        return pact
    units = df.loc[miss, "Standard Units"].astype(str).str.lower()
    vals = pd.to_numeric(df.loc[miss, "Standard Value"], errors="coerce")
    derived = pd.Series(np.nan, index=df.index, dtype=float)
    nm = miss & units.str.contains("nm", na=False) & vals.notna() & (vals > 0)
    derived.loc[nm] = 9.0 - np.log10(vals.loc[nm])
    um = miss & units.str.contains(r"um|µm|microm", na=False, regex=True) & vals.notna() & (vals > 0)
    derived.loc[um] = 6.0 - np.log10(vals.loc[um])
    return pact.fillna(derived)


def curate_urat1_raw(
    path: str | Path,
    *,
    pactivity_min: float = 4.0,
    pactivity_max: float = 10.0,
    max_std: float = 0.5,
    max_range: float = 1.0,
) -> pd.DataFrame:
    """Molecule-level URAT1 curation for regression."""
    df = pd.read_csv(path, low_memory=False)
    df = df[_clean_relation(df["Standard Relation"]) == "="]
    df = df[df["Smiles"].notna()]
    df["pActivity"] = fill_pactivity_column(df)
    df = df[df["pActivity"].between(pactivity_min, pactivity_max)]

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
                "molecule_name": g["Molecule Name"].iloc[0] if "Molecule Name" in g else "",
                "canonical_smiles": canon,
                "pActivity": float(np.median(vals)),
                "n_measurements": int(len(vals)),
            }
        )
    out = pd.DataFrame(rows).drop_duplicates("canonical_smiles").reset_index(drop=True)
    out["scaffold"] = out["canonical_smiles"].map(murcko_scaffold)
    return out


def curate_nlrp3_records(
    path: str | Path,
    *,
    pactivity_min: float = 4.0,
    pactivity_max: float = 10.0,
    min_assay_compounds: int = 5,
    active_threshold: float = 6.0,
) -> pd.DataFrame:
    """Record-level NLRP3 curation for assay-conditioned classification."""
    df = pd.read_csv(path, low_memory=False)
    df = df[df["Assay Description"].str.contains(r"IL-1", case=False, na=False)]
    df = df[df["Assay Type"] == "B"]
    df = df[_clean_relation(df["Standard Relation"]) == "="]
    df = df[df["Smiles"].notna()]
    df["pActivity"] = pd.to_numeric(df["pChEMBL Value"], errors="coerce")
    df = df[df["pActivity"].between(pactivity_min, pactivity_max)]

    if min_assay_compounds > 0:
        assay_counts = df.groupby("Assay ChEMBL ID").size()
        keep_assays = set(assay_counts[assay_counts >= min_assay_compounds].index)
        df = df[df["Assay ChEMBL ID"].isin(keep_assays)]

    df["canonical_smiles"] = df["Smiles"].map(canonicalize)
    df = df[df["canonical_smiles"].notna()].copy()
    df["active"] = (df["pActivity"] >= active_threshold).astype(int)
    df["scaffold"] = df["canonical_smiles"].map(murcko_scaffold)
    df["assay_cell_type"] = df["Assay Cell Type"].fillna("unknown").astype(str)
    return df.reset_index(drop=True)


def feature_names(morgan_bits: int = 2048) -> list[str]:
    return [f"Bit_{i}" for i in range(morgan_bits)] + RDKIT_DESC_NAMES


def featurize_smiles(smiles_list: list[str], morgan_bits: int = 2048) -> np.ndarray:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    desc_funcs = _get_desc_funcs()
    n_desc = len(desc_funcs)
    x = np.zeros((len(smiles_list), morgan_bits + n_desc), dtype=np.float32)
    for i, smi in enumerate(smiles_list):
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=morgan_bits)
        DataStructs.ConvertToNumpyArray(fp, x[i, :morgan_bits])
        x[i, morgan_bits:] = [f(mol) for f in desc_funcs]
    return x


def max_tanimoto_to_library(query_smiles: str, library_smiles: list[str], morgan_bits: int = 2048) -> float:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    qmol = Chem.MolFromSmiles(query_smiles)
    if qmol is None:
        return 0.0
    qfp = AllChem.GetMorganFingerprintAsBitVect(qmol, 2, nBits=morgan_bits)
    best = 0.0
    for smi in library_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=morgan_bits)
        best = max(best, DataStructs.TanimotoSimilarity(qfp, fp))
    return float(best)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
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


def classification_metrics(y_true: np.ndarray, y_prob: np.ndarray) -> dict:
    mask = ~np.isnan(y_true) & ~np.isnan(y_prob)
    if mask.sum() == 0 or len(np.unique(y_true[mask])) < 2:
        return {"auroc": np.nan, "auprc": np.nan, "n": 0}
    yt, yp = y_true[mask].astype(int), y_prob[mask]
    return {
        "auroc": float(roc_auc_score(yt, yp)),
        "auprc": float(average_precision_score(yt, yp)),
        "n": int(mask.sum()),
    }


def enrichment_factor(y_true: np.ndarray, y_score: np.ndarray, fraction: float = 0.1) -> float:
    """EF@fraction for binary actives (1) vs inactives (0). Higher score = more likely active."""
    mask = ~np.isnan(y_true) & ~np.isnan(y_score)
    if mask.sum() == 0:
        return np.nan
    yt = y_true[mask].astype(int)
    ys = y_score[mask]
    n = max(1, int(len(yt) * fraction))
    order = np.argsort(-ys)
    top = yt[order[:n]]
    overall = yt.mean()
    if overall == 0:
        return np.nan
    return float(top.mean() / overall)


def regression_enrichment_factor(
    y_true: np.ndarray, y_pred: np.ndarray, threshold: float = 6.0, fraction: float = 0.1
) -> float:
    y_bin = (y_true >= threshold).astype(int)
    return enrichment_factor(y_bin, y_pred, fraction=fraction)


def roc_auc_binary(y_true: np.ndarray, y_score: np.ndarray, threshold: float = 6.0) -> float:
    from sklearn.metrics import roc_auc_score

    y_bin = (y_true >= threshold).astype(int)
    if len(np.unique(y_bin)) < 2:
        return float("nan")
    return float(roc_auc_score(y_bin, y_score))


def scaffold_cv_indices(smiles: list[str], n_splits: int = 5) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = [murcko_scaffold(s) for s in smiles]
    gkf = GroupKFold(n_splits=n_splits)
    return list(gkf.split(np.arange(len(smiles)), groups=groups))


def scaffold_holdout_split(smiles: list[str], test_frac: float = 0.1, val_frac: float = 0.1, seed: int = 42):
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


class SplitConformalRegressor:
    """Split conformal prediction intervals for regression."""

    def __init__(self, alpha: float = 0.1):
        self.alpha = alpha
        self.q: float | None = None

    def fit(self, y_true: np.ndarray, y_pred: np.ndarray) -> "SplitConformalRegressor":
        residuals = np.abs(y_true - y_pred)
        self.q = float(np.quantile(residuals, 1 - self.alpha))
        return self

    def interval(self, y_pred: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.q is None:
            raise RuntimeError("Call fit() before interval()")
        return y_pred - self.q, y_pred + self.q

    @property
    def half_width(self) -> float:
        return float(self.q) if self.q is not None else np.nan


def assay_one_hot_matrix(assay_ids: pd.Series, top_assays: list[str]) -> np.ndarray:
    cols = {a: i for i, a in enumerate(top_assays)}
    mat = np.zeros((len(assay_ids), len(top_assays) + 1), dtype=np.float32)
    for i, aid in enumerate(assay_ids.astype(str)):
        if aid in cols:
            mat[i, cols[aid]] = 1.0
        else:
            mat[i, -1] = 1.0
    return mat


def save_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
