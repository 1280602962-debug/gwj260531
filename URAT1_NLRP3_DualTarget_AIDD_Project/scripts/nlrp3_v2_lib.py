"""Model and calibrator classes for the NLRP3 v2 classifier.

Kept in an importable module (not in the training entry point) so that the
pickled bundles in results/nlrp3_model_v2/ can be loaded by any script.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from utils_ml import assay_one_hot_matrix, featurize_smiles

N_TOP_ASSAYS = 25
N_DEPLOY_ASSAYS = 5


def xgb_classifier(seed: int) -> xgb.XGBClassifier:
    """Same hyper-parameters as the historical production classifier."""
    return xgb.XGBClassifier(
        n_estimators=400,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=seed,
        n_jobs=-1,
        eval_metric="logloss",
    )


def logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


class PlattCalibrator:
    def fit(self, p: np.ndarray, y: np.ndarray) -> "PlattCalibrator":
        self.lr = LogisticRegression(C=1e6, solver="lbfgs")
        self.lr.fit(logit(p).reshape(-1, 1), y)
        return self

    def predict(self, p: np.ndarray) -> np.ndarray:
        return self.lr.predict_proba(logit(p).reshape(-1, 1))[:, 1]


class IsotonicCalibrator:
    def fit(self, p: np.ndarray, y: np.ndarray) -> "IsotonicCalibrator":
        self.iso = IsotonicRegression(out_of_bounds="clip")
        self.iso.fit(p, y)
        return self

    def predict(self, p: np.ndarray) -> np.ndarray:
        return np.asarray(self.iso.predict(p), dtype=float)


def fit_calibrators(p_cal: np.ndarray, y_cal: np.ndarray, min_n: int = 30) -> dict:
    """Calibrators are fitted only on molecules the base model never saw."""
    out: dict = {"none": None}
    if len(np.unique(y_cal)) == 2 and len(y_cal) >= min_n:
        out["platt"] = PlattCalibrator().fit(p_cal, y_cal)
        out["isotonic"] = IsotonicCalibrator().fit(p_cal, y_cal)
    return out


def apply_calibrator(cal, p: np.ndarray) -> np.ndarray:
    return np.asarray(p, dtype=float) if cal is None else cal.predict(p)


class StructureOnlyModel:
    """Molecule-level classifier on Morgan2048 + 12 RDKit descriptors."""

    name = "structure_only"

    def fit(self, mol_fit: pd.DataFrame, seed: int) -> "StructureOnlyModel":
        x = featurize_smiles(mol_fit["canonical_smiles"].tolist())
        self.scaler = StandardScaler().fit(x)
        self.clf = xgb_classifier(seed)
        self.clf.fit(self.scaler.transform(x), mol_fit["active"].values.astype(int))
        return self

    def predict_molecules(self, smiles: list[str]) -> np.ndarray:
        x = self.scaler.transform(featurize_smiles(list(smiles)))
        return self.clf.predict_proba(x)[:, 1]


class AssayConditionedModel:
    """Record-level classifier with assay one-hot, as in the historical pipeline.

    The top-25 assay vocabulary and the five deployment contexts are derived only
    from the records available for fitting, never from held-out molecules.
    """

    name = "assay_conditioned"

    def fit(self, rec_fit: pd.DataFrame, seed: int) -> "AssayConditionedModel":
        counts = rec_fit["Assay ChEMBL ID"].value_counts()
        self.top_assay_ids = counts.head(N_TOP_ASSAYS).index.astype(str).tolist()
        self.deploy_assay_ids = counts.head(N_DEPLOY_ASSAYS).index.astype(str).tolist()

        x_mol = featurize_smiles(rec_fit["canonical_smiles"].tolist())
        x = np.hstack([x_mol, assay_one_hot_matrix(rec_fit["Assay ChEMBL ID"], self.top_assay_ids)])
        self.scaler = StandardScaler().fit(x)
        w = rec_fit["Assay ChEMBL ID"].map(lambda a: 1.0 / np.sqrt(counts.get(a, 1))).values
        self.clf = xgb_classifier(seed)
        self.clf.fit(self.scaler.transform(x), rec_fit["active"].values.astype(int), sample_weight=w)
        return self

    def _predict(self, smiles: list[str], assay_ids: pd.Series) -> np.ndarray:
        x = np.hstack(
            [featurize_smiles(list(smiles)), assay_one_hot_matrix(assay_ids, self.top_assay_ids)]
        )
        return self.clf.predict_proba(self.scaler.transform(x))[:, 1]

    def predict_molecules(self, smiles: list[str]) -> np.ndarray:
        """Deployment rule: max over the five fixed assay contexts."""
        smiles = list(smiles)
        probs = [
            self._predict(smiles, pd.Series([aid] * len(smiles)))
            for aid in self.deploy_assay_ids
        ]
        return np.max(np.vstack(probs), axis=0)

    def predict_molecules_legacy(self, records: pd.DataFrame, smiles: list[str]) -> np.ndarray:
        """Historical rule: score each record with its true assay ID, take max per molecule."""
        tmp = records.copy()
        tmp["prob"] = self._predict(tmp["canonical_smiles"].tolist(), tmp["Assay ChEMBL ID"])
        return tmp.groupby("canonical_smiles")["prob"].max().reindex(list(smiles)).values
