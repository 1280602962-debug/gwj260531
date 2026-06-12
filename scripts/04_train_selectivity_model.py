#!/usr/bin/env python3
"""
Phase 3: Train multitask + selectivity models for JNK1/2/3.

Pipeline:
  Stage A: XGBoost single-target baselines (JNK1/2/3)
  Stage B: XGBoost multi-output regression (MTL-style)
  Stage C: Selectivity regression (delta_min) + classification
  Stage D: Stacking ensemble → saved model bundle

Validation: scaffold split (Bemis-Murcko)

Usage:
    python scripts/04_train_selectivity_model.py --input data/processed --output models
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.multioutput import MultiOutputRegressor

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "targets.yaml"


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


def smiles_to_fp_matrix(smiles_list: list[str], radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    X = np.zeros((len(smiles_list), n_bits), dtype=np.int8)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            DataStructs.ConvertToNumpyArray(fp, X[i])
    return X


def scaffold_split(smiles: list[str], test_frac: float = 0.1, val_frac: float = 0.1, seed: int = 42):
    scaffolds = [murcko_scaffold(s) for s in smiles]
    groups = np.array(scaffolds)
    idx = np.arange(len(smiles))

    gss_test = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    train_val_idx, test_idx = next(gss_test.split(idx, groups=groups))

    sub_groups = groups[train_val_idx]
    val_size = val_frac / (1 - test_frac)
    gss_val = GroupShuffleSplit(n_splits=1, test_size=val_size, random_state=seed)
    train_idx, val_idx = next(gss_val.split(train_val_idx, groups=sub_groups))

    return train_val_idx[train_idx], train_val_idx[val_idx], test_idx


def regression_metrics(y_true, y_pred) -> dict:
    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
    if mask.sum() == 0:
        return {"rmse": np.nan, "mae": np.nan, "r2": np.nan, "n": 0}
    yt, yp = y_true[mask], y_pred[mask]
    rmse = float(np.sqrt(mean_squared_error(yt, yp)))
    mae = float(mean_absolute_error(yt, yp))
    ss_res = np.sum((yt - yp) ** 2)
    ss_tot = np.sum((yt - yt.mean()) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot else np.nan
    return {"rmse": rmse, "mae": mae, "r2": r2, "n": int(mask.sum())}


def train_xgb_regressor(X_train, y_train, X_val, y_val, **kwargs):
    import xgboost as xgb

    model = xgb.XGBRegressor(
        n_estimators=kwargs.get("n_estimators", 500),
        max_depth=kwargs.get("max_depth", 6),
        learning_rate=kwargs.get("learning_rate", 0.05),
        subsample=0.8,
        colsample_bytree=0.8,
        reg_alpha=0.1,
        reg_lambda=1.0,
        early_stopping_rounds=30,
        eval_metric="rmse",
        random_state=42,
        n_jobs=-1,
    )
    model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)
    return model


def train_single_target_baselines(datasets: dict[str, pd.DataFrame], output_dir: Path) -> dict:
    """Stage A: independent XGBoost per isoform."""
    results = {}
    models = {}
    for iso, df in datasets.items():
        df = df.dropna(subset=["canonical_smiles", "pActivity"])
        smiles = df["canonical_smiles"].tolist()
        y = df["pActivity"].values
        X = smiles_to_fp_matrix(smiles)

        tr, va, te = scaffold_split(smiles)
        model = train_xgb_regressor(X[tr], y[tr], X[va], y[va])
        pred = model.predict(X[te])
        metrics = regression_metrics(y[te], pred)
        results[iso] = metrics
        models[iso] = model
        joblib.dump(model, output_dir / f"baseline_{iso.lower()}.joblib")
        logger.info("%s baseline scaffold-test: R²=%.3f, RMSE=%.3f", iso, metrics["r2"], metrics["rmse"])
    return {"metrics": results, "models": models}


def train_mtl_regressor(paired: pd.DataFrame, output_dir: Path) -> dict:
    """Stage B: multi-output regression on paired set."""
    paired = paired.dropna(subset=["canonical_smiles"])
    smiles = paired["canonical_smiles"].tolist()
    targets = ["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]
    Y = paired[targets].values.astype(float)
    X = smiles_to_fp_matrix(smiles)

    tr, va, te = scaffold_split(smiles)
    base = train_xgb_regressor(
        X[tr], Y[tr, 0], X[va], Y[va, 0]
    )  # template params
    import xgboost as xgb

    mtl = MultiOutputRegressor(
        xgb.XGBRegressor(
            n_estimators=500,
            max_depth=6,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=42,
            n_jobs=-1,
        )
    )
    mtl.fit(X[tr], Y[tr])
    pred = mtl.predict(X[te])

    metrics = {}
    for i, t in enumerate(targets):
        metrics[t] = regression_metrics(Y[te, i], pred[:, i])
        logger.info("MTL %s: R²=%.3f", t, metrics[t]["r2"])

    joblib.dump(mtl, output_dir / "mtl_regressor.joblib")
    return {"metrics": metrics, "model": mtl}


def train_selectivity_models(paired: pd.DataFrame, config: dict, output_dir: Path) -> dict:
    """Stage C: delta_min regression + JNK1-selective classification."""
    import xgboost as xgb

    paired = paired.dropna(subset=["canonical_smiles"]).copy()
    paired = paired[paired["n_isoforms"] >= 2] if "n_isoforms" in paired.columns else paired

    smiles = paired["canonical_smiles"].tolist()
    X = smiles_to_fp_matrix(smiles)

    # Regression target
    y_delta = paired["delta_min"].values.astype(float)
    valid_delta = ~np.isnan(y_delta)

    tr, va, te = scaffold_split(smiles)
    tr_d = np.array([i for i in tr if valid_delta[i]])
    va_d = np.array([i for i in va if valid_delta[i]])
    te_d = np.array([i for i in te if valid_delta[i]])

    delta_model = train_xgb_regressor(X[tr_d], y_delta[tr_d], X[va_d], y_delta[va_d])
    pred_delta = delta_model.predict(X[te_d])
    delta_metrics = regression_metrics(y_delta[te_d], pred_delta)
    joblib.dump(delta_model, output_dir / "selectivity_delta_regressor.joblib")
    logger.info("Selectivity Δ regression: R²=%.3f", delta_metrics["r2"])

    # Classification
    y_class = (paired["sel_class"] == "JNK1-selective").astype(int).values
    if y_class.sum() >= 5:
        clf = xgb.XGBClassifier(
            n_estimators=300,
            max_depth=5,
            learning_rate=0.05,
            scale_pos_weight=max(1, (len(y_class) - y_class.sum()) / max(y_class.sum(), 1)),
            random_state=42,
            eval_metric="logloss",
            n_jobs=-1,
        )
        clf.fit(X[tr], y_class[tr], eval_set=[(X[va], y_class[va])], verbose=False)
        prob = clf.predict_proba(X[te])[:, 1]
        pred_c = (prob >= 0.5).astype(int)
        clf_metrics = {
            "auc": float(roc_auc_score(y_class[te], prob)) if len(np.unique(y_class[te])) > 1 else np.nan,
            "accuracy": float(accuracy_score(y_class[te], pred_c)),
            "f1": float(f1_score(y_class[te], pred_c, zero_division=0)),
        }
        joblib.dump(clf, output_dir / "selectivity_classifier.joblib")
        logger.info("Selectivity classifier AUC=%.3f", clf_metrics["auc"])
    else:
        clf_metrics = {"note": "Insufficient JNK1-selective samples for classifier"}

    return {"delta_regression": delta_metrics, "classification": clf_metrics}


def build_ensemble_bundle(output_dir: Path, config: dict):
    """Stage D: save unified model bundle with scoring function metadata."""
    bundle = {
        "version": "1.0",
        "targets": config["targets"],
        "scoring_weights": {
            "w_pAct_JNK1": 0.35,
            "w_delta_min": 0.30,
            "w_neg_JNK23": 0.20,
            "w_qed": 0.10,
            "w_sa": -0.05,
        },
        "thresholds": config["selectivity"],
        "models": {
            "mtl": "mtl_regressor.joblib",
            "delta": "selectivity_delta_regressor.joblib",
            "classifier": "selectivity_classifier.joblib",
        },
    }
    with open(output_dir / "model_bundle.json", "w") as f:
        json.dump(bundle, f, indent=2)
    logger.info("Saved model bundle metadata")


def main():
    parser = argparse.ArgumentParser(description="Train JNK selectivity models")
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "processed")
    parser.add_argument("--output", type=Path, default=ROOT / "models")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config = load_config()
    datasets = {}
    for iso in ["JNK1", "JNK2", "JNK3"]:
        path = args.input / f"{iso.lower()}_curated.csv"
        if path.exists():
            datasets[iso] = pd.read_csv(path)

    report = {}

    # Stage A
    if datasets:
        stage_a = train_single_target_baselines(datasets, args.output)
        report["stage_a_baselines"] = stage_a["metrics"]

    # Stage B & C require paired set
    paired_path = args.input / "paired_set.csv"
    if paired_path.exists():
        paired = pd.read_csv(paired_path)
        report["stage_b_mtl"] = train_mtl_regressor(paired, args.output)
        report["stage_c_selectivity"] = train_selectivity_models(paired, config, args.output)
    else:
        logger.warning("paired_set.csv not found. Run 01_download with --build-paired")

    build_ensemble_bundle(args.output, config)

    # Save full metrics report
    def serialize(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, dict):
            return {k: serialize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [serialize(v) for v in obj]
        return obj

    with open(args.output / "training_report.json", "w") as f:
        json.dump(serialize(report), f, indent=2)

    # Convenience alias for downstream scripts
    if (args.output / "mtl_regressor.joblib").exists():
        joblib.dump(
            {
                "mtl": joblib.load(args.output / "mtl_regressor.joblib"),
                "delta": joblib.load(args.output / "selectivity_delta_regressor.joblib")
                if (args.output / "selectivity_delta_regressor.joblib").exists()
                else None,
                "classifier": joblib.load(args.output / "selectivity_classifier.joblib")
                if (args.output / "selectivity_classifier.joblib").exists()
                else None,
            },
            args.output / "best_model.joblib",
        )
        logger.info("Saved best_model.joblib (ensemble bundle)")

    logger.info("Training complete → %s", args.output)


if __name__ == "__main__":
    main()
