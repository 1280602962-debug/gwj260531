#!/usr/bin/env python3
"""
Phase 3: Train multitask + selectivity models for JNK1/2/3 (v2 features).

Pipeline:
  Stage A: Load v2 single-target XGBoost models (Morgan + RDKit descriptors)
  Stage B: Multi-output regression (MTL) on paired set
  Stage C: Selectivity regression (delta_min) + JNK1-selective classification
  Stage D: Unified model bundle for screening / SHAP

Validation: Bemis-Murcko scaffold split

Usage:
    python scripts/04_train_selectivity_model.py --input data/processed --output models
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    roc_auc_score,
)
from sklearn.multioutput import MultiOutputRegressor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import FIGSIZE_DOUBLE, FIGSIZE_SQUARE, apply_journal_style, save_figure  # noqa: E402
from utils_ml import (  # noqa: E402
    featurize_smiles,
    regression_metrics,
    scaffold_holdout_split,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CONFIG_PATH = ROOT / "config" / "targets.yaml"
MORGAN_BITS = 2048
V2_MODEL_DIR = ROOT / "models" / "xgboost"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def xgb_params(config: dict) -> dict:
    cfg = config.get("training", {}).get("xgboost", {})
    return {
        "n_estimators": cfg.get("n_estimators", 2500),
        "max_depth": cfg.get("max_depth", 7),
        "learning_rate": cfg.get("learning_rate", 0.015),
        "subsample": cfg.get("subsample", 0.85),
        "colsample_bytree": cfg.get("colsample_bytree", 0.65),
        "min_child_weight": cfg.get("min_child_weight", 3),
        "reg_alpha": cfg.get("reg_alpha", 0.8),
        "reg_lambda": cfg.get("reg_lambda", 2.5),
        "early_stopping_rounds": cfg.get("early_stopping_rounds", 100),
    }


def train_xgb_regressor(X_train, y_train, X_val, y_val, params: dict):
    import xgboost as xgb

    model = xgb.XGBRegressor(
        **{k: v for k, v in params.items() if k != "early_stopping_rounds"},
        eval_metric="rmse",
        random_state=42,
        n_jobs=-1,
    )
    model.fit(
        X_train,
        y_train,
        eval_set=[(X_val, y_val)],
        verbose=False,
    )
    return model


def load_or_train_single_target(
    iso: str,
    df: pd.DataFrame,
    output_dir: Path,
    params: dict,
) -> tuple[object, dict]:
    """Stage A: reuse v2 per-isoform models when available."""
    v2_path = V2_MODEL_DIR / f"xgboost_{iso.lower()}.joblib"
    df = df.dropna(subset=["canonical_smiles", "pActivity"])
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    X = featurize_smiles(smiles, morgan_bits=MORGAN_BITS)
    tr, va, te = scaffold_holdout_split(smiles)

    if v2_path.exists():
        model = joblib.load(v2_path)
        logger.info("Loaded existing v2 model for %s from %s", iso, v2_path)
    else:
        model = train_xgb_regressor(X[tr], y[tr], X[va], y[va], params)
        v2_path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(model, v2_path)
        logger.info("Trained and saved v2 model for %s", iso)

    pred = model.predict(X[te])
    metrics = regression_metrics(y[te], pred)
    joblib.dump(model, output_dir / f"baseline_{iso.lower()}.joblib")
    logger.info("%s baseline scaffold-test: R²=%.3f, RMSE=%.3f", iso, metrics["r2"], metrics["rmse"])
    return model, metrics


def plot_parity(y_true, y_pred, title: str, path: Path) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
    yt, yp = y_true[mask], y_pred[mask]
    ax.scatter(yt, yp, s=12, alpha=0.6, edgecolors="none")
    lo, hi = min(yt.min(), yp.min()), max(yt.max(), yp.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xlabel("Observed pActivity")
    ax.set_ylabel("Predicted pActivity")
    ax.set_title(title)
    save_figure(path, fig)


def plot_delta_parity(y_true, y_pred, path: Path) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
    yt, yp = y_true[mask], y_pred[mask]
    ax.scatter(yt, yp, s=12, alpha=0.6, edgecolors="none")
    lo, hi = min(yt.min(), yp.min()), max(yt.max(), yp.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xlabel("Observed Δmin (log units)")
    ax.set_ylabel("Predicted Δmin (log units)")
    ax.set_title("Selectivity Regression (Scaffold Holdout)")
    save_figure(path, fig)


def plot_sel_class_distribution(paired: pd.DataFrame, path: Path) -> None:
    if "sel_class" not in paired.columns:
        return
    apply_journal_style()
    counts = paired["sel_class"].value_counts()
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    counts.plot(kind="bar", ax=ax, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Selectivity Class")
    ax.set_ylabel("Compound Count")
    ax.set_title("Paired Set Selectivity Labels")
    ax.tick_params(axis="x", rotation=35)
    save_figure(path, fig)


def train_mtl_regressor(paired: pd.DataFrame, output_dir: Path, params: dict, plot_dir: Path) -> dict:
    """Stage B: multi-output regression on paired set (median-imputed missing targets)."""
    import xgboost as xgb

    paired = paired.dropna(subset=["canonical_smiles"]).copy()
    if "n_isoforms" in paired.columns:
        paired = paired[paired["n_isoforms"] >= 2]

    smiles = paired["canonical_smiles"].tolist()
    targets = ["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]
    Y_raw = paired[targets].values.astype(float)
    Y = Y_raw.copy()
    for j in range(Y.shape[1]):
        col = Y[:, j]
        mask = ~np.isnan(col)
        if mask.any():
            Y[np.isnan(col), j] = np.median(col[mask])

    X = featurize_smiles(smiles, morgan_bits=MORGAN_BITS)

    tr, va, te = scaffold_holdout_split(smiles)
    base_params = {k: v for k, v in params.items() if k != "early_stopping_rounds"}
    mtl = MultiOutputRegressor(
        xgb.XGBRegressor(**base_params, random_state=42, n_jobs=-1)
    )
    mtl.fit(X[tr], Y[tr])
    pred = mtl.predict(X[te])

    metrics = {}
    for i, t in enumerate(targets):
        metrics[t] = regression_metrics(Y_raw[te, i], pred[:, i])
        logger.info("MTL %s: R²=%.3f", t, metrics[t]["r2"])
        plot_parity(
            Y_raw[te, i],
            pred[:, i],
            f"MTL {t} (Scaffold Holdout)",
            plot_dir / f"mtl_parity_{t.lower()}.png",
        )

    joblib.dump(mtl, output_dir / "mtl_regressor.joblib")
    return {"metrics": metrics}


def train_selectivity_models(
    paired: pd.DataFrame,
    config: dict,
    output_dir: Path,
    params: dict,
    plot_dir: Path,
) -> dict:
    """Stage C: delta_min regression + JNK1-selective classification."""
    import xgboost as xgb

    paired = paired.dropna(subset=["canonical_smiles"]).copy()
    if "n_isoforms" in paired.columns:
        paired = paired[paired["n_isoforms"] >= 2]

    plot_sel_class_distribution(paired, plot_dir / "sel_class_distribution.png")

    smiles = paired["canonical_smiles"].tolist()
    X = featurize_smiles(smiles, morgan_bits=MORGAN_BITS)

    y_delta = paired["delta_min"].values.astype(float)
    valid_delta = ~np.isnan(y_delta)

    tr, va, te = scaffold_holdout_split(smiles)
    tr_d = np.array([i for i in tr if valid_delta[i]])
    va_d = np.array([i for i in va if valid_delta[i]])
    te_d = np.array([i for i in te if valid_delta[i]])

    delta_model = train_xgb_regressor(X[tr_d], y_delta[tr_d], X[va_d], y_delta[va_d], params)
    pred_delta = delta_model.predict(X[te_d])
    delta_metrics = regression_metrics(y_delta[te_d], pred_delta)
    plot_delta_parity(y_delta[te_d], pred_delta, plot_dir / "delta_parity_holdout.png")
    joblib.dump(delta_model, output_dir / "selectivity_delta_regressor.joblib")
    logger.info("Selectivity Δ regression: R²=%.3f", delta_metrics["r2"])

    y_class = (paired["sel_class"] == "JNK1-selective").astype(int).values
    clf = None
    clf_metrics: dict = {}
    if y_class.sum() >= 5:
        clf = xgb.XGBClassifier(
            n_estimators=500,
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
            "n_positive_train": int(y_class[tr].sum()),
            "n_positive_test": int(y_class[te].sum()),
            "auc": float(roc_auc_score(y_class[te], prob)) if len(np.unique(y_class[te])) > 1 else np.nan,
            "accuracy": float(accuracy_score(y_class[te], pred_c)),
            "f1": float(f1_score(y_class[te], pred_c, zero_division=0)),
        }
        joblib.dump(clf, output_dir / "selectivity_classifier.joblib")
        logger.info("Selectivity classifier AUC=%.3f", clf_metrics["auc"])
    else:
        clf_metrics = {
            "note": "Insufficient JNK1-selective samples for classifier",
            "n_jnk1_selective": int(y_class.sum()),
        }

    return {"delta_regression": delta_metrics, "classification": clf_metrics, "classifier": clf}


def build_ensemble_bundle(
    output_dir: Path,
    config: dict,
    single_target: dict[str, str],
) -> None:
    """Stage D: save unified model bundle metadata."""
    bundle = {
        "version": "2.0",
        "feature_spec": {
            "morgan_bits": MORGAN_BITS,
            "n_descriptors": 12,
            "total_features": MORGAN_BITS + 12,
        },
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
            "single_target": single_target,
        },
    }
    with open(output_dir / "model_bundle.json", "w") as f:
        json.dump(bundle, f, indent=2)
    logger.info("Saved model bundle metadata")


def serialize(obj):
    if isinstance(obj, (np.integer, np.floating)):
        val = float(obj)
        return val if np.isfinite(val) else None
    if isinstance(obj, dict):
        return {k: serialize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [serialize(v) for v in obj]
    return obj


def main():
    parser = argparse.ArgumentParser(description="Train JNK selectivity models (v2)")
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "processed")
    parser.add_argument("--output", type=Path, default=ROOT / "models")
    parser.add_argument("--plots", type=Path, default=ROOT / "results" / "training")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.plots.mkdir(parents=True, exist_ok=True)

    config = load_config()
    params = xgb_params(config)
    report: dict = {}
    single_target_models: dict[str, object] = {}
    single_target_paths: dict[str, str] = {}

    datasets = {}
    for iso in ["JNK1", "JNK2", "JNK3"]:
        path = args.input / f"{iso.lower()}_curated.csv"
        if path.exists():
            datasets[iso] = pd.read_csv(path)

    if datasets:
        stage_a_metrics = {}
        for iso, df in datasets.items():
            model, metrics = load_or_train_single_target(iso, df, args.output, params)
            stage_a_metrics[iso] = metrics
            single_target_models[iso] = model
            single_target_paths[iso] = str(V2_MODEL_DIR / f"xgboost_{iso.lower()}.joblib")
        report["stage_a_baselines"] = stage_a_metrics

    paired_path = args.input / "paired_set.csv"
    mtl_model = delta_model = classifier = None
    if paired_path.exists():
        paired = pd.read_csv(paired_path)
        report["stage_b_mtl"] = train_mtl_regressor(paired, args.output, params, args.plots)
        mtl_model = joblib.load(args.output / "mtl_regressor.joblib")
        report["stage_b_mtl"] = {"metrics": report["stage_b_mtl"]["metrics"]}
        stage_c = train_selectivity_models(paired, config, args.output, params, args.plots)
        report["stage_c_selectivity"] = {
            "delta_regression": stage_c["delta_regression"],
            "classification": stage_c["classification"],
        }
        delta_model = joblib.load(args.output / "selectivity_delta_regressor.joblib")
        clf_path = args.output / "selectivity_classifier.joblib"
        classifier = joblib.load(clf_path) if clf_path.exists() else stage_c.get("classifier")
    else:
        logger.warning("paired_set.csv not found. Run 00_prepare_user_data.py first.")

    build_ensemble_bundle(args.output, config, single_target_paths)

    with open(args.output / "training_report.json", "w") as f:
        json.dump(serialize(report), f, indent=2)

    results_report = args.plots / "training_report.json"
    with open(results_report, "w") as f:
        json.dump(serialize(report), f, indent=2)

    if mtl_model is not None and delta_model is not None:
        joblib.dump(
            {
                "version": "2.0",
                "feature_spec": {
                    "morgan_bits": MORGAN_BITS,
                    "n_descriptors": 12,
                    "total_features": MORGAN_BITS + 12,
                },
                "mtl": mtl_model,
                "delta": delta_model,
                "classifier": classifier,
                "single_target": single_target_models,
                "thresholds": config["selectivity"],
                "scoring_weights": {
                    "w_pAct_JNK1": 0.35,
                    "w_delta_min": 0.30,
                    "w_neg_JNK23": 0.20,
                    "w_qed": 0.10,
                    "w_sa": -0.05,
                },
            },
            args.output / "best_model.joblib",
        )
        logger.info("Saved best_model.joblib (v2 ensemble bundle)")

    logger.info("Training complete → %s", args.output)


if __name__ == "__main__":
    main()
