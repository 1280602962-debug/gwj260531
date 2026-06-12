#!/usr/bin/env python3
"""
Compare Chemprop 2.0 MTL vs XGBoost MTL on identical scaffold splits.

Runs:
  1. Data preparation (if needed)
  2. Dataset similarity analysis
  3. XGBoost MTL training + test evaluation
  4. Chemprop MTL training + test evaluation
  5. Model selection report

Usage:
    python scripts/07_compare_models.py
    python scripts/07_compare_models.py --skip-chemprop  # XGBoost only
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.multioutput import MultiOutputRegressor
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
TARGET_COLS = ["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]
N_BITS = 2048


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


def smiles_to_fp_matrix(smiles_list: list[str]) -> np.ndarray:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    X = np.zeros((len(smiles_list), N_BITS), dtype=np.int8)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=N_BITS)
            DataStructs.ConvertToNumpyArray(fp, X[i])
    return X


def train_xgboost_mtl(train_df, val_df, test_df, output_dir: Path) -> dict:
    import xgboost as xgb

    output_dir.mkdir(parents=True, exist_ok=True)
    smiles_train = train_df["smiles"].tolist()
    Y_train = train_df[TARGET_COLS].values.astype(float)
    X_train = smiles_to_fp_matrix(smiles_train)

    smiles_val = val_df["smiles"].tolist()
    Y_val = val_df[TARGET_COLS].values.astype(float)
    X_val = smiles_to_fp_matrix(smiles_val)

    mtl = MultiOutputRegressor(
        xgb.XGBRegressor(
            n_estimators=400,
            max_depth=6,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=42,
            n_jobs=-1,
            early_stopping_rounds=30,
        )
    )
    # MultiOutputRegressor doesn't support eval_set directly; train per-task with shared params
    models = []
    for i, col in enumerate(TARGET_COLS):
        mask_tr = ~np.isnan(Y_train[:, i])
        mask_va = ~np.isnan(Y_val[:, i])
        m = xgb.XGBRegressor(
            n_estimators=400,
            max_depth=6,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=42,
            n_jobs=-1,
            early_stopping_rounds=30,
        )
        m.fit(
            X_train[mask_tr],
            Y_train[mask_tr, i],
            eval_set=[(X_val[mask_va], Y_val[mask_va, i])],
            verbose=False,
        )
        models.append(m)

    joblib.dump(models, output_dir / "xgboost_mtl_models.joblib")

    metrics = {}
    smiles_test = test_df["smiles"].tolist()
    X_test = smiles_to_fp_matrix(smiles_test)
    Y_test = test_df[TARGET_COLS].values.astype(float)

    for i, col in enumerate(TARGET_COLS):
        mask = ~np.isnan(Y_test[:, i])
        if mask.sum() == 0:
            continue
        pred = models[i].predict(X_test[mask])
        metrics[col] = regression_metrics(Y_test[mask, i], pred)

    return metrics


def aggregate_score(metrics: dict) -> dict:
    """Compute mean metrics across tasks (only tasks with data)."""
    r2s, rmses, spearmans, ns = [], [], [], []
    for v in metrics.values():
        if v["n"] > 0 and not np.isnan(v["r2"]):
            r2s.append(v["r2"])
            rmses.append(v["rmse"])
            spearmans.append(v["spearman"])
            ns.append(v["n"])
    return {
        "mean_r2": float(np.mean(r2s)) if r2s else np.nan,
        "mean_rmse": float(np.mean(rmses)) if rmses else np.nan,
        "mean_spearman": float(np.mean(spearmans)) if spearmans else np.nan,
        "total_eval_points": int(sum(ns)),
    }


def select_winner(xgb_metrics, chemprop_metrics) -> dict:
    xgb_agg = aggregate_score(xgb_metrics)
    cp_agg = aggregate_score(chemprop_metrics)

    # Primary: mean R²; tiebreaker: mean Spearman, then lower RMSE
    def score(agg):
        return (
            agg["mean_r2"] if not np.isnan(agg["mean_r2"]) else -999,
            agg["mean_spearman"] if not np.isnan(agg["mean_spearman"]) else -999,
            -agg["mean_rmse"] if not np.isnan(agg["mean_rmse"]) else -999,
        )

    xgb_s, cp_s = score(xgb_agg), score(cp_agg)
    if cp_s > xgb_s:
        winner = "Chemprop 2.0 MTL"
        reason = "Higher mean R² and/or Spearman on scaffold-test set"
    elif xgb_s > cp_s:
        winner = "XGBoost MTL"
        reason = "Higher mean R² and/or Spearman on scaffold-test set"
    else:
        winner = "Tie"
        reason = "Equivalent aggregate performance"

    return {
        "winner": winner,
        "reason": reason,
        "xgboost_aggregate": xgb_agg,
        "chemprop_aggregate": cp_agg,
        "recommendation": (
            f"Use **{winner}** as primary activity predictor for virtual screening."
            if winner != "Tie"
            else "Use ensemble of both models."
        ),
    }


def write_report(
    output_path: Path,
    data_summary: dict,
    xgb_metrics: dict,
    chemprop_metrics: dict,
    selection: dict,
):
    lines = [
        "# JNK1/2/3 Model Comparison Report",
        "",
        "## Data Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
    ]
    for k, v in data_summary.items():
        lines.append(f"| {k} | {v} |")

    lines += ["", "## Test Set Performance (Scaffold Split)", ""]

    for model_name, metrics in [("XGBoost MTL", xgb_metrics), ("Chemprop 2.0 MTL", chemprop_metrics)]:
        lines += [f"### {model_name}", "", "| Target | R² | RMSE | MAE | Spearman | n |", "|--------|-----|------|-----|----------|---|"]
        for col, m in metrics.items():
            lines.append(
                f"| {col} | {m['r2']:.3f} | {m['rmse']:.3f} | {m['mae']:.3f} | {m['spearman']:.3f} | {m['n']} |"
            )
        agg = aggregate_score(metrics)
        lines.append(
            f"| **Mean** | **{agg['mean_r2']:.3f}** | **{agg['mean_rmse']:.3f}** | — | **{agg['mean_spearman']:.3f}** | {agg['total_eval_points']} |"
        )
        lines.append("")

    lines += [
        "## Model Selection",
        "",
        f"**Winner: {selection['winner']}**",
        "",
        f"Reason: {selection['reason']}",
        "",
        f"Recommendation: {selection['recommendation']}",
        "",
        "## Notes",
        "",
        "- Both models trained on identical scaffold-based train/val/test splits.",
        "- Missing JNK isoform labels handled natively (Chemprop mask; XGBoost per-task training).",
        "- JNK1 has fewer data points; compare JNK1 task performance carefully.",
        "- For selectivity modeling + SHAP, continue using XGBoost selective models (script 04/05).",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def run_script(script: str, extra_args: list[str] | None = None):
    cmd = [sys.executable, str(ROOT / "scripts" / script)] + (extra_args or [])
    logger.info("Running %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Compare Chemprop vs XGBoost")
    parser.add_argument("--skip-prepare", action="store_true")
    parser.add_argument("--skip-similarity", action="store_true")
    parser.add_argument("--skip-chemprop", action="store_true")
    parser.add_argument("--epochs", type=int, default=40)
    args = parser.parse_args()

    processed = ROOT / "data" / "processed"
    results = ROOT / "results"
    results.mkdir(exist_ok=True)

    if not args.skip_prepare:
        run_script("00_prepare_user_data.py")

    if not args.skip_similarity:
        run_script("02_dataset_similarity.py", ["--input", str(processed), "--output", str(results / "similarity")])

    splits_dir = processed / "splits"
    train_df = pd.read_csv(splits_dir / "train.csv")
    val_df = pd.read_csv(splits_dir / "val.csv")
    test_df = pd.read_csv(splits_dir / "test.csv")

    logger.info("Training XGBoost MTL ...")
    xgb_metrics = train_xgboost_mtl(train_df, val_df, test_df, ROOT / "models" / "xgboost")
    logger.info("XGBoost metrics: %s", xgb_metrics)

    chemprop_metrics = {}
    if not args.skip_chemprop:
        run_script(
            "04b_train_chemprop_mtl.py",
            ["--splits-dir", str(splits_dir), "--output", str(ROOT / "models" / "chemprop"), "--epochs", str(args.epochs)],
        )
        metrics_path = ROOT / "models" / "chemprop" / "chemprop_metrics.json"
        if metrics_path.exists():
            with open(metrics_path) as f:
                chemprop_metrics = json.load(f)["test_metrics"]

    selection = select_winner(xgb_metrics, chemprop_metrics)

    data_summary = {}
    summary_path = processed / "data_summary.json"
    if summary_path.exists():
        with open(summary_path) as f:
            data_summary = json.load(f)

    comparison = {
        "xgboost": xgb_metrics,
        "chemprop": chemprop_metrics,
        "selection": selection,
        "data_summary": data_summary,
    }
    comparison_dir = results / "model_comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    with open(comparison_dir / "comparison.json", "w") as f:
        json.dump(comparison, f, indent=2)

    write_report(comparison_dir / "MODEL_COMPARISON_REPORT.md", data_summary, xgb_metrics, chemprop_metrics, selection)

    logger.info("=" * 60)
    logger.info("WINNER: %s", selection["winner"])
    logger.info("Report: %s", comparison_dir / "MODEL_COMPARISON_REPORT.md")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
