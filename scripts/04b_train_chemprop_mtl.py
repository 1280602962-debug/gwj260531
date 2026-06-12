#!/usr/bin/env python3
"""
Phase 3b: Train Chemprop 2.0 native multitask model on JNK1/2/3.

Uses pre-generated scaffold splits from 00_prepare_user_data.py.

Usage:
    python scripts/04b_train_chemprop_mtl.py \
        --splits-dir data/processed/splits \
        --output models/chemprop
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
TARGET_COLS = ["pAct_JNK1", "pAct_JNK2", "pAct_JNK3"]


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
    from scipy.stats import spearmanr

    rho, _ = spearmanr(yt, yp)
    return {"rmse": rmse, "mae": mae, "r2": r2, "spearman": float(rho), "n": int(mask.sum())}


def train_via_cli(train_path, val_path, test_path, output_dir, epochs=40, batch_size=64):
    output_dir.mkdir(parents=True, exist_ok=True)
    chemprop_bin = shutil.which("chemprop") or "chemprop"
    cmd = [
        chemprop_bin,
        "train",
        "-i",
        str(train_path),
        str(val_path),
        str(test_path),
        "-s",
        "smiles",
        "--target-columns",
        *TARGET_COLS,
        "--metrics",
        "rmse",
        "mae",
        "r2",
        "--epochs",
        str(epochs),
        "--batch-size",
        str(batch_size),
        "--accelerator",
        "cpu",
        "--num-workers",
        "0",
        "-o",
        str(output_dir),
    ]
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("Chemprop CLI stderr:\n%s", result.stderr[-3000:])
        raise RuntimeError(f"chemprop train failed (exit {result.returncode})")
    logger.info("Chemprop training completed")
    return result.stdout


def predict_via_cli(model_dir: Path, test_path: Path, pred_path: Path):
    chemprop_bin = shutil.which("chemprop") or "chemprop"
    # Chemprop saves checkpoints in output dir; use the directory or .pt files
    model_paths = list(model_dir.rglob("*.pt"))
    if not model_paths:
        model_paths = [model_dir]
    else:
        model_paths = [model_paths[0]]
    cmd = [
        chemprop_bin,
        "predict",
        "-i",
        str(test_path),
        "-s",
        "smiles",
        "--model-paths",
        *[str(p) for p in model_paths],
        "-o",
        str(pred_path),
    ]
    logger.info("Predicting: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("Chemprop predict stderr:\n%s", result.stderr[-3000:])
        raise RuntimeError(f"chemprop predict failed (exit {result.returncode})")


def find_model_dir(output_dir: Path) -> Path:
    for p in sorted(output_dir.rglob("model_0")):
        return p.parent
    checkpoints = list(output_dir.rglob("*.pt"))
    if checkpoints:
        return checkpoints[0].parent
    return output_dir


def evaluate_predictions(test_df: pd.DataFrame, pred_df: pd.DataFrame) -> dict:
    merged = test_df.merge(pred_df, on="smiles", suffixes=("_true", "_pred"), how="inner")
    metrics = {}
    for col in TARGET_COLS:
        true_col = col
        pred_col = f"{col}_pred" if f"{col}_pred" in merged.columns else col
        # chemprop may name predictions same as targets
        if pred_col not in merged.columns:
            for c in merged.columns:
                if col in c and c != col and "pred" in c.lower():
                    pred_col = c
                    break
        if true_col not in merged.columns:
            continue
        # After merge with suffixes
        tcol = true_col if true_col in merged.columns else f"{true_col}_true"
        pcol = None
        for candidate in [f"{col}_pred", col, f"{col}_pred_pred"]:
            if candidate in merged.columns:
                pcol = candidate
                break
        if pcol is None:
            # chemprop output columns often match target names exactly in separate file
            pred_only = pred_df.set_index("smiles")
            y_true = merged[tcol].values.astype(float)
            y_pred = pred_only.loc[merged["smiles"], col].values.astype(float)
        else:
            y_true = merged[tcol].values.astype(float)
            y_pred = merged[pcol].values.astype(float)
        metrics[col] = regression_metrics(y_true, y_pred)
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Train Chemprop 2.0 MTL model")
    parser.add_argument("--splits-dir", type=Path, default=ROOT / "data" / "processed" / "splits")
    parser.add_argument("--output", type=Path, default=ROOT / "models" / "chemprop")
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()

    train_path = args.splits_dir / "train.csv"
    val_path = args.splits_dir / "val.csv"
    test_path = args.splits_dir / "test.csv"
    for p in [train_path, val_path, test_path]:
        if not p.exists():
            raise SystemExit(f"Missing split file: {p}. Run 00_prepare_user_data.py first.")

    if args.output.exists():
        shutil.rmtree(args.output)
    train_via_cli(train_path, val_path, test_path, args.output, args.epochs, args.batch_size)

    model_dir = find_model_dir(args.output)
    pred_path = args.output / "test_predictions.csv"
    predict_via_cli(model_dir, test_path, pred_path)

    test_df = pd.read_csv(test_path)
    pred_df = pd.read_csv(pred_path)

    # Align predictions: chemprop output has smiles + target columns
    metrics = {}
    pred_indexed = pred_df.set_index("smiles")
    for col in TARGET_COLS:
        if col not in pred_indexed.columns:
            continue
        rows = test_df[test_df[col].notna()].copy()
        if len(rows) == 0:
            continue
        y_true = rows[col].values.astype(float)
        y_pred = pred_indexed.loc[rows["smiles"], col].values.astype(float)
        metrics[col] = regression_metrics(y_true, y_pred)

    report = {"model": "Chemprop 2.0 MTL", "test_metrics": metrics}
    with open(args.output / "chemprop_metrics.json", "w") as f:
        json.dump(report, f, indent=2)

    logger.info("Chemprop test metrics:")
    for k, v in metrics.items():
        logger.info("  %s: R²=%.3f, RMSE=%.3f, Spearman=%.3f (n=%d)", k, v["r2"], v["rmse"], v["spearman"], v["n"])


if __name__ == "__main__":
    main()
