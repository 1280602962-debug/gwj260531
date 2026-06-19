#!/usr/bin/env python3
"""
Train Chemprop 2.0 single-target models for JNK1, JNK2, JNK3.

Uses per-isoform scaffold splits from data/processed/splits/{jnk1,jnk2,jnk3}/.

Usage:
    python3 scripts/04b_train_chemprop_mtl.py
    python3 scripts/04b_train_chemprop_mtl.py --isoform JNK1 --epochs 100
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from utils_ml import regression_metrics  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ISOFORMS = ["JNK1", "JNK2", "JNK3"]


def load_config() -> dict:
    with open(ROOT / "config" / "targets.yaml") as f:
        return yaml.safe_load(f)


def train_chemprop(train_path, val_path, test_path, target_col, output_dir, cfg):
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
        target_col,
        "--metrics",
        "rmse",
        "mae",
        "r2",
        "--epochs",
        str(cfg.get("epochs", 100)),
        "--batch-size",
        str(cfg.get("batch_size", 50)),
        "--patience",
        str(cfg.get("patience", 20)),
        "--ffn-hidden-dim",
        str(cfg.get("ffn_hidden_dim", 500)),
        "--depth",
        str(cfg.get("depth", 4)),
        "--accelerator",
        "cpu",
        "--num-workers",
        "0",
        "-o",
        str(output_dir),
    ]
    logger.info("Chemprop train [%s]: %s", target_col, " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(result.stderr[-4000:])
        raise RuntimeError(f"chemprop train failed for {target_col}")
    return output_dir


def predict_chemprop(model_pt, test_path, pred_path):
    chemprop_bin = shutil.which("chemprop") or "chemprop"
    cmd = [
        chemprop_bin,
        "predict",
        "-i",
        str(test_path),
        "-s",
        "smiles",
        "--model-paths",
        str(model_pt),
        "-o",
        str(pred_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def find_best_pt(output_dir: Path) -> Path:
    pts = list(output_dir.rglob("best.pt"))
    if not pts:
        pts = list(output_dir.rglob("*.pt"))
    if not pts:
        raise FileNotFoundError(f"No checkpoint in {output_dir}")
    return pts[0]


def train_isoform(isoform: str, splits_dir: Path, output_root: Path, cfg: dict) -> dict:
    iso_dir = splits_dir / isoform.lower()
    target_col = f"pAct_{isoform}"
    out_dir = output_root / isoform.lower()
    if out_dir.exists():
        shutil.rmtree(out_dir)

    train_chemprop(
        iso_dir / "train.csv",
        iso_dir / "val.csv",
        iso_dir / "test.csv",
        target_col,
        out_dir,
        cfg,
    )
    model_pt = find_best_pt(out_dir)
    pred_path = out_dir / "test_predictions.csv"
    predict_chemprop(model_pt, iso_dir / "test.csv", pred_path)

    test_df = pd.read_csv(iso_dir / "test.csv")
    pred_df = pd.read_csv(pred_path)

    merged = test_df.merge(pred_df, on="smiles", suffixes=("_true", "_pred"))
    true_col = f"{target_col}_true" if f"{target_col}_true" in merged.columns else target_col
    pred_col = f"{target_col}_pred" if f"{target_col}_pred" in merged.columns else target_col
    if true_col not in merged.columns or pred_col not in merged.columns:
        # fallback: second pAct column after merge
        pcols = [c for c in merged.columns if "pAct" in c]
        true_col, pred_col = pcols[0], pcols[1]

    y_true = merged[true_col].values.astype(float)
    y_pred = merged[pred_col].values.astype(float)
    metrics = regression_metrics(y_true, y_pred)
    metrics_path = out_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump({"isoform": isoform, "holdout_test": metrics}, f, indent=2)
    logger.info("%s Chemprop holdout R²=%.3f Spearman=%.3f", isoform, metrics["r2"], metrics["spearman"])
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Train Chemprop per isoform")
    parser.add_argument("--splits-dir", type=Path, default=ROOT / "data" / "processed" / "splits")
    parser.add_argument("--output", type=Path, default=ROOT / "models" / "chemprop")
    parser.add_argument("--isoform", choices=ISOFORMS, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    args = parser.parse_args()

    config = load_config()
    cp_cfg = dict(config.get("training", {}).get("chemprop", {}))
    if args.epochs:
        cp_cfg["epochs"] = args.epochs

    isoforms = [args.isoform] if args.isoform else ISOFORMS
    all_metrics = {}
    for iso in isoforms:
        all_metrics[iso] = train_isoform(iso, args.splits_dir, args.output, cp_cfg)

    with open(args.output / "chemprop_metrics.json", "w") as f:
        json.dump({"model": "Chemprop 2.0 single-target", "holdout_test": all_metrics}, f, indent=2)


if __name__ == "__main__":
    main()
