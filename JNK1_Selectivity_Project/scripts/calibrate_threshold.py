#!/usr/bin/env python3
"""
Benchmark-calibrated F1 threshold for JNK family ML pre-filter.

Usage:
    python3 scripts/calibrate_threshold.py \
        --benchmarks data/benchmarks/literature_benchmarks.csv \
        --output results/calibration
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from utils_ml import featurize_smiles  # noqa: E402


def predict_family(smiles: list[str], model_dir: Path) -> np.ndarray:
    X = featurize_smiles(smiles)
    preds = []
    for iso in ["JNK1", "JNK2", "JNK3"]:
        path = model_dir / f"xgboost_{iso.lower()}.joblib"
        if not path.exists():
            raise FileNotFoundError(f"Missing model: {path}. Run scripts/07_compare_models.py first.")
        model = joblib.load(path)
        preds.append(model.predict(X))
    return np.nanmax(np.vstack(preds), axis=0)


def main() -> None:
    parser = argparse.ArgumentParser(description="Calibrate F1 p_family threshold using benchmarks")
    parser.add_argument("--benchmarks", type=Path, default=ROOT / "data" / "benchmarks" / "literature_benchmarks.csv")
    parser.add_argument("--models", type=Path, default=ROOT / "models" / "xgboost")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "calibration")
    parser.add_argument("--thresholds", type=str, default="5.0,5.5,6.0,6.5,7.0")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.benchmarks)
    df["must_pass_F1"] = df["must_pass_F1"].astype(str).str.lower().eq("yes")
    p_family = predict_family(df["smiles"].tolist(), args.models)
    df["p_family_pred"] = p_family

    thresholds = [float(x) for x in args.thresholds.split(",")]
    scan_rows = []
    must = df[df["must_pass_F1"]]
    for t in thresholds:
        passed = df[p_family >= t]
        recall = float((must["p_family_pred"] >= t).mean()) if len(must) else np.nan
        scan_rows.append(
            {
                "threshold": t,
                "n_passed_benchmarks": int((p_family >= t).sum()),
                "n_total_benchmarks": len(df),
                "must_pass_recall": recall,
                "n_must_pass_total": int(len(must)),
            }
        )

    scan = pd.DataFrame(scan_rows)
    scan.to_csv(args.output / "threshold_scan.csv", index=False)

    # pick: recall >= 0.9 with lowest threshold; else max recall
    ok = scan[scan["must_pass_recall"] >= 0.9]
    if len(ok):
        best = ok.sort_values("threshold").iloc[0]
    else:
        best = scan.sort_values(["must_pass_recall", "threshold"], ascending=[False, True]).iloc[0]

    df.to_csv(args.output / "benchmark_predictions.csv", index=False)
    rec = {
        "recommended_threshold": float(best["threshold"]),
        "must_pass_recall": float(best["must_pass_recall"]),
        "note": "p_family = max(pred_JNK1, pred_JNK2, pred_JNK3)",
    }
    with open(args.output / "threshold_recommendation.json", "w") as f:
        json.dump(rec, f, indent=2)

    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
