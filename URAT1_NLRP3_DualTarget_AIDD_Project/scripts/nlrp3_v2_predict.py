#!/usr/bin/env python3
"""Score molecules with an NLRP3 v2 bundle produced by nlrp3_model_v2_validate.py.

This is the prediction entry point for the next stage (re-scoring the clinical
library). It is deliberately kept separate from the historical
scripts/screen_repurposing_library.py so that the frozen 1588-molecule pool and
everything downstream of it stay untouched.

    python3 scripts/nlrp3_v2_predict.py \
        --input data/repurposing/repurposing_manifest.csv \
        --model results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib \
        --calibration none --threshold 0.5 \
        --output results/nlrp3_model_v2/predictions/<name>.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from utils_ml import canonicalize, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def score(bundle: dict, smiles: list[str], calibration: str) -> np.ndarray:
    raw = bundle["model"].predict_molecules(smiles)
    if calibration == "none":
        return raw
    cal = bundle["calibrators"].get(calibration)
    if cal is None:
        raise ValueError(
            f"calibration '{calibration}' not available in bundle "
            f"(available: {sorted(bundle['calibrators'])})"
        )
    return cal.predict(raw)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True, help="CSV with a SMILES column")
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--smiles-column", default=None)
    ap.add_argument("--calibration", default="none", choices=["none", "platt", "isotonic"])
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    df = pd.read_csv(args.input, low_memory=False)
    col = args.smiles_column or ("canonical_smiles" if "canonical_smiles" in df.columns else "smiles")
    df = df.dropna(subset=[col]).copy()
    df["canonical_smiles"] = df[col].map(canonicalize)
    df = df[df["canonical_smiles"].notna()].drop_duplicates("canonical_smiles").reset_index(drop=True)

    bundle = joblib.load(args.model)
    p = score(bundle, df["canonical_smiles"].tolist(), args.calibration)

    df["p_active_nlrp3_v2"] = p
    df["nlrp3_v2_pred_active"] = df["p_active_nlrp3_v2"] >= args.threshold
    df["nlrp3_v2_percentile"] = df["p_active_nlrp3_v2"].rank(pct=True) * 100.0
    df["scaffold"] = df["canonical_smiles"].map(murcko_scaffold)
    df = df.sort_values("p_active_nlrp3_v2", ascending=False).reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    summary = {
        "input": str(args.input),
        "model": str(args.model),
        "prediction_rule": bundle["prediction_rule"],
        "calibration": args.calibration,
        "threshold": args.threshold,
        "n_scored": int(len(df)),
        "n_pred_active": int(df["nlrp3_v2_pred_active"].sum()),
        "output": str(args.output),
    }
    with open(args.output.with_suffix(".summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
