#!/usr/bin/env python3
"""
TAPE-GATE Stage 2: Train asymmetric dual-evidence models.

URAT1: XGBoost regression + split conformal UQ (5-fold scaffold CV)
NLRP3: Assay-conditioned XGBoost classifier (5-fold scaffold CV by molecule)

Outputs models, CV metrics, and screening suitability verdict under results/training/.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler

from utils_ml import (
    SplitConformalRegressor,
    assay_one_hot_matrix,
    classification_metrics,
    enrichment_factor,
    featurize_smiles,
    murcko_scaffold,
    regression_enrichment_factor,
    regression_metrics,
    save_json,
    scaffold_cv_indices,
    scaffold_holdout_split,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results" / "training"

# Screening suitability thresholds (conservative but achievable for ~500-800 compound QSAR)
URAT1_THRESHOLDS = {
    "spearman": 0.50,
    "r2": 0.25,
    "ef_10pct": 1.5,
    "rmse": 1.0,
}
NLRP3_THRESHOLDS = {
    "auroc": 0.65,
    "auprc": 0.60,
    "ef_10pct": 1.5,
}


def _xgb_regressor(seed: int = 42) -> xgb.XGBRegressor:
    return xgb.XGBRegressor(
        n_estimators=400,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=seed,
        n_jobs=-1,
        objective="reg:squarederror",
    )


def _xgb_classifier(seed: int = 42) -> xgb.XGBClassifier:
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


def train_urat1_cv(df: pd.DataFrame, n_splits: int = 5, seed: int = 42) -> dict:
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    x_mol = featurize_smiles(smiles)

    fold_metrics = []
    oof_pred = np.full(len(y), np.nan)
    oof_lo = np.full(len(y), np.nan)
    oof_hi = np.full(len(y), np.nan)

    for fold_i, (tr_idx, te_idx) in enumerate(scaffold_cv_indices(smiles, n_splits=n_splits)):
        x_tr, x_te = x_mol[tr_idx], x_mol[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]

        # calibration split inside training fold (15% scaffolds)
        n_tr = len(tr_idx)
        if n_tr < 20:
            fit_idx, cal_idx = tr_idx, tr_idx
        else:
            rel_tr = np.arange(n_tr)
            sub_groups = [murcko_scaffold(smiles[i]) for i in tr_idx]
            gss_cal = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=seed + fold_i)
            fit_rel, cal_rel = next(gss_cal.split(rel_tr, groups=sub_groups))
            fit_idx = tr_idx[fit_rel]
            cal_idx = tr_idx[cal_rel]

        scaler = StandardScaler()
        x_fit = scaler.fit_transform(x_mol[fit_idx])
        x_cal = scaler.transform(x_mol[cal_idx])
        x_te_s = scaler.transform(x_te)

        model = _xgb_regressor(seed=seed + fold_i)
        model.fit(x_fit, y[fit_idx])

        y_cal_pred = model.predict(x_cal)
        conformal = SplitConformalRegressor(alpha=0.1).fit(y[cal_idx], y_cal_pred)

        y_te_pred = model.predict(x_te_s)
        lo, hi = conformal.interval(y_te_pred)
        oof_pred[te_idx] = y_te_pred
        oof_lo[te_idx] = lo
        oof_hi[te_idx] = hi

        m = regression_metrics(y_te, y_te_pred)
        m["fold"] = fold_i
        m["ef_10pct"] = regression_enrichment_factor(y_te, y_te_pred, threshold=6.0, fraction=0.1)
        m["mean_interval_width"] = float(np.mean(hi - lo))
        fold_metrics.append(m)

    agg = {
        "rmse": float(np.mean([m["rmse"] for m in fold_metrics])),
        "mae": float(np.mean([m["mae"] for m in fold_metrics])),
        "r2": float(np.mean([m["r2"] for m in fold_metrics])),
        "spearman": float(np.mean([m["spearman"] for m in fold_metrics])),
        "ef_10pct": float(np.mean([m["ef_10pct"] for m in fold_metrics])),
        "mean_interval_width": float(np.mean([m["mean_interval_width"] for m in fold_metrics])),
        "n_compounds": int(len(df)),
        "n_folds": n_splits,
        "fold_metrics": fold_metrics,
        "oof_predictions": {
            "canonical_smiles": smiles,
            "y_true": y.tolist(),
            "y_pred": oof_pred.tolist(),
            "y_lo": oof_lo.tolist(),
            "y_hi": oof_hi.tolist(),
        },
    }
    return agg


def fit_urat1_final(df: pd.DataFrame, seed: int = 42) -> dict:
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    x_mol = featurize_smiles(smiles)

    tr_idx, cal_idx, _ = scaffold_holdout_split(smiles, test_frac=0.1, val_frac=0.1, seed=seed)
    scaler = StandardScaler()
    x_tr = scaler.fit_transform(x_mol[tr_idx])
    x_cal = scaler.transform(x_mol[cal_idx])

    model = _xgb_regressor(seed=seed)
    model.fit(x_tr, y[tr_idx])

    y_cal_pred = model.predict(x_cal)
    conformal = SplitConformalRegressor(alpha=0.1).fit(y[cal_idx], y_cal_pred)

    return {"model": model, "scaler": scaler, "conformal": conformal, "feature_type": "morgan2048+rdkit"}


def train_nlrp3_cv(df: pd.DataFrame, n_splits: int = 5, seed: int = 42, top_assays: int = 25) -> dict:
  # group CV by molecule scaffold (all assay records for a molecule stay in same fold)
    mol_df = (
        df.groupby("canonical_smiles", as_index=False)
        .agg(active=("active", "max"), scaffold=("scaffold", "first"))
        .reset_index(drop=True)
    )
    smiles = mol_df["canonical_smiles"].tolist()
    y_mol = mol_df["active"].values.astype(float)

    assay_counts = df["Assay ChEMBL ID"].value_counts()
    top_assay_ids = assay_counts.head(top_assays).index.astype(str).tolist()

    x_mol = featurize_smiles(smiles)
    fold_metrics = []
    oof_prob = np.full(len(smiles), np.nan)

    for fold_i, (tr_mol_idx, te_mol_idx) in enumerate(scaffold_cv_indices(smiles, n_splits=n_splits)):
        tr_smiles = set(mol_df.iloc[tr_mol_idx]["canonical_smiles"])
        te_smiles = set(mol_df.iloc[te_mol_idx]["canonical_smiles"])

        tr_records = df[df["canonical_smiles"].isin(tr_smiles)].copy()
        te_records = df[df["canonical_smiles"].isin(te_smiles)].copy()

        x_tr_mol = featurize_smiles(tr_records["canonical_smiles"].tolist())
        x_tr_assay = assay_one_hot_matrix(tr_records["Assay ChEMBL ID"], top_assay_ids)
        x_tr = np.hstack([x_tr_mol, x_tr_assay])

        x_te_mol = featurize_smiles(te_records["canonical_smiles"].tolist())
        x_te_assay = assay_one_hot_matrix(te_records["Assay ChEMBL ID"], top_assay_ids)
        x_te = np.hstack([x_te_mol, x_te_assay])

        y_tr = tr_records["active"].values.astype(int)
        y_te_rec = te_records["active"].values.astype(int)

        assay_w = tr_records["Assay ChEMBL ID"].value_counts().to_dict()
        sample_w = tr_records["Assay ChEMBL ID"].map(lambda a: 1.0 / np.sqrt(assay_w.get(a, 1))).values

        scaler = StandardScaler()
        x_tr_s = scaler.fit_transform(x_tr)

        model = _xgb_classifier(seed=seed + fold_i)
        model.fit(x_tr_s, y_tr, sample_weight=sample_w)

        x_te_s = scaler.transform(x_te)
        prob_te_rec = model.predict_proba(x_te_s)[:, 1]

        # molecule-level OOF: max prob across assay records
        te_records = te_records.copy()
        te_records["prob"] = prob_te_rec
        mol_prob = te_records.groupby("canonical_smiles")["prob"].max()
        for j, smi in enumerate(mol_df.iloc[te_mol_idx]["canonical_smiles"]):
            idx = te_mol_idx[j]
            oof_prob[idx] = mol_prob.get(smi, np.nan)

        y_te_mol = mol_df.iloc[te_mol_idx]["active"].values.astype(int)
        prob_te_mol = mol_prob.reindex(mol_df.iloc[te_mol_idx]["canonical_smiles"]).values

        m = classification_metrics(y_te_mol, prob_te_mol)
        m["fold"] = fold_i
        m["ef_10pct"] = enrichment_factor(y_te_mol, prob_te_mol, fraction=0.1)
        fold_metrics.append(m)

    agg = {
        "auroc": float(np.mean([m["auroc"] for m in fold_metrics])),
        "auprc": float(np.mean([m["auprc"] for m in fold_metrics])),
        "ef_10pct": float(np.mean([m["ef_10pct"] for m in fold_metrics])),
        "n_compounds": int(mol_df["canonical_smiles"].nunique()),
        "n_records": int(len(df)),
        "n_assays_used": int(len(top_assay_ids)),
        "n_folds": n_splits,
        "fold_metrics": fold_metrics,
        "top_assay_ids": top_assay_ids,
        "oof_predictions": {
            "canonical_smiles": smiles,
            "y_true": y_mol.tolist(),
            "y_prob": oof_prob.tolist(),
        },
    }
    return agg


def fit_nlrp3_final(df: pd.DataFrame, top_assay_ids: list[str], seed: int = 42) -> dict:
    smiles = df["canonical_smiles"].unique().tolist()
    tr_idx, cal_idx, _ = scaffold_holdout_split(smiles, test_frac=0.1, val_frac=0.1, seed=seed)
    cal_smiles = set(np.array(smiles)[cal_idx])

    x_mol = featurize_smiles(df["canonical_smiles"].tolist())
    x_assay = assay_one_hot_matrix(df["Assay ChEMBL ID"], top_assay_ids)
    x = np.hstack([x_mol, x_assay])
    y = df["active"].values.astype(int)

    assay_w = df["Assay ChEMBL ID"].value_counts().to_dict()
    sample_w = df["Assay ChEMBL ID"].map(lambda a: 1.0 / np.sqrt(assay_w.get(a, 1))).values

    scaler = StandardScaler()
    x_s = scaler.fit_transform(x)

    base = _xgb_classifier(seed=seed)
    base.fit(x_s, y, sample_weight=sample_w)

  # manual isotonic calibration on held-out scaffolds
    cal_mask = df["canonical_smiles"].isin(cal_smiles).values
    calibrator = None
    if cal_mask.sum() >= 30 and len(np.unique(y[cal_mask])) == 2:
        raw_prob = base.predict_proba(x_s[cal_mask])[:, 1]
        calibrator = IsotonicRegression(out_of_bounds="clip")
        calibrator.fit(raw_prob, y[cal_mask])

    return {
        "model": base,
        "calibrator": calibrator,
        "scaler": scaler,
        "top_assay_ids": top_assay_ids,
        "feature_type": "morgan2048+rdkit+assay_onehot",
        "calibrated": calibrator is not None,
    }


def assess_screening_suitability(target: str, metrics: dict) -> dict:
    if target == "urat1":
        checks = {
            "spearman_ok": metrics["spearman"] >= URAT1_THRESHOLDS["spearman"],
            "r2_ok": metrics["r2"] >= URAT1_THRESHOLDS["r2"],
            "ef_ok": metrics["ef_10pct"] >= URAT1_THRESHOLDS["ef_10pct"],
            "rmse_ok": metrics["rmse"] <= URAT1_THRESHOLDS["rmse"],
        }
        passed = sum(checks.values())
        suitable = passed >= 3  # at least 3/4 criteria
    else:
        checks = {
            "auroc_ok": metrics["auroc"] >= NLRP3_THRESHOLDS["auroc"],
            "auprc_ok": metrics["auprc"] >= NLRP3_THRESHOLDS["auprc"],
            "ef_ok": metrics["ef_10pct"] >= NLRP3_THRESHOLDS["ef_10pct"],
        }
        passed = sum(checks.values())
        suitable = passed >= 2  # at least 2/3 criteria

    return {
        "suitable_for_screening": bool(suitable),
        "checks": checks,
        "criteria_passed": int(passed),
        "thresholds": URAT1_THRESHOLDS if target == "urat1" else NLRP3_THRESHOLDS,
        "recommendation": (
            "Model meets minimum CV criteria; usable as primary ML filter before docking."
            if suitable
            else "Model below screening thresholds; use with caution or rely more on structure evidence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=PROCESSED)
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    urat1 = pd.read_csv(args.data_dir / "urat1_curated.csv")
    nlrp3 = pd.read_csv(args.data_dir / "nlrp3_records.csv")

    print("Training URAT1 regression model (scaffold CV + conformal UQ)...")
    urat1_cv = train_urat1_cv(urat1, n_splits=args.n_splits, seed=args.seed)
    urat1_suit = assess_screening_suitability("urat1", urat1_cv)
    urat1_final = fit_urat1_final(urat1, seed=args.seed)

    print("Training NLRP3 assay-conditioned classifier (scaffold CV)...")
    nlrp3_cv = train_nlrp3_cv(nlrp3, n_splits=args.n_splits, seed=args.seed)
    nlrp3_suit = assess_screening_suitability("nlrp3", nlrp3_cv)
    nlrp3_final = fit_nlrp3_final(nlrp3, nlrp3_cv["top_assay_ids"], seed=args.seed)

    joblib.dump(urat1_final, args.output / "urat1_model.joblib")
    joblib.dump(nlrp3_final, args.output / "nlrp3_model.joblib")

    pd.DataFrame(urat1_cv["oof_predictions"]).to_csv(args.output / "urat1_oof_predictions.csv", index=False)
    pd.DataFrame(nlrp3_cv["oof_predictions"]).to_csv(args.output / "nlrp3_oof_predictions.csv", index=False)

    report = {
        "framework": "TAPE-GATE",
        "urat1": {
            "model": "XGBoost regression + split conformal (alpha=0.1)",
            "cv_metrics": {k: urat1_cv[k] for k in ["rmse", "mae", "r2", "spearman", "ef_10pct", "mean_interval_width", "n_compounds"]},
            "fold_metrics": urat1_cv["fold_metrics"],
            "screening_assessment": urat1_suit,
        },
        "nlrp3": {
            "model": "Assay-conditioned XGBoost classifier (top-25 assays one-hot)",
            "cv_metrics": {k: nlrp3_cv[k] for k in ["auroc", "auprc", "ef_10pct", "n_compounds", "n_records", "n_assays_used"]},
            "fold_metrics": nlrp3_cv["fold_metrics"],
            "screening_assessment": nlrp3_suit,
        },
        "overall_screening_ready": bool(urat1_suit["suitable_for_screening"] and nlrp3_suit["suitable_for_screening"]),
    }

    # remove nested oof from json (saved as CSV)
    save_json(args.output / "training_report.json", report)

    print("\n=== URAT1 CV ===")
    print(f"  RMSE={urat1_cv['rmse']:.3f}  R2={urat1_cv['r2']:.3f}  Spearman={urat1_cv['spearman']:.3f}  EF@10%={urat1_cv['ef_10pct']:.2f}")
    print(f"  Screening suitable: {urat1_suit['suitable_for_screening']} ({urat1_suit['recommendation']})")

    print("\n=== NLRP3 CV ===")
    print(f"  AUROC={nlrp3_cv['auroc']:.3f}  AUPRC={nlrp3_cv['auprc']:.3f}  EF@10%={nlrp3_cv['ef_10pct']:.2f}")
    print(f"  Screening suitable: {nlrp3_suit['suitable_for_screening']} ({nlrp3_suit['recommendation']})")

    print(f"\nOverall ML screening ready: {report['overall_screening_ready']}")
    print(f"Models saved to {args.output}")


if __name__ == "__main__":
    main()
