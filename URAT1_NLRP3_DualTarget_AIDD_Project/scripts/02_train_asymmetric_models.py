#!/usr/bin/env python3
"""Train asymmetric dual-evidence models.

URAT1: XGBoost regression + split conformal UQ (SI / contrast only; not clinical ranking)
NLRP3: Assay-conditioned XGBoost classifier (library shrink, P≥0.5)
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
    roc_auc_binary,
    save_json,
    scaffold_cv_indices,
    scaffold_holdout_split,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
AUXILIARY = PROJECT_ROOT / "data" / "auxiliary"
RESULTS = PROJECT_ROOT / "results" / "training"

# URAT1: strict QSAR / virtual-screening criteria
# NOTE: EF@10% at p>=6 is misleading when base active rate ~57% (theoretical max EF ~1.75).
# Use strong-active enrichment (p>=7) and OOF-level metrics instead.
URAT1_THRESHOLDS = {
    "spearman_oof": 0.65,
    "r2_oof": 0.45,
    "rmse_oof": 0.70,
    "roc_auc_p7": 0.80,
    "ef_5pct_p7": 2.5,
    "fold_spearman_min": 0.45,  # worst fold must not collapse
}
NLRP3_THRESHOLDS = {
    "auroc": 0.70,
    "auprc": 0.65,
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


def load_oat_transfer_data(path: Path | None = None) -> pd.DataFrame:
    oat_path = path or (AUXILIARY / "oat_combined_transfer.csv")
    if not oat_path.exists():
        raise FileNotFoundError(
            f"OAT transfer file not found: {oat_path}. Run scripts/00b_prepare_auxiliary_data.py first."
        )
    df = pd.read_csv(oat_path)
    required = {"canonical_smiles", "pActivity"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"OAT transfer CSV missing columns: {sorted(missing)}")
    return df.dropna(subset=["canonical_smiles", "pActivity"]).reset_index(drop=True)


def _pretrain_xgb_on_oat(
    oat_df: pd.DataFrame,
    scaler: StandardScaler,
    exclude_smiles: set[str],
    seed: int,
) -> xgb.XGBRegressor | None:
    mask = ~oat_df["canonical_smiles"].isin(exclude_smiles)
    oat_fold = oat_df.loc[mask]
    if len(oat_fold) < 10:
        return None

    x_oat = scaler.transform(featurize_smiles(oat_fold["canonical_smiles"].tolist()))
    y_oat = oat_fold["pActivity"].values.astype(float)
    model = _xgb_regressor(seed=seed)
    model.fit(x_oat, y_oat)
    return model


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


def train_urat1_cv(
    df: pd.DataFrame,
    n_splits: int = 5,
    seed: int = 42,
    oat_df: pd.DataFrame | None = None,
    use_oat_transfer: bool = False,
) -> dict:
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

        test_smiles = {smiles[i] for i in te_idx}
        oat_pretrain = None
        if use_oat_transfer and oat_df is not None:
            oat_pretrain = _pretrain_xgb_on_oat(oat_df, scaler, exclude_smiles=test_smiles, seed=seed + fold_i)

        model = _xgb_regressor(seed=seed + fold_i)
        if oat_pretrain is not None:
            model.fit(x_fit, y[fit_idx], xgb_model=oat_pretrain.get_booster())
        else:
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
        m["ef_10pct_p6"] = regression_enrichment_factor(y_te, y_te_pred, threshold=6.0, fraction=0.1)
        m["ef_5pct_p7"] = regression_enrichment_factor(y_te, y_te_pred, threshold=7.0, fraction=0.05)
        m["roc_auc_p7"] = roc_auc_binary(y_te, y_te_pred, threshold=7.0)
        m["mean_interval_width"] = float(np.mean(hi - lo))
        m["oat_pretrain_used"] = oat_pretrain is not None
        fold_metrics.append(m)

    oof_metrics = regression_metrics(y, oof_pred)
    oof_metrics["ef_10pct_p6"] = regression_enrichment_factor(y, oof_pred, threshold=6.0, fraction=0.1)
    oof_metrics["ef_5pct_p7"] = regression_enrichment_factor(y, oof_pred, threshold=7.0, fraction=0.05)
    oof_metrics["roc_auc_p7"] = roc_auc_binary(y, oof_pred, threshold=7.0)
    oof_metrics["active_rate_p6"] = float((y >= 6).mean())
    oof_metrics["active_rate_p7"] = float((y >= 7).mean())
    oof_metrics["ef_p6_theoretical_max"] = float(1.0 / max((y >= 6).mean(), 1e-6))

    agg = {
        "rmse": oof_metrics["rmse"],
        "mae": oof_metrics["mae"],
        "r2": oof_metrics["r2"],
        "spearman": oof_metrics["spearman"],
        "ef_10pct_p6": oof_metrics["ef_10pct_p6"],
        "ef_5pct_p7": oof_metrics["ef_5pct_p7"],
        "roc_auc_p7": oof_metrics["roc_auc_p7"],
        "active_rate_p6": oof_metrics["active_rate_p6"],
        "active_rate_p7": oof_metrics["active_rate_p7"],
        "ef_p6_theoretical_max": oof_metrics["ef_p6_theoretical_max"],
        "mean_interval_width": float(np.mean(oof_hi - oof_lo)),
        "n_compounds": int(len(df)),
        "n_folds": n_splits,
        "oat_transfer": bool(use_oat_transfer and oat_df is not None),
        "oat_pretrain_compounds": int(len(oat_df)) if oat_df is not None else 0,
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


def fit_urat1_final(
    df: pd.DataFrame,
    seed: int = 42,
    oat_df: pd.DataFrame | None = None,
    use_oat_transfer: bool = False,
) -> dict:
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    x_mol = featurize_smiles(smiles)

    tr_idx, cal_idx, _ = scaffold_holdout_split(smiles, test_frac=0.1, val_frac=0.1, seed=seed)
    scaler = StandardScaler()
    x_tr = scaler.fit_transform(x_mol[tr_idx])
    x_cal = scaler.transform(x_mol[cal_idx])

    oat_pretrain = None
    if use_oat_transfer and oat_df is not None:
        oat_pretrain = _pretrain_xgb_on_oat(oat_df, scaler, exclude_smiles=set(), seed=seed)

    model = _xgb_regressor(seed=seed)
    if oat_pretrain is not None:
        model.fit(x_tr, y[tr_idx], xgb_model=oat_pretrain.get_booster())
    else:
        model.fit(x_tr, y[tr_idx])

    y_cal_pred = model.predict(x_cal)
    conformal = SplitConformalRegressor(alpha=0.1).fit(y[cal_idx], y_cal_pred)

    return {
        "model": model,
        "scaler": scaler,
        "conformal": conformal,
        "feature_type": "morgan2048+rdkit",
        "oat_transfer": bool(oat_pretrain is not None),
    }


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
        fold_spearmans = [f["spearman"] for f in metrics.get("fold_metrics", [])]
        worst_fold_spearman = min(fold_spearmans) if fold_spearmans else float("nan")
        checks = {
            "spearman_oof_ok": metrics["spearman"] >= URAT1_THRESHOLDS["spearman_oof"],
            "r2_oof_ok": metrics["r2"] >= URAT1_THRESHOLDS["r2_oof"],
            "rmse_oof_ok": metrics["rmse"] <= URAT1_THRESHOLDS["rmse_oof"],
            "roc_auc_p7_ok": metrics["roc_auc_p7"] >= URAT1_THRESHOLDS["roc_auc_p7"],
            "ef_5pct_p7_ok": metrics["ef_5pct_p7"] >= URAT1_THRESHOLDS["ef_5pct_p7"],
            "fold_stability_ok": worst_fold_spearman >= URAT1_THRESHOLDS["fold_spearman_min"],
        }
        passed = sum(checks.values())
        suitable = passed >= 5  # require 5/6 strict checks
        note = (
            "EF@10% at p>=6 is NOT used: with ~57% actives the theoretical maximum EF is only ~1.75, "
            "so values near 1.8 do not indicate good virtual screening performance."
        )
    else:
        checks = {
            "auroc_ok": metrics["auroc"] >= NLRP3_THRESHOLDS["auroc"],
            "auprc_ok": metrics["auprc"] >= NLRP3_THRESHOLDS["auprc"],
            "ef_ok": metrics["ef_10pct"] >= NLRP3_THRESHOLDS["ef_10pct"],
        }
        passed = sum(checks.values())
        suitable = passed >= 2
        note = None

    return {
        "suitable_for_screening": bool(suitable),
        "checks": checks,
        "criteria_passed": int(passed),
        "criteria_total": len(checks),
        "thresholds": URAT1_THRESHOLDS if target == "urat1" else NLRP3_THRESHOLDS,
        "metric_note": note,
        "recommendation": (
            "Passes strict OOF scaffold-CV criteria; still requires benchmark backtest before library screening."
            if suitable
            else "Does NOT meet strict screening criteria. Do not use as primary ML filter; rely on structure-based scoring."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=PROCESSED)
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--oat-transfer",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="OAT1/OAT3 sequential pretrain before URAT1 fine-tune (default: enabled)",
    )
    parser.add_argument("--oat-csv", type=Path, default=AUXILIARY / "oat_combined_transfer.csv")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    urat1 = pd.read_csv(args.data_dir / "urat1_curated.csv")
    nlrp3 = pd.read_csv(args.data_dir / "nlrp3_records.csv")

    oat_df = None
    if args.oat_transfer:
        oat_df = load_oat_transfer_data(args.oat_csv)
        print(f"OAT transfer enabled: {len(oat_df)} auxiliary compounds from {args.oat_csv.name}")
    else:
        print("OAT transfer disabled (URAT1 baseline from scratch)")

    print("Training URAT1 regression model (scaffold CV + conformal UQ)...")
    urat1_cv = train_urat1_cv(
        urat1,
        n_splits=args.n_splits,
        seed=args.seed,
        oat_df=oat_df,
        use_oat_transfer=args.oat_transfer,
    )
    urat1_suit = assess_screening_suitability("urat1", urat1_cv)
    urat1_final = fit_urat1_final(
        urat1,
        seed=args.seed,
        oat_df=oat_df,
        use_oat_transfer=args.oat_transfer,
    )

    print("Training NLRP3 assay-conditioned classifier (scaffold CV)...")
    nlrp3_cv = train_nlrp3_cv(nlrp3, n_splits=args.n_splits, seed=args.seed)
    nlrp3_suit = assess_screening_suitability("nlrp3", nlrp3_cv)
    nlrp3_final = fit_nlrp3_final(nlrp3, nlrp3_cv["top_assay_ids"], seed=args.seed)

    joblib.dump(urat1_final, args.output / "urat1_model.joblib")
    joblib.dump(nlrp3_final, args.output / "nlrp3_model.joblib")

    pd.DataFrame(urat1_cv["oof_predictions"]).to_csv(args.output / "urat1_oof_predictions.csv", index=False)
    pd.DataFrame(nlrp3_cv["oof_predictions"]).to_csv(args.output / "nlrp3_oof_predictions.csv", index=False)

    report = {
        "framework": "asymmetric_dual_evidence",
        "urat1": {
            "model": "XGBoost regression + split conformal (alpha=0.1)",
            "transfer_learning": {
                "enabled": bool(args.oat_transfer),
                "method": "sequential_finetune",
                "auxiliary": "slc22_oat_transfer",
                "oat_pretrain_compounds": int(len(oat_df)) if oat_df is not None else 0,
                "final_model_oat_transfer": bool(urat1_final.get("oat_transfer", False)),
            },
            "cv_metrics": {
                k: urat1_cv[k]
                for k in [
                    "rmse", "mae", "r2", "spearman",
                    "ef_10pct_p6", "ef_5pct_p7", "roc_auc_p7",
                    "active_rate_p6", "active_rate_p7", "ef_p6_theoretical_max",
                    "mean_interval_width", "n_compounds",
                    "oat_transfer", "oat_pretrain_compounds",
                ]
            },
            "fold_metrics": urat1_cv["fold_metrics"],
            "screening_assessment": urat1_suit,
        },
        "nlrp3": {
            "model": "Assay-conditioned XGBoost classifier (top-25 assays one-hot)",
            "cv_metrics": {k: nlrp3_cv[k] for k in ["auroc", "auprc", "ef_10pct", "n_compounds", "n_records", "n_assays_used"]},
            "fold_metrics": nlrp3_cv["fold_metrics"],
            "screening_assessment": nlrp3_suit,
        },
        "overall_screening_ready": bool(
            urat1_suit["suitable_for_screening"]
            and nlrp3_suit["suitable_for_screening"]
        ),
        "warning": (
            "CV pass alone is insufficient for URAT1. Run 07_benchmark_backtest.py; "
            "URAT1 requires benchmark recovery before library screening."
        ),
    }

    # remove nested oof from json (saved as CSV)
    save_json(args.output / "training_report.json", report)
    docs_dir = PROJECT_ROOT / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    save_json(docs_dir / "MODEL_TRAINING_SUMMARY.json", report)

    print("\n=== URAT1 CV (OOF, scaffold split) ===")
    print(f"  RMSE={urat1_cv['rmse']:.3f}  R2={urat1_cv['r2']:.3f}  Spearman={urat1_cv['spearman']:.3f}")
    print(f"  ROC-AUC(p>=7)={urat1_cv['roc_auc_p7']:.3f}  EF@5%(p>=7)={urat1_cv['ef_5pct_p7']:.2f}")
    print(f"  EF@10%(p>=6)={urat1_cv['ef_10pct_p6']:.2f}  [misleading: theoretical max={urat1_cv['ef_p6_theoretical_max']:.2f} at 57% actives]")
    print(f"  OAT transfer: {urat1_cv.get('oat_transfer', False)}")
    print(f"  Strict screening suitable: {urat1_suit['suitable_for_screening']} ({urat1_suit['criteria_passed']}/{urat1_suit['criteria_total']} checks)")

    print("\n=== NLRP3 CV ===")
    print(f"  AUROC={nlrp3_cv['auroc']:.3f}  AUPRC={nlrp3_cv['auprc']:.3f}  EF@10%={nlrp3_cv['ef_10pct']:.2f}")
    print(f"  Strict screening suitable: {nlrp3_suit['suitable_for_screening']} ({nlrp3_suit['criteria_passed']}/{nlrp3_suit['criteria_total']} checks)")

    print(f"\nOverall ML screening ready: {report['overall_screening_ready']}")
    print(f"Models saved to {args.output}")


if __name__ == "__main__":
    main()
