#!/usr/bin/env python3
"""
Non-docking computational module A — ML rigor / trustworthiness battery.

Adds the validation evidence a QSAR reviewer expects, WITHOUT retraining or
overwriting the production models / screening scores. It recomputes scaffold-CV
out-of-fold predictions purely for validation and reports:

  1. Y-scrambling (label-permutation) test  -> proves signal is not fingerprint
     artifact; empirical p-value of the real OOF metric vs permuted null.
  2. Applicability domain (AD)              -> nearest-neighbour Tanimoto to the
     training set for training compounds, benchmarks and the Pareto shortlist;
     flags out-of-domain predictions.
  3. NLRP3 probability calibration          -> reliability curve + Brier score on
     OOF predictions.

Targets:
  URAT1  : XGBoost regression (Morgan+RDKit), scaffold 5-fold CV, Spearman metric.
  NLRP3  : XGBoost classifier (Morgan+RDKit), molecule-grouped scaffold CV, AUROC.

Inputs (read-only):
  data/processed/urat1_curated.csv
  data/processed/nlrp3_records.csv
  data/benchmarks/literature_benchmarks.csv
  data/repurposing/pareto/pareto_shortlist.csv

Outputs:
  results/model_validation/yscramble_urat1.json
  results/model_validation/yscramble_nlrp3.json
  results/model_validation/applicability_domain.csv
  results/model_validation/nlrp3_calibration.csv
  results/model_validation/ml_rigor_summary.json

Usage:
  python3 scripts/12_ml_rigor_validation.py --n-permutations 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

from utils_ml import (
    classification_metrics,
    featurize_smiles,
    regression_metrics,
    scaffold_cv_indices,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
BENCHMARKS = PROJECT_ROOT / "data" / "benchmarks"
PARETO_DIR = PROJECT_ROOT / "data" / "repurposing" / "pareto"
OUT_DIR = PROJECT_ROOT / "results" / "model_validation"

SMILES_COL = "canonical_smiles"


def _xgb_reg(seed: int = 42) -> xgb.XGBRegressor:
    return xgb.XGBRegressor(
        n_estimators=400, max_depth=5, learning_rate=0.05, subsample=0.8,
        colsample_bytree=0.8, reg_lambda=1.0, random_state=seed, n_jobs=-1,
        objective="reg:squarederror",
    )


def _xgb_clf(seed: int = 42) -> xgb.XGBClassifier:
    return xgb.XGBClassifier(
        n_estimators=400, max_depth=4, learning_rate=0.05, subsample=0.8,
        colsample_bytree=0.8, reg_lambda=1.0, random_state=seed, n_jobs=-1,
        objective="binary:logistic", eval_metric="logloss",
    )


def oof_regression(X, y, folds) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for tr, te in folds:
        m = _xgb_reg()
        m.fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return pred


def oof_classification(X, y, folds) -> np.ndarray:
    prob = np.full(len(y), np.nan)
    for tr, te in folds:
        if len(np.unique(y[tr])) < 2:
            continue
        m = _xgb_clf()
        m.fit(X[tr], y[tr])
        prob[te] = m.predict_proba(X[te])[:, 1]
    return prob


def yscramble_regression(X, y, folds, n_perm: int, seed: int) -> dict:
    real_pred = oof_regression(X, y, folds)
    real = regression_metrics(y, real_pred)["spearman"]
    rng = np.random.default_rng(seed)
    perm = []
    for _ in range(n_perm):
        yp = y.copy()
        rng.shuffle(yp)
        pr = oof_regression(X, yp, folds)
        perm.append(regression_metrics(yp, pr)["spearman"])
    perm = np.array(perm)
    p_val = float((np.sum(perm >= real) + 1) / (n_perm + 1))
    return {
        "metric": "spearman_oof",
        "real": round(float(real), 4),
        "permuted_mean": round(float(np.nanmean(perm)), 4),
        "permuted_std": round(float(np.nanstd(perm)), 4),
        "permuted_max": round(float(np.nanmax(perm)), 4),
        "n_permutations": n_perm,
        "empirical_p_value": p_val,
        "passes": bool(p_val < 0.05 and real > np.nanmax(perm)),
    }


def yscramble_classification(X, y, folds, n_perm: int, seed: int) -> dict:
    real_prob = oof_classification(X, y, folds)
    real = classification_metrics(y, real_prob)["auroc"]
    rng = np.random.default_rng(seed)
    perm = []
    for _ in range(n_perm):
        yp = y.copy()
        rng.shuffle(yp)
        pr = oof_classification(X, yp, folds)
        perm.append(classification_metrics(yp, pr)["auroc"])
    perm = np.array(perm, dtype=float)
    p_val = float((np.sum(perm >= real) + 1) / (n_perm + 1))
    return {
        "metric": "auroc_oof",
        "real": round(float(real), 4),
        "permuted_mean": round(float(np.nanmean(perm)), 4),
        "permuted_std": round(float(np.nanstd(perm)), 4),
        "permuted_max": round(float(np.nanmax(perm)), 4),
        "n_permutations": n_perm,
        "empirical_p_value": p_val,
        "passes": bool(p_val < 0.05 and real > np.nanmax(perm)),
    }


def fps_bitvects(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        out.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None)
    return out


def nn_tanimoto_to_ref(query_fps, ref_fps, exclude_self: bool = False) -> list[float]:
    from rdkit import DataStructs

    ref = [f for f in ref_fps if f is not None]
    out = []
    for i, q in enumerate(query_fps):
        if q is None:
            out.append(float("nan"))
            continue
        sims = DataStructs.BulkTanimotoSimilarity(q, ref)
        if exclude_self and sims:
            sims = sorted(sims, reverse=True)
            out.append(float(sims[1]) if len(sims) > 1 else float(sims[0]))
        else:
            out.append(float(max(sims)) if sims else float("nan"))
    return out


def applicability_domain(urat1, nlrp3_actives_smiles, benchmarks, shortlist, ad_pct=5) -> pd.DataFrame:
    """AD threshold = ad_pct-th percentile of train intra-NN Tanimoto (leave-one-out).

    A query is 'in domain' if its nearest-neighbour Tanimoto to train >= threshold.
    Computed against URAT1 train set (the reference chemistry for the docking-led arm).
    """
    train_smiles = urat1[SMILES_COL].tolist()
    train_fps = fps_bitvects(train_smiles)
    intra = nn_tanimoto_to_ref(train_fps, train_fps, exclude_self=True)
    threshold = float(np.nanpercentile(intra, ad_pct))

    rows = []
    # Benchmarks
    if "canonical_smiles" in benchmarks.columns:
        bsmi = benchmarks["canonical_smiles"].astype(str).tolist()
        bname = benchmarks.get("compound_name", benchmarks.get("compound_id", pd.Series([""] * len(bsmi)))).astype(str).tolist()
    else:
        bsmi, bname = [], []
    bfps = fps_bitvects(bsmi)
    b_nn = nn_tanimoto_to_ref(bfps, train_fps)
    for name, smi, nn in zip(bname, bsmi, b_nn):
        rows.append({"set": "benchmark", "name": name, "nn_tanimoto_urat1_train": round(nn, 3) if nn == nn else None,
                     "in_domain": bool(nn >= threshold) if nn == nn else None})
    # Shortlist
    ssmi = shortlist[SMILES_COL].astype(str).tolist()
    sname = shortlist.get("name", pd.Series([""] * len(ssmi))).astype(str).tolist()
    sfps = fps_bitvects(ssmi)
    s_nn = nn_tanimoto_to_ref(sfps, train_fps)
    for name, smi, nn in zip(sname, ssmi, s_nn):
        rows.append({"set": "pareto_shortlist", "name": name, "nn_tanimoto_urat1_train": round(nn, 3) if nn == nn else None,
                     "in_domain": bool(nn >= threshold) if nn == nn else None})
    df = pd.DataFrame(rows)
    df.attrs["threshold"] = threshold
    return df, threshold


def calibration_curve(y_true, y_prob, n_bins=10) -> pd.DataFrame:
    mask = ~np.isnan(y_true) & ~np.isnan(y_prob)
    yt, yp = y_true[mask].astype(int), y_prob[mask]
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.digitize(yp, bins) - 1
    idx = np.clip(idx, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        sel = idx == b
        if sel.sum() == 0:
            continue
        rows.append({
            "bin_lower": round(float(bins[b]), 3),
            "bin_upper": round(float(bins[b + 1]), 3),
            "n": int(sel.sum()),
            "mean_predicted": round(float(yp[sel].mean()), 4),
            "observed_freq": round(float(yt[sel].mean()), 4),
        })
    brier = float(np.mean((yp - yt) ** 2))
    df = pd.DataFrame(rows)
    df.attrs["brier"] = brier
    return df, brier


def main() -> None:
    parser = argparse.ArgumentParser(description="ML rigor battery (non-docking module A)")
    parser.add_argument("--n-permutations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = {"module": "A_ml_rigor_validation", "n_permutations": args.n_permutations}

    # ---------- URAT1 regression ----------
    urat1 = pd.read_csv(PROCESSED / "urat1_curated.csv")
    Xu = featurize_smiles(urat1[SMILES_COL].tolist())
    yu = urat1["pActivity"].to_numpy(dtype=float)
    folds_u = scaffold_cv_indices(urat1[SMILES_COL].tolist(), n_splits=5)
    print("URAT1 y-scrambling...")
    ys_u = yscramble_regression(Xu, yu, folds_u, args.n_permutations, args.seed)
    with open(args.output_dir / "yscramble_urat1.json", "w") as f:
        json.dump(ys_u, f, indent=2)
    summary["urat1_yscramble"] = ys_u
    print(f"  real Spearman={ys_u['real']} perm_max={ys_u['permuted_max']} p={ys_u['empirical_p_value']} pass={ys_u['passes']}")

    # ---------- NLRP3 classification ----------
    nlrp3 = pd.read_csv(PROCESSED / "nlrp3_records.csv")
    nlrp3_mol = nlrp3.drop_duplicates(SMILES_COL).reset_index(drop=True)
    Xn = featurize_smiles(nlrp3_mol[SMILES_COL].tolist())
    yn = nlrp3_mol["active"].to_numpy(dtype=float)
    folds_n = scaffold_cv_indices(nlrp3_mol[SMILES_COL].tolist(), n_splits=5)
    print("NLRP3 y-scrambling...")
    ys_n = yscramble_classification(Xn, yn, folds_n, args.n_permutations, args.seed)
    with open(args.output_dir / "yscramble_nlrp3.json", "w") as f:
        json.dump(ys_n, f, indent=2)
    summary["nlrp3_yscramble"] = ys_n
    print(f"  real AUROC={ys_n['real']} perm_max={ys_n['permuted_max']} p={ys_n['empirical_p_value']} pass={ys_n['passes']}")

    # ---------- NLRP3 calibration ----------
    print("NLRP3 calibration...")
    prob_n = oof_classification(Xn, yn, folds_n)
    cal, brier = calibration_curve(yn, prob_n)
    cal.to_csv(args.output_dir / "nlrp3_calibration.csv", index=False)
    summary["nlrp3_brier_score"] = round(brier, 4)
    print(f"  Brier score={round(brier, 4)}")

    # ---------- Applicability domain ----------
    print("Applicability domain...")
    benchmarks = pd.read_csv(BENCHMARKS / "literature_benchmarks.csv")
    shortlist = pd.read_csv(PARETO_DIR / "pareto_shortlist.csv")
    nlrp3_actives = nlrp3_mol[nlrp3_mol["active"] == 1][SMILES_COL].tolist()
    ad, threshold = applicability_domain(urat1, nlrp3_actives, benchmarks, shortlist)
    ad.to_csv(args.output_dir / "applicability_domain.csv", index=False)
    summary["ad_threshold_urat1_train"] = round(threshold, 4)

    with open(args.output_dir / "ml_rigor_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== ML rigor summary ===")
    print(f"  URAT1 y-scramble p={ys_u['empirical_p_value']} (real {ys_u['real']} vs perm max {ys_u['permuted_max']})")
    print(f"  NLRP3 y-scramble p={ys_n['empirical_p_value']} (real {ys_n['real']} vs perm max {ys_n['permuted_max']})")
    print(f"  NLRP3 Brier={summary['nlrp3_brier_score']}  AD threshold(URAT1)={summary['ad_threshold_urat1_train']}")
    print("\nApplicability domain (shortlist):")
    print(ad[ad["set"] == "pareto_shortlist"].to_string(index=False))


if __name__ == "__main__":
    main()
