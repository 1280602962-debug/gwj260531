#!/usr/bin/env python3
"""Compare URAT1 ML/DL models: identical scaffold-CV + benchmark recovery."""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from utils_ml import (
    assay_one_hot_matrix,
    canonicalize,
    featurize_smiles,
    max_tanimoto_to_library,
    regression_enrichment_factor,
    regression_metrics,
    roc_auc_binary,
    save_json,
    scaffold_cv_indices,
)

warnings.filterwarnings("ignore")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
OUT = PROJECT_ROOT / "results" / "urat1_model_comparison"

BENCHMARKS = [
    ("lesinurad", "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12", 7.0),
    ("benzbromarone", "CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1", 7.5),
    ("verinurad", "CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O", 8.0),
    ("dotinurad", "O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21", 8.2),
]


def eval_model(name: str, smiles: list[str], y: np.ndarray, train_fn, pred_fn) -> dict:
    oof_pred = np.full(len(y), np.nan)
    for tr, te in scaffold_cv_indices(smiles, 5):
        bundle = train_fn([smiles[i] for i in tr], y[tr])
        for i in te:
            oof_pred[i] = pred_fn(bundle, smiles[i])
    met = regression_metrics(y, oof_pred)
    bench = []
    full_bundle = train_fn(smiles, y)
    train_set = set(smiles)
    for bname, smi, lit in BENCHMARKS:
        c = canonicalize(smi)
        pred = float(pred_fn(full_bundle, c))
        bench.append({
            "compound": bname,
            "in_train": c in train_set,
            "pred": round(pred, 3),
            "lit": lit,
            "pass": pred >= 6.0,
            "max_tc": round(max_tanimoto_to_library(c, smiles), 3),
        })
    return {
        "model": name,
        "rmse": met["rmse"],
        "r2": met["r2"],
        "spearman": met["spearman"],
        "ef_5pct_p7": regression_enrichment_factor(y, oof_pred, 7.0, 0.05),
        "roc_auc_p7": roc_auc_binary(y, oof_pred, 7.0),
        "benchmark_pass": f"{sum(b['pass'] for b in bench)}/4",
        "benchmark_detail": bench,
    }


def run_chemprop(smiles: list[str], y: np.ndarray) -> dict:
    oof_pred = np.full(len(y), np.nan)
    bench_fold_preds = {n: [] for n, _, _ in BENCHMARKS}
    for tr, te in scaffold_cv_indices(smiles, 5):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            pd.DataFrame({"smiles": [smiles[i] for i in tr], "target": y[tr]}).to_csv(tdp / "train.csv", index=False)
            pd.DataFrame({"smiles": [smiles[i] for i in te], "target": y[te]}).to_csv(tdp / "test.csv", index=False)
            subprocess.run([
                "chemprop", "train", "--data-path", str(tdp / "train.csv"),
                "--task-type", "regression", "--smiles-columns", "smiles",
                "--target-columns", "target", "--epochs", "30",
                "--batch-size", "64", "--save-dir", str(tdp / "model"),
                "--num-workers", "0", "--accelerator", "cpu",
            ], check=True, capture_output=True)
            subprocess.run([
                "chemprop", "predict", "--test-path", str(tdp / "test.csv"),
                "--model-path", str(tdp / "model"), "--preds-path", str(tdp / "preds.csv"),
                "--num-workers", "0", "--accelerator", "cpu",
            ], check=True, capture_output=True)
            preds = pd.read_csv(tdp / "preds.csv")
            for j, i in enumerate(te):
                oof_pred[i] = preds["target"].iloc[j]
            for bname, smi, _ in BENCHMARKS:
                c = canonicalize(smi)
                pd.DataFrame({"smiles": [c], "target": [0.0]}).to_csv(tdp / "bench.csv", index=False)
                subprocess.run([
                    "chemprop", "predict", "--test-path", str(tdp / "bench.csv"),
                    "--model-path", str(tdp / "model"), "--preds-path", str(tdp / "bench_p.csv"),
                    "--num-workers", "0", "--accelerator", "cpu",
                ], check=True, capture_output=True)
                bench_fold_preds[bname].append(float(pd.read_csv(tdp / "bench_p.csv")["target"].iloc[0]))
    met = regression_metrics(y, oof_pred)
    bench = []
    for bname, _, lit in BENCHMARKS:
        p = float(np.mean(bench_fold_preds[bname]))
        bench.append({"compound": bname, "pred": round(p, 3), "lit": lit, "pass": p >= 6.0})
    return {
        "model": "Chemprop_D-MPNN",
        "rmse": met["rmse"], "r2": met["r2"], "spearman": met["spearman"],
        "ef_5pct_p7": regression_enrichment_factor(y, oof_pred, 7.0, 0.05),
        "roc_auc_p7": roc_auc_binary(y, oof_pred, 7.0),
        "benchmark_pass": f"{sum(b['pass'] for b in bench)}/4",
        "benchmark_detail": bench,
    }


def run_assay_conditioned(df: pd.DataFrame, raw: pd.DataFrame) -> dict:
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    assay_top = raw["Assay ChEMBL ID"].value_counts().head(20).index.astype(str).tolist()
    oof_pred = np.full(len(y), np.nan)
    for tr, te in scaffold_cv_indices(smiles, 5):
        tr_set = {smiles[i] for i in tr}
        tr_rec = raw[raw["canonical_smiles"].isin(tr_set)]
        x_tr = np.hstack([
            featurize_smiles(tr_rec["canonical_smiles"].tolist()),
            assay_one_hot_matrix(tr_rec["Assay ChEMBL ID"], assay_top),
        ])
        sc = StandardScaler().fit(x_tr)
        m = xgb.XGBRegressor(n_estimators=400, max_depth=5, learning_rate=0.05, random_state=42, n_jobs=-1)
        m.fit(sc.transform(x_tr), tr_rec["pActivity"].values)
        for i in te:
            x_te = np.hstack([
                featurize_smiles([smiles[i]]),
                assay_one_hot_matrix(pd.Series([assay_top[0]]), assay_top),
            ])
            oof_pred[i] = m.predict(sc.transform(x_te))[0]
    met = regression_metrics(y, oof_pred)
    tr_rec = raw
    x_tr = np.hstack([
        featurize_smiles(tr_rec["canonical_smiles"].tolist()),
        assay_one_hot_matrix(tr_rec["Assay ChEMBL ID"], assay_top),
    ])
    sc = StandardScaler().fit(x_tr)
    m = xgb.XGBRegressor(n_estimators=400, max_depth=5, learning_rate=0.05, random_state=42, n_jobs=-1)
    m.fit(sc.transform(x_tr), tr_rec["pActivity"].values)
    bench = []
    for bname, smi, lit in BENCHMARKS:
        c = canonicalize(smi)
        x = np.hstack([featurize_smiles([c]), assay_one_hot_matrix(pd.Series([assay_top[0]]), assay_top)])
        p = float(m.predict(sc.transform(x))[0])
        bench.append({"compound": bname, "pred": round(p, 3), "lit": lit, "pass": p >= 6.0})
    return {
        "model": "XGBoost_assay_conditioned",
        "rmse": met["rmse"], "r2": met["r2"], "spearman": met["spearman"],
        "ef_5pct_p7": regression_enrichment_factor(y, oof_pred, 7.0, 0.05),
        "roc_auc_p7": roc_auc_binary(y, oof_pred, 7.0),
        "benchmark_pass": f"{sum(b['pass'] for b in bench)}/4",
        "benchmark_detail": bench,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PROCESSED / "urat1_curated.csv")
    smiles = df["canonical_smiles"].tolist()
    y = df["pActivity"].values.astype(float)
    models = []

    def train_xgb(s, ya):
        sc = StandardScaler()
        Xs = sc.fit_transform(featurize_smiles(s))
        m = xgb.XGBRegressor(n_estimators=400, max_depth=5, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, random_state=42, n_jobs=-1)
        m.fit(Xs, ya)
        return {"m": m, "sc": sc}

    def pred_xgb(b, smi):
        return b["m"].predict(b["sc"].transform(featurize_smiles([smi])))[0]

    models.append(eval_model("XGBoost_Morgan+RDKit (current)", smiles, y, train_xgb, pred_xgb))

    def train_xgb_fp(s, ya):
        sc = StandardScaler()
        Xs = sc.fit_transform(featurize_smiles(s)[:, :2048])
        m = xgb.XGBRegressor(n_estimators=400, max_depth=5, learning_rate=0.05, random_state=42, n_jobs=-1)
        m.fit(Xs, ya)
        return {"m": m, "sc": sc, "fp": True}

    def pred_xgb_fp(b, smi):
        return b["m"].predict(b["sc"].transform(featurize_smiles([smi])[:, :2048]))[0]

    models.append(eval_model("XGBoost_Morgan2048_only", smiles, y, train_xgb_fp, pred_xgb_fp))

    def train_rf(s, ya):
        m = RandomForestRegressor(500, max_depth=12, min_samples_leaf=2, random_state=42, n_jobs=-1)
        m.fit(featurize_smiles(s), ya)
        return {"m": m}

    def pred_rf(b, smi):
        return b["m"].predict(featurize_smiles([smi]))[0]

    models.append(eval_model("RandomForest_Morgan+RDKit", smiles, y, train_rf, pred_rf))

    def train_svr(s, ya):
        sc = StandardScaler()
        Xs = sc.fit_transform(featurize_smiles(s))
        pca = PCA(n_components=min(50, len(s) - 1, Xs.shape[1]))
        m = SVR(C=10, gamma="scale", epsilon=0.1)
        m.fit(pca.fit_transform(Xs), ya)
        return {"m": m, "sc": sc, "pca": pca}

    def pred_svr(b, smi):
        x = b["pca"].transform(b["sc"].transform(featurize_smiles([smi])))
        return b["m"].predict(x)[0]

    models.append(eval_model("SVR_PCA50 (PLK1-style)", smiles, y, train_svr, pred_svr))

    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    def fps(slist):
        return [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 2048) for s in slist]

    def train_knn(s, ya):
        return {"fps": fps(s), "y": np.array(ya), "k": 5}

    def pred_knn(b, smi):
        q = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2, 2048)
        sims = np.array([DataStructs.TanimotoSimilarity(q, f) for f in b["fps"]])
        order = np.argsort(-sims)[: b["k"]]
        w = sims[order]
        return float(np.average(b["y"][order], weights=w) if w.sum() else b["y"].mean())

    models.append(eval_model("kNN_Tanimoto_k5", smiles, y, train_knn, pred_knn))

    def train_krr(s, ya):
        sc = StandardScaler()
        m = KernelRidge(alpha=1.0, kernel="rbf", gamma=0.001)
        m.fit(sc.fit_transform(featurize_smiles(s)), ya)
        return {"m": m, "sc": sc}

    def pred_krr(b, smi):
        return b["m"].predict(b["sc"].transform(featurize_smiles([smi])))[0]

    models.append(eval_model("KernelRidge_RBF", smiles, y, train_krr, pred_krr))

    raw = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "URAT1_CHEMBL_cf12.csv", low_memory=False)
    raw = raw[raw["Standard Relation"].astype(str).str.strip("'\"") == "="]
    raw["pActivity"] = pd.to_numeric(raw["pChEMBL Value"], errors="coerce")
    raw = raw[raw["pActivity"].between(4, 10)]
    raw["canonical_smiles"] = raw["Smiles"].map(canonicalize)
    raw = raw[raw["canonical_smiles"].notna()]
    models.append(run_assay_conditioned(df, raw))

    try:
        print("Running Chemprop (5-fold, ~5-10 min)...")
        models.append(run_chemprop(smiles, y))
    except Exception as e:
        models.append({"model": "Chemprop_D-MPNN", "error": str(e)})

    summary = [{k: m.get(k) for k in ["model", "rmse", "r2", "spearman", "ef_5pct_p7", "roc_auc_p7", "benchmark_pass", "error"]} for m in models]
    df_out = pd.DataFrame(summary)
    print(df_out.to_string(index=False))
    save_json(OUT / "comparison_full.json", {"models": models})
    df_out.to_csv(OUT / "comparison_summary.csv", index=False)
    print(f"\nSaved to {OUT}")


if __name__ == "__main__":
    main()
