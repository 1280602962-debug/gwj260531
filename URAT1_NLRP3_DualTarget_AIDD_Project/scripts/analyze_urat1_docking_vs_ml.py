#!/usr/bin/env python3
"""
URAT1 retrospective: 8973 docking vs ML on shared manifest.

Compares:
  - Subset A (labeled actives) vs D (unlabeled decoys): enrichment by S_U (9DKB XP)
  - Benchmark 4-drug recovery: ML rank vs docking rank
  - Optional hybrid ML (Morgan+RDKit + glide_xp feature) ablation on subset A

Outputs:
  results/docking/urat1_docking_vs_ml_summary.json
  results/docking/urat1_benchmark_rankings_docking.csv
  results/docking/urat1_enrichment_subset_a_vs_d.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from utils_ml import (
    featurize_smiles,
    regression_metrics,
    roc_auc_binary,
    scaffold_cv_indices,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MERGED = PROJECT_ROOT / "data" / "docking" / "8973_9DKB_with_manifest.csv"
BENCHMARKS = PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks_summary.csv"
URAT1_MODEL = PROJECT_ROOT / "results" / "training" / "urat1_model.joblib"
URAT1_CURATED = PROJECT_ROOT / "data" / "processed" / "urat1_curated.csv"
OUT = PROJECT_ROOT / "results" / "docking"

URAT1_BENCH = ["lesinurad", "benzbromarone", "verinurad", "dotinurad"]


def _percentile_rank(series: pd.Series, values: pd.Series) -> np.ndarray:
    return values.map(lambda v: (series <= v).mean() * 100.0 if pd.notna(v) else np.nan).values


def enrichment_a_vs_d(merged: pd.DataFrame) -> dict:
    a = merged[(merged["subset"] == "A") & merged["docked"]].copy()
    d = merged[(merged["subset"] == "D") & merged["docked"]].copy()
    if a.empty or d.empty:
        return {"error": "subset A or D empty after docking merge"}

    # Active if pActivity >= 6 (project convention)
    a["active"] = pd.to_numeric(a["pActivity"], errors="coerce") >= 6.0
    scores = pd.concat([a, d], ignore_index=True)
    scores["is_active_arm"] = scores["subset"] == "A"
    scores["score"] = scores["s_u_raw"]

    y = scores["is_active_arm"].astype(int).values
    s = scores["score"].values
    auc = roc_auc_binary(y, s)

    # EF@5% and EF@10% on full docked pool (A positives vs D negatives)
    order = np.argsort(-s)
    y_sorted = y[order]
    n = len(y_sorted)
    ef = {}
    for pct in (5, 10):
        k = max(1, int(np.ceil(n * pct / 100.0)))
        ef[f"ef_{pct}pct_a_vs_d"] = float(y_sorted[:k].sum() / max(y.sum(), 1) * 100.0 / pct)

    return {
        "n_subset_a_docked": int(len(a)),
        "n_subset_d_docked": int(len(d)),
        "n_active_ge6_in_a": int(a["active"].sum()),
        "roc_auc_active_a_vs_decoy_d": round(float(auc), 4),
        **{k: round(v, 3) for k, v in ef.items()},
        "spearman_pactivity_vs_su_on_a": round(
            float(a["pActivity"].corr(a["s_u_raw"], method="spearman")), 4
        )
        if a["pActivity"].notna().sum() >= 10
        else None,
    }


def benchmark_rankings(merged: pd.DataFrame, model_path: Path) -> pd.DataFrame:
    docked = merged[merged["docked"]].copy()
    su = docked["s_u_raw"]

    panel_path = PROJECT_ROOT / "data" / "distill" / "teacher_gate_qc_panel_b_direction.csv"
    panel = pd.read_csv(panel_path) if panel_path.exists() else pd.DataFrame()
    bench_lit = pd.read_csv(BENCHMARKS)
    bench_lit = bench_lit[bench_lit["target"] == "URAT1"]

    predict_map: dict[str, float] = {}
    if model_path.exists():
        bundle = joblib.load(model_path)
        model = bundle["model"]
        scaler = bundle["scaler"]
        for smi in docked["canonical_smiles"].unique():
            x = scaler.transform(featurize_smiles([smi]))
            predict_map[smi] = float(model.predict(x)[0])

    rows = []
    for drug in URAT1_BENCH:
        smi = None
        if not panel.empty:
            hit = panel[panel["compound_name"].str.lower() == drug]
            if len(hit):
                smi = hit["canonical_smiles"].iloc[0]
        if smi is None:
            rows.append({"compound": drug, "docked": False})
            continue

        row = docked[docked["canonical_smiles"] == smi]
        if not len(row):
            rows.append({"compound": drug, "canonical_smiles": smi, "docked": False})
            continue

        r = row.iloc[0]
        su_pct = float((su <= r["s_u_raw"]).mean() * 100.0)
        ml = predict_map.get(smi)
        ml_pct = None
        if ml is not None and predict_map:
            all_ml = np.array(list(predict_map.values()))
            ml_pct = float((all_ml <= ml).mean() * 100.0)
        lit = bench_lit[bench_lit["compound_name"].str.lower() == drug]["primary_pactivity"]
        lit_p = float(lit.iloc[0]) if len(lit) else None
        rows.append(
            {
                "compound": drug,
                "canonical_smiles": smi,
                "docked": True,
                "glide_score_xp": float(r["glide_score_xp"]),
                "s_u_percentile": float(r["s_u_percentile"]),
                "ml_pred_pactivity": ml,
                "ml_percentile_vs_8973": ml_pct,
                "literature_pactivity": lit_p,
                "docking_pass_top10pct": bool(su_pct >= 90.0),
                "ml_pass_pred_ge6": bool(ml >= 6.0) if ml is not None else None,
            }
        )
    return pd.DataFrame(rows)


def hybrid_cv_ablation(merged: pd.DataFrame, n_splits: int = 5, seed: int = 42) -> dict:
    """Subset A only: does adding glide_xp improve scaffold CV? (SI ablation)"""
    import xgboost as xgb
    from sklearn.preprocessing import StandardScaler

    a = merged[(merged["subset"] == "A") & merged["docked"]].copy()
    a = a.dropna(subset=["pActivity", "s_u_raw"])
    if len(a) < 50:
        return {"skipped": True, "reason": "too few docked labeled compounds"}

    smiles = a["canonical_smiles"].tolist()
    y = a["pActivity"].astype(float).values
    x_fp = featurize_smiles(smiles)
    x_hybrid = np.column_stack([x_fp, a["s_u_raw"].values])
    groups = a["scaffold"].fillna("unknown").values

    def _run(X: np.ndarray) -> dict:
        oof = np.zeros(len(y))
        for tr, te in scaffold_cv_indices(smiles, n_splits=n_splits):
            scaler = StandardScaler()
            x_tr = scaler.fit_transform(X[tr])
            x_te = scaler.transform(X[te])
            model = xgb.XGBRegressor(
                n_estimators=300,
                max_depth=5,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                random_state=seed,
                n_jobs=-1,
            )
            model.fit(x_tr, y[tr])
            oof[te] = model.predict(x_te)
        m = regression_metrics(y, oof)
        return {"spearman": round(float(m["spearman"]), 4), "r2": round(float(m["r2"]), 4)}

    return {
        "n_compounds": int(len(a)),
        "fp_only": _run(x_fp),
        "fp_plus_glide_xp": _run(x_hybrid),
        "note": "Hybrid requires docking at inference; use for ablation/SI only.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED)
    parser.add_argument("--model", type=Path, default=URAT1_MODEL)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--hybrid-cv", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    merged = pd.read_csv(args.merged)
    enrich = enrichment_a_vs_d(merged)
    bench_df = benchmark_rankings(merged, args.model)
    bench_df.to_csv(args.output_dir / "urat1_benchmark_rankings_docking.csv", index=False)

    dock_pass = int(bench_df.get("docking_pass_top10pct", pd.Series(dtype=bool)).sum())
    ml_pass = int(bench_df.get("ml_pass_pred_ge6", pd.Series(dtype=bool)).sum())

    summary = {
        "enrichment": enrich,
        "benchmark_recovery": {
            "docking_top10pct_pass": f"{dock_pass}/{len(bench_df)}",
            "ml_pred_ge6_pass": f"{ml_pass}/{len(bench_df)}",
        },
        "recommendation": (
            "Use 9DKB XP (S_U) as primary URAT1 evidence on 8973; "
            "keep ML as auxiliary unless hybrid ablation shows clear gain."
        ),
    }
    if args.hybrid_cv:
        summary["hybrid_ablation"] = hybrid_cv_ablation(merged)

    out_path = args.output_dir / "urat1_docking_vs_ml_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
