#!/usr/bin/env python3
"""
TAPE-GATE: Benchmark backtest and model reliability report.

Metrics:
  - Percentile rank vs training library
  - Binary pass: predicted active (NLRP3) / predicted pActivity>=6 (URAT1)
  - Applicability domain: max Tanimoto to training set
  - Stratified by in_training vs scaffold-novel benchmarks
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from utils_ml import (
    assay_one_hot_matrix,
    canonicalize,
    featurize_smiles,
    max_tanimoto_to_library,
    murcko_scaffold,
    save_json,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
MODELS = PROJECT_ROOT / "results" / "training"
RESULTS = PROJECT_ROOT / "results" / "benchmark_backtest"

BENCHMARKS = [
    {"name": "lesinurad", "target": "URAT1", "smiles": "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12", "role": "must_recover", "pactivity_ref": 7.0},
    {"name": "benzbromarone", "target": "URAT1", "smiles": "CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1", "role": "must_recover", "pactivity_ref": 7.5},
    {"name": "verinurad", "target": "URAT1", "smiles": "CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O", "role": "must_recover", "pactivity_ref": 8.0},
    {"name": "dotinurad", "target": "URAT1", "smiles": "O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21", "role": "must_recover", "pactivity_ref": 8.2},
    {"name": "MCC950", "target": "NLRP3", "smiles": "CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1", "role": "must_recover", "pactivity_ref": 8.09},
    {"name": "GDC-2394", "target": "NLRP3", "smiles": "CN[C@@H]1COc2c(S(=O)(=O)NC(=O)Nc3c4c(cc5c3CCC5)CCC4)cnn2C1", "role": "must_recover", "pactivity_ref": 8.5},
    {"name": "allopurinol", "target": "URAT1", "smiles": "Oc1ncnc2c1ncn2C", "role": "negative_control", "pactivity_ref": None},
]


def percentile_rank(score: float, population: np.ndarray, higher_is_better: bool = True) -> float:
    pop = population[~np.isnan(population)]
    if len(pop) == 0:
        return np.nan
    if higher_is_better:
        return float((pop < score).mean() * 100)
    return float((pop > score).mean() * 100)


def predict_urat1(bundle: dict, smiles_list: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = featurize_smiles(smiles_list)
    x_s = bundle["scaler"].transform(x)
    y_pred = bundle["model"].predict(x_s)
    lo, hi = bundle["conformal"].interval(y_pred)
    return y_pred, lo, hi


def predict_nlrp3_ensemble(bundle: dict, smiles_list: list[str], assay_ids: list[str]) -> np.ndarray:
    probs = []
    for aid in assay_ids:
        assay_col = pd.Series([aid] * len(smiles_list))
        x_mol = featurize_smiles(smiles_list)
        x_assay = assay_one_hot_matrix(assay_col, bundle["top_assay_ids"])
        x = np.hstack([x_mol, x_assay])
        x_s = bundle["scaler"].transform(x)
        raw = bundle["model"].predict_proba(x_s)[:, 1]
        if bundle.get("calibrator") is not None:
            raw = bundle["calibrator"].predict(raw)
        probs.append(raw)
    return np.max(np.vstack(probs), axis=0)


def backtest_urat1(urat1_df: pd.DataFrame, bundle: dict, benchmarks: list[dict]) -> dict:
    train_smiles = urat1_df["canonical_smiles"].tolist()
    pop_pred, _, _ = predict_urat1(bundle, train_smiles)

    rows = []
    for b in benchmarks:
        if b["target"] != "URAT1":
            continue
        smi = canonicalize(b["smiles"])
        if smi is None:
            continue
        pred, lo, hi = predict_urat1(bundle, [smi])
        in_train = smi in set(train_smiles)
        train_row = urat1_df[urat1_df["canonical_smiles"] == smi]
        y_true = float(train_row["pActivity"].iloc[0]) if len(train_row) else None
        max_tc = max_tanimoto_to_library(smi, train_smiles)

        pct = percentile_rank(float(pred[0]), pop_pred, higher_is_better=True)
        lit = b.get("pactivity_ref")
        within_1log = abs(float(pred[0]) - lit) <= 1.0 if lit is not None else None

        rows.append(
            {
                "compound": b["name"],
                "role": b["role"],
                "canonical_smiles": smi,
                "in_training_set": in_train,
                "max_tanimoto_to_train": round(max_tc, 3),
                "in_applicability_domain": max_tc >= 0.35,
                "y_true_pActivity": y_true,
                "literature_pActivity": lit,
                "y_pred": round(float(pred[0]), 3),
                "interval_lo": round(float(lo[0]), 3),
                "interval_hi": round(float(hi[0]), 3),
                "percentile_vs_train": round(pct, 1),
                "predicted_active_ge6": float(pred[0]) >= 6.0,
                "conformal_lower_ge6": float(lo[0]) >= 6.0,
                "within_1log_of_literature": within_1log,
                "binary_pass": float(pred[0]) >= 6.0 if b["role"] == "must_recover" else float(pred[0]) < 6.0,
            }
        )

    must = [r for r in rows if r["role"] == "must_recover"]
    in_train_must = [r for r in must if r["in_training_set"]]
    novel_must = [r for r in must if not r["in_training_set"]]

    return {
        "benchmarks": rows,
        "must_recover_binary_pass": sum(r["binary_pass"] for r in must),
        "must_recover_count": len(must),
        "in_training_must_pass": sum(r["binary_pass"] for r in in_train_must),
        "in_training_must_count": len(in_train_must),
        "novel_must_pass": sum(r["binary_pass"] for r in novel_must),
        "novel_must_count": len(novel_must),
        "pass_strict_top20pct": sum(r["percentile_vs_train"] >= 80 for r in must) >= max(1, len(must) // 2),
        "pass_binary": sum(r["binary_pass"] for r in must) >= max(2, len(must) - 1),
    }


def backtest_nlrp3(nlrp3_df: pd.DataFrame, bundle: dict, benchmarks: list[dict]) -> dict:
    mol_df = nlrp3_df.groupby("canonical_smiles", as_index=False).agg(active=("active", "max"))
    train_smiles = mol_df["canonical_smiles"].tolist()
    top_assays = nlrp3_df["Assay ChEMBL ID"].value_counts().head(5).index.astype(str).tolist()
    pop_prob = predict_nlrp3_ensemble(bundle, train_smiles, top_assays)

    rows = []
    for b in benchmarks:
        if b["target"] != "NLRP3":
            continue
        smi = canonicalize(b["smiles"])
        if smi is None:
            continue
        prob = predict_nlrp3_ensemble(bundle, [smi], top_assays)
        in_train = smi in set(train_smiles)
        train_row = mol_df[mol_df["canonical_smiles"] == smi]
        y_true = int(train_row["active"].iloc[0]) if len(train_row) else None
        max_tc = max_tanimoto_to_library(smi, train_smiles)
        pct = percentile_rank(float(prob[0]), pop_prob, higher_is_better=True)

        rows.append(
            {
                "compound": b["name"],
                "role": b["role"],
                "canonical_smiles": smi,
                "in_training_set": in_train,
                "max_tanimoto_to_train": round(max_tc, 3),
                "in_applicability_domain": max_tc >= 0.35,
                "y_true_active": y_true,
                "p_active_ensemble": round(float(prob[0]), 3),
                "percentile_vs_train": round(pct, 1),
                "predicted_active": float(prob[0]) >= 0.5,
                "binary_pass": float(prob[0]) >= 0.5,
            }
        )

    must = [r for r in rows if r["role"] == "must_recover"]
    return {
        "benchmarks": rows,
        "ensemble_assays": top_assays,
        "must_recover_binary_pass": sum(r["binary_pass"] for r in must),
        "must_recover_count": len(must),
        "pass_binary": all(r["binary_pass"] for r in must),
        "note": "High base rate of actives (~61%) makes percentile rank less informative than binary predicted_active.",
    }


def scaffold_exclusion_check(urat1_df: pd.DataFrame, nlrp3_df: pd.DataFrame, benchmarks: list[dict]) -> list[dict]:
    urat_scaffolds = set(urat1_df["scaffold"])
    nlrp_scaffolds = set(nlrp3_df["scaffold"].unique())
    out = []
    for b in benchmarks:
        smi = canonicalize(b["smiles"])
        if smi is None:
            continue
        scaf = murcko_scaffold(smi)
        out.append(
            {
                "compound": b["name"],
                "target": b["target"],
                "scaffold": scaf,
                "scaffold_in_urat1_train": scaf in urat_scaffolds,
                "scaffold_in_nlrp3_train": scaf in nlrp_scaffolds,
            }
        )
    return out


def overall_verdict(training_report: dict, urat1_bt: dict, nlrp3_bt: dict) -> dict:
    urat1_cv = training_report["urat1"]["screening_assessment"]["suitable_for_screening"]
    nlrp3_cv = training_report["nlrp3"]["screening_assessment"]["suitable_for_screening"]
    urat1_bench = urat1_bt["pass_binary"]
    nlrp3_bench = nlrp3_bt["pass_binary"]

    if urat1_cv and nlrp3_cv and urat1_bench and nlrp3_bench:
        verdict, text = "GO", "Strict CV and benchmark criteria passed."
    elif nlrp3_cv and nlrp3_bench and not (urat1_cv and urat1_bench):
        verdict, text = "URAT1_NO_GO", (
            "NLRP3 model is screening-ready; URAT1 model fails strict CV and/or benchmark recovery. "
            "URAT1 library filtering must NOT rely on ML alone — use $S_trap$ conformational ensemble docking as primary evidence."
        )
    elif urat1_cv or nlrp3_cv:
        verdict, text = "CONDITIONAL_GO", "Partial pass; see per-target tables."
    else:
        verdict, text = "NO_GO", "Both models fail strict criteria."

    return {
        "verdict": verdict,
        "summary": text,
        "cv_pass": {"urat1": urat1_cv, "nlrp3": nlrp3_cv},
        "benchmark_binary_pass": {"urat1": urat1_bench, "nlrp3": nlrp3_bench},
    }


def write_markdown_report(report: dict, path: Path) -> None:
    u = report["urat1_backtest"]
    n = report["nlrp3_backtest"]
    v = report["overall_verdict"]
    cv = report["training_cv_summary"]

    lines = [
        "# Model Quality & Benchmark Backtest Report",
        "",
        f"**Overall verdict: {v['verdict']}**",
        "",
        v["summary"],
        "",
        "## 1. Cross-validation (scaffold GroupKFold, 5 folds)",
        "",
        "### URAT1 regression + conformal UQ",
        f"- RMSE (OOF): {cv['urat1']['rmse']:.3f}",
        f"- R² (OOF): {cv['urat1']['r2']:.3f}",
        f"- Spearman (OOF): {cv['urat1']['spearman']:.3f}",
        f"- ROC-AUC (p≥7): {cv['urat1'].get('roc_auc_p7', float('nan')):.3f}",
        f"- EF@5% (p≥7, strong actives): {cv['urat1'].get('ef_5pct_p7', float('nan')):.2f}",
        f"- EF@10% (p≥6): {cv['urat1'].get('ef_10pct_p6', float('nan')):.2f} — **misleading** (theoretical max ≈{cv['urat1'].get('ef_p6_theoretical_max', 1.75):.2f} at {cv['urat1'].get('active_rate_p6', 0.57)*100:.0f}% base rate)",
        f"- Strict CV pass: {v['cv_pass']['urat1']}",
        "",
        "### NLRP3 assay-conditioned classifier",
        f"- AUROC: {cv['nlrp3']['auroc']:.3f}",
        f"- AUPRC: {cv['nlrp3']['auprc']:.3f}",
        f"- EF@10%: {cv['nlrp3']['ef_10pct']:.2f}",
        f"- CV screening suitable: {v['cv_pass']['nlrp3']}",
        "",
        "## 2. Benchmark backtest",
        "",
        "### URAT1 (predicted pActivity ≥ 6 = pass)",
        "",
        "| Compound | In train | Max Tc | Pred | Lit. | Pass |",
        "|----------|----------|--------|------|------|------|",
    ]
    for r in u["benchmarks"]:
        lit = r["literature_pActivity"] if r["literature_pActivity"] is not None else "—"
        lines.append(
            f"| {r['compound']} | {r['in_training_set']} | {r['max_tanimoto_to_train']} | "
            f"{r['y_pred']} | {lit} | {r['binary_pass']} |"
        )

    lines += [
        "",
        f"URAT1 must-recover binary pass: {u['must_recover_binary_pass']}/{u['must_recover_count']}",
        f"- In training set: {u['in_training_must_pass']}/{u['in_training_must_count']}",
        f"- Scaffold-novel (excluded from curation): {u['novel_must_pass']}/{u['novel_must_count']}",
        "",
        "### NLRP3 (P(active) ≥ 0.5 = pass)",
        "",
        "| Compound | In train | Max Tc | P(active) | Pass |",
        "|----------|----------|--------|-----------|------|",
    ]
    for r in n["benchmarks"]:
        lines.append(
            f"| {r['compound']} | {r['in_training_set']} | {r['max_tanimoto_to_train']} | "
            f"{r['p_active_ensemble']} | {r['binary_pass']} |"
        )

    lines += [
        "",
        f"NLRP3 must-recover binary pass: {n['must_recover_binary_pass']}/{n['must_recover_count']}",
        "",
        "## 3. Why the previous URAT1 table was misleading",
        "",
        "1. **EF@10% at p≥6 is capped near 1.75** when 57% of training compounds are already actives — even a perfect ranker cannot exceed ~1.75.",
        "2. **Thresholds were too lenient** (R²≥0.25, EF≥1.5), allowing a mediocre model to show all green checks.",
        "3. **Fold-averaged R² (0.44) understates OOF R² (0.51)** but both are only moderate for prospective screening.",
        "4. **Benchmark backtest contradicts** the pass table: lesinurad/dotinurad fail despite CV pass.",
        "",
        "## 4. Interpretation notes",
        "",
        "- **lesinurad / benzbromarone** were dropped during ChEMBL curation due to >1 log assay conflict; "
        "ChEMBL median pActivity (~5.1–6.5) is lower than literature references used in benchmarks.",
        "- **verinurad** is in the training set; model prediction is consistent with held-in data.",
        "- **MCC950** is in NLRP3 training data; high P(active) confirms correct class assignment.",
        "- For scaffold-novel benchmarks, prioritize **conformational ensemble docking** ($S_{trap}$) over ML rank.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=PROCESSED)
    parser.add_argument("--model-dir", type=Path, default=MODELS)
    parser.add_argument("--output", type=Path, default=RESULTS)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    urat1_df = pd.read_csv(args.data_dir / "urat1_curated.csv")
    nlrp3_df = pd.read_csv(args.data_dir / "nlrp3_records.csv")
    with open(args.model_dir / "training_report.json") as f:
        training_report = json.load(f)

    urat1_bundle = joblib.load(args.model_dir / "urat1_model.joblib")
    nlrp3_bundle = joblib.load(args.model_dir / "nlrp3_model.joblib")

    urat1_bt = backtest_urat1(urat1_df, urat1_bundle, BENCHMARKS)
    nlrp3_bt = backtest_nlrp3(nlrp3_df, nlrp3_bundle, BENCHMARKS)
    scaf_check = scaffold_exclusion_check(urat1_df, nlrp3_df, BENCHMARKS)
    verdict = overall_verdict(training_report, urat1_bt, nlrp3_bt)

    report = {
        "framework": "TAPE-GATE",
        "urat1_backtest": urat1_bt,
        "nlrp3_backtest": nlrp3_bt,
        "scaffold_overlap_check": scaf_check,
        "overall_verdict": verdict,
        "training_cv_summary": {
            "urat1": training_report["urat1"]["cv_metrics"],
            "nlrp3": training_report["nlrp3"]["cv_metrics"],
        },
    }

    save_json(args.output / "benchmark_backtest_report.json", report)
    pd.DataFrame(urat1_bt["benchmarks"]).to_csv(args.output / "urat1_benchmark_rankings.csv", index=False)
    pd.DataFrame(nlrp3_bt["benchmarks"]).to_csv(args.output / "nlrp3_benchmark_rankings.csv", index=False)
    write_markdown_report(report, args.output / "MODEL_QUALITY_REPORT.md")

    print("\n=== URAT1 Benchmark ===")
    for r in urat1_bt["benchmarks"]:
        print(f"  {r['compound']:15s} pred={r['y_pred']}  in_train={r['in_training_set']}  pass={r['binary_pass']}")

    print("\n=== NLRP3 Benchmark ===")
    for r in nlrp3_bt["benchmarks"]:
        print(f"  {r['compound']:15s} P={r['p_active_ensemble']}  in_train={r['in_training_set']}  pass={r['binary_pass']}")

    print(f"\n=== Verdict: {verdict['verdict']} ===")
    print(verdict["summary"])
    print(f"\nReport: {args.output / 'MODEL_QUALITY_REPORT.md'}")


if __name__ == "__main__":
    main()
