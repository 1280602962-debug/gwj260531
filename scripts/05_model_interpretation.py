#!/usr/bin/env python3
"""
Phase 5: Model interpretability — SHAP analysis for selectivity model.

Supports:
  - TreeExplainer for XGBoost delta / MTL models
  - Global summary (beeswarm) + dependence plots
  - Per-compound waterfall for benchmark molecules
  - Morgan bit → substructure mapping

Usage:
    python scripts/05_model_interpretation.py \
        --model models/best_model.joblib \
        --data data/processed/paired_set.csv \
        --output results/shap
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import FIGSIZE_DOUBLE, FIGSIZE_SINGLE, apply_journal_style, save_figure  # noqa: E402
from utils_ml import featurize_smiles, feature_names  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def get_feature_spec(bundle: dict) -> tuple[int, list[str]]:
    spec = bundle.get("feature_spec", {}) if isinstance(bundle, dict) else {}
    morgan_bits = int(spec.get("morgan_bits", 2048))
    return morgan_bits, feature_names(morgan_bits)


def bit_to_substructure(smiles: str, bit_id: int, morgan_bits: int = 2048, radius: int = 2) -> str | None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    info = {}
    AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=morgan_bits, bitInfo=info)
    if bit_id not in info:
        return None
    atom_id, r = info[bit_id][0]
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, r, atom_id)
    submol = Chem.PathToSubmol(mol, env)
    if submol:
        return Chem.MolToSmiles(submol)
    return None


def plot_global_summary(shap_values: np.ndarray, X: np.ndarray, names: list[str], output_dir: Path, title: str):
    import shap

    apply_journal_style()
    plt.figure(figsize=FIGSIZE_DOUBLE)
    shap.summary_plot(shap_values, X, feature_names=names, show=False, max_display=25)
    plt.title(title)
    save_figure(output_dir / "shap_summary_beeswarm.png")


def plot_dependence(shap_values: np.ndarray, X: np.ndarray, top_features: list[str], names: list[str], output_dir: Path):
    import shap

    for feat in top_features[:5]:
        apply_journal_style()
        plt.figure(figsize=FIGSIZE_SINGLE)
        shap.dependence_plot(feat, shap_values, X, feature_names=names, show=False)
        plt.title(f"SHAP Dependence: {feat}")
        safe = feat.replace("/", "_")
        save_figure(output_dir / f"shap_dependence_{safe}.png")


def plot_waterfall(explainer, X_row: np.ndarray, names: list[str], output_dir: Path, name: str):
    import shap

    sv = explainer.shap_values(X_row.reshape(1, -1))
    if isinstance(sv, list):
        sv = sv[0]
    base = explainer.expected_value
    if isinstance(base, (list, np.ndarray)):
        base = base[0]
    apply_journal_style()
    plt.figure(figsize=FIGSIZE_DOUBLE)
    shap.waterfall_plot(
        shap.Explanation(values=sv[0], base_values=base, data=X_row, feature_names=names),
        show=False,
        max_display=15,
    )
    plt.title(f"SHAP Waterfall: {name}")
    save_figure(output_dir / f"shap_waterfall_{name.replace(' ', '_')}.png")


def top_substructure_report(
    shap_values: np.ndarray,
    X: np.ndarray,
    smiles_list: list[str],
    output_dir: Path,
    morgan_bits: int,
    top_k: int = 20,
):
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    top_idx = np.argsort(mean_abs_shap)[::-1][:top_k]

    rows = []
    top_feature_names = []
    for idx in top_idx:
        fname = feature_names(morgan_bits)[idx]
        top_feature_names.append(fname)
        if not fname.startswith("Bit_"):
            rows.append(
                {
                    "feature": fname,
                    "mean_abs_shap": float(mean_abs_shap[idx]),
                    "mean_shap": float(shap_values[:, idx].mean()),
                    "direction": "promotes" if shap_values[:, idx].mean() > 0 else "inhibits",
                    "example_substructure": None,
                    "example_smiles": smiles_list[0],
                }
            )
            continue
        bit = int(fname.split("_")[1])
        on_idx = np.where(X[:, idx] == 1)[0]
        example_smi = smiles_list[on_idx[0]] if len(on_idx) else smiles_list[0]
        sub = bit_to_substructure(example_smi, bit, morgan_bits)
        rows.append(
            {
                "feature": fname,
                "mean_abs_shap": float(mean_abs_shap[idx]),
                "mean_shap": float(shap_values[:, idx].mean()),
                "direction": "promotes" if shap_values[:, idx].mean() > 0 else "inhibits",
                "example_substructure": sub,
                "example_smiles": example_smi,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(output_dir / "top_shap_features.csv", index=False)
    return df, top_feature_names


def compare_classes_shap(shap_values, y_class, output_dir: Path, morgan_bits: int):
    sel_mask = y_class == 1
    if sel_mask.sum() < 2 or (~sel_mask).sum() < 2:
        return
    mean_sel = shap_values[sel_mask].mean(axis=0)
    mean_non = shap_values[~sel_mask].mean(axis=0)
    diff = mean_sel - mean_non
    top_diff = np.argsort(np.abs(diff))[::-1][:20]
    names = feature_names(morgan_bits)
    rows = [{"feature": names[b], "shap_diff_sel_vs_non": float(diff[b])} for b in top_diff]
    pd.DataFrame(rows).to_csv(output_dir / "shap_class_comparison.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description="SHAP interpretability analysis")
    parser.add_argument("--model", type=Path, default=ROOT / "models" / "best_model.joblib")
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "processed" / "paired_set.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "shap")
    parser.add_argument("--sample-n", type=int, default=500, help="Background sample size for SHAP")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    if not args.model.exists():
        raise SystemExit(f"Model not found: {args.model}. Run 04_train_selectivity_model.py first.")
    if not args.data.exists():
        raise SystemExit(f"Data not found: {args.data}")

    bundle = joblib.load(args.model)
    delta_model = bundle.get("delta") if isinstance(bundle, dict) else bundle
    if delta_model is None:
        delta_model = bundle.get("mtl") if isinstance(bundle, dict) else bundle

    morgan_bits, names = get_feature_spec(bundle)

    df = pd.read_csv(args.data).dropna(subset=["canonical_smiles"])
    if "n_isoforms" in df.columns:
        df = df[df["n_isoforms"] >= 2]
    if len(df) > args.sample_n:
        df = df.sample(args.sample_n, random_state=42)

    smiles = df["canonical_smiles"].tolist()
    X = featurize_smiles(smiles, morgan_bits=morgan_bits)

    import shap

    logger.info("Computing SHAP values for %d compounds ...", len(smiles))
    explainer = shap.TreeExplainer(delta_model)
    shap_values = explainer.shap_values(X)
    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    np.save(args.output / "shap_values.npy", shap_values)
    np.save(args.output / "feature_matrix.npy", X)

    plot_global_summary(
        shap_values,
        X,
        names,
        args.output,
        "SHAP: JNK1 Selectivity (Δmin Model)",
    )
    shap_df, top_features = top_substructure_report(shap_values, X, smiles, args.output, morgan_bits)
    plot_dependence(shap_values, X, top_features, names, args.output)

    if "sel_class" in df.columns:
        y_class = (df["sel_class"] == "JNK1-selective").astype(int).values
        compare_classes_shap(shap_values, y_class, args.output, morgan_bits)

        sel = df[df["sel_class"] == "JNK1-selective"]
        if len(sel):
            idx = df.index.get_loc(sel.index[0])
            plot_waterfall(explainer, X[idx], names, args.output, "JNK1_selective_example")

    report = {
        "n_compounds": len(smiles),
        "top_10_features": shap_df.head(10).to_dict(orient="records"),
        "interpretation_guide": (
            "Positive mean SHAP → feature increases predicted JNK1 selectivity (Δmin). "
            "Cross-reference Morgan bits with MMP rules and isoform structural differences."
        ),
    }
    with open(args.output / "shap_report.json", "w") as f:
        json.dump(report, f, indent=2)

    logger.info("SHAP analysis complete → %s", args.output)


if __name__ == "__main__":
    main()
