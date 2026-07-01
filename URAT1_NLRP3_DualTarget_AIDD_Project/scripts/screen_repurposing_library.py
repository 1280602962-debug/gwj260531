#!/usr/bin/env python3
"""
NLRP3 (and optional URAT1) ML screening on ChEMBL repurposing manifest.

Typical workflow:
  1. NLRP3 assay-conditioned ensemble on full clinical library (fast)
  2. Export top-N / top-percentile for dual-target docking follow-up

Example:
  python3 scripts/screen_repurposing_library.py \\
    --input /path/to/repurposing_manifest.csv \\
    --panel clinical_all \\
    --top-n 150 --top-pct 5
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from utils_ml import canonicalize, featurize_smiles, max_tanimoto_to_library, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS = PROJECT_ROOT / "results" / "training"
PROCESSED = PROJECT_ROOT / "data" / "processed"
DEFAULT_OUT = PROJECT_ROOT / "results" / "repurposing"


def _load_assay_ids(bundle: dict, nlrp3_df: pd.DataFrame, n_assays: int) -> list[str]:
    if "top_assay_ids" in bundle and bundle["top_assay_ids"]:
        return list(bundle["top_assay_ids"][:n_assays])
    return nlrp3_df["Assay ChEMBL ID"].value_counts().head(n_assays).index.astype(str).tolist()


def predict_nlrp3_ensemble(bundle: dict, smiles_list: list[str], assay_ids: list[str]) -> np.ndarray:
    from utils_ml import assay_one_hot_matrix

    x_mol = featurize_smiles(smiles_list)
    probs = []
    for aid in assay_ids:
        assay_col = pd.Series([aid] * len(smiles_list))
        x_assay = assay_one_hot_matrix(assay_col, bundle["top_assay_ids"])
        x = np.hstack([x_mol, x_assay])
        x_s = bundle["scaler"].transform(x)
        raw = bundle["model"].predict_proba(x_s)[:, 1]
        if bundle.get("calibrator") is not None:
            raw = bundle["calibrator"].predict(raw)
        probs.append(raw)
    return np.max(np.vstack(probs), axis=0)


def predict_urat1(bundle: dict, smiles_list: list[str]) -> np.ndarray:
    x = featurize_smiles(smiles_list)
    x_s = bundle["scaler"].transform(x)
    return bundle["model"].predict(x_s)


def select_panel(df: pd.DataFrame, panel: str) -> pd.DataFrame:
    if panel == "clinical_all":
        return df.copy()
    if panel == "phase_ge3":
        return df[pd.to_numeric(df["max_phase"], errors="coerce") >= 3].copy()
    if panel == "primary_atc_phase":
        return df[df["library_panel"] == "primary_atc_phase"].copy()
    if panel == "phase_only":
        return df[df["library_panel"] == "phase_only"].copy()
    raise ValueError(f"Unknown panel: {panel}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Screen repurposing library with NLRP3 ML")
    parser.add_argument("--input", type=Path, required=True, help="repurposing_manifest.csv")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--panel",
        choices=["clinical_all", "phase_ge3", "primary_atc_phase", "phase_only"],
        default="clinical_all",
        help="Which library subset to score (default: full clinical manifest)",
    )
    parser.add_argument("--top-n", type=int, default=150, help="Top N for docking shortlist")
    parser.add_argument("--top-pct", type=float, default=5.0, help="Also include top X percentile")
    parser.add_argument("--nlrp3-threshold", type=float, default=0.5, help="Binary active cutoff")
    parser.add_argument("--with-urat1-ml", action="store_true", help="Also score URAT1 ML (auxiliary)")
    parser.add_argument("--n-ensemble-assays", type=int, default=5)
    parser.add_argument("--skip-tanimoto", action="store_true", help="Skip slow applicability-domain Tanimoto")
    parser.add_argument(
        "--export-p05-pool",
        action="store_true",
        help="Export all compounds with P(active)>=threshold as docking_pool_p05.csv (main workflow)",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.input, low_memory=False)
    smi_col = "canonical_smiles" if "canonical_smiles" in raw.columns else "smiles"
    panel_df = select_panel(raw, args.panel)
    panel_df = panel_df.dropna(subset=[smi_col]).copy()
    panel_df["canonical_smiles"] = panel_df[smi_col].map(
        lambda s: canonicalize(s) if pd.notna(s) else None
    )
    panel_df = panel_df[panel_df["canonical_smiles"].notna()].drop_duplicates("canonical_smiles")

    nlrp3_bundle = joblib.load(MODELS / "nlrp3_model.joblib")
    nlrp3_train = pd.read_csv(PROCESSED / "nlrp3_records.csv")
    train_smiles = nlrp3_train["canonical_smiles"].unique().tolist()
    assay_ids = _load_assay_ids(nlrp3_bundle, nlrp3_train, args.n_ensemble_assays)

    smiles = panel_df["canonical_smiles"].tolist()
    print(f"Scoring {len(smiles)} compounds (panel={args.panel})...")
    p_active = predict_nlrp3_ensemble(nlrp3_bundle, smiles, assay_ids)

    out = panel_df.copy()
    out["p_active_nlrp3"] = p_active
    out["nlrp3_percentile"] = out["p_active_nlrp3"].rank(method="average", pct=True) * 100.0
    out["nlrp3_pred_active"] = out["p_active_nlrp3"] >= args.nlrp3_threshold
    if not args.skip_tanimoto:
        out["max_tanimoto_nlrp3_train"] = [
            round(max_tanimoto_to_library(s, train_smiles), 3) for s in smiles
        ]
    else:
        out["max_tanimoto_nlrp3_train"] = np.nan
    out["scaffold"] = out["canonical_smiles"].map(murcko_scaffold)

    if args.with_urat1_ml:
        urat1_bundle = joblib.load(MODELS / "urat1_model.joblib")
        out["urat1_ml_pactivity"] = predict_urat1(urat1_bundle, smiles)
        out["urat1_ml_percentile"] = out["urat1_ml_pactivity"].rank(method="average", pct=True) * 100.0

    out = out.sort_values("p_active_nlrp3", ascending=False).reset_index(drop=True)
    out["rank_nlrp3"] = np.arange(1, len(out) + 1)

    scores_path = args.output_dir / f"nlrp3_ml_scores_{args.panel}.csv"
    out.to_csv(scores_path, index=False)

    n_top = max(args.top_n, int(np.ceil(len(out) * args.top_pct / 100.0)))
    docking_pool = out.head(n_top).copy()
    docking_pool["selection_reason"] = "nlrp3_ml_top"
    dock_path = args.output_dir / f"nlrp3_top_for_dual_docking_{args.panel}.csv"
    docking_pool.to_csv(dock_path, index=False)

    p05_path = None
    if args.export_p05_pool:
        p05_pool = out[out["nlrp3_pred_active"]].copy()
        p05_pool["selection_reason"] = "nlrp3_p_active_ge_threshold"
        p05_path = args.output_dir / (
            "docking_pool_p05.csv"
            if args.panel == "clinical_all"
            else f"docking_pool_p05_{args.panel}.csv"
        )
        p05_pool.to_csv(p05_path, index=False)

    # Known controls in library
    controls = ["lesinurad", "benzbromarone", "verinurad", "dotinurad", "colchicine", "allopurinol", "febuxostat"]
    name_col = "name" if "name" in out.columns else "pref_name"
    ctrl_rows = []
    for c in controls:
        if name_col not in out.columns:
            break
        hit = out[out[name_col].astype(str).str.upper() == c.upper()]
        if len(hit):
            ctrl_rows.append(
                {
                    "name": c,
                    "rank_nlrp3": int(hit["rank_nlrp3"].iloc[0]),
                    "p_active_nlrp3": float(hit["p_active_nlrp3"].iloc[0]),
                    "nlrp3_percentile": float(hit["nlrp3_percentile"].iloc[0]),
                }
            )

    summary = {
        "input": str(args.input),
        "panel": args.panel,
        "n_scored": int(len(out)),
        "n_pred_active_ge_threshold": int(out["nlrp3_pred_active"].sum()),
        "threshold": args.nlrp3_threshold,
        "ensemble_assays": assay_ids,
        "top_n_export": int(len(docking_pool)),
        "top_pct": args.top_pct,
        "known_controls": ctrl_rows,
        "outputs": {
            "full_scores": str(scores_path),
            "docking_shortlist": str(dock_path),
            **({"docking_pool_p05": str(p05_path)} if p05_path else {}),
        },
        "next_step": (
            "Dock docking_pool_p05.csv (P>=0.5) at URAT1 9DKB XP and NLRP3 8ETR XP; "
            "integrate with scripts/merge_docking_pareto.py."
        ),
    }
    summary_path = args.output_dir / f"nlrp3_screening_summary_{args.panel}.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"\nTop 10 NLRP3 ML hits:")
    cols = [name_col, "chembl_id", "p_active_nlrp3", "nlrp3_percentile", "max_phase", "library_panel"]
    cols = [c for c in cols if c in out.columns]
    print(out[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
