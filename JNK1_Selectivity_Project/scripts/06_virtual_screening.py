#!/usr/bin/env python3
"""
Phase 4: Virtual screening for JNK1-selective hits.

Funnel:
  F1 preprocess → F2 drug-like → F3 JNK1 activity → F4 selectivity
  → F5 SA/QED → F6 diversity selection

Usage:
    python scripts/06_virtual_screening.py \
        --model models/best_model.joblib \
        --library data/libraries/screening_demo.smi \
        --output results/screening
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import FIGSIZE_DOUBLE, apply_journal_style, save_figure  # noqa: E402
from utils_ml import featurize_smiles  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CONFIG_PATH = ROOT / "config" / "targets.yaml"
DEFAULT_LIBRARY = ROOT / "data" / "libraries" / "screening_demo.smi"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def read_smiles_library(path: Path) -> list[str]:
    smiles = []
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi = line.split()[0] if " " in line else line.split(",")[0]
            smiles.append(smi)
    return smiles


def preprocess_smiles(smiles: str) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def lipinski_filter(smiles: str) -> bool:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    return 200 <= mw <= 600 and -1 <= logp <= 5 and hbd <= 5 and hba <= 10


def compute_sa_qed(smiles: str) -> tuple[float, float]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 10.0, 0.0
    qed = QED.qed(mol)
    try:
        from rdkit.Contrib.SA_Score import sascorer

        sa = float(sascorer.calculateScore(mol))
    except Exception:
        # Fallback: map Bertz complexity to ~1–10 (lower = easier to synthesize)
        sa = min(10.0, max(1.0, Descriptors.BertzCT(mol) / 150.0))
    return sa, qed


def get_morgan_bits(bundle: dict) -> int:
    if isinstance(bundle, dict) and "feature_spec" in bundle:
        return int(bundle["feature_spec"].get("morgan_bits", 2048))
    return 2048


def predict_batch(bundle: dict, X: np.ndarray) -> dict[str, np.ndarray]:
    mtl = bundle.get("mtl")
    delta = bundle.get("delta")
    clf = bundle.get("classifier")
    single = bundle.get("single_target") or {}

    out: dict[str, np.ndarray] = {}

    if single and all(k in single for k in ("JNK1", "JNK2", "JNK3")):
        out["pred_pAct_JNK1"] = single["JNK1"].predict(X)
        out["pred_pAct_JNK2"] = single["JNK2"].predict(X)
        out["pred_pAct_JNK3"] = single["JNK3"].predict(X)
        out["pred_delta_min_computed"] = out["pred_pAct_JNK1"] - np.nanmax(
            np.vstack([out["pred_pAct_JNK2"], out["pred_pAct_JNK3"]]),
            axis=0,
        )
    elif mtl is not None:
        preds = mtl.predict(X)
        out["pred_pAct_JNK1"] = preds[:, 0]
        out["pred_pAct_JNK2"] = preds[:, 1]
        out["pred_pAct_JNK3"] = preds[:, 2]
        out["pred_delta_min_computed"] = preds[:, 0] - np.nanmax(preds[:, 1:3], axis=1)

    if delta is not None:
        out["pred_delta_min"] = delta.predict(X)
    if clf is not None:
        out["prob_JNK1_selective"] = clf.predict_proba(X)[:, 1]
    return out


def compute_final_score(row: pd.Series, weights: dict) -> float:
    score = 0.0
    if "pred_pAct_JNK1" in row and not pd.isna(row["pred_pAct_JNK1"]):
        score += weights["w_pAct_JNK1"] * row["pred_pAct_JNK1"] / 10.0
    delta = row.get("pred_delta_min", row.get("pred_delta_min_computed", 0))
    if not pd.isna(delta):
        score += weights["w_delta_min"] * max(0, delta) / 3.0
    jnk23 = max(row.get("pred_pAct_JNK2", 0), row.get("pred_pAct_JNK3", 0))
    if not pd.isna(jnk23):
        score -= weights["w_neg_JNK23"] * max(0, jnk23 - 5.0) / 5.0
    if "qed" in row:
        score += weights["w_qed"] * row["qed"]
    if "sa" in row:
        score += weights["w_sa"] * row["sa"] / 10.0
    return score


def butina_diverse_selection(df: pd.DataFrame, n: int = 100, cutoff: float = 0.7) -> pd.DataFrame:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.ML.Cluster import Butina

    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols if m]
    dists = []
    nf = len(fps)
    for i in range(1, nf):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    clusters = Butina.ClusterData(dists, nf, 1 - cutoff, isDistData=True)
    selected_idx = [cluster[0] for cluster in clusters]
    ranked = df.iloc[selected_idx].sort_values("final_score", ascending=False).head(n)
    return ranked


def plot_funnel(stats: dict, output_dir: Path) -> None:
    apply_journal_style()
    stages = ["Input", "Preprocessed", "Drug-like", "JNK1 Active", "Selective"]
    keys = ["input", "preprocessed", "druglike", "active", "selective"]
    values = [stats[k] for k in keys]
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    bars = ax.bar(stages, values, color="#55A868", edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Compound Count")
    ax.set_title("Virtual Screening Funnel")
    ax.tick_params(axis="x", rotation=20)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{val:,}", ha="center", va="bottom", fontsize=7)
    save_figure(output_dir / "screening_funnel.png", fig)


def plot_score_distribution(hits: pd.DataFrame, output_dir: Path) -> None:
    if hits.empty or "final_score" not in hits.columns:
        return
    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    ax.hist(hits["final_score"], bins=30, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Composite Screening Score")
    ax.set_ylabel("Hit Count")
    ax.set_title("Distribution of Final Screening Scores")
    save_figure(output_dir / "score_distribution.png", fig)


def annotate_benchmarks(hits: pd.DataFrame, config: dict, processed_dir: Path) -> pd.DataFrame:
    """Tag known benchmark compounds if present in curated data."""
    if hits.empty:
        return hits
    bench_rows = []
    for bench in config.get("benchmarks", []):
        chembl_id = bench.get("chembl_id")
        for iso in ["jnk1", "jnk2", "jnk3"]:
            path = processed_dir / f"{iso}_curated.csv"
            if not path.exists():
                continue
            sub = pd.read_csv(path)
            row = sub[sub["molecule_chembl_id"] == chembl_id]
            if len(row):
                smi = row["canonical_smiles"].iloc[0]
                match = hits[hits["smiles"] == smi]
                if len(match):
                    bench_rows.append(
                        {
                            "name": bench["name"],
                            "role": bench["role"],
                            "smiles": smi,
                            "final_score": float(match["final_score"].iloc[0]),
                            "pred_pAct_JNK1": float(match.get("pred_pAct_JNK1", pd.Series([np.nan])).iloc[0]),
                            "pred_delta_min": float(
                                match.get("pred_delta_min", match.get("pred_delta_min_computed", pd.Series([np.nan]))).iloc[0]
                            ),
                        }
                    )
                break
    if bench_rows:
        pd.DataFrame(bench_rows).to_csv(process_dir.parent / "screening" / "benchmark_validation.csv", index=False)
    return hits


def main():
    parser = argparse.ArgumentParser(description="Virtual screening for JNK1-selective inhibitors")
    parser.add_argument("--model", type=Path, default=ROOT / "models" / "best_model.joblib")
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY, help="SMILES file (.smi/.csv/.smi.gz)")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "screening")
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument("--top-n", type=int, default=500)
    parser.add_argument("--diverse-n", type=int, default=100)
    parser.add_argument("--processed", type=Path, default=ROOT / "data" / "processed")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config = load_config()
    thresholds = config["selectivity"]
    if not args.model.exists():
        raise SystemExit(f"Model not found: {args.model}")

    bundle = joblib.load(args.model)
    weights = bundle.get(
        "scoring_weights",
        {
            "w_pAct_JNK1": 0.35,
            "w_delta_min": 0.30,
            "w_neg_JNK23": 0.20,
            "w_qed": 0.10,
            "w_sa": -0.05,
        },
    )
    morgan_bits = get_morgan_bits(bundle)

    if not args.library.exists():
        raise SystemExit(f"Library not found: {args.library}. Run run_selectivity_pipeline.py first.")

    logger.info("Reading library: %s", args.library)
    all_smiles = read_smiles_library(args.library)
    logger.info("Total input SMILES: %d", len(all_smiles))

    hits = []
    stats = {"input": len(all_smiles), "preprocessed": 0, "druglike": 0, "active": 0, "selective": 0}

    for start in tqdm(range(0, len(all_smiles), args.batch_size), desc="Screening batches"):
        batch_smi = all_smiles[start : start + args.batch_size]

        clean = []
        for s in batch_smi:
            cs = preprocess_smiles(s)
            if cs:
                clean.append(cs)
        stats["preprocessed"] += len(clean)

        druglike = [s for s in clean if lipinski_filter(s)]
        stats["druglike"] += len(druglike)
        if not druglike:
            continue

        X = featurize_smiles(druglike, morgan_bits=morgan_bits)
        preds = predict_batch(bundle, X)
        batch_df = pd.DataFrame({"smiles": druglike, **preds})

        jnk1_col = "pred_pAct_JNK1"
        if jnk1_col in batch_df.columns:
            batch_df = batch_df[batch_df[jnk1_col] >= 7.0]
        stats["active"] += len(batch_df)

        delta_col = "pred_delta_min" if "pred_delta_min" in batch_df.columns else "pred_delta_min_computed"
        if delta_col in batch_df.columns:
            batch_df = batch_df[batch_df[delta_col] >= thresholds["delta_log_threshold"]]
        if "pred_pAct_JNK2" in batch_df.columns:
            batch_df = batch_df[
                batch_df["pred_pAct_JNK2"] < batch_df.get("pred_pAct_JNK1", 99) - thresholds["delta_log_threshold"] * 0.5
            ]
        if "pred_pAct_JNK3" in batch_df.columns:
            batch_df = batch_df[
                batch_df["pred_pAct_JNK3"] < batch_df.get("pred_pAct_JNK1", 99) - thresholds["delta_log_threshold"] * 0.5
            ]
        stats["selective"] += len(batch_df)

        sa_qed = [compute_sa_qed(s) for s in batch_df["smiles"]]
        batch_df["sa"] = [x[0] for x in sa_qed]
        batch_df["qed"] = [x[1] for x in sa_qed]
        batch_df = batch_df[(batch_df["sa"] <= 6.0) & (batch_df["qed"] >= 0.35)]

        batch_df["final_score"] = batch_df.apply(lambda r: compute_final_score(r, weights), axis=1)
        hits.append(batch_df)

    if not hits:
        logger.warning("No hits passed the funnel. Check model and thresholds.")
        plot_funnel(stats, args.output)
        with open(args.output / "screening_report.json", "w") as f:
            json.dump({"funnel_stats": stats, "thresholds": thresholds}, f, indent=2)
        return

    all_hits = pd.concat(hits, ignore_index=True).sort_values("final_score", ascending=False)
    all_hits.to_csv(args.output / "all_hits.csv", index=False)

    top = all_hits.head(args.top_n)
    top.to_csv(args.output / f"top{args.top_n}.csv", index=False)

    diverse = butina_diverse_selection(top, n=args.diverse_n)
    diverse.to_csv(args.output / f"top{args.diverse_n}_diverse.csv", index=False)

    plot_funnel(stats, args.output)
    plot_score_distribution(all_hits, args.output)

    bench_rows = []
    for bench in config.get("benchmarks", []):
        chembl_id = bench.get("chembl_id")
        for iso in ["jnk1", "jnk2", "jnk3"]:
            path = args.processed / f"{iso}_curated.csv"
            if not path.exists():
                continue
            sub = pd.read_csv(path)
            row = sub[sub["molecule_chembl_id"] == chembl_id]
            if len(row):
                smi = row["canonical_smiles"].iloc[0]
                match = all_hits[all_hits["smiles"] == smi]
                bench_rows.append(
                    {
                        "name": bench["name"],
                        "role": bench["role"],
                        "smiles": smi,
                        "in_hits": bool(len(match)),
                        "final_score": float(match["final_score"].iloc[0]) if len(match) else np.nan,
                        "pred_pAct_JNK1": float(match["pred_pAct_JNK1"].iloc[0]) if len(match) else np.nan,
                        "pred_delta_min": float(
                            match.get("pred_delta_min", match.get("pred_delta_min_computed", pd.Series([np.nan]))).iloc[0]
                        )
                        if len(match)
                        else np.nan,
                    }
                )
                break
    if bench_rows:
        pd.DataFrame(bench_rows).to_csv(args.output / "benchmark_validation.csv", index=False)

    report = {
        "funnel_stats": stats,
        "thresholds": thresholds,
        "weights": weights,
        "top10": top.head(10)[["smiles", "final_score"]].to_dict(orient="records"),
        "benchmarks": bench_rows,
    }
    with open(args.output / "screening_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info("Screening complete: %d hits → Top %d → Diverse %d", len(all_hits), len(top), len(diverse))
    logger.info("Results → %s", args.output)


if __name__ == "__main__":
    main()
