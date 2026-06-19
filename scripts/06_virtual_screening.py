#!/usr/bin/env python3
"""
Phase 4: Virtual screening — F1 ML pre-filter + F2 drug-like + ranking.

Funnel (v2, benchmark-calibrated):
  F0 preprocess → F2 drug-like → F1 p_family ≥ threshold → F5 SA/QED → rank

Isoform selectivity is NOT filtered by ML delta (use F3 docking downstream).

Usage:
    python3 scripts/06_virtual_screening.py \
        --models-dir models/xgboost \
        --library data/libraries/screening_demo.smi \
        --output results/screening_v2
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
DEFAULT_MODELS_DIR = ROOT / "models" / "xgboost"
DEFAULT_BENCHMARKS = ROOT / "data" / "benchmarks" / "literature_benchmarks.csv"
ISOFORMS = ["JNK1", "JNK2", "JNK3"]


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
        sa = min(10.0, max(1.0, Descriptors.BertzCT(mol) / 150.0))
    return sa, qed


def load_xgboost_models(model_dir: Path) -> dict[str, object]:
    models = {}
    for iso in ISOFORMS:
        path = model_dir / f"xgboost_{iso.lower()}.joblib"
        if not path.exists():
            raise FileNotFoundError(f"Missing model: {path}. Run scripts/07_compare_models.py first.")
        models[iso] = joblib.load(path)
    return models


def predict_batch(models: dict[str, object], smiles: list[str]) -> pd.DataFrame:
    X = featurize_smiles(smiles)
    out = {"smiles": smiles}
    for iso in ISOFORMS:
        out[f"pred_pAct_{iso}"] = models[iso].predict(X)
    df = pd.DataFrame(out)
    df["p_family"] = df[[f"pred_pAct_{iso}" for iso in ISOFORMS]].max(axis=1)
    df["pred_delta_min_computed"] = df["pred_pAct_JNK1"] - df[["pred_pAct_JNK2", "pred_pAct_JNK3"]].max(axis=1)
    return df


def compute_final_score(row: pd.Series, weights: dict) -> float:
    score = weights.get("w_p_family", 0.55) * row["p_family"] / 10.0
    if "pred_pAct_JNK1" in row and not pd.isna(row["pred_pAct_JNK1"]):
        score += weights.get("w_pAct_JNK1", 0.15) * row["pred_pAct_JNK1"] / 10.0
    if "qed" in row:
        score += weights.get("w_qed", 0.20) * row["qed"]
    if "sa" in row:
        score += weights.get("w_sa", 0.10) * (10.0 - row["sa"]) / 10.0
    return float(score)


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
    return df.iloc[selected_idx].sort_values("final_score", ascending=False).head(n)


def plot_funnel(stats: dict, output_dir: Path) -> None:
    apply_journal_style()
    stages = ["Input", "Preprocessed", "Drug-like", "F1 p_family", "SA/QED"]
    keys = ["input", "preprocessed", "druglike", "f1_pass", "sa_qed_pass"]
    values = [stats[k] for k in keys]
    fig, ax = plt.subplots(figsize=FIGSIZE_DOUBLE)
    bars = ax.bar(stages, values, color="#55A868", edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Compound Count")
    ax.set_title("Virtual Screening Funnel (F1 p_family + F2 drug-like)")
    ax.tick_params(axis="x", rotation=25)
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


def validate_literature_benchmarks(
    all_smiles: set[str],
    hits: pd.DataFrame,
    benchmarks_path: Path,
    models: dict[str, object],
    p_family_threshold: float,
) -> list[dict]:
    if not benchmarks_path.exists():
        return []

    bench = pd.read_csv(benchmarks_path)
    preds = predict_batch(models, bench["smiles"].tolist())
    rows = []
    for i, b in bench.iterrows():
        smi = bench.loc[i, "smiles"]
        in_library = smi in all_smiles
        p_family = float(preds.loc[i, "p_family"])
        f1_pass = p_family >= p_family_threshold
        in_hits = smi in set(hits["smiles"]) if not hits.empty else False
        rows.append(
            {
                "name": b["name"],
                "expected_profile": b.get("expected_profile", ""),
                "in_demo_library": in_library,
                "p_family_pred": p_family,
                "f1_pass": f1_pass,
                "in_final_hits": in_hits,
                "pred_pAct_JNK1": float(preds.loc[i, "pred_pAct_JNK1"]),
                "pred_pAct_JNK2": float(preds.loc[i, "pred_pAct_JNK2"]),
                "pred_pAct_JNK3": float(preds.loc[i, "pred_pAct_JNK3"]),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Virtual screening — F1 p_family pre-filter (v2)")
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS_DIR, help="Per-isoform XGBoost models")
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY, help="SMILES file (.smi/.csv/.smi.gz)")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "screening_v2")
    parser.add_argument("--benchmarks", type=Path, default=DEFAULT_BENCHMARKS)
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument("--top-n", type=int, default=500)
    parser.add_argument("--diverse-n", type=int, default=100)
    parser.add_argument("--p-family-threshold", type=float, default=None, help="Override config screening.p_family_threshold")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config = load_config()
    screening = config.get("screening", {})
    p_family_threshold = args.p_family_threshold or screening.get("p_family_threshold", 6.0)
    sa_max = screening.get("sa_max", 6.0)
    qed_min = screening.get("qed_min", 0.35)
    weights = screening.get(
        "scoring_weights",
        {"w_p_family": 0.55, "w_pAct_JNK1": 0.15, "w_qed": 0.20, "w_sa": 0.10},
    )

    models = load_xgboost_models(args.models_dir)
    logger.info("Loaded XGBoost models from %s", args.models_dir)
    logger.info("F1 threshold: p_family >= %.1f", p_family_threshold)

    if not args.library.exists():
        raise SystemExit(f"Library not found: {args.library}. Run scripts/build_demo_library.py first.")

    all_smiles_raw = read_smiles_library(args.library)
    logger.info("Reading library: %s (%d SMILES)", args.library, len(all_smiles_raw))

    hits: list[pd.DataFrame] = []
    stats = {"input": len(all_smiles_raw), "preprocessed": 0, "druglike": 0, "f1_pass": 0, "sa_qed_pass": 0}
    canonical_seen: set[str] = set()

    for start in tqdm(range(0, len(all_smiles_raw), args.batch_size), desc="Screening batches"):
        batch_smi = all_smiles_raw[start : start + args.batch_size]

        clean = []
        for s in batch_smi:
            cs = preprocess_smiles(s)
            if cs and cs not in canonical_seen:
                canonical_seen.add(cs)
                clean.append(cs)
        stats["preprocessed"] += len(clean)

        druglike = [s for s in clean if lipinski_filter(s)]
        stats["druglike"] += len(druglike)
        if not druglike:
            continue

        batch_df = predict_batch(models, druglike)
        batch_df = batch_df[batch_df["p_family"] >= p_family_threshold]
        stats["f1_pass"] += len(batch_df)
        if batch_df.empty:
            continue

        sa_qed = [compute_sa_qed(s) for s in batch_df["smiles"]]
        batch_df = batch_df.copy()
        batch_df["sa"] = [x[0] for x in sa_qed]
        batch_df["qed"] = [x[1] for x in sa_qed]
        batch_df = batch_df[(batch_df["sa"] <= sa_max) & (batch_df["qed"] >= qed_min)]
        stats["sa_qed_pass"] += len(batch_df)
        if batch_df.empty:
            continue

        batch_df["final_score"] = batch_df.apply(lambda r: compute_final_score(r, weights), axis=1)
        hits.append(batch_df)

    bench_rows = validate_literature_benchmarks(
        canonical_seen, pd.DataFrame(), args.benchmarks, models, p_family_threshold
    )

    if not hits:
        logger.warning("No hits passed the funnel.")
        plot_funnel(stats, args.output)
        report = {
            "funnel_stats": stats,
            "p_family_threshold": p_family_threshold,
            "sa_max": sa_max,
            "qed_min": qed_min,
            "benchmarks": bench_rows,
        }
        with open(args.output / "screening_report.json", "w") as f:
            json.dump(report, f, indent=2)
        if bench_rows:
            pd.DataFrame(bench_rows).to_csv(args.output / "benchmark_validation.csv", index=False)
        return

    all_hits = pd.concat(hits, ignore_index=True).sort_values("final_score", ascending=False)
    all_hits.to_csv(args.output / "all_hits.csv", index=False)

    top = all_hits.head(args.top_n)
    top.to_csv(args.output / f"top{args.top_n}.csv", index=False)

    diverse = butina_diverse_selection(top, n=args.diverse_n)
    diverse.to_csv(args.output / f"top{args.diverse_n}_diverse.csv", index=False)

    plot_funnel(stats, args.output)
    plot_score_distribution(all_hits, args.output)

    bench_rows = validate_literature_benchmarks(
        canonical_seen, all_hits, args.benchmarks, models, p_family_threshold
    )
    if bench_rows:
        pd.DataFrame(bench_rows).to_csv(args.output / "benchmark_validation.csv", index=False)

    report = {
        "funnel_stats": stats,
        "p_family_threshold": p_family_threshold,
        "sa_max": sa_max,
        "qed_min": qed_min,
        "weights": weights,
        "top10": top.head(10)[["smiles", "p_family", "final_score"]].to_dict(orient="records"),
        "benchmarks": bench_rows,
    }
    with open(args.output / "screening_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info(
        "Screening complete: %d hits → Top %d → Diverse %d",
        len(all_hits),
        len(top),
        len(diverse),
    )
    logger.info("Funnel: %s", stats)
    logger.info("Results → %s", args.output)


if __name__ == "__main__":
    main()
