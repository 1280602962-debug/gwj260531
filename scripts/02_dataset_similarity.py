#!/usr/bin/env python3
"""
Phase 1: Compare similarity among JNK1/2/3 datasets.

Analyses:
  - Morgan fingerprint Tanimoto cross-similarity
  - UMAP chemical space visualization
  - Bemis-Murcko scaffold overlap
  - pActivity distribution comparison (KS test)
  - Paired compound selectivity distribution

Usage:
    python scripts/02_dataset_similarity.py --input data/processed --output results/similarity
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp
from sklearn.decomposition import PCA

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]


def morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def fp_array(fp) -> np.ndarray:
    arr = np.zeros((1,), dtype=int)
    from rdkit.DataStructs import ConvertToNumpyArray

    ConvertToNumpyArray(fp, arr)
    return arr


def tanimoto(a: np.ndarray, b: np.ndarray) -> float:
    inter = np.sum(a & b)
    union = np.sum(a | b)
    return float(inter / union) if union else 0.0


def murcko_scaffold(smiles: str) -> str | None:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def load_datasets(input_dir: Path) -> dict[str, pd.DataFrame]:
    datasets = {}
    for isoform in ["jnk1", "jnk2", "jnk3"]:
        path = input_dir / f"{isoform}_curated.csv"
        if path.exists():
            datasets[isoform.upper()] = pd.read_csv(path)
            logger.info("Loaded %s: %d compounds", isoform.upper(), len(datasets[isoform.upper()]))
        else:
            logger.warning("Missing %s", path)
    return datasets


def cross_similarity(datasets: dict[str, pd.DataFrame], sample_n: int = 2000) -> pd.DataFrame:
    """Mean cross-dataset Tanimoto similarity."""
    isoforms = list(datasets.keys())
    fps = {}
    for iso, df in datasets.items():
        smiles = df["canonical_smiles"].dropna().unique()
        if len(smiles) > sample_n:
            smiles = np.random.choice(smiles, sample_n, replace=False)
        fps[iso] = {s: morgan_fp(s) for s in smiles}
        fps[iso] = {s: fp for s, fp in fps[iso].items() if fp is not None}

    matrix = pd.DataFrame(index=isoforms, columns=isoforms, dtype=float)
    for a in isoforms:
        for b in isoforms:
            if a == b:
                sims = []
                arr_a = list(fps[a].values())
                idx = np.random.choice(len(arr_a), min(500, len(arr_a)), replace=False)
                for i in idx:
                    for j in idx:
                        if i < j:
                            sims.append(tanimoto(fp_array(arr_a[i]), fp_array(arr_a[j])))
                matrix.loc[a, b] = np.mean(sims) if sims else 1.0
            else:
                sims = []
                keys_a = list(fps[a].keys())[:300]
                keys_b = list(fps[b].keys())[:300]
                for sa in keys_a:
                    fa = fp_array(fps[a][sa])
                    for sb in keys_b:
                        sims.append(tanimoto(fa, fp_array(fps[b][sb])))
                matrix.loc[a, b] = np.mean(sims) if sims else np.nan
    return matrix


def scaffold_overlap(datasets: dict[str, pd.DataFrame]) -> dict:
    scaffolds = {}
    for iso, df in datasets.items():
        scaffolds[iso] = set(
            filter(None, [murcko_scaffold(s) for s in df["canonical_smiles"].dropna()])
        )

    all_iso = list(scaffolds.keys())
    shared_all = set.intersection(*scaffolds.values()) if scaffolds else set()
    result = {
        "unique_scaffolds": {iso: len(s) for iso, s in scaffolds.items()},
        "shared_all_three": len(shared_all),
        "pairwise_jaccard": {},
    }
    for i, a in enumerate(all_iso):
        for b in all_iso[i + 1 :]:
            inter = len(scaffolds[a] & scaffolds[b])
            union = len(scaffolds[a] | scaffolds[b])
            result["pairwise_jaccard"][f"{a}-{b}"] = inter / union if union else 0.0
    return result


def activity_distribution(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for iso, df in datasets.items():
        p = df["pActivity"].dropna()
        rows.append(
            {
                "isoform": iso,
                "n": len(p),
                "median": p.median(),
                "mean": p.mean(),
                "std": p.std(),
                "active_frac": (p >= 6.5).mean(),
            }
        )
    return pd.DataFrame(rows)


def ks_tests(datasets: dict[str, pd.DataFrame]) -> dict:
    results = {}
    if "JNK1" in datasets and "JNK2" in datasets:
        stat, p = ks_2samp(
            datasets["JNK1"]["pActivity"].dropna(),
            datasets["JNK2"]["pActivity"].dropna(),
        )
        results["JNK1_vs_JNK2"] = {"statistic": stat, "p_value": p}
    if "JNK1" in datasets and "JNK3" in datasets:
        stat, p = ks_2samp(
            datasets["JNK1"]["pActivity"].dropna(),
            datasets["JNK3"]["pActivity"].dropna(),
        )
        results["JNK1_vs_JNK3"] = {"statistic": stat, "p_value": p}
    return results


def plot_umap(datasets: dict[str, pd.DataFrame], output_dir: Path, sample_n: int = 3000):
    frames = []
    for iso, df in datasets.items():
        sub = df[["canonical_smiles", "pActivity"]].dropna().copy()
        sub["isoform"] = iso
        if len(sub) > sample_n:
            sub = sub.sample(sample_n, random_state=42)
        frames.append(sub)
    combined = pd.concat(frames, ignore_index=True)

    fps = []
    valid_idx = []
    for i, smi in enumerate(combined["canonical_smiles"]):
        fp = morgan_fp(smi)
        if fp is not None:
            fps.append(fp_array(fp))
            valid_idx.append(i)
    X = np.vstack(fps)
    combined = combined.iloc[valid_idx].reset_index(drop=True)

    try:
        import umap

        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
        coords = reducer.fit_transform(X)
        method = "UMAP"
    except ImportError:
        coords = PCA(n_components=2, random_state=42).fit_transform(X)
        method = "PCA"

    combined["x"] = coords[:, 0]
    combined["y"] = coords[:, 1]

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(data=combined, x="x", y="y", hue="isoform", alpha=0.5, s=15, ax=ax)
    ax.set_title(f"JNK1/2/3 Chemical Space ({method})")
    ax.set_xlabel(f"{method} 1")
    ax.set_ylabel(f"{method} 2")
    fig.tight_layout()
    fig.savefig(output_dir / "chemical_space.png", dpi=200)
    plt.close(fig)
    logger.info("Saved chemical space plot (%s)", method)


def plot_similarity_heatmap(matrix: pd.DataFrame, output_dir: Path):
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(matrix.astype(float), annot=True, fmt=".3f", cmap="YlOrRd", ax=ax)
    ax.set_title("Cross-Dataset Mean Tanimoto Similarity")
    fig.tight_layout()
    fig.savefig(output_dir / "cross_similarity_heatmap.png", dpi=200)
    plt.close(fig)


def analyze_paired(paired_path: Path, output_dir: Path):
    if not paired_path.exists():
        logger.warning("Paired set not found: %s", paired_path)
        return
    paired = pd.read_csv(paired_path)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, col in zip(axes, ["delta_12", "delta_13"]):
        if col in paired.columns:
            paired[col].dropna().hist(bins=40, ax=ax, edgecolor="black")
            ax.axvline(1.0, color="red", linestyle="--", label="Δ=1 (10× SI)")
            ax.set_xlabel(col)
            ax.set_ylabel("Count")
            ax.legend()
    fig.suptitle("Selectivity Distribution (Paired Compounds)")
    fig.tight_layout()
    fig.savefig(output_dir / "selectivity_distribution.png", dpi=200)
    plt.close(fig)

    if "sel_class" in paired.columns:
        counts = paired["sel_class"].value_counts()
        counts.to_csv(output_dir / "sel_class_counts.csv")
        logger.info("Selectivity class counts:\n%s", counts)


def main():
    parser = argparse.ArgumentParser(description="JNK1/2/3 dataset similarity analysis")
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "processed")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "similarity")
    parser.add_argument("--sample-n", type=int, default=2000)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    datasets = load_datasets(args.input)
    if len(datasets) < 2:
        raise SystemExit("Need at least 2 curated datasets. Run 01_download_chembl_data.py first.")

    # 3.2 Cross-similarity
    sim_matrix = cross_similarity(datasets, sample_n=args.sample_n)
    sim_matrix.to_csv(args.output / "cross_similarity_matrix.csv")
    plot_similarity_heatmap(sim_matrix, args.output)

    # 3.3 Scaffold overlap
    scaffold_stats = scaffold_overlap(datasets)
    with open(args.output / "scaffold_overlap.json", "w") as f:
        json.dump(scaffold_stats, f, indent=2)

    # 3.4 Activity distribution
    act_stats = activity_distribution(datasets)
    act_stats.to_csv(args.output / "activity_distribution.csv", index=False)
    ks = ks_tests(datasets)
    with open(args.output / "ks_tests.json", "w") as f:
        json.dump(ks, f, indent=2)

    # 3.2 UMAP
    plot_umap(datasets, args.output)

    # 3.5 Paired analysis
    analyze_paired(args.input / "paired_set.csv", args.output)

    # Summary report
    report = {
        "cross_similarity": sim_matrix.to_dict(),
        "scaffold_overlap": scaffold_stats,
        "activity_distribution": act_stats.to_dict(orient="records"),
        "ks_tests": ks,
        "recommendations": [],
    }
    mean_cross = sim_matrix.values[np.triu_indices_from(sim_matrix.values, k=1)].mean()
    if mean_cross > 0.25:
        report["recommendations"].append("High chemical space overlap → use MTL / transfer learning")
    if scaffold_stats["shared_all_three"] > 50:
        report["recommendations"].append("Many shared scaffolds → strict scaffold split required")

    with open(args.output / "similarity_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info("Analysis complete. Results in %s", args.output)


if __name__ == "__main__":
    main()
