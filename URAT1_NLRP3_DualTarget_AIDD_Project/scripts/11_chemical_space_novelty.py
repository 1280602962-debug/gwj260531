#!/usr/bin/env python3
"""
Non-docking computational module D — chemical space & novelty analysis.

(1) Chemical space: Morgan(ECFP4) fingerprints for the clinical library, the
    dual-docked pool, the Pareto shortlist, and the URAT1/NLRP3 training actives,
    projected to 2D by PCA (UMAP optionally if umap-learn is available). Exports
    coordinates so the manuscript can show where hits sit relative to knowns.

(2) Novelty: nearest-neighbour Tanimoto (ECFP4) from every shortlisted / pooled
    compound to the URAT1 and NLRP3 known-active sets, quantifying how novel each
    hit is versus prior chemistry.

Pure downstream analysis of EXISTING data: no docking, no score changes.

Inputs (read-only):
  data/repurposing/pareto/pareto_shortlist.csv
  data/repurposing/pareto/pareto_merged_scores.csv
  data/processed/urat1_curated.csv        (URAT1 actives; canonical_smiles)
  data/processed/nlrp3_records.csv         (NLRP3 records; canonical_smiles, active)

Outputs:
  results/cheminformatics/chemical_space_pca.csv
  results/cheminformatics/novelty_shortlist.csv
  results/cheminformatics/novelty_pool.csv
  results/cheminformatics/chemspace_novelty_summary.json

Usage:
  python3 scripts/11_chemical_space_novelty.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARETO_DIR = PROJECT_ROOT / "data" / "repurposing" / "pareto"
PROCESSED = PROJECT_ROOT / "data" / "processed"
OUT_DIR = PROJECT_ROOT / "results" / "cheminformatics"

SMILES_COL = "canonical_smiles"
NAME_COL = "name"
N_BITS = 2048
RADIUS = 2


def fp_from_smiles(smiles_list):
    """Return (valid_indices, list_of_ExplicitBitVect)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    idx, fps = [], []
    for i, s in enumerate(smiles_list):
        if not isinstance(s, str):
            continue
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(m, RADIUS, nBits=N_BITS))
        idx.append(i)
    return idx, fps


def fps_to_matrix(fps) -> np.ndarray:
    from rdkit import DataStructs

    arr = np.zeros((len(fps), N_BITS), dtype=np.int8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, arr[i])
    return arr


def nearest_tanimoto(query_fps, ref_fps) -> list[float]:
    from rdkit import DataStructs

    out = []
    for q in query_fps:
        if not ref_fps:
            out.append(float("nan"))
            continue
        sims = DataStructs.BulkTanimotoSimilarity(q, ref_fps)
        out.append(float(max(sims)))
    return out


def load_smiles(path: Path, smiles_col: str = SMILES_COL, active_only: bool = False) -> list[str]:
    df = pd.read_csv(path, low_memory=False)
    if active_only and "active" in df.columns:
        df = df[df["active"] == 1]
    if smiles_col not in df.columns:
        raise ValueError(f"{smiles_col} not in {path} (has {list(df.columns)[:8]}...)")
    return df[smiles_col].dropna().astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description="Chemical space & novelty (non-docking module D)")
    parser.add_argument("--shortlist", type=Path, default=PARETO_DIR / "pareto_shortlist.csv")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    from sklearn.decomposition import PCA

    # ---- Reference active sets ----
    urat1_smiles = load_smiles(PROCESSED / "urat1_curated.csv")
    nlrp3_smiles = load_smiles(PROCESSED / "nlrp3_records.csv", active_only=True)
    u_idx, u_fps = fp_from_smiles(urat1_smiles)
    n_idx, n_fps = fp_from_smiles(nlrp3_smiles)

    # ---- Query sets ----
    shortlist = pd.read_csv(args.shortlist)
    pool = pd.read_csv(args.pool)
    s_idx, s_fps = fp_from_smiles(shortlist[SMILES_COL].tolist())
    p_idx, p_fps = fp_from_smiles(pool[SMILES_COL].tolist())

    # ---- Novelty: nearest Tanimoto to URAT1 and NLRP3 actives ----
    def novelty_frame(df, idx, fps):
        sub = df.iloc[idx].copy().reset_index(drop=True)
        sub_meta = sub[[c for c in (NAME_COL, SMILES_COL, "max_phase",
                                     "s_u_percentile", "s_n_percentile", "pareto_front")
                        if c in sub.columns]].copy()
        sub_meta["nn_tanimoto_urat1_active"] = np.round(nearest_tanimoto(fps, u_fps), 3)
        sub_meta["nn_tanimoto_nlrp3_active"] = np.round(nearest_tanimoto(fps, n_fps), 3)
        sub_meta["novel_vs_urat1_tc_lt_0.4"] = sub_meta["nn_tanimoto_urat1_active"] < 0.4
        sub_meta["novel_vs_nlrp3_tc_lt_0.4"] = sub_meta["nn_tanimoto_nlrp3_active"] < 0.4
        return sub_meta

    nov_short = novelty_frame(shortlist, s_idx, s_fps)
    nov_pool = novelty_frame(pool, p_idx, p_fps)
    nov_short.to_csv(args.output_dir / "novelty_shortlist.csv", index=False)
    nov_pool.to_csv(args.output_dir / "novelty_pool.csv", index=False)

    # ---- Chemical space PCA (fit on combined fps) ----
    groups, all_fps, labels = [], [], []
    for label, fps in (
        ("urat1_active", u_fps),
        ("nlrp3_active", n_fps),
        ("dual_docked_pool", p_fps),
        ("pareto_shortlist", s_fps),
    ):
        all_fps.extend(fps)
        labels.extend([label] * len(fps))
    X = fps_to_matrix(all_fps).astype(np.float32)
    pca = PCA(n_components=2, random_state=args.seed)
    coords = pca.fit_transform(X)
    space = pd.DataFrame({
        "set": labels,
        "pc1": np.round(coords[:, 0], 4),
        "pc2": np.round(coords[:, 1], 4),
    })
    # attach names for the shortlist rows (last block)
    space_names = [None] * len(space)
    start = len(u_fps) + len(n_fps) + len(p_fps)
    if NAME_COL in shortlist.columns:
        for j, i in enumerate(s_idx):
            space_names[start + j] = shortlist.iloc[i][NAME_COL]
    space["name"] = space_names
    space.to_csv(args.output_dir / "chemical_space_pca.csv", index=False)

    summary = {
        "module": "D_chemical_space_novelty",
        "n_urat1_active": len(u_fps),
        "n_nlrp3_active": len(n_fps),
        "n_pool": len(p_fps),
        "n_shortlist": len(s_fps),
        "pca_explained_variance_ratio": [round(float(v), 4) for v in pca.explained_variance_ratio_],
        "shortlist_novelty": nov_short.to_dict(orient="records"),
        "pool_median_nn_tc_urat1": round(float(nov_pool["nn_tanimoto_urat1_active"].median()), 3),
        "pool_median_nn_tc_nlrp3": round(float(nov_pool["nn_tanimoto_nlrp3_active"].median()), 3),
    }
    with open(args.output_dir / "chemspace_novelty_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=== Chemical space & novelty ===")
    print(f"  URAT1 actives={len(u_fps)}  NLRP3 actives={len(n_fps)}  pool={len(p_fps)}  shortlist={len(s_fps)}")
    print(f"  PCA explained variance: {summary['pca_explained_variance_ratio']}")
    print("\nShortlist novelty (nearest Tanimoto to known actives):")
    cols = [NAME_COL, "nn_tanimoto_urat1_active", "nn_tanimoto_nlrp3_active"]
    print(nov_short[[c for c in cols if c in nov_short.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
