#!/usr/bin/env python3
"""SI audit: scaffold/similarity leakage between URAT1 actives, weak-actives, and decoys.

Answers three reviewer-facing questions without rebuilding or re-docking the
TrueDecoy/RandomDecoy benchmark:

  1. Do any of the 80 experimental weak-actives share a Murcko scaffold, or
     have high Tanimoto similarity, with the 469 actives? (property-matched
     decoys are already capped at max Tanimoto-to-active <= 0.5 by
     construction; weak-actives were not put through that filter.)
  2. RandomDecoy negatives have no active-similarity exclusion by design
     (Gu-style random draw). What fraction are, in fact, close analogs of an
     active (TC > 0.5)?
  3. Sensitivity check: recompute each protocol's RandomDecoy EF@1%/EF@5%/AUC
     after excluding those close-analog RandomDecoy negatives (TC <= 0.5
     only), using the already-archived P0-P5 scores. Does Pi* selection
     change?

Does not rebuild the benchmark, does not redock, does not change the locked
Pi* = P2 decision.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRUE_BENCH = PROJECT_ROOT / "data" / "benchmarks" / "urat1_true_decoy" / "true_decoy_benchmark.csv"
RANDOM_BENCH = PROJECT_ROOT / "data" / "benchmarks" / "urat1_true_decoy" / "random_decoy_benchmark.csv"
MOL_SCORES = PROJECT_ROOT / "data" / "benchmarks" / "protocol_selection" / "mol_protocol_scores.csv"
OUT_DIR = PROJECT_ROOT / "data" / "si" / "decoy_leakage_audit"

PROTOCOLS = [
    ("P0", "P0_CNNscore", True),
    ("P1", "P1_vina_affinity", False),
    ("P2", "P2_CNNaffinity", True),
    ("P3", "P3_gnina_affinity", False),
    ("P4", "P4_RTM_vina", True),
    ("P5", "P5_RTM_gnina", True),
]


def morgan_fps(smiles_list: list[str]):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol else None)
    return fps


def max_tc_to_set(query_fps, ref_fps) -> np.ndarray:
    from rdkit import DataStructs

    ref_fps = [f for f in ref_fps if f is not None]
    out = np.full(len(query_fps), np.nan)
    for i, qfp in enumerate(query_fps):
        if qfp is None or not ref_fps:
            continue
        out[i] = float(max(DataStructs.BulkTanimotoSimilarity(qfp, ref_fps)))
    return out


def enrichment_factor(y: np.ndarray, score: np.ndarray, fraction: float, higher_is_better: bool) -> tuple[float, str]:
    mask = ~np.isnan(score)
    y, score = y[mask], score[mask]
    n = max(1, int(np.floor(fraction * len(y))))
    order = np.argsort(-score if higher_is_better else score)
    top = y[order[:n]]
    base = y.mean()
    ef = float(top.mean() / base) if base > 0 else np.nan
    return ef, f"{int(top.sum())}/{n}"


def roc_auc(y: np.ndarray, score: np.ndarray, higher_is_better: bool) -> float:
    from sklearn.metrics import roc_auc_score

    mask = ~np.isnan(score)
    y, score = y[mask], score[mask]
    s = score if higher_is_better else -score
    return float(roc_auc_score(y, s))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    true_bench = pd.read_csv(TRUE_BENCH, low_memory=False)
    random_bench = pd.read_csv(RANDOM_BENCH, low_memory=False)

    actives = true_bench[true_bench["decoy_source"] == "active"].reset_index(drop=True)
    weak = true_bench[true_bench["decoy_source"] == "experimental_inactive"].reset_index(drop=True)
    matched = true_bench[true_bench["decoy_source"] == "property_matched"].reset_index(drop=True)
    random_decoys = random_bench[random_bench["decoy_source"] == "library_random"].reset_index(drop=True)

    active_scaffolds = set(actives["scaffold"].dropna())
    active_fps = morgan_fps(actives["canonical_smiles"].tolist())

    # 1. Weak-active vs active: scaffold overlap + similarity
    weak_scaffold_overlap = set(weak["scaffold"].dropna()) & active_scaffolds
    weak_tc = max_tc_to_set(morgan_fps(weak["canonical_smiles"].tolist()), active_fps)

    # 2. RandomDecoy vs active: scaffold overlap + similarity (no filter applied at build time)
    rand_scaffold_overlap = set(random_decoys["scaffold"].dropna()) & active_scaffolds
    rand_tc = max_tc_to_set(morgan_fps(random_decoys["canonical_smiles"].tolist()), active_fps)
    random_decoys = random_decoys.copy()
    random_decoys["max_tc_active"] = rand_tc

    # property-matched decoys already have max_tc_active <= 0.5 by construction; report for completeness
    matched_tc = matched["max_tc_active"].astype(float).values

    summary = {
        "n_actives": int(len(actives)),
        "n_weak_actives": int(len(weak)),
        "n_property_matched_decoys": int(len(matched)),
        "n_random_decoys": int(len(random_decoys)),
        "weak_active": {
            "n_scaffold_overlap_with_actives": len(weak_scaffold_overlap),
            "frac_scaffold_overlap": len(weak_scaffold_overlap) / max(len(set(weak["scaffold"].dropna())), 1),
            "max_tc_to_active_mean": float(np.nanmean(weak_tc)),
            "max_tc_to_active_median": float(np.nanmedian(weak_tc)),
            "max_tc_to_active_p90": float(np.nanpercentile(weak_tc, 90)),
            "n_tc_gt_0.5": int(np.nansum(weak_tc > 0.5)),
            "n_tc_gt_0.85": int(np.nansum(weak_tc > 0.85)),
        },
        "property_matched_decoy": {
            "max_tc_to_active_mean": float(np.nanmean(matched_tc)),
            "max_tc_to_active_max": float(np.nanmax(matched_tc)),
            "n_tc_gt_0.5": int(np.nansum(matched_tc > 0.5)),
        },
        "random_decoy": {
            "n_scaffold_overlap_with_actives": len(rand_scaffold_overlap),
            "frac_scaffold_overlap": len(rand_scaffold_overlap) / max(len(set(random_decoys["scaffold"].dropna())), 1),
            "max_tc_to_active_mean": float(np.nanmean(rand_tc)),
            "max_tc_to_active_median": float(np.nanmedian(rand_tc)),
            "max_tc_to_active_p90": float(np.nanpercentile(rand_tc, 90)),
            "n_tc_gt_0.5": int(np.nansum(rand_tc > 0.5)),
            "frac_tc_gt_0.5": float(np.nanmean(rand_tc > 0.5)),
            "n_tc_gt_0.85": int(np.nansum(rand_tc > 0.85)),
        },
    }

    random_decoys[["canonical_smiles", "scaffold", "max_tc_active"]].to_csv(
        OUT_DIR / "random_decoy_max_tc_to_active.csv", index=False
    )
    weak_out = weak[["canonical_smiles", "scaffold"]].copy()
    weak_out["max_tc_active"] = weak_tc
    weak_out.to_csv(OUT_DIR / "weak_active_max_tc_to_active.csv", index=False)

    # 3. Sensitivity: RandomDecoy EF/AUC recomputed after dropping close-analog decoys (TC > 0.5)
    mol_scores = pd.read_csv(MOL_SCORES, low_memory=False)
    smi_col = "canonical_smiles" if "canonical_smiles" in mol_scores.columns else "Smiles"
    id_col = None
    for cand in ("mol_id", "canonical_smiles", "id"):
        if cand in mol_scores.columns:
            id_col = cand
            break

    rand_full = random_bench.merge(mol_scores, left_on="canonical_smiles", right_on=smi_col, how="left")
    keep_smiles = set(random_decoys.loc[random_decoys["max_tc_active"] <= 0.5, "canonical_smiles"])
    active_smiles = set(actives["canonical_smiles"])
    rand_filtered = rand_full[
        rand_full["canonical_smiles"].isin(active_smiles) | rand_full["canonical_smiles"].isin(keep_smiles)
    ]

    rows = []
    for pid, col, higher in PROTOCOLS:
        if col not in rand_full.columns:
            continue
        for label, frame in (("random_full", rand_full), ("random_tc_le_0.5", rand_filtered)):
            y = frame["label"].values.astype(int)
            score = pd.to_numeric(frame[col], errors="coerce").values
            ef1, hits1 = enrichment_factor(y, score, 0.01, higher)
            ef5, hits5 = enrichment_factor(y, score, 0.05, higher)
            try:
                auc = roc_auc(y, score, higher)
            except ValueError:
                auc = float("nan")
            rows.append({
                "protocol": pid, "benchmark": label, "n": int((~pd.isna(score)).sum()),
                "EF1pct": ef1, "hits1pct": hits1, "EF5pct": ef5, "AUC": auc,
            })
    pd.DataFrame(rows).to_csv(OUT_DIR / "random_decoy_ef_before_after_tc_filter.csv", index=False)
    summary["random_decoy_n_dropped_tc_gt_0.5"] = int(len(random_decoys) - len(keep_smiles))
    summary["random_decoy_n_kept_after_filter"] = int(len(keep_smiles))

    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
