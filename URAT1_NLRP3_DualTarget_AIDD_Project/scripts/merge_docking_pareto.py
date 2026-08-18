#!/usr/bin/env python3
"""Merge URAT1 + NLRP3 docking scores on the P(active)>=0.5 pool.

Production engine is gnina P2. Do not mix historical Glide scores into a P2 ranking table.
Raw Pareto is audit-only; follow-up uses scripts/14_candidate_nomination.py.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dock_score_utils import (
    DOCK_SCORE_ALIASES,
    NAME_ALIASES,
    SMILES_ALIASES,
    STATUS_ALIASES,
    best_pose_per_compound,
    docking_status_from_score,
    ensure_dock_score_column,
    percentile_rank,
    pick_col,
)
from utils_ml import canonicalize

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = PROJECT_ROOT / "results" / "repurposing"


def _load_dock_table(path: Path, pdb_label: str) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    smi_col = pick_col(raw.columns, SMILES_ALIASES)
    score_col = pick_col(raw.columns, DOCK_SCORE_ALIASES)
    status_col = pick_col(raw.columns, STATUS_ALIASES)
    pdb_col = pick_col(raw.columns, ["pdb_id", "pdb", "structure_id"])
    if smi_col is None:
        raise ValueError(f"{path}: no SMILES column found (tried {SMILES_ALIASES})")
    if score_col is None:
        raise ValueError(f"{path}: no docking score column found (tried {DOCK_SCORE_ALIASES})")

    out = raw.copy()
    out["canonical_smiles"] = out[smi_col].map(lambda s: canonicalize(s) if pd.notna(s) else None)
    out = out[out["canonical_smiles"].notna()].copy()
    out = ensure_dock_score_column(out, score_col)
    out["docking_status"] = docking_status_from_score(
        out["dock_score"],
        out[status_col] if status_col else None,
    )
    if pdb_col:
        out["pdb_id"] = out[pdb_col].astype(str).str.upper()
    else:
        out["pdb_id"] = pdb_label

    out = best_pose_per_compound(out, score_col="dock_score")
    return out[["canonical_smiles", "dock_score", "glide_score_xp", "docking_status", "pdb_id"]]


def pareto_front(su: np.ndarray, sn: np.ndarray) -> np.ndarray:
    """Return boolean mask of non-dominated points (maximize both axes)."""
    n = len(su)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        if not is_pareto[i]:
            continue
        better = (su >= su[i]) & (sn >= sn[i]) & ((su > su[i]) | (sn > sn[i]))
        better[i] = False
        if better.any():
            is_pareto[i] = False
    return is_pareto


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge dual docking + NLRP3 ML for Pareto shortlist")
    parser.add_argument("--ml-scores", type=Path, required=True, help="nlrp3_ml_scores_*.csv")
    parser.add_argument("--urat1-dock", type=Path, required=True, help="URAT1 9DKB docking export")
    parser.add_argument("--nlrp3-dock", type=Path, required=True, help="NLRP3 7ALV docking export")
    parser.add_argument("--nlrp3-pdb", type=str, default=None, help="NLRP3 PDB label (default: from dock file or 7ALV)")
    parser.add_argument("--pool", type=Path, default=None, help="docking_pool_p05.csv (optional filter)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--sn-mode", choices=["ml", "dock", "both"], default="both", help="S_N axis source")
    parser.add_argument("--min-su", type=float, default=0.0, help="Min S_U percentile for shortlist")
    parser.add_argument("--min-sn", type=float, default=0.0, help="Min S_N percentile for shortlist")
    parser.add_argument(
        "--mw-min", type=float, default=200.0,
        help="Oral-ish MW lower bound (Da) for druglike Pareto shortlist",
    )
    parser.add_argument(
        "--mw-max", type=float, default=550.0,
        help="Oral-ish MW upper bound (Da) for druglike Pareto shortlist; demotes macrolides",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ml = pd.read_csv(args.ml_scores, low_memory=False)
    if "canonical_smiles" not in ml.columns:
        smi = pick_col(ml.columns, SMILES_ALIASES)
        if smi is None:
            raise ValueError("--ml-scores: missing canonical_smiles")
        ml["canonical_smiles"] = ml[smi].map(lambda s: canonicalize(s) if pd.notna(s) else None)

    pool_keys = None
    if args.pool:
        pool = pd.read_csv(args.pool, low_memory=False)
        pool_keys = set(pool["canonical_smiles"].dropna().astype(str))

    urat1 = _load_dock_table(args.urat1_dock, "9DKB")
    nlrp3_pdb = args.nlrp3_pdb or "7ALV"
    nlrp3 = _load_dock_table(args.nlrp3_dock, nlrp3_pdb)

    merged = ml.merge(urat1, on="canonical_smiles", how="inner", suffixes=("", "_urat1"))
    merged = merged.merge(
        nlrp3.rename(
            columns={
                "dock_score": "nlrp3_dock_score",
                "glide_score_xp": "nlrp3_glide_score_xp",
                "docking_status": "nlrp3_docking_status",
                "pdb_id": "nlrp3_pdb_id",
            }
        ),
        on="canonical_smiles",
        how="inner",
    )
    if pool_keys is not None:
        merged = merged[merged["canonical_smiles"].astype(str).isin(pool_keys)].copy()

    # Prefer dock_score; glide_score_xp kept as identical legacy alias
    su_col = "dock_score" if "dock_score" in merged.columns else "glide_score_xp"
    sn_col = "nlrp3_dock_score" if "nlrp3_dock_score" in merged.columns else "nlrp3_glide_score_xp"
    merged["s_u_percentile"] = percentile_rank(merged[su_col], higher_is_better=False)
    merged["s_n_ml_percentile"] = merged["p_active_nlrp3"].rank(method="average", pct=True) * 100.0
    merged["s_n_dock_percentile"] = percentile_rank(merged[sn_col], higher_is_better=False)

    if args.sn_mode == "ml":
        merged["s_n_percentile"] = merged["s_n_ml_percentile"]
    elif args.sn_mode == "dock":
        merged["s_n_percentile"] = merged["s_n_dock_percentile"]
    else:
        merged["s_n_percentile"] = merged[["s_n_ml_percentile", "s_n_dock_percentile"]].max(axis=1)

    merged = merged.sort_values(["s_u_percentile", "s_n_percentile"], ascending=False).reset_index(drop=True)
    su = merged["s_u_percentile"].to_numpy(dtype=float)
    sn = merged["s_n_percentile"].to_numpy(dtype=float)
    merged["pareto_front"] = pareto_front(su, sn)

    shortlist = merged[
        merged["pareto_front"]
        & (merged["s_u_percentile"] >= args.min_su)
        & (merged["s_n_percentile"] >= args.min_sn)
    ].copy()
    shortlist_path = args.output_dir / "pareto_shortlist.csv"
    shortlist.to_csv(shortlist_path, index=False)

    # Chemistry-aware Pareto slice: same non-dominated set, but drop oversized molecules
    # that inflate contact-based docking ranks (macrolides / polyketides). Follow-up
    # story molecules still come from Module F nomination, not this file alone.
    if "mw" in merged.columns:
        mw = pd.to_numeric(merged["mw"], errors="coerce")
        merged["mw_oral_ok"] = mw.between(args.mw_min, args.mw_max, inclusive="both")
    else:
        merged["mw_oral_ok"] = True
        print("WARNING: no mw column in merge; druglike shortlist cannot apply MW window")

    shortlist_druglike = shortlist.copy()
    if "mw" in shortlist_druglike.columns:
        shortlist_druglike = shortlist_druglike[
            pd.to_numeric(shortlist_druglike["mw"], errors="coerce")
            .between(args.mw_min, args.mw_max, inclusive="both")
        ].copy()
    shortlist_druglike_path = args.output_dir / "pareto_shortlist_druglike.csv"
    shortlist_druglike.to_csv(shortlist_druglike_path, index=False)

    merged_path = args.output_dir / "pareto_merged_scores.csv"
    merged.to_csv(merged_path, index=False)

    name_col = pick_col(merged.columns, NAME_ALIASES) or "name"
    controls = ["lesinurad", "benzbromarone", "verinurad", "dotinurad", "colchicine", "allopurinol", "febuxostat"]
    ctrl_hits = []
    if name_col in merged.columns:
        for drug in controls:
            hit = merged[merged[name_col].astype(str).str.upper() == drug.upper()]
            if len(hit):
                row = hit.iloc[0]
                ctrl_hits.append(
                    {
                        "name": drug,
                        "s_u_percentile": float(row["s_u_percentile"]),
                        "s_n_percentile": float(row["s_n_percentile"]),
                        "pareto_front": bool(row["pareto_front"]),
                        "mw_oral_ok": bool(row["mw_oral_ok"]) if "mw_oral_ok" in row.index else None,
                    }
                )

    summary = {
        "n_ml_scores": int(len(ml)),
        "n_merged_dual_dock": int(len(merged)),
        "n_pareto_front": int(merged["pareto_front"].sum()),
        "n_shortlist": int(len(shortlist)),
        "n_shortlist_druglike": int(len(shortlist_druglike)),
        "mw_window_da": [args.mw_min, args.mw_max],
        "sn_mode": args.sn_mode,
        "note": (
            "pareto_shortlist.csv is raw docking Pareto (audit only). "
            "Use pareto_shortlist_druglike.csv + Module F nomination for follow-up."
        ),
        "known_controls": ctrl_hits,
        "outputs": {
            "merged": str(merged_path),
            "shortlist": str(shortlist_path),
            "shortlist_druglike": str(shortlist_druglike_path),
        },
    }
    summary_path = args.output_dir / "pareto_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
