#!/usr/bin/env python3
"""Normalize Maestro Canvas docking CSV (XP_OUT) for merge_docking_pareto.py.

DEPRECATED: Project now uses open-source AutoDock Vina. Prefer:
  scripts/normalize_docking_export.py
See docs/OPEN_SOURCE_DOCKING.md

Maestro exports use s_canvas_repurposing_id + r_glide_XP_GScore, not canonical_smiles.
This script joins docking_pool_p05.csv to attach canonical_smiles.

Example:
  python3 scripts/normalize_canvas_docking_export.py \\
    --input ~/Maestro/dual_dockingvsw_1-XP_OUT_2_853b.csv \\
    --pdb 9DKB \\
    --output results/repurposing/docking_raw/urat1_9dkb_p05.csv

  python3 scripts/normalize_canvas_docking_export.py \\
    --input ~/Maestro/dual_dockingvsw_1-XP_OUT_1_caad.csv \\
    --pdb 7ALV \\
    --output results/repurposing/docking_raw/nlrp3_7alv_p05.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POOL = PROJECT_ROOT / "data" / "repurposing" / "screening" / "docking_pool_p05.csv"

REPURPOSING_ALIASES = ["repurposing_id", "s_canvas_repurposing_id", "s_canvas_repurposing\\_id"]
XP_SCORE_ALIASES = ["r_glide_XP_GScore", "r_glide_xp_gscore", "glide_score_xp", "r_i_glide_gscore"]
GRID_ALIASES = ["s_i_glide_gridfile", "gridfile", "grid"]


def _pick_col(df: pd.DataFrame, aliases: list[str]) -> str | None:
    norm = {c.lower().strip(): c for c in df.columns}
    for a in aliases:
        key = a.lower().strip()
        if key in norm:
            return norm[key]
    return None


def normalize_export(
    input_csv: Path,
    pool_csv: Path,
    pdb_id: str,
    output_csv: Path,
) -> dict:
    raw = pd.read_csv(input_csv, low_memory=False)
    pool = pd.read_csv(pool_csv, low_memory=False)

    rid_col = _pick_col(raw, REPURPOSING_ALIASES)
    score_col = _pick_col(raw, XP_SCORE_ALIASES)
    grid_col = _pick_col(raw, GRID_ALIASES)
    if rid_col is None:
        raise ValueError(f"No repurposing id column in {input_csv}")
    if score_col is None:
        raise ValueError(f"No XP score column in {input_csv}")
    if "repurposing_id" not in pool.columns or "canonical_smiles" not in pool.columns:
        raise ValueError(f"{pool_csv} must contain repurposing_id and canonical_smiles")

    merged = raw.merge(
        pool[["repurposing_id", "canonical_smiles", "name", "chembl_id"]],
        left_on=rid_col,
        right_on="repurposing_id",
        how="left",
    )
    merged["glide_score_xp"] = pd.to_numeric(merged[score_col], errors="coerce")
    merged["docking_status"] = merged["glide_score_xp"].notna().map({True: "docked", False: "missing"})
    merged["pdb_id"] = pdb_id.upper()
    if grid_col:
        merged["gridfile"] = merged[grid_col].astype(str)

    out_cols = [
        "repurposing_id",
        "name",
        "chembl_id",
        "canonical_smiles",
        "glide_score_xp",
        "docking_status",
        "pdb_id",
    ]
    if grid_col:
        out_cols.append("gridfile")
    out = merged[out_cols].copy()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_csv, index=False)

    pool_ids = set(pool["repurposing_id"].astype(str))
    dock_ids = set(out["repurposing_id"].astype(str))
    summary = {
        "input": str(input_csv),
        "output": str(output_csv),
        "pdb_id": pdb_id.upper(),
        "n_input_rows": int(len(raw)),
        "n_output_rows": int(len(out)),
        "n_with_canonical_smiles": int(out["canonical_smiles"].notna().sum()),
        "n_in_pool": int(len(dock_ids & pool_ids)),
        "n_pool_missing_dock": int(len(pool_ids - dock_ids)),
        "xp_score_column": score_col,
        "grid_sample": out["gridfile"].iloc[0] if grid_col and len(out) else None,
        "glide_score_xp_min": float(out["glide_score_xp"].min()) if out["glide_score_xp"].notna().any() else None,
        "glide_score_xp_max": float(out["glide_score_xp"].max()) if out["glide_score_xp"].notna().any() else None,
    }
    qc_path = output_csv.with_suffix(".qc.json")
    qc_path.write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize Maestro Canvas XP CSV for Pareto merge")
    parser.add_argument("--input", type=Path, required=True, help="Maestro XP_OUT *.csv")
    parser.add_argument("--pdb", type=str, required=True, help="9DKB or 7ALV or 8ETR")
    parser.add_argument("--pool", type=Path, default=DEFAULT_POOL, help="docking_pool_p05.csv")
    parser.add_argument("--output", type=Path, required=True, help="Normalized output CSV")
    args = parser.parse_args()
    summary = normalize_export(args.input, args.pool, args.pdb, args.output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
