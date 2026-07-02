#!/usr/bin/env python3
"""Normalize docking exports from Vina/smina batch runs or legacy Glide Canvas.

Output schema (tool-agnostic):
  repurposing_id, name, chembl_id, canonical_smiles,
  dock_score, docking_status, pdb_id, docking_engine

Legacy alias `glide_score_xp` is duplicated from dock_score for downstream scripts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from dock_score_utils import (
    DOCK_SCORE_ALIASES,
    REPURPOSING_ALIASES,
    SMILES_ALIASES,
    docking_status_from_score,
    ensure_dock_score_column,
    pick_col,
)
from utils_ml import canonicalize

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POOL = PROJECT_ROOT / "data" / "repurposing" / "screening" / "docking_pool_p05.csv"


def normalize_vina_or_glide(
    input_csv: Path,
    pdb_id: str,
    output_csv: Path,
    pool_csv: Path | None = None,
    engine: str = "auto",
) -> dict:
    raw = pd.read_csv(input_csv, low_memory=False)
    rid_col = pick_col(raw.columns, REPURPOSING_ALIASES)
    smi_col = pick_col(raw.columns, SMILES_ALIASES)
    score_col = pick_col(raw.columns, DOCK_SCORE_ALIASES)
    engine_col = pick_col(raw.columns, ["docking_engine", "engine", "tool"])

    if score_col is None:
        raise ValueError(f"No score column in {input_csv}")

    out = raw.copy()
    out = ensure_dock_score_column(out, score_col)

    if smi_col:
        out["canonical_smiles"] = out[smi_col].map(lambda s: canonicalize(s) if pd.notna(s) else None)
    elif pool_csv and rid_col:
        pool = pd.read_csv(pool_csv, low_memory=False)
        out = out.merge(
            pool[["repurposing_id", "canonical_smiles", "name", "chembl_id"]],
            left_on=rid_col,
            right_on="repurposing_id",
            how="left",
        )
    else:
        raise ValueError("Need SMILES column or pool CSV with repurposing_id join")

    if rid_col and "repurposing_id" not in out.columns:
        out["repurposing_id"] = out[rid_col]

    status_col = pick_col(raw.columns, ["docking_status", "status", "pose_status"])
    out["docking_status"] = docking_status_from_score(
        out["dock_score"],
        out[status_col] if status_col else None,
    )
    out["pdb_id"] = pdb_id.upper()
    if engine_col:
        out["docking_engine"] = out[engine_col].astype(str)
    elif engine != "auto":
        out["docking_engine"] = engine
    else:
        out["docking_engine"] = "vina" if "vina" in input_csv.name.lower() else "glide_legacy"

    out = out.sort_values("dock_score", ascending=True, na_position="last")
    out = out.groupby("canonical_smiles", as_index=False).first()

    cols = [
        "repurposing_id",
        "name",
        "chembl_id",
        "canonical_smiles",
        "dock_score",
        "glide_score_xp",
        "docking_status",
        "pdb_id",
        "docking_engine",
    ]
    cols = [c for c in cols if c in out.columns]
    result = out[cols].copy()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)

    summary = {
        "input": str(input_csv),
        "output": str(output_csv),
        "pdb_id": pdb_id.upper(),
        "docking_engine": str(result["docking_engine"].iloc[0]) if len(result) else engine,
        "n_input_rows": int(len(raw)),
        "n_output_rows": int(len(result)),
        "n_docked": int((result["docking_status"] == "docked").sum()),
        "dock_score_min": float(result["dock_score"].min()) if result["dock_score"].notna().any() else None,
        "dock_score_max": float(result["dock_score"].max()) if result["dock_score"].notna().any() else None,
        "score_column_used": score_col,
    }
    output_csv.with_suffix(".qc.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize Vina/smina/Glide docking CSV")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--pdb", type=str, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pool", type=Path, default=DEFAULT_POOL)
    parser.add_argument("--engine", type=str, default="auto", help="vina | smina | glide_legacy | auto")
    args = parser.parse_args()
    summary = normalize_vina_or_glide(args.input, args.pdb, args.output, args.pool, args.engine)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
