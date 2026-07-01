#!/usr/bin/env python3
"""
Merge Maestro/Glide 8973 @ 9DKB XP exports with distill_manifest.csv.

Inputs (flexible Maestro CSV):
  - One or more Glide XP score tables (per-ligand best pose kept)

Outputs:
  data/docking/8973_9DKB_merged.csv
  data/docking/8973_9DKB_with_manifest.csv
  data/docking/8973_docking_qc_summary.json

Example:
  python3 scripts/merge_8973_docking_results.py \\
    --glide-csv results/docking/raw/9DKB_xp_scores.csv \\
    --pdb 9DKB
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from utils_ml import canonicalize

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT_ROOT / "data" / "distill" / "distill_manifest.csv"
BENCHMARK_PANEL = PROJECT_ROOT / "data" / "distill" / "teacher_gate_qc_panel_b_direction.csv"
DEFAULT_OUT = PROJECT_ROOT / "data" / "docking"

SMILES_ALIASES = [
    "smiles",
    "canonical_smiles",
    "ligprep_smiles",
    "r_m_chemaxon_smiles",
    "s_m_entry_name",
]
SCORE_ALIASES = [
    "r_i_glide_xp",
    "r_i_glide xp",
    "glide xp gscore",
    "glide_score_xp",
    "glide_xp",
    "xp gscore",
    "docking score",
    "r_i_docking_score",
]
STATUS_ALIASES = [
    "pose",
    "pose_status",
    "glide pose",
    "r_i_glide_pose",
    "status",
]
NAME_ALIASES = ["title", "s_m_title", "name", "compound_name", "ligand"]


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip().lower())


def _pick_col(columns: list[str], aliases: list[str]) -> str | None:
    norm = {_norm(c): c for c in columns}
    for a in aliases:
        if a in norm:
            return norm[a]
    for c in columns:
        cn = _norm(c)
        for a in aliases:
            if a.replace(" ", "") in cn.replace(" ", ""):
                return c
    return None


def read_glide_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xls", ".xlsx"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, low_memory=False)

    cols = list(df.columns)
    smi_col = _pick_col(cols, SMILES_ALIASES)
    score_col = _pick_col(cols, SCORE_ALIASES)
    if smi_col is None or score_col is None:
        raise ValueError(
            f"Cannot map SMILES/score in {path.name}. Columns: {cols[:25]}"
        )

    status_col = _pick_col(cols, STATUS_ALIASES)
    name_col = _pick_col(cols, NAME_ALIASES)

    out = pd.DataFrame()
    out["smiles_raw"] = df[smi_col].astype(str)
    out["glide_score_xp"] = pd.to_numeric(df[score_col], errors="coerce")
    out["pose_status"] = df[status_col].astype(str) if status_col else "unknown"
    out["compound_name"] = df[name_col].astype(str) if name_col else out["smiles_raw"]
    out["source_file"] = path.name
    return out


def best_pose_per_smiles(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["canonical_smiles"] = df["smiles_raw"].map(
        lambda s: canonicalize(s) if pd.notna(s) and str(s).lower() not in {"nan", ""} else None
    )
    df = df[df["canonical_smiles"].notna() & df["glide_score_xp"].notna()].copy()

    failed_mask = df["pose_status"].str.lower().str.contains(
        "fail|no pose|unconverged|skip", na=False
    )
    df.loc[failed_mask, "glide_score_xp"] = np.nan
    df = df[df["glide_score_xp"].notna()].copy()

    # Glide: more negative = better binding
    df = df.sort_values("glide_score_xp", ascending=True)
    df = df.drop_duplicates("canonical_smiles", keep="first")
    return df.reset_index(drop=True)


def merge_manifest(dock: pd.DataFrame, manifest: pd.DataFrame, pdb: str) -> pd.DataFrame:
    m = manifest.copy()
    merged = m.merge(
        dock[
            [
                "canonical_smiles",
                "glide_score_xp",
                "pose_status",
                "compound_name",
                "source_file",
            ]
        ],
        on="canonical_smiles",
        how="left",
    )
    merged["pdb_id"] = pdb
    merged["docked"] = merged["glide_score_xp"].notna()
    merged["s_u_raw"] = -merged["glide_score_xp"]  # higher = better for ranking
    docked = merged[merged["docked"]].copy()
    merged["s_u_percentile"] = np.nan
    if len(docked):
        ranks = docked["s_u_raw"].rank(method="average", pct=True) * 100.0
        merged.loc[docked.index, "s_u_percentile"] = ranks
    return merged


def qc_summary(merged: pd.DataFrame, dock: pd.DataFrame) -> dict:
    panel = pd.read_csv(BENCHMARK_PANEL) if BENCHMARK_PANEL.exists() else pd.DataFrame()
    bench_rows = []
    if not panel.empty and "canonical_smiles" in panel.columns:
        for _, r in panel.iterrows():
            smi = r["canonical_smiles"]
            hit = merged[merged["canonical_smiles"] == smi]
            bench_rows.append(
                {
                    "compound_name": r.get("compound_name", r.get("compound_id")),
                    "canonical_smiles": smi,
                    "docked": bool(hit["docked"].any()) if len(hit) else False,
                    "glide_score_xp": float(hit["glide_score_xp"].iloc[0]) if len(hit) and hit["docked"].any() else None,
                    "s_u_percentile": float(hit["s_u_percentile"].iloc[0]) if len(hit) and hit["docked"].any() else None,
                }
            )

    subset_cov = {}
    for sub in sorted(merged["subset"].dropna().unique()):
        sub_df = merged[merged["subset"] == sub]
        subset_cov[str(sub)] = {
            "n_manifest": int(len(sub_df)),
            "n_docked": int(sub_df["docked"].sum()),
            "coverage_pct": round(100.0 * sub_df["docked"].mean(), 2),
        }

    return {
        "n_manifest": int(len(merged)),
        "n_unique_docked": int(len(dock)),
        "n_merged_docked": int(merged["docked"].sum()),
        "coverage_pct": round(100.0 * merged["docked"].mean(), 2),
        "subset_coverage": subset_cov,
        "benchmark_panel_9dkb": bench_rows,
        "missing_smiles_count": int((~merged["docked"]).sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge 8973 Glide XP with distill manifest")
    parser.add_argument(
        "--glide-csv",
        type=Path,
        action="append",
        required=True,
        help="Maestro/Glide XP export CSV or Excel (repeatable)",
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--pdb", type=str, default="9DKB")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tables = [read_glide_table(p) for p in args.glide_csv]
    dock_raw = pd.concat(tables, ignore_index=True)
    dock = best_pose_per_smiles(dock_raw)

    manifest = pd.read_csv(args.manifest)
    merged = merge_manifest(dock, manifest, args.pdb)

    merged_path = args.output_dir / "8973_9DKB_with_manifest.csv"
    dock_path = args.output_dir / "8973_9DKB_merged.csv"
    merged.to_csv(merged_path, index=False)
    dock.to_csv(dock_path, index=False)

    summary = qc_summary(merged, dock)
    summary["inputs"] = [str(p) for p in args.glide_csv]
    summary["outputs"] = {"merged": str(merged_path), "dock_best_pose": str(dock_path)}
    summary_path = args.output_dir / "8973_docking_qc_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {merged_path} ({summary['n_merged_docked']}/{summary['n_manifest']} docked)")


if __name__ == "__main__":
    main()
