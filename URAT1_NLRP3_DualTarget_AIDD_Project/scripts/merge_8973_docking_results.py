#!/usr/bin/env python3
"""
Merge Vina/smina (or legacy Glide) 8973 @ 9DKB exports with distill_manifest.csv.

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

from dock_score_utils import DOCK_SCORE_ALIASES, pick_col
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
SCORE_ALIASES = DOCK_SCORE_ALIASES
STATUS_ALIASES = [
    "docking_status",
    "pose",
    "pose_status",
    "glide pose",
    "r_i_glide_pose",
    "status",
]
NAME_ALIASES = ["title", "s_m_title", "name", "compound_name", "ligand"]
TITLE_ALIASES = ["title", "s_m_title", "name"]


def inchikey_block1(smiles: str) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return None


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
    title_col = _pick_col(cols, TITLE_ALIASES)

    out = pd.DataFrame()
    out["smiles_raw"] = df[smi_col].astype(str)
    out["dock_score"] = pd.to_numeric(df[score_col], errors="coerce")
    out["glide_score_xp"] = out["dock_score"]
    out["pose_status"] = df[status_col].astype(str) if status_col else "unknown"
    out["compound_name"] = df[name_col].astype(str) if name_col else out["smiles_raw"]
    out["dock_title"] = df[title_col].astype(str) if title_col else None
    out["source_file"] = path.name
    return out


def prepare_dock_table(df: pd.DataFrame) -> pd.DataFrame:
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
    df["inchikey_block1"] = df["canonical_smiles"].map(inchikey_block1)
    return df.reset_index(drop=True)


def best_pose_per_smiles(df: pd.DataFrame) -> pd.DataFrame:
    df = prepare_dock_table(df)
    # Glide: more negative = better binding
    df = df.sort_values("glide_score_xp", ascending=True)
    df = df.drop_duplicates("canonical_smiles", keep="first")
    return df.reset_index(drop=True)


def _parse_dock_title(title: str) -> tuple[str, int] | None:
    m = re.match(r"^([ADEC])(\d+)$", str(title).strip())
    if not m:
        return None
    return m.group(1), int(m.group(2))


def build_title_index_map(manifest: pd.DataFrame) -> dict[str, str]:
    """Map dock title (e.g. A1) -> manifest canonical_smiles by subset row order."""
    mapping: dict[str, str] = {}
    for sub in manifest["subset"].dropna().unique():
        sub_df = manifest[manifest["subset"] == sub].reset_index(drop=True)
        for i, row in sub_df.iterrows():
            mapping[f"{sub}{i + 1}"] = row["canonical_smiles"]
    return mapping


def merge_manifest(dock_raw: pd.DataFrame, manifest: pd.DataFrame, pdb: str) -> pd.DataFrame:
    dock = prepare_dock_table(dock_raw)
    title_map = build_title_index_map(manifest)

    # Best pose per title (Maestro may repeat titles only once)
    dock_by_title = (
        dock[dock["dock_title"].notna()]
        .sort_values("glide_score_xp", ascending=True)
        .drop_duplicates("dock_title", keep="first")
    )
    title_scores = {
        str(r["dock_title"]): float(r["glide_score_xp"])
        for _, r in dock_by_title.iterrows()
        if _parse_dock_title(r["dock_title"])
    }

    ik_scores = (
        dock.sort_values("glide_score_xp", ascending=True)
        .drop_duplicates("inchikey_block1", keep="first")
        .set_index("inchikey_block1")["glide_score_xp"]
        .to_dict()
    )
    canon_scores = (
        dock.sort_values("glide_score_xp", ascending=True)
        .drop_duplicates("canonical_smiles", keep="first")
        .set_index("canonical_smiles")["glide_score_xp"]
        .to_dict()
    )

    m = manifest.copy()
    m["inchikey_block1"] = m["canonical_smiles"].map(inchikey_block1)
    glide = []
    join_method = []
    dock_title_out = []

    for sub_idx, row in m.iterrows():
        sub = row["subset"]
        # position within subset for title lookup
        sub_rows = m[m["subset"] == sub]
        pos = list(sub_rows.index).index(sub_idx) + 1
        title_key = f"{sub}{pos}"
        score = None
        method = None

        if title_key in title_scores:
            score = title_scores[title_key]
            method = "title"
        elif pd.notna(row["inchikey_block1"]) and row["inchikey_block1"] in ik_scores:
            score = ik_scores[row["inchikey_block1"]]
            method = "inchikey"
        elif row["canonical_smiles"] in canon_scores:
            score = canon_scores[row["canonical_smiles"]]
            method = "canonical_smiles"

        glide.append(score)
        join_method.append(method)
        dock_title_out.append(title_key if method == "title" else None)

    merged = m.copy()
    merged["dock_title"] = dock_title_out
    merged["join_method"] = join_method
    merged["glide_score_xp"] = glide
    merged["pdb_id"] = pdb
    merged["docked"] = merged["glide_score_xp"].notna()
    merged["s_u_raw"] = -merged["glide_score_xp"]
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
            if not len(hit) or not hit["docked"].any():
                # fallback: inchikey match for benchmark
                ik = inchikey_block1(smi)
                if ik:
                    hit = merged[(merged["inchikey_block1"] == ik) & merged["docked"]]
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
        "join_method_counts": merged["join_method"].value_counts(dropna=False).to_dict(),
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
    merged = merge_manifest(dock_raw, manifest, args.pdb)

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
