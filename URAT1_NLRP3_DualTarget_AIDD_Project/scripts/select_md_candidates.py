#!/usr/bin/env python3
"""
Select MD-candidate molecules from the chemistry-aware nomination outputs.

Selection is NOT a re-run of protocol choice; it ranks already-nominated
candidates (results/candidates/nominated_candidates.csv), preferring
preferred_candidate (clean + oral MW window + absorption) and scaffold
diversity, then recovers repurposing_id / canonical_smiles / docking status
by joining back to the Pareto merged-scores table.

Do NOT pick MD leads from raw pareto_shortlist.csv (docking-only; often macrolides).
Do NOT use the committed Glide-era CSV in data/md_candidates/ as the current list.
Current intended MD set is documented in docs/MANUSCRIPT.md (GSK-3008348, Vecabrutinib + controls).

Output: a compact CSV that scripts/export_md_ready_candidates.py consumes to
produce receptor + ligand files for external MD (no MD is run here).

Example:
  python3 scripts/select_md_candidates.py \\
    --n-novel 4 --n-controls 2 \\
    --output data/md_candidates/md_candidate_selection.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NOMINATED = PROJECT_ROOT / "results" / "candidates" / "nominated_candidates.csv"
DEFAULT_DIVERSE = PROJECT_ROOT / "results" / "candidates" / "nominated_shortlist_diverse.csv"
DEFAULT_PARETO = PROJECT_ROOT / "data" / "repurposing" / "pareto" / "pareto_merged_scores.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "md_candidates" / "md_candidate_selection.csv"

# Preferred known controls if present among is_known_reference rows (URAT1 / NLRP3 tool compounds).
CONTROL_PRIORITY = [
    "LESINURAD",
    "MCC950",
    "DOTINURAD",
    "BENZBROMARONE",
    "VERINURAD",
    "GDC-2394",
]


def _norm_name(s: object) -> str:
    return str(s).strip().upper()


def load_and_join(nominated_path: Path, pareto_path: Path) -> pd.DataFrame:
    nom = pd.read_csv(nominated_path)
    pareto = pd.read_csv(pareto_path, low_memory=False)

    keep_cols = [
        "repurposing_id",
        "chembl_id",
        "name",
        "canonical_smiles",
        "scaffold",
        "docking_status",
        "pdb_id",
        "nlrp3_docking_status",
        "nlrp3_pdb_id",
    ]
    keep_cols = [c for c in keep_cols if c in pareto.columns]
    pareto_sub = pareto[keep_cols].drop_duplicates(subset=["chembl_id", "name"], keep="first")

    # nominated may already carry canonical_smiles/scaffold; suffix pareto copies
    overlap = [c for c in ("canonical_smiles", "scaffold") if c in nom.columns and c in pareto_sub.columns]
    merged = nom.merge(
        pareto_sub,
        on=["chembl_id", "name"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_pareto"),
    )
    for c in overlap:
        if f"{c}_pareto" in merged.columns:
            merged[c] = merged[c].fillna(merged[f"{c}_pareto"])
            merged = merged.drop(columns=[f"{c}_pareto"])

    n_missing = merged["repurposing_id"].isna().sum()
    if n_missing:
        print(
            f"WARNING: {n_missing}/{len(merged)} nominated candidates could not be joined to "
            f"repurposing_id via (chembl_id, name); they will be excluded from MD selection."
        )
    merged = merged.dropna(subset=["repurposing_id"]).copy()

    merged["has_9dkb_pose"] = merged.get("docking_status", pd.Series(dtype=object)) == "docked"
    merged["has_7alv_pose"] = merged.get("nlrp3_docking_status", pd.Series(dtype=object)) == "docked"
    return merged


def _diversify(df: pd.DataFrame, n: int) -> pd.DataFrame:
    if df.empty or n <= 0:
        return df.iloc[0:0].copy()
    seen: set[str] = set()
    keep: list[int] = []
    for idx, row in df.iterrows():
        sc = row.get("scaffold")
        if pd.isna(sc) or sc is None or str(sc).strip() == "":
            sc = f"__name__:{row.get('name', idx)}"
        sc = str(sc)
        if sc in seen:
            continue
        seen.add(sc)
        keep.append(idx)
        if len(keep) >= n:
            break
    return df.loc[keep].copy()


def select_candidates(
    df: pd.DataFrame,
    n_novel: int,
    n_controls: int,
    prefer_diverse_file: pd.DataFrame | None = None,
) -> pd.DataFrame:
    df = df[df["has_9dkb_pose"] | df["has_7alv_pose"]].copy()

    # Prefer Module-F preferred_candidate; fall back to clean_candidate + mw_oral_ok
    if "preferred_candidate" in df.columns:
        novel_mask = (df["is_known_reference"] == False) & (df["preferred_candidate"] == True)  # noqa: E712
    else:
        novel_mask = (df["is_known_reference"] == False) & (df["clean_candidate"] == True)  # noqa: E712
        if "mw_oral_ok" in df.columns:
            novel_mask = novel_mask & (df["mw_oral_ok"] == True)  # noqa: E712

    sort_cols = [c for c in (
        "nlrp3_structure_supported",
        "max_phase",
        "dual_structure_balance",
        "chemistry_rank_score",
        "qed",
    ) if c in df.columns]
    ascending = [False] * len(sort_cols)
    novel_pool = df[novel_mask].sort_values(sort_cols, ascending=ascending) if sort_cols else df[novel_mask]

    # If a precomputed diverse shortlist exists, honor its order for names present here
    if prefer_diverse_file is not None and len(prefer_diverse_file) and len(novel_pool):
        order = { _norm_name(n): i for i, n in enumerate(prefer_diverse_file["name"].tolist()) }
        novel_pool = novel_pool.copy()
        novel_pool["_diverse_order"] = novel_pool["name"].map(_norm_name).map(
            lambda n: order.get(n, 10_000)
        )
        novel_pool = novel_pool.sort_values(
            ["_diverse_order"] + sort_cols,
            ascending=[True] + ascending,
        )

    novel = _diversify(novel_pool, n_novel)
    novel["md_category"] = "novel_candidate"

    control_pool = df[df["is_known_reference"] == True].copy()  # noqa: E712
    control_pool["_priority"] = control_pool["name"].map(_norm_name).apply(
        lambda n: CONTROL_PRIORITY.index(n) if n in CONTROL_PRIORITY else len(CONTROL_PRIORITY)
    )
    ctrl_sort = ["_priority"] + ([c for c in ("dual_structure_balance",) if c in control_pool.columns])
    ctrl_asc = [True] + [False] * (len(ctrl_sort) - 1)
    control_pool = control_pool.sort_values(ctrl_sort, ascending=ctrl_asc)
    controls = control_pool.head(n_controls).copy()
    controls["md_category"] = "known_control"

    out = pd.concat([novel, controls], ignore_index=True)
    out = out.drop(columns=["_priority", "_diverse_order"], errors="ignore")
    out.insert(0, "md_rank", range(1, len(out) + 1))

    def rationale(row: pd.Series) -> str:
        if row["md_category"] == "known_control":
            return "Known reference compound; positive-control MD system for calibration."
        bits = []
        if row.get("preferred_candidate") is True:
            bits.append("preferred_candidate (clean + oral MW + absorption)")
        elif row.get("clean_candidate") is True:
            bits.append("clean_candidate (PAINS/Brenk clear, Lipinski+Veber)")
        if "mw" in row.index and pd.notna(row["mw"]):
            bits.append(f"MW={row['mw']:.0f}")
        if "dual_structure_balance" in row.index and pd.notna(row["dual_structure_balance"]):
            bits.append(f"dual_structure_balance={row['dual_structure_balance']:.1f}")
        if "chemistry_rank_score" in row.index and pd.notna(row["chemistry_rank_score"]):
            bits.append(f"chemistry_rank_score={row['chemistry_rank_score']:.3f}")
        bits.append("scaffold-diversified")
        return ", ".join(bits)

    out["selection_rationale"] = out.apply(rationale, axis=1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Select MD-candidate molecules (no MD is run here)")
    parser.add_argument("--nominated", type=Path, default=DEFAULT_NOMINATED)
    parser.add_argument("--diverse", type=Path, default=DEFAULT_DIVERSE)
    parser.add_argument("--pareto", type=Path, default=DEFAULT_PARETO)
    parser.add_argument("--n-novel", type=int, default=4)
    parser.add_argument("--n-controls", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if not args.nominated.exists():
        raise FileNotFoundError(
            f"{args.nominated} not found. Run scripts/14_candidate_nomination.py first "
            "(after the P2 dual-target docking + Pareto merge)."
        )
    if not args.pareto.exists():
        raise FileNotFoundError(
            f"{args.pareto} not found. Run scripts/run_funnel_p2.sh first."
        )

    diverse = pd.read_csv(args.diverse) if args.diverse.exists() else None
    joined = load_and_join(args.nominated, args.pareto)
    selection = select_candidates(joined, args.n_novel, args.n_controls, prefer_diverse_file=diverse)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "md_rank",
        "md_category",
        "repurposing_id",
        "chembl_id",
        "name",
        "canonical_smiles",
        "mw",
        "scaffold",
        "s_u_percentile",
        "s_n_percentile",
        "dual_structure_balance",
        "chemistry_rank_score",
        "preferred_candidate",
        "clean_candidate",
        "is_known_reference",
        "has_9dkb_pose",
        "has_7alv_pose",
        "selection_rationale",
    ]
    cols = [c for c in cols if c in selection.columns]
    selection[cols].to_csv(args.output, index=False)

    print(f"Selected {len(selection)} MD candidates -> {args.output}")
    show = [c for c in ("md_rank", "md_category", "name", "mw", "has_9dkb_pose", "has_7alv_pose") if c in selection.columns]
    print(selection[show].to_string(index=False))


if __name__ == "__main__":
    main()
