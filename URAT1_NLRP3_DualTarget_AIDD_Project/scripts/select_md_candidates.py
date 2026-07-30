#!/usr/bin/env python3
"""
Select MD-candidate molecules from the P2 dual-target funnel outputs.

Selection is NOT a re-run of protocol choice; it only ranks already-nominated
candidates (results/candidates/nominated_candidates.csv) plus known reference
compounds, and recovers repurposing_id / canonical_smiles / docking status by
joining back to the Pareto merged-scores table (which carries those columns).

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
        "docking_status",
        "pdb_id",
        "nlrp3_docking_status",
        "nlrp3_pdb_id",
    ]
    keep_cols = [c for c in keep_cols if c in pareto.columns]
    pareto_sub = pareto[keep_cols].drop_duplicates(subset=["chembl_id", "name"], keep="first")

    merged = nom.merge(pareto_sub, on=["chembl_id", "name"], how="left", validate="one_to_one")
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


def select_candidates(
    df: pd.DataFrame,
    n_novel: int,
    n_controls: int,
) -> pd.DataFrame:
    df = df[df["has_9dkb_pose"] | df["has_7alv_pose"]].copy()

    novel_pool = df[
        (df["is_known_reference"] == False)  # noqa: E712
        & (df["clean_candidate"] == True)  # noqa: E712
    ].sort_values("dual_structure_balance", ascending=False)
    novel = novel_pool.head(n_novel).copy()
    novel["md_category"] = "novel_candidate"

    control_pool = df[df["is_known_reference"] == True].copy()  # noqa: E712
    control_pool["_priority"] = control_pool["name"].map(_norm_name).apply(
        lambda n: CONTROL_PRIORITY.index(n) if n in CONTROL_PRIORITY else len(CONTROL_PRIORITY)
    )
    control_pool = control_pool.sort_values(["_priority", "dual_structure_balance"], ascending=[True, False])
    controls = control_pool.head(n_controls).copy()
    controls["md_category"] = "known_control"

    out = pd.concat([novel, controls], ignore_index=True)
    out = out.drop(columns=["_priority"], errors="ignore")
    out.insert(0, "md_rank", range(1, len(out) + 1))

    def rationale(row: pd.Series) -> str:
        if row["md_category"] == "known_control":
            return "Known reference compound; positive-control MD system for calibration."
        return (
            f"clean_candidate=True, no PAINS/Brenk, Lipinski+Veber pass, "
            f"dual_structure_balance={row['dual_structure_balance']:.1f}"
        )

    out["selection_rationale"] = out.apply(rationale, axis=1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Select MD-candidate molecules (no MD is run here)")
    parser.add_argument("--nominated", type=Path, default=DEFAULT_NOMINATED)
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

    joined = load_and_join(args.nominated, args.pareto)
    selection = select_candidates(joined, args.n_novel, args.n_controls)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "md_rank",
        "md_category",
        "repurposing_id",
        "chembl_id",
        "name",
        "canonical_smiles",
        "s_u_percentile",
        "s_n_percentile",
        "dual_structure_balance",
        "clean_candidate",
        "is_known_reference",
        "has_9dkb_pose",
        "has_7alv_pose",
        "selection_rationale",
    ]
    cols = [c for c in cols if c in selection.columns]
    selection[cols].to_csv(args.output, index=False)

    print(f"Selected {len(selection)} MD candidates -> {args.output}")
    print(selection[["md_rank", "md_category", "name", "has_9dkb_pose", "has_7alv_pose"]].to_string(index=False))


if __name__ == "__main__":
    main()
