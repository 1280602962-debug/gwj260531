#!/usr/bin/env python3
"""
Non-docking computational module F — chemistry-aware candidate nomination.

Answers "besides EGCG, what other candidates are there, and how do we find them?"
without re-docking. It relaxes the thin Pareto front to a top-k% dual percentile
gate, then applies transparent, reviewer-defensible chemistry filters and ranks
survivors so cleaner dual-node candidates surface — not docking-score giants
(macrolides / polyketides) that inflate contact-based scores.

Selection funnel (on EXISTING dual-dock + ML percentile data):
  1. Dual gate        : S_U percentile >= tau AND S_N percentile >= tau
  2. Structural-alert : no PAINS, no Brenk                         (module B)
  3. Drug-likeness    : Lipinski pass AND Veber pass               (module C)
  4. Oral MW window   : mw_min <= MW <= mw_max (default 200–550 Da)
                        — hard demotion for macrolide-scale actives
  5. Evidence type    : flag NLRP3 structure-supported vs ML-only
  6. Rank             : preferred chemistry -> dual-structure balance ->
                        phase -> QED; then scaffold-diversify top-N

Known benchmark uricosurics/tools are labelled (positive controls), so genuinely
NEW repurposing candidates can be read off separately.

Inputs (read-only):
  data/repurposing/pareto/pareto_merged_scores.csv
  results/cheminformatics/filters_pool.csv
  results/cheminformatics/admet_pool.csv
  results/cheminformatics/novelty_pool.csv   (optional)

Outputs:
  results/candidates/nominated_candidates.csv
  results/candidates/nominated_shortlist_diverse.csv
  results/candidates/candidate_nomination_summary.json

Usage:
  python3 scripts/14_candidate_nomination.py --tau 90
  python3 scripts/14_candidate_nomination.py --tau 90 --mw-max 550 --top-diverse 12
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARETO_DIR = PROJECT_ROOT / "data" / "repurposing" / "pareto"
CHEM_DIR = PROJECT_ROOT / "results" / "cheminformatics"
OUT_DIR = PROJECT_ROOT / "results" / "candidates"

# Known benchmark / reference drugs (positive controls, not novel hits)
KNOWN_REFERENCE = {
    "LESINURAD", "VERINURAD", "DOTINURAD", "BENZBROMARONE", "COLCHICINE",
    "ALLOPURINOL", "FEBUXOSTAT", "PROBENECID", "MCC950", "GDC-2394", "SULFINPYRAZONE",
}

# Soft URAT1-relevant carboxylic acid cue (ranking only; not a hard filter).
_ACID_RE = re.compile(r"C\(=O\)O|C\(=O\)\[O-\]|C\(O\)=O")


def _has_carboxylic_acid(smi: object) -> bool:
    if not isinstance(smi, str) or not smi:
        return False
    return bool(_ACID_RE.search(smi))


def diversify_by_scaffold(df: pd.DataFrame, n: int, scaffold_col: str = "scaffold") -> pd.DataFrame:
    """Greedy keep first row per Murcko scaffold (df must already be rank-sorted)."""
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    seen: set[str] = set()
    keep_idx: list[int] = []
    for idx, row in df.iterrows():
        sc = row.get(scaffold_col)
        if pd.isna(sc) or sc is None or str(sc).strip() == "":
            sc = f"__name__:{row.get('name', idx)}"
        sc = str(sc)
        if sc in seen:
            continue
        seen.add(sc)
        keep_idx.append(idx)
        if len(keep_idx) >= n:
            break
    out = df.loc[keep_idx].copy().reset_index(drop=True)
    out.insert(0, "diverse_rank", range(1, len(out) + 1))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Candidate nomination (chemistry-aware module F)")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--filters", type=Path, default=CHEM_DIR / "filters_pool.csv")
    parser.add_argument("--admet", type=Path, default=CHEM_DIR / "admet_pool.csv")
    parser.add_argument("--novelty", type=Path, default=CHEM_DIR / "novelty_pool.csv")
    parser.add_argument("--tau", type=float, default=90.0, help="Dual percentile gate")
    parser.add_argument("--mw-min", type=float, default=200.0, help="Oral-ish MW lower bound (Da)")
    parser.add_argument("--mw-max", type=float, default=550.0, help="Oral-ish MW upper bound (Da)")
    parser.add_argument("--top-diverse", type=int, default=12, help="Scaffold-diverse shortlist size")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pool)
    f = pd.read_csv(args.filters)
    a = pd.read_csv(args.admet)

    # Prefer ADMET mw/qed when present; keep pool scaffold / smiles for diversity + acid cue.
    m = df.merge(
        f[["name", "pains_any", "brenk", "nih", "aggregation_risk_heuristic"]],
        on="name", how="left",
    ).merge(
        a[["name", "qed", "mw", "clogp", "tpsa", "hbd", "hba", "heavy_atoms",
           "lipinski_pass", "veber_pass", "ghose_pass", "oral_absorption_ok"]],
        on="name", how="left", suffixes=("_pool", ""),
    )
    if "mw" not in m.columns and "mw_pool" in m.columns:
        m["mw"] = m["mw_pool"]
    elif "mw_pool" in m.columns:
        m["mw"] = m["mw"].fillna(m["mw_pool"])

    if args.novelty.exists():
        nov = pd.read_csv(args.novelty)
        keep = [c for c in ("name", "nn_tanimoto_urat1_active", "nn_tanimoto_nlrp3_active") if c in nov.columns]
        m = m.merge(nov[keep], on="name", how="left")

    # 1. Dual gate
    gate = m[(m["s_u_percentile"] >= args.tau) & (m["s_n_percentile"] >= args.tau)].copy()

    # 2+3. Clean structural + classic drug-likeness
    for col in ("pains_any", "brenk"):
        gate[col] = gate[col].fillna(False).astype(bool)
    gate["structurally_clean"] = ~gate["pains_any"] & ~gate["brenk"]
    gate["druglike"] = (gate["lipinski_pass"] == True) & (gate["veber_pass"] == True)  # noqa: E712
    gate["clean_candidate"] = gate["structurally_clean"] & gate["druglike"]

    # 4. Oral MW window — demotes macrolides / oversized polyketides that dominate docking Pareto
    mw = pd.to_numeric(gate["mw"], errors="coerce")
    gate["mw"] = mw
    gate["mw_oral_ok"] = mw.between(args.mw_min, args.mw_max, inclusive="both")
    oral_ok = gate["oral_absorption_ok"].fillna(False).astype(bool) if "oral_absorption_ok" in gate.columns else False
    gate["oral_absorption_ok"] = oral_ok
    gate["preferred_candidate"] = (
        gate["clean_candidate"] & gate["mw_oral_ok"].fillna(False) & gate["oral_absorption_ok"]
    )

    # Soft URAT1-relevant acid cue (ranking only)
    smi_col = "canonical_smiles" if "canonical_smiles" in gate.columns else None
    if smi_col:
        gate["has_carboxylic_acid"] = gate[smi_col].map(_has_carboxylic_acid)
    else:
        gate["has_carboxylic_acid"] = False

    # 5. Evidence type for NLRP3 axis
    gate["nlrp3_structure_supported"] = gate["s_n_dock_percentile"] >= args.tau
    gate["nlrp3_ml_only_risk"] = (gate["s_n_ml_percentile"] >= args.tau) & (gate["s_n_dock_percentile"] < args.tau)

    # Reference vs novel
    gate["is_known_reference"] = gate["name"].astype(str).str.upper().isin(KNOWN_REFERENCE)

    # 6. Ranking: chemistry first, then dual-structure balance — NOT raw docking score
    gate["dual_structure_balance"] = gate[["s_u_percentile", "s_n_dock_percentile"]].min(axis=1)
    # Mild size preference toward typical oral small-molecule MW (~350 Da)
    gate["mw_centrality"] = 1.0 - ((gate["mw"] - 350.0).abs() / 200.0).clip(lower=0.0, upper=1.0)
    gate["chemistry_rank_score"] = (
        0.40 * gate["qed"].fillna(0.0)
        + 0.25 * gate["mw_centrality"].fillna(0.0)
        + 0.20 * gate["ghose_pass"].fillna(False).astype(float)
        + 0.15 * gate["has_carboxylic_acid"].astype(float)
    )

    gate = gate.sort_values(
        [
            "preferred_candidate",
            "clean_candidate",
            "mw_oral_ok",
            "is_known_reference",
            "nlrp3_structure_supported",
            "max_phase",
            "dual_structure_balance",
            "chemistry_rank_score",
            "qed",
        ],
        ascending=[False, False, False, True, False, False, False, False, False],
    ).reset_index(drop=True)

    cols = [
        "name", "chembl_id", "max_phase", "canonical_smiles", "scaffold", "mw",
        "s_u_percentile", "s_n_percentile",
        "s_n_ml_percentile", "s_n_dock_percentile", "p_active_nlrp3",
        "pains_any", "brenk", "qed", "lipinski_pass", "veber_pass", "ghose_pass",
        "oral_absorption_ok", "mw_oral_ok",
        "structurally_clean", "druglike", "clean_candidate", "preferred_candidate",
        "has_carboxylic_acid",
        "nlrp3_structure_supported", "nlrp3_ml_only_risk", "is_known_reference",
        "dual_structure_balance", "chemistry_rank_score",
    ]
    cols = [c for c in cols if c in gate.columns]
    out = gate[cols]
    out.to_csv(args.output_dir / "nominated_candidates.csv", index=False)

    preferred_novel = gate[(gate["preferred_candidate"]) & (~gate["is_known_reference"])]
    clean_novel = gate[(gate["clean_candidate"]) & (~gate["is_known_reference"])]
    clean_novel_struct = clean_novel[clean_novel["nlrp3_structure_supported"]]
    preferred_novel_struct = preferred_novel[preferred_novel["nlrp3_structure_supported"]]

    # Scaffold-diverse follow-up shortlist (what MD / story picking should read)
    diverse_pool = preferred_novel if len(preferred_novel) else clean_novel
    # Prefer structure-supported, later clinical phase, then dual-structure / chemistry scores
    diverse_pool = diverse_pool.sort_values(
        [
            "nlrp3_structure_supported",
            "max_phase",
            "dual_structure_balance",
            "chemistry_rank_score",
            "qed",
        ],
        ascending=[False, False, False, False, False],
    )
    diverse = diversify_by_scaffold(diverse_pool, n=args.top_diverse, scaffold_col="scaffold")
    diverse_cols = [c for c in (["diverse_rank"] + cols) if c in diverse.columns]
    diverse[diverse_cols].to_csv(args.output_dir / "nominated_shortlist_diverse.csv", index=False)

    summary = {
        "module": "F_candidate_nomination",
        "tau": args.tau,
        "mw_min": args.mw_min,
        "mw_max": args.mw_max,
        "ranking_note": (
            "preferred_candidate = clean (no PAINS/Brenk, Lipinski+Veber) "
            "+ oral MW window + oral_absorption_ok; "
            "scaffold-diverse shortlist is the follow-up list — not raw Pareto / docking score."
        ),
        "n_dual_gate": int(len(gate)),
        "n_clean_candidate": int(gate["clean_candidate"].sum()),
        "n_preferred_candidate": int(gate["preferred_candidate"].sum()),
        "n_clean_novel": int(len(clean_novel)),
        "n_preferred_novel": int(len(preferred_novel)),
        "n_clean_novel_structure_supported": int(len(clean_novel_struct)),
        "n_preferred_novel_structure_supported": int(len(preferred_novel_struct)),
        "n_diverse_shortlist": int(len(diverse)),
        "top_preferred_novel": preferred_novel[[
            "name", "max_phase", "mw", "s_u_percentile", "s_n_dock_percentile",
            "p_active_nlrp3", "qed", "chemistry_rank_score", "nlrp3_structure_supported",
        ]].head(12).to_dict(orient="records") if len(preferred_novel) else [],
        "top_diverse_shortlist": diverse[[
            c for c in (
                "diverse_rank", "name", "max_phase", "mw", "scaffold",
                "dual_structure_balance", "qed", "chemistry_rank_score",
            ) if c in diverse.columns
        ]].to_dict(orient="records") if len(diverse) else [],
        "known_reference_in_gate": gate[gate["is_known_reference"]][[
            "name", "max_phase", "mw", "s_u_percentile", "s_n_percentile", "preferred_candidate",
        ]].to_dict(orient="records"),
        "demoted_high_mw_in_gate": gate[
            (~gate["mw_oral_ok"].fillna(False)) & gate["clean_candidate"]
        ][["name", "mw", "s_u_percentile", "s_n_percentile"]].head(20).to_dict(orient="records"),
    }
    with open(args.output_dir / "candidate_nomination_summary.json", "w") as f_:
        json.dump(summary, f_, indent=2)

    print(f"=== Candidate nomination (tau={args.tau}, MW {args.mw_min}-{args.mw_max}) ===")
    print(f"  dual gate: {len(gate)}")
    print(f"  clean candidates: {int(gate['clean_candidate'].sum())}")
    print(f"  preferred (clean + oral MW + absorption): {int(gate['preferred_candidate'].sum())}")
    print(f"  preferred & novel: {len(preferred_novel)}")
    print(f"  preferred & novel & NLRP3 structure-supported: {len(preferred_novel_struct)}")
    print(f"  scaffold-diverse shortlist: {len(diverse)}")
    print("\nTop preferred novel candidates (chemistry-first, not docking-only):")
    show = ["name", "max_phase", "mw", "s_u_percentile", "s_n_dock_percentile",
            "qed", "chemistry_rank_score", "nlrp3_structure_supported"]
    show = [c for c in show if c in preferred_novel.columns]
    if len(preferred_novel):
        print(preferred_novel[show].head(12).to_string(index=False))
    else:
        print("  (none — relax --mw-max or inspect clean_candidate demotions)")
    print("\nScaffold-diverse follow-up shortlist:")
    dshow = [c for c in ("diverse_rank", "name", "mw", "dual_structure_balance", "qed") if c in diverse.columns]
    if len(diverse):
        print(diverse[dshow].to_string(index=False))
    print("\nKnown reference drugs in gate (positive controls):")
    print(gate[gate["is_known_reference"]][[
        "name", "max_phase", "mw", "s_u_percentile", "s_n_percentile", "preferred_candidate",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
