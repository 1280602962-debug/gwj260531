#!/usr/bin/env python3
"""
Non-docking computational module F — chemistry-aware candidate nomination.

Answers how chemistry-aware nomination is applied after P2 dual-dock percentiles,
without re-docking. Relaxes the thin Pareto front to a dual-dock percentile gate,

  S_U percentile >= tau  AND  S_N_dock percentile >= tau

ML (p_active_nlrp3) remains an annotation / soft rank signal, not the hard gate.

Selection funnel (on EXISTING dual-dock + ML percentile data):
  1. Dual-DOCK gate   : S_U >= tau AND S_N_dock >= tau   (default; --gate-mode)
  2. Structural-alert : no PAINS, no Brenk
  3. Drug-likeness    : Veber + Ro5 HBD/HBA/logP (MW via oral window, not Lipinski MW≤500)
  4. Oral MW window   : mw_min <= MW <= mw_max (default 200–550 Da)
  5. Evidence type    : flag NLRP3 structure-supported vs ML-only
  6. Rank             : preferred chemistry -> dual-structure balance ->
                        phase -> QED; then scaffold-diversify top-N

Known benchmark uricosurics/tools are labelled (positive controls).

Inputs (read-only):
  data/repurposing/p2/pareto_merged_scores.csv
  data/repurposing/p2/filters_pool.csv
  data/repurposing/p2/admet_pool.csv
  data/repurposing/p2/novelty_pool.csv   (optional)

Outputs:
  data/repurposing/p2/nominated_candidates.csv
  data/repurposing/p2/nominated_shortlist_diverse.csv
  data/repurposing/p2/candidate_nomination_summary.json

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
PARETO_DIR = PROJECT_ROOT / "data" / "repurposing" / "p2"
CHEM_DIR = PROJECT_ROOT / "data" / "repurposing" / "p2"
OUT_DIR = PROJECT_ROOT / "data" / "repurposing" / "p2"

# Known benchmark / reference drugs (positive controls, not novel hits)
KNOWN_REFERENCE = {
    "LESINURAD", "VERINURAD", "DOTINURAD", "BENZBROMARONE", "COLCHICINE",
    "ALLOPURINOL", "FEBUXOSTAT", "PROBENECID", "MCC950", "GDC-2394", "SULFINPYRAZONE",
}

# Soft URAT1-relevant carboxylic acid cue (ranking only; not a hard filter).
_ACID_RE = re.compile(r"C\(=O\)O|C\(=O\)\[O-\]|C\(O\)=O")

# Macrolide / polyketide Murcko scaffolds seen inflating the raw Pareto front.
_MACRO_SCAFFOLD_MARKERS = (
    "O=C1CC(OC2CCCCO2)CC(OC2CCCCO2)C2CC=C(CCCCO1)O2",  # erythromycin-like
    "O=C1CCC(OCC=Cc2cnc3ccccc3c2)C(OC2CCCCO2)CC(=O)CC(=O)OCC2OC(=O)NC2C1",  # cethromycin-like
    "O=C1CCCC=CCC=CCC(C=Cc2",  # epothilone-like (KOS)
)


def _has_carboxylic_acid(smi: object) -> bool:
    if not isinstance(smi, str) or not smi:
        return False
    return bool(_ACID_RE.search(smi))


def _is_macro_scaffold(scaffold: object) -> bool:
    if not isinstance(scaffold, str) or not scaffold:
        return False
    return any(m in scaffold for m in _MACRO_SCAFFOLD_MARKERS)


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
    parser.add_argument(
        "--gate-mode",
        choices=["dual_dock", "max_sn"],
        default="dual_dock",
        help="dual_dock: S_U & S_N_dock >= tau (recommended). max_sn: original max(ML,dock) S_N.",
    )
    parser.add_argument("--mw-min", type=float, default=200.0, help="Oral-ish MW lower bound (Da)")
    parser.add_argument("--mw-max", type=float, default=550.0, help="Oral-ish MW upper bound (Da)")
    parser.add_argument("--top-diverse", type=int, default=12, help="Scaffold-diverse shortlist size")
    parser.add_argument(
        "--require-structure-for-preferred",
        action="store_true",
        default=True,
        help="Preferred candidates must have NLRP3 docking percentile >= tau (default on)",
    )
    parser.add_argument(
        "--allow-ml-only-preferred",
        action="store_true",
        help="Override: allow preferred without NLRP3 docking support",
    )
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    if args.allow_ml_only_preferred:
        args.require_structure_for_preferred = False
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.filters.exists() or not args.admet.exists():
        raise SystemExit(
            f"Missing cheminformatics tables.\n"
            f"  expected: {args.filters}\n"
            f"            {args.admet}\n"
            f"Run: python3 scripts/build_cheminformatics_pool.py"
        )

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

    # 1. Dual-dock gate (S_U and S_N,dock). ML is annotation, not the hard gate.
    if args.gate_mode == "dual_dock":
        gate = m[(m["s_u_percentile"] >= args.tau) & (m["s_n_dock_percentile"] >= args.tau)].copy()
    else:
        gate = m[(m["s_u_percentile"] >= args.tau) & (m["s_n_percentile"] >= args.tau)].copy()

    # 2+3. Clean structural + oral drug-likeness
    # Note: classic Lipinski MW<=500 conflicts with the 200–550 oral window used to
    # demote macrolides while still keeping mid-size clinical compounds. We therefore
    # treat full lipinski_pass as annotation, and define druglike as Veber + Ro5
    # HBD/HBA/logP (MW handled separately by mw_oral_ok).
    for col in ("pains_any", "brenk"):
        gate[col] = gate[col].fillna(False).astype(bool)
    gate["structurally_clean"] = ~gate["pains_any"] & ~gate["brenk"]
    clogp = pd.to_numeric(gate.get("clogp"), errors="coerce")
    hbd = pd.to_numeric(gate.get("hbd"), errors="coerce")
    hba = pd.to_numeric(gate.get("hba"), errors="coerce")
    gate["ro5_props_ok"] = (
        (gate["veber_pass"] == True)  # noqa: E712
        & clogp.fillna(99).le(5.5)
        & hbd.fillna(99).le(5)
        & hba.fillna(99).le(10)
    )
    gate["druglike"] = gate["ro5_props_ok"]
    gate["lipinski_pass"] = gate["lipinski_pass"].fillna(False).astype(bool)
    gate["clean_candidate"] = gate["structurally_clean"] & gate["druglike"]

    # 4. Oral MW window — demotes macrolides / oversized polyketides
    mw = pd.to_numeric(gate["mw"], errors="coerce")
    gate["mw"] = mw
    gate["mw_oral_ok"] = mw.between(args.mw_min, args.mw_max, inclusive="both")
    oral_ok = gate["oral_absorption_ok"].fillna(False).astype(bool) if "oral_absorption_ok" in gate.columns else False
    gate["oral_absorption_ok"] = oral_ok
    gate["macro_scaffold_risk"] = gate["scaffold"].map(_is_macro_scaffold) if "scaffold" in gate.columns else False

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

    preferred = (
        gate["clean_candidate"]
        & gate["mw_oral_ok"].fillna(False)
        & gate["oral_absorption_ok"]
        & ~gate["macro_scaffold_risk"].fillna(False)
    )
    if args.require_structure_for_preferred:
        preferred = preferred & gate["nlrp3_structure_supported"]
    gate["preferred_candidate"] = preferred

    # 6. Ranking: chemistry first, then dual-structure balance — NOT raw docking score
    gate["dual_structure_balance"] = gate[["s_u_percentile", "s_n_dock_percentile"]].min(axis=1)
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
        "dock_score", "nlrp3_dock_score",
        "pains_any", "brenk", "qed", "lipinski_pass", "veber_pass", "ghose_pass",
        "ro5_props_ok", "oral_absorption_ok", "mw_oral_ok", "macro_scaffold_risk",
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

    # Scaffold-diverse follow-up shortlist
    diverse_pool = preferred_novel if len(preferred_novel) else clean_novel
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
        "gate_mode": args.gate_mode,
        "mw_min": args.mw_min,
        "mw_max": args.mw_max,
        "require_structure_for_preferred": args.require_structure_for_preferred,
        "ranking_note": (
            "gate_mode=dual_dock uses S_U & S_N_dock (not max with ML). "
            "preferred_candidate = clean (no PAINS/Brenk, Veber + Ro5 HBD/HBA/logP) "
            "+ oral MW window + oral_absorption_ok + not macro scaffold"
            + (" + NLRP3 dock >= tau" if args.require_structure_for_preferred else "")
            + "; scaffold-diverse shortlist is the follow-up list — not raw Pareto."
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
            c for c in (
                "name", "max_phase", "mw", "s_u_percentile", "s_n_dock_percentile",
                "p_active_nlrp3", "qed", "chemistry_rank_score", "nlrp3_structure_supported",
            ) if c in preferred_novel.columns
        ]].head(12).to_dict(orient="records") if len(preferred_novel) else [],
        "top_diverse_shortlist": diverse[[
            c for c in (
                "diverse_rank", "name", "max_phase", "mw", "scaffold",
                "dual_structure_balance", "qed", "chemistry_rank_score",
                "s_u_percentile", "s_n_dock_percentile",
            ) if c in diverse.columns
        ]].to_dict(orient="records") if len(diverse) else [],
        "known_reference_in_gate": gate[gate["is_known_reference"]][[
            c for c in (
                "name", "max_phase", "mw", "s_u_percentile", "s_n_percentile",
                "s_n_dock_percentile", "preferred_candidate",
            ) if c in gate.columns
        ]].to_dict(orient="records"),
        "demoted_high_mw_in_gate": gate[
            (~gate["mw_oral_ok"].fillna(False)) & gate["clean_candidate"]
        ][["name", "mw", "s_u_percentile", "s_n_dock_percentile"]].head(20).to_dict(orient="records")
        if "s_n_dock_percentile" in gate.columns else [],
        "raw_pareto_front_note": (
            "Raw docking Pareto (4 hits) is dominated by macrolide/erythromycin scaffolds "
            "and is NOT the nomination list."
        ),
    }
    with open(args.output_dir / "candidate_nomination_summary.json", "w") as f_:
        json.dump(summary, f_, indent=2)

    print(f"=== Candidate nomination (gate={args.gate_mode}, tau={args.tau}, MW {args.mw_min}-{args.mw_max}) ===")
    print(f"  dual gate: {len(gate)}")
    print(f"  clean candidates: {int(gate['clean_candidate'].sum())}")
    print(f"  preferred (clean + oral MW + absorption + structure): {int(gate['preferred_candidate'].sum())}")
    print(f"  preferred & novel: {len(preferred_novel)}")
    print(f"  preferred & novel & NLRP3 structure-supported: {len(preferred_novel_struct)}")
    print(f"  scaffold-diverse shortlist: {len(diverse)}")
    print("\nTop preferred novel candidates (chemistry-first, dual-dock gated):")
    show = ["name", "max_phase", "mw", "s_u_percentile", "s_n_dock_percentile",
            "qed", "chemistry_rank_score", "nlrp3_structure_supported", "scaffold"]
    show = [c for c in show if c in preferred_novel.columns]
    if len(preferred_novel):
        print(preferred_novel[show].head(12).to_string(index=False))
    else:
        print("  (none — inspect clean_candidate / MW demotions; try --tau 85)")
    print("\nScaffold-diverse follow-up shortlist:")
    dshow = [c for c in (
        "diverse_rank", "name", "mw", "s_u_percentile", "s_n_dock_percentile",
        "dual_structure_balance", "qed",
    ) if c in diverse.columns]
    if len(diverse):
        print(diverse[dshow].to_string(index=False))
    print("\nKnown reference drugs in gate (positive controls):")
    ref = gate[gate["is_known_reference"]]
    if len(ref):
        print(ref[[
            "name", "max_phase", "mw", "s_u_percentile", "s_n_dock_percentile", "preferred_candidate",
        ]].to_string(index=False))
    else:
        print("  (none in this gate)")


if __name__ == "__main__":
    main()
