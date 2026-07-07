#!/usr/bin/env python3
"""
Non-docking computational module F — data-driven candidate nomination.

Answers "besides EGCG, what other candidates are there, and how do we find them?"
without re-docking. It relaxes the thin 6-point Pareto front to a top-k% dual
percentile gate, then applies transparent, reviewer-defensible filters and ranks
the survivors so cleaner dual-node candidates surface.

Selection funnel (all on EXISTING Glide XP data):
  1. Dual gate      : S_U percentile >= tau AND S_N percentile >= tau
  2. Structural-alert clean : no PAINS, no Brenk           (module B)
  3. Drug-likeness  : Lipinski pass AND Veber pass          (module C)
  4. Evidence type  : flag whether the NLRP3 axis is structure-supported
                      (s_n_dock_percentile high) or ML-only (confounding risk)
  5. Rank           : clinical phase -> dual-structure balance -> QED

Known benchmark uricosurics/tools are labelled (positive controls), so genuinely
NEW repurposing candidates can be read off separately.

Inputs (read-only):
  data/repurposing/pareto/pareto_merged_scores.csv
  results/cheminformatics/filters_pool.csv
  results/cheminformatics/admet_pool.csv
  results/cheminformatics/novelty_pool.csv   (optional)

Outputs:
  results/candidates/nominated_candidates.csv
  results/candidates/candidate_nomination_summary.json

Usage:
  python3 scripts/14_candidate_nomination.py --tau 90
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Candidate nomination (non-docking module F)")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--filters", type=Path, default=CHEM_DIR / "filters_pool.csv")
    parser.add_argument("--admet", type=Path, default=CHEM_DIR / "admet_pool.csv")
    parser.add_argument("--novelty", type=Path, default=CHEM_DIR / "novelty_pool.csv")
    parser.add_argument("--tau", type=float, default=90.0, help="Dual percentile gate")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pool)
    f = pd.read_csv(args.filters)
    a = pd.read_csv(args.admet)

    m = df.merge(
        f[["name", "pains_any", "brenk", "nih", "aggregation_risk_heuristic"]],
        on="name", how="left",
    ).merge(
        a[["name", "qed", "mw", "clogp", "tpsa", "hbd", "hba",
           "lipinski_pass", "veber_pass", "oral_absorption_ok"]],
        on="name", how="left",
    )
    if args.novelty.exists():
        nov = pd.read_csv(args.novelty)
        keep = [c for c in ("name", "nn_tanimoto_urat1_active", "nn_tanimoto_nlrp3_active") if c in nov.columns]
        m = m.merge(nov[keep], on="name", how="left")

    # 1. Dual gate
    gate = m[(m["s_u_percentile"] >= args.tau) & (m["s_n_percentile"] >= args.tau)].copy()

    # 2+3. Clean structural + drug-likeness
    for col in ("pains_any", "brenk"):
        gate[col] = gate[col].fillna(False).astype(bool)
    gate["structurally_clean"] = ~gate["pains_any"] & ~gate["brenk"]
    gate["druglike"] = (gate["lipinski_pass"] == True) & (gate["veber_pass"] == True)  # noqa: E712
    gate["clean_candidate"] = gate["structurally_clean"] & gate["druglike"]

    # 4. Evidence type for NLRP3 axis
    gate["nlrp3_structure_supported"] = gate["s_n_dock_percentile"] >= args.tau
    gate["nlrp3_ml_only_risk"] = (gate["s_n_ml_percentile"] >= args.tau) & (gate["s_n_dock_percentile"] < args.tau)

    # Reference vs novel
    gate["is_known_reference"] = gate["name"].astype(str).str.upper().isin(KNOWN_REFERENCE)

    # 5. Ranking: dual-structure balance rewards both axes being structure-high
    gate["dual_structure_balance"] = gate[["s_u_percentile", "s_n_dock_percentile"]].min(axis=1)
    gate = gate.sort_values(
        ["clean_candidate", "is_known_reference", "max_phase", "dual_structure_balance", "qed"],
        ascending=[False, True, False, False, False],
    ).reset_index(drop=True)

    cols = [
        "name", "chembl_id", "max_phase", "s_u_percentile", "s_n_percentile",
        "s_n_ml_percentile", "s_n_dock_percentile", "p_active_nlrp3",
        "pains_any", "brenk", "qed", "lipinski_pass", "veber_pass",
        "structurally_clean", "druglike", "clean_candidate",
        "nlrp3_structure_supported", "nlrp3_ml_only_risk", "is_known_reference",
        "dual_structure_balance",
    ]
    cols = [c for c in cols if c in gate.columns]
    out = gate[cols]
    out.to_csv(args.output_dir / "nominated_candidates.csv", index=False)

    clean_novel = gate[(gate["clean_candidate"]) & (~gate["is_known_reference"])]
    clean_novel_struct = clean_novel[clean_novel["nlrp3_structure_supported"]]

    summary = {
        "module": "F_candidate_nomination",
        "tau": args.tau,
        "n_dual_gate": int(len(gate)),
        "n_clean_candidate": int(gate["clean_candidate"].sum()),
        "n_clean_novel": int(len(clean_novel)),
        "n_clean_novel_structure_supported": int(len(clean_novel_struct)),
        "top_clean_novel": clean_novel[["name", "max_phase", "s_u_percentile",
                                        "s_n_dock_percentile", "p_active_nlrp3", "qed",
                                        "nlrp3_structure_supported"]].head(12).to_dict(orient="records"),
        "known_reference_in_gate": gate[gate["is_known_reference"]][["name", "max_phase",
                                        "s_u_percentile", "s_n_percentile"]].to_dict(orient="records"),
    }
    with open(args.output_dir / "candidate_nomination_summary.json", "w") as f_:
        json.dump(summary, f_, indent=2)

    print(f"=== Candidate nomination (tau={args.tau}) ===")
    print(f"  dual gate: {len(gate)}   clean candidates: {gate['clean_candidate'].sum()}")
    print(f"  clean & novel (excl. known drugs): {len(clean_novel)}")
    print(f"  clean & novel & NLRP3 structure-supported: {len(clean_novel_struct)}")
    print("\nTop clean novel candidates:")
    show = ["name", "max_phase", "s_u_percentile", "s_n_dock_percentile",
            "p_active_nlrp3", "qed", "nlrp3_structure_supported"]
    print(clean_novel[show].head(12).to_string(index=False))
    print("\nKnown reference drugs in gate (positive controls):")
    print(gate[gate["is_known_reference"]][["name", "max_phase", "s_u_percentile", "s_n_percentile"]].to_string(index=False))


if __name__ == "__main__":
    main()
