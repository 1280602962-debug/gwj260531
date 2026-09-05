#!/usr/bin/env python3
"""C5 tier assignment — purely mechanical, no hardcoded candidate names.

Replaces the hardcoded PRIMARY_TIER1/PRIMARY_TIER2/BACKUP_TIER lists in
build_c1_acid_shortlist_a2.py (which is why PF-04620110 ended up "primary"
despite failing the NLRP3 structural gate on all 3 seeds: it was never
selected by a gate at all, it was a hardcoded name).

tier1 = A1(seed42) dual-geometry pass  AND  A2 dual-structural in >=2/3 seeds
         AND  chemistry-audited-eligible
tier2 = (A2 dual-structural >=2/3 seeds AND chemistry-audited-eligible) - tier1

No re-docking. Reads only frozen C1 outputs.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
C1 = PROJECT_ROOT / "data/campaigns/c1"
OUT = PROJECT_ROOT / "data/campaigns/c5/03_tiering"

# Beta-lactam / cephalosporin core: 4-membered ring with N-C(=O), fused or not.
BETA_LACTAM_SMARTS = "[#6]1[#6][#7][#6]1=O"
STRUCTURAL_CONTROL_NOT_CANDIDATE = {"GSK-3008348 FREE BASE"}


def load_a1_dual(seed: int = 42) -> set[str]:
    df = pd.read_csv(C1 / f"07_clinical_dock/acid_dual_a1_frozen/acid_dual_keep_seed{seed}.csv")
    return set(df.loc[df["keep_dual_acid_geometry"] == True, "ligand_id"])  # noqa: E712


def load_a2_structural_ge2of3(seeds: list[int]) -> set[str]:
    cnt: Counter[str] = Counter()
    for s in seeds:
        df = pd.read_csv(C1 / f"07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{s}.csv")
        passing = set(df.loc[df["keep_dual_acid_structural"] == True, "ligand_id"])  # noqa: E712
        for lid in passing:
            cnt[lid] += 1
    return {lid for lid, c in cnt.items() if c >= 2}


def flag_beta_lactam(smiles_series: pd.Series) -> pd.Series:
    try:
        from rdkit import Chem
    except ImportError:
        return pd.Series([None] * len(smiles_series), index=smiles_series.index)
    patt = Chem.MolFromSmarts(BETA_LACTAM_SMARTS)
    out = []
    for smi in smiles_series:
        m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
        out.append(bool(m and patt and m.HasSubstructMatch(patt)))
    return pd.Series(out, index=smiles_series.index)


def main() -> None:
    a1_dual = load_a1_dual(42)
    a2_ge2 = load_a2_structural_ge2of3([42, 43, 44])
    audited = pd.read_csv(C1 / "08_nomination/acid_a2_eligible_audited.csv")
    audited_ids = set(audited["ligand_id"])

    eligible = a2_ge2 & audited_ids
    tier1_ids = a1_dual & eligible
    tier2_ids = eligible - tier1_ids

    cols = [
        "ligand_id", "name", "chembl_id", "max_phase", "p_active_nlrp3", "nlrp3_percentile",
        "qed", "mw", "canonical_smiles", "n_seed_pass", "pains_any", "brenk", "soft_excluded",
    ]
    cols = [c for c in cols if c in audited.columns]

    def build_table(ids: set[str], tier_label: str) -> pd.DataFrame:
        t = audited[audited["ligand_id"].isin(ids)][cols].copy()
        t["tier"] = tier_label
        if "canonical_smiles" in t.columns:
            t["beta_lactam_flag"] = flag_beta_lactam(t["canonical_smiles"])
        t["is_structural_control_not_candidate"] = t["name"].isin(STRUCTURAL_CONTROL_NOT_CANDIDATE)
        return t.sort_values("ligand_id")

    t1 = build_table(tier1_ids, "tier1")
    t2 = build_table(tier2_ids, "tier2")

    OUT.mkdir(parents=True, exist_ok=True)
    t1.to_csv(OUT / "tier1_candidates.csv", index=False)
    t2.to_csv(OUT / "tier2_candidates.csv", index=False)

    summary = {
        "method": "mechanical_gate_intersection_no_hardcoded_names",
        "a1_seed42_dual_pass": len(a1_dual),
        "a2_dual_structural_ge2of3": len(a2_ge2),
        "chemistry_audited_eligible": len(audited_ids),
        "tier1_n": len(t1),
        "tier2_n": len(t2),
        "tier1_beta_lactam_flagged": int(t1.get("beta_lactam_flag", pd.Series(dtype=bool)).sum()) if "beta_lactam_flag" in t1 else None,
        "tier2_beta_lactam_flagged": int(t2.get("beta_lactam_flag", pd.Series(dtype=bool)).sum()) if "beta_lactam_flag" in t2 else None,
        "tier1_structural_control_not_candidate": t1.loc[t1["is_structural_control_not_candidate"], "name"].tolist(),
        "superseded_hand_curated_file": "data/campaigns/c1/08_nomination/acid_shortlist_a2_competition.csv",
        "note": (
            "The superseded file hardcoded PRIMARY_TIER1=[PF-04620110] etc. in "
            "build_c1_acid_shortlist_a2.py; PF-04620110 was never selected by any "
            "gate (it fails the NLRP3 structural gate on 3/3 seeds) and is absent "
            "from tier1/tier2 here."
        ),
    }
    import json
    (OUT / "tier_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
