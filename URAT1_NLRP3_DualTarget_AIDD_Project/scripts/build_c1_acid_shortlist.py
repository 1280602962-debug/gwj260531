#!/usr/bin/env python3
"""Build C1 Acid-track nomination tables after dual-geometry keep.

No docking percentiles. Gates:
  - dual keep from acid_dual_keep_seed42.csv (A1 Arg ≤ 7.7027 Å + NLRP3 pose)
  - annotate classic Arg ≤ 4.0 Å tier
  - PAINS/Brenk from filters_pool; chemistry soft already upstream
  - draft shortlist: scaffold-diverse top from classic+clean pool
    target 2 primary + ≤3 backup (campaign §6 Acid)

Outputs under data/campaigns/c1/08_nomination/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]

KNOWN_REFERENCE = {
    "LESINURAD",
    "VERINURAD",
    "DOTINURAD",
    "BENZBROMARONE",
    "PROBENECID",
    "PULIGINURAD",
    "SHR-4640",
    "COLCHICINE",
    "MCC950",
}


def diversify_by_scaffold(df: pd.DataFrame, n: int) -> pd.DataFrame:
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    seen: set[str] = set()
    keep_idx: list = []
    for idx, row in df.iterrows():
        sc = row.get("scaffold")
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
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--keep",
        type=Path,
        default=PROJECT_ROOT
        / "data/campaigns/c1/07_clinical_dock/acid_dual/acid_dual_keep_seed42.csv",
    )
    ap.add_argument(
        "--pool",
        type=Path,
        default=PROJECT_ROOT / "data/repurposing/screening/docking_pool_p05.csv",
    )
    ap.add_argument(
        "--filters",
        type=Path,
        default=PROJECT_ROOT / "data/repurposing/p2/filters_pool.csv",
    )
    ap.add_argument(
        "--admet",
        type=Path,
        default=PROJECT_ROOT / "data/repurposing/p2/admet_pool.csv",
    )
    ap.add_argument(
        "--acid-chem",
        type=Path,
        default=PROJECT_ROOT
        / "data/campaigns/c1/07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/08_nomination",
    )
    ap.add_argument("--arg-a1", type=float, default=7.7027)
    ap.add_argument("--arg-classic", type=float, default=4.0)
    ap.add_argument("--n-primary", type=int, default=2)
    ap.add_argument("--n-backup", type=int, default=3)
    args = ap.parse_args()

    keep = pd.read_csv(args.keep)
    dual = keep[keep["keep_dual_acid_geometry"] == True].copy()  # noqa: E712
    dual = dual.rename(columns={"ligand_id": "repurposing_id"})

    pool = pd.read_csv(args.pool)
    filt = pd.read_csv(args.filters)
    admet = pd.read_csv(args.admet)
    acid = pd.read_csv(args.acid_chem)

    d = dual.merge(
        pool[
            [
                "repurposing_id",
                "name",
                "chembl_id",
                "canonical_smiles",
                "inchi_key",
                "max_phase",
                "atc_code",
                "scaffold",
                "p_active_nlrp3",
                "nlrp3_percentile",
                "selection_reason",
            ]
        ],
        on="repurposing_id",
        how="left",
    )
    d = d.merge(filt, on="name", how="left")
    d = d.merge(admet, on="name", how="left", suffixes=("", "_admet"))
    chem_cols = [
        "repurposing_id",
        "pass_veber",
        "pass_ro5_hb",
        "pass_mw_200_550",
        "pass_chemistry_soft",
        "has_carboxylate",
        "has_tetrazole",
        "has_acyl_sulfonamide",
        "mw",
        "hbd",
        "hba",
        "logp",
        "tpsa",
        "rotbonds",
    ]
    d = d.merge(acid[chem_cols], on="repurposing_id", how="left", suffixes=("", "_acid"))

    d["arg_pass_a1"] = d["acid_arg477_min_A"] <= args.arg_a1
    d["arg_pass_classic_4A"] = d["acid_arg477_min_A"] <= args.arg_classic
    d["pains_any"] = d["pains_any"].fillna(False).astype(bool)
    d["brenk"] = d["brenk"].fillna(False).astype(bool)
    d["alert_clean"] = (~d["pains_any"]) & (~d["brenk"])
    d["is_known_reference"] = d["name"].astype(str).str.upper().isin(KNOWN_REFERENCE)
    d["is_cephalosporin_like"] = d["name"].astype(str).str.upper().str.startswith("CEF")
    # Soft demotion: antibacterial quinolones are weak gout dual-node stories
    name_u = d["name"].astype(str).str.upper()
    d["is_antibacterial_soft"] = (
        name_u.str.contains("FLOXACIN")
        | name_u.str.contains("PIROMIDIC")
        | name_u.str.startswith("CEF")
    )

    # Ranking for draft shortlist: classic Arg + alert-clean first
    eligible = d[
        d["arg_pass_classic_4A"]
        & d["alert_clean"]
        & d["pass_chemistry_soft"].fillna(False)
        & (~d["is_known_reference"])
        & (~d["is_antibacterial_soft"])
    ].copy()
    eligible = eligible.sort_values(
        by=["max_phase", "acid_arg477_min_A", "qed", "u_CNNaffinity", "n_CNNaffinity"],
        ascending=[False, True, False, False, False],
    )

    n_take = args.n_primary + args.n_backup
    draft = diversify_by_scaffold(eligible, n_take)
    roles = (["primary"] * args.n_primary) + (["backup"] * max(0, len(draft) - args.n_primary))
    draft["shortlist_role"] = roles[: len(draft)]
    draft["nomination_status"] = "draft_frozen_for_review"
    draft["md_authorized"] = False

    args.output_dir.mkdir(parents=True, exist_ok=True)
    keep_annot = d.sort_values("acid_arg477_min_A").reset_index(drop=True)
    keep_path = args.output_dir / "acid_dual_geometry_keep_annotated.csv"
    draft_path = args.output_dir / "acid_shortlist_draft.csv"
    eligible_path = args.output_dir / "acid_classic4A_alert_clean_pool.csv"
    keep_annot.to_csv(keep_path, index=False)
    draft.to_csv(draft_path, index=False)
    eligible.to_csv(eligible_path, index=False)

    summary = {
        "seed": 42,
        "arg_threshold_a1_A": args.arg_a1,
        "arg_threshold_classic_A": args.arg_classic,
        "n_dual_geometry_keep": int(len(d)),
        "n_arg_classic_4A": int(d["arg_pass_classic_4A"].sum()),
        "n_alert_clean": int(d["alert_clean"].sum()),
        "n_brenk": int(d["brenk"].sum()),
        "n_pains": int(d["pains_any"].sum()),
        "n_classic4A_alert_clean": int(len(eligible)),
        "n_shortlist_draft": int(len(draft)),
        "n_primary": int((draft["shortlist_role"] == "primary").sum()) if len(draft) else 0,
        "n_backup": int((draft["shortlist_role"] == "backup").sum()) if len(draft) else 0,
        "primary_names": draft.loc[draft["shortlist_role"] == "primary", "name"].tolist()
        if len(draft)
        else [],
        "backup_names": draft.loc[draft["shortlist_role"] == "backup", "name"].tolist()
        if len(draft)
        else [],
        "md_authorized": False,
        "note": (
            "Acid-pose dual-node hypotheses only; no percentile ranking. "
            "Draft shortlist is scaffold-diverse from Arg≤4 Å + alert-clean. "
            "MD (L7) stays closed until explicit authorization."
        ),
        "outputs": {
            "annotated_keep": str(keep_path.relative_to(PROJECT_ROOT)),
            "eligible_pool": str(eligible_path.relative_to(PROJECT_ROOT)),
            "shortlist_draft": str(draft_path.relative_to(PROJECT_ROOT)),
        },
    }
    sum_path = args.output_dir / "acid_nomination_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
