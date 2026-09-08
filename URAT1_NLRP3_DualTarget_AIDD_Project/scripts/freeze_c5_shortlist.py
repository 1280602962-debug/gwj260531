#!/usr/bin/env python3
"""Freeze the C5 reportable shortlist (0 new docks).

Membership stays the mechanical W3 intersection (A1 ∩ A2-structural ≥2/3 ∩
chemistry audit). W2 IFP is annotated with frozen thresholds; it does not
reshuffle names. Beta-lactams are dropped from reportable tier-2. GSK-3008348
stays a structural control, not a candidate.

Requires existing clinical 9DKB SDFs (acid_dual / acid_dual_a2).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_c5_w2_urat1_ifp_gate import (  # noqa: E402
    ARG_JSON,
    KEY_JSON,
    REC_PDBQT,
    build_crystal_anchors,
    evaluate_sdf,
    load_key_map,
    load_receptor_heavy,
)

TIER1 = PROJECT_ROOT / "data/campaigns/c5/03_tiering/tier1_candidates.csv"
TIER2 = PROJECT_ROOT / "data/campaigns/c5/03_tiering/tier2_candidates.csv"
OUT = PROJECT_ROOT / "data/campaigns/c5/04_shortlist_frozen"
SDF = {
    42: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual/urat1_9dkb/seed42",
    43: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/urat1_9dkb/seed43",
    44: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/urat1_9dkb/seed44",
}


def exclusion_reason(row: pd.Series) -> str:
    reasons = []
    if bool(row["is_structural_control_not_candidate"]):
        reasons.append("structural_control_not_candidate")
    if bool(row["beta_lactam_flag"]):
        reasons.append("beta_lactam")
    return ";".join(reasons)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "per_seed").mkdir(exist_ok=True)

    t1 = pd.read_csv(TIER1)
    t2 = pd.read_csv(TIER2)
    pool = pd.concat([t1, t2], ignore_index=True)

    key_map = load_key_map()
    receptor_heavy = load_receptor_heavy(REC_PDBQT)
    arg_atoms = json.loads(ARG_JSON.read_text())["atoms"]
    anchors = build_crystal_anchors(key_map, receptor_heavy, arg_atoms)
    thresholds = anchors["thresholds"]

    rows = []
    missing_sdfs = []
    for _, rec in pool.iterrows():
        lid = rec["ligand_id"]
        seed_pass = []
        seed_rows = []
        for seed, root in SDF.items():
            sdf = root / f"{lid}_out.sdf"
            exists = sdf.exists() and sdf.stat().st_size > 0
            if not exists:
                missing_sdfs.append({"ligand_id": lid, "seed": seed, "sdf": str(sdf)})
            ev = evaluate_sdf(
                sdf,
                lid,
                key_map,
                receptor_heavy,
                anchors["ref_heavy"],
                anchors["ref_ifp"],
                arg_atoms,
                thresholds,
            )
            ev["seed"] = seed
            ev["sdf_exists"] = exists
            seed_rows.append(ev)
            seed_pass.append(bool(ev.get("keep_urat1_ifp")))
        pd.DataFrame(seed_rows).to_csv(OUT / "per_seed" / f"{lid}_w2.csv", index=False)
        rows.append(
            {
                **{k: rec[k] for k in rec.index},
                "w2_ifp_seed42": seed_pass[0],
                "w2_ifp_seed43": seed_pass[1],
                "w2_ifp_seed44": seed_pass[2],
                "w2_ifp_n_seeds": int(sum(seed_pass)),
                "w2_ifp_ge_2of3": int(sum(seed_pass)) >= 2,
                "acid_arg477_min_A_seed42": seed_rows[0].get("acid_arg477_min_A"),
                "ifp_jaccard_seed42": seed_rows[0].get("ifp_jaccard_vs_crystal_union"),
                "n_key_contacts_seed42": seed_rows[0].get("n_key_contacts"),
            }
        )

    annot = pd.DataFrame(rows)
    annot["excluded_reason"] = annot.apply(exclusion_reason, axis=1)
    annot["reportable"] = annot["excluded_reason"].eq("")

    primary = annot[(annot.tier == "tier1") & annot.reportable].copy()
    backup = annot[(annot.tier == "tier2") & annot.reportable].copy()

    annot.to_csv(OUT / "tier_pool_w2_annotated.csv", index=False)
    primary.to_csv(OUT / "primary_candidates.csv", index=False)
    backup.to_csv(OUT / "backup_candidates.csv", index=False)
    annot[annot.reportable].to_csv(OUT / "reportable_candidates.csv", index=False)

    excluded = annot.loc[~annot.reportable, ["ligand_id", "name", "tier", "excluded_reason"]]
    summary = {
        "method": "W3_mechanical_tiers_plus_W2_IFP_annotation_no_name_shuffle",
        "w2_thresholds_frozen": thresholds,
        "tier1_n": int((annot.tier == "tier1").sum()),
        "tier2_n": int((annot.tier == "tier2").sum()),
        "primary_n": int(len(primary)),
        "backup_n": int(len(backup)),
        "reportable_n": int(annot.reportable.sum()),
        "excluded": excluded.to_dict("records"),
        "primary_names": primary["name"].tolist(),
        "backup_names": backup["name"].tolist(),
        "primary_w2_ifp_ge_2of3_n": int(primary["w2_ifp_ge_2of3"].sum()),
        "backup_w2_ifp_ge_2of3_n": int(backup["w2_ifp_ge_2of3"].sum()),
        "missing_sdfs": missing_sdfs,
        "note": (
            "W2 IFP is an annotation on the frozen W3 intersection. "
            "It does not promote or demote names. GSK-3008348 is a control. "
            "Three cephalosporins are excluded from reportable backup."
        ),
    }
    (OUT / "shortlist_freeze_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    printable = {k: summary[k] for k in summary if k != "w2_thresholds_frozen"}
    print(json.dumps(printable, indent=2))


if __name__ == "__main__":
    main()
