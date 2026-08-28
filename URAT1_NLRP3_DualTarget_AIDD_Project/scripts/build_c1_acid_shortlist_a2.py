#!/usr/bin/env python3
"""Freeze C1 Acid shortlist under Amendment A2 (competition pool + clinical audit).

Requires A2 dual-geometry tables. Multi-seed stability applied when seed43/44 exist.
MD remains closed (md_authorized=false) until explicit authorization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# User-prioritized competition tier (pre-A2 clinical audit; not docking ranks)
PRIMARY_TIER1 = ["PF-04620110"]
PRIMARY_TIER2 = ["ADMILPARANT", "RUNCACIGUAT", "LANIFIBRANOR"]
BACKUP_TIER = ["PSI-697", "PF-03882845"]
STRUCTURAL_CONTROL = ["GSK-3008348"]

SOFT_EXCLUDE_SUBSTR = [
    "FLOXACIN",
    "PIROMIDIC",
    "TERBOGREL",
    "KIO-100",
    "MK-5108",
    "CEFCANEL",
    "CEFAZAFLUR",
    "CEFOVECIN",
    "MARBOFLOXACIN",
]


def load_a2_seed(seed: int, a2_dir: Path) -> pd.DataFrame:
    for name in (f"acid_dual_keep_a2_seed{seed}.csv", f"acid_dual_keep_seed{seed}.csv"):
        p = a2_dir / name
        if p.exists():
            return pd.read_csv(p)
    return pd.DataFrame()


def stability_table(a2_dir: Path, seeds: list[int]) -> pd.DataFrame:
    frames = []
    for s in seeds:
        df = load_a2_seed(s, a2_dir)
        if df.empty:
            continue
        df = df[["ligand_id", "keep_dual_acid_geometry", "acid_arg477_min_A"]].copy()
        df = df.rename(
            columns={
                "keep_dual_acid_geometry": f"keep_seed{s}",
                "acid_arg477_min_A": f"arg_seed{s}",
            }
        )
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = frames[0]
    for df in frames[1:]:
        out = out.merge(df, on="ligand_id", how="outer")
    keep_cols = [c for c in out.columns if c.startswith("keep_seed")]
    out["n_seed_pass"] = out[keep_cols].fillna(False).sum(axis=1).astype(int)
    out["stable_ge_2of3"] = out["n_seed_pass"] >= 2
    return out


def name_match(name: str, query: str) -> bool:
    u = str(name).upper()
    q = query.upper().replace("-", "")
    return q in u.replace("-", "") or u == query.upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--a2-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2",
    )
    ap.add_argument("--seeds", nargs="*", type=int, default=[42, 43, 44])
    ap.add_argument("--min-seeds", type=int, default=2)
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/08_nomination",
    )
    args = ap.parse_args()

    pool = pd.read_csv(PROJECT_ROOT / "data/repurposing/screening/docking_pool_p05.csv")
    filt = pd.read_csv(PROJECT_ROOT / "data/repurposing/p2/filters_pool.csv")
    admet = pd.read_csv(PROJECT_ROOT / "data/repurposing/p2/admet_pool.csv")

    stab = stability_table(args.a2_dir, args.seeds)
    if stab.empty:
        raise SystemExit("No A2 seed tables found")

    s42 = load_a2_seed(42, args.a2_dir)
    base = s42.merge(pool, left_on="ligand_id", right_on="repurposing_id", how="left")
    base = base.merge(stab, on="ligand_id", how="left")
    base = base.merge(filt, on="name", how="left")
    base = base.merge(admet, on="name", how="left", suffixes=("", "_admet"))
    # NLRP3 structural metrics (Amendment A2b)
    n_path = args.a2_dir / "nlrp3_structural_metrics_seed42.csv"
    if n_path.exists():
        nmet = pd.read_csv(n_path)[
            [
                "ligand_id",
                "keep_nlrp3_pose",
                "keep_nlrp3_structural",
                "pocket_overlap_frac",
                "ifp_jaccard_vs_np3146",
                "n_key_contacts",
                "key_recovery_frac",
                "CNNscore",
            ]
        ].rename(columns={"CNNscore": "n_CNNscore_structural"})
        drop_cols = [c for c in ("keep_nlrp3_pose", "keep_nlrp3_structural") if c in base.columns]
        if drop_cols:
            base = base.drop(columns=drop_cols)
        base = base.merge(nmet, on="ligand_id", how="left")
    else:
        base["keep_nlrp3_structural"] = False
        base["pocket_overlap_frac"] = None
        base["ifp_jaccard_vs_np3146"] = None
        base["n_key_contacts"] = None

    base["alert_clean"] = (~base.pains_any.fillna(False)) & (~base.brenk.fillna(False))
    base["arg_classic_4A"] = base.acid_arg477_min_A <= 4.0
    base["soft_excluded"] = base.name.astype(str).str.upper().apply(
        lambda n: any(x in n for x in SOFT_EXCLUDE_SUBSTR)
    )

    # Dual keep: URAT1 A2 geometry + NLRP3 loose pose (sensitivity) vs structural (nomination)
    if "keep_dual_acid_geometry" not in base.columns:
        base["keep_dual_acid_geometry"] = base.get("keep_urat1_acid", False) & base.get(
            "keep_nlrp3_pose", False
        )

    eligible = base[
        (base.keep_dual_acid_geometry == True)  # noqa: E712
        & base.alert_clean
        & (~base.soft_excluded)
    ].copy()

    # If multi-seed available, prefer stable; else seed42-only with flag
    available_seeds = [
        s
        for s in args.seeds
        if (args.a2_dir / f"acid_dual_keep_a2_seed{s}.csv").exists()
        or (args.a2_dir / f"acid_dual_keep_seed{s}.csv").exists()
    ]
    if len(available_seeds) >= 2:
        eligible = eligible[eligible.stable_ge_2of3.fillna(False)].copy()
        stability_note = f">= {args.min_seeds}/{len(available_seeds)} seeds dual-pass"
    else:
        stability_note = "seed42 only; seeds 43/44 pending — stability rule not yet applied"

    rows = []
    role_map = {}

    def _flag(series_row, col: str) -> bool:
        if col not in series_row.index:
            return False
        v = series_row[col]
        if pd.isna(v):
            return False
        return bool(v)

    def pick(name_query: str, role: str, tier: str) -> None:
        hit = eligible[eligible.name.apply(lambda n: name_match(n, name_query))]
        if hit.empty:
            return
        r = hit.sort_values(["arg_classic_4A", "acid_arg477_min_A"], ascending=[False, True]).iloc[0]
        if r.ligand_id in role_map:
            return
        role_map[r.ligand_id] = role
        pose_ok = _flag(r, "keep_nlrp3_pose")
        struct_ok = _flag(r, "keep_nlrp3_structural")
        rows.append(
            {
                "shortlist_role": role,
                "competition_tier": tier,
                "ligand_id": r.ligand_id,
                "name": r["name"],
                "chembl_id": r.chembl_id,
                "acid_arg477_min_A": r.acid_arg477_min_A,
                "arg_classic_4A": bool(r.arg_classic_4A),
                "max_phase": r.max_phase,
                "p_active_nlrp3": r.p_active_nlrp3,
                "qed": r.qed,
                "n_seed_pass": int(r["n_seed_pass"]) if "n_seed_pass" in r.index and pd.notna(r["n_seed_pass"]) else 1,
                "stable_ge_2of3": _flag(r, "stable_ge_2of3"),
                "keep_nlrp3_pose": pose_ok,
                "keep_nlrp3_structural": struct_ok,
                "pocket_overlap_frac": r["pocket_overlap_frac"] if "pocket_overlap_frac" in r.index else None,
                "ifp_jaccard_vs_np3146": r["ifp_jaccard_vs_np3146"] if "ifp_jaccard_vs_np3146" in r.index else None,
                "n_key_contacts": r["n_key_contacts"] if "n_key_contacts" in r.index else None,
                "nlrp3_claim": (
                    "NP3-146-compatible pocket pose"
                    if struct_ok
                    else "loose pocket only — not IFP-matched; pathway evidence ≠ direct NACHT binding"
                ),
            }
        )

    for nm in PRIMARY_TIER1:
        pick(nm, "primary", "tier1_mechanistic_anchor")
    if len([x for x in rows if x["shortlist_role"] == "primary"]) < 2:
        for nm in PRIMARY_TIER2:
            pick(nm, "primary", "tier2_translational")
            if len([x for x in rows if x["shortlist_role"] == "primary"]) >= 2:
                break
    for nm in BACKUP_TIER + ["LANIFIBRANOR", "RUNCACIGUAT"]:
        if len([x for x in rows if x["shortlist_role"] == "backup"]) >= 3:
            break
        pick(nm, "backup", "backup_translational")
    for nm in STRUCTURAL_CONTROL:
        hit = base[base.name.apply(lambda n: name_match(n, nm))]
        if not hit.empty:
            r = hit.iloc[0]
            rows.append(
                {
                    "shortlist_role": "structural_control",
                    "competition_tier": "urat1_geometry_positive_non_claim",
                    "ligand_id": r.ligand_id,
                    "name": r["name"],
                    "chembl_id": r.chembl_id,
                    "acid_arg477_min_A": r.acid_arg477_min_A,
                    "arg_classic_4A": bool(r.acid_arg477_min_A <= 4.0) if pd.notna(r.acid_arg477_min_A) else False,
                    "max_phase": r.max_phase,
                    "p_active_nlrp3": r.p_active_nlrp3,
                    "qed": r["qed"] if "qed" in r.index else None,
                    "n_seed_pass": int(r["n_seed_pass"]) if "n_seed_pass" in r.index and pd.notna(r["n_seed_pass"]) else 0,
                    "stable_ge_2of3": _flag(r, "stable_ge_2of3"),
                    "keep_nlrp3_pose": _flag(r, "keep_nlrp3_pose"),
                    "keep_nlrp3_structural": _flag(r, "keep_nlrp3_structural"),
                    "pocket_overlap_frac": r["pocket_overlap_frac"] if "pocket_overlap_frac" in r.index else None,
                    "ifp_jaccard_vs_np3146": r["ifp_jaccard_vs_np3146"] if "ifp_jaccard_vs_np3146" in r.index else None,
                    "n_key_contacts": r["n_key_contacts"] if "n_key_contacts" in r.index else None,
                    "nlrp3_claim": "structural control only",
                }
            )

    shortlist = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    shortlist.to_csv(args.output_dir / "acid_shortlist_a2_competition.csv", index=False)
    eligible.to_csv(args.output_dir / "acid_a2_eligible_audited.csv", index=False)

    summary = {
        "amendment": "A2",
        "md_authorized": False,
        "stability_rule": stability_note,
        "available_seeds": available_seeds,
        "n_a2_dual_seed42": int((base.keep_dual_acid_geometry == True).sum()),  # noqa: E712
        "n_eligible_after_audit": int(len(eligible)),
        "n_shortlist_rows": int(len(shortlist)),
        "primary_names": shortlist[shortlist.shortlist_role == "primary"].name.tolist(),
        "backup_names": shortlist[shortlist.shortlist_role == "backup"].name.tolist(),
        "structural_controls": shortlist[shortlist.shortlist_role == "structural_control"].name.tolist(),
        "note": (
            "Competition shortlist from A2 geometry + clinical audit. "
            "Not activity-ranked. MD closed until authorization."
        ),
    }
    (args.output_dir / "acid_nomination_a2_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
