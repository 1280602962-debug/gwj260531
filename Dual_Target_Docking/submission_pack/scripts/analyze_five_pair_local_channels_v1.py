#!/usr/bin/env python3
"""K=4-comparable stats on locally uploaded five-pair score channels.

Does not dock or rescore. Reads existing CSVs. Destination identity:
PROJECT_IDENTITY_LOCK_V1.md. Does not restock Table 2.

Primary production numbers remain Vina seed 20260727.
CNN is joined by ligand ID (not pair_hint). Primary CNN readout is
cnn_affinity; cnn_score is a sensitivity column only.
Holdout is unused-pool confirmation, not external validation.
Independent GNINA is JAK1/TYK2 only.

Production Vina bootstrap seeds match analyze_five_pair_stack_v1.py so
the already-frozen table2_comparable_theta6_v1.csv numbers reproduce.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from analyze_five_pair_stack_v1 import (  # noqa: E402
    AN,
    LOCAL,
    PAIRS,
    SEED,
    TAB,
    boot_paired,
    boot_pm_ci,
    boot_single,
    caliper_table,
    directional,
    largest_fragment,
    r4,
    stable_offset,
    write_csv,
)

OUT = LOCAL / "tables" / "five_pair_local_channels_v1"
HO_DIR = LOCAL / "tables" / "five_pair_dump_gated_v1"
HO_SCORES = LOCAL / "allpairs_stack" / "holdout" / "tables" / "holdout_scores_vina_mode1_v1.csv"
FROZEN_T2 = LOCAL / "tables" / "five_pair_stack_v1" / "table2_comparable_theta6_v1.csv"
FROZEN_WP = LOCAL / "tables" / "five_pair_stack_v1" / "wrong_pocket_paired_delta_v1.csv"
SEEDS = (20260727, 20260811, 20260812, 20260813, 20260814)
HOLDOUT_PREFIX = {
    "F2/F10": "HOF2F10",
    "JAK1/TYK2": "HOJ1TYK2",
    "JAK1/JAK2": "HOJ1J2",
    "PPARG/PPARA": "HOPGPA",
    "PPARA/PPARD": "HOPAPD",
}


def load_status_scores(path: Path, pair_field: str, score_field: str, ok=("success",)):
    by = defaultdict(lambda: defaultdict(dict))
    status = []
    if not path.exists():
        return by, status
    for r in csv.DictReader(path.open()):
        pair = r.get(pair_field) or r.get("pair")
        st = r.get("status") or "success"
        rec = {
            "pair": pair,
            "ligand": r["ligand"],
            "target": r["target"],
            "status": st,
            "reason": r.get("reason") or "",
        }
        status.append(rec)
        if st not in ok:
            continue
        try:
            by[pair][r["ligand"]][r["target"]] = float(r[score_field])
        except (TypeError, ValueError):
            continue
    return by, status


def load_cnn_by_ligand(path: Path, field: str):
    """Join by ligand ID. Ignore pair_hint (shared 6N7A / 6LXA are concatenated)."""
    by_lig = defaultdict(dict)
    n = 0
    seen = set()
    dups = 0
    for r in csv.DictReader(path.open()):
        key = (r["ligand"], r["target"])
        if key in seen:
            dups += 1
        seen.add(key)
        try:
            by_lig[r["ligand"]][r["target"]] = float(r[field])
            n += 1
        except (TypeError, ValueError):
            continue
    return by_lig, n, dups


def props_from_smiles(smiles: str, lig: str) -> dict | None:
    mol = largest_fragment(smiles)
    if mol is None:
        return None
    ha = mol.GetNumHeavyAtoms()
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        scaf = ""
    return {
        "heavy": float(ha),
        "mw": float(Descriptors.MolWt(mol)),
        "clogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(Descriptors.TPSA(mol)),
        "charge": float(sum(a.GetFormalCharge() for a in mol.GetAtoms())),
        "rotatable": float(Lipinski.NumRotatableBonds(mol)),
        "scaffold": scaf or f"NONE:{lig}",
        "fp": AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
    }


def recs_from_panel(spec, score_map, cls_field="theta6_class"):
    panel = list(csv.DictReader(spec["panel"].open()))
    recs, missing = [], []
    for row in panel:
        lig = row["panel_id"]
        props = props_from_smiles(row["canonical_smiles"], lig)
        if props is None:
            continue
        m = score_map.get(lig, {})
        sa, sb = m.get(spec["target_a"]), m.get(spec["target_b"])
        rec = {
            "pair": spec["pair"],
            "ligand": lig,
            "cls": row[cls_field],
            "pA": float(row["pchembl_A"]),
            "pB": float(row["pchembl_B"]),
            "vina_A": sa,
            "vina_B": sb,
            **props,
        }
        if sa is not None and sb is not None:
            rec["vina_mean"] = (sa + sb) / 2.0
            rec["vina_worst"] = min(sa, sb)
            recs.append(rec)
        else:
            missing.append(lig)
    return recs, missing


def recs_from_holdout(spec, score_map):
    pref = HOLDOUT_PREFIX[spec["pair"]]
    panel = list(csv.DictReader((HO_DIR / f"holdout_panel_{pref}_v1.csv").open()))
    recs, missing = [], []
    for row in panel:
        lig = row["holdout_id"]
        props = props_from_smiles(row["canonical_smiles"], lig)
        if props is None:
            continue
        m = score_map.get(lig, {})
        sa, sb = m.get(spec["target_a"]), m.get(spec["target_b"])
        rec = {
            "pair": spec["pair"],
            "ligand": lig,
            "cls": row["class"],
            "pA": float(row["pchembl_A"]),
            "pB": float(row["pchembl_B"]),
            "vina_A": sa,
            "vina_B": sb,
            **props,
        }
        if sa is not None and sb is not None:
            rec["vina_mean"] = (sa + sb) / 2.0
            rec["vina_worst"] = min(sa, sb)
            recs.append(rec)
        else:
            missing.append(lig)
    return recs, missing


def _t2_seed(pair: str, channel: str, match_stack: bool) -> int:
    if match_stack:
        return SEED + stable_offset(pair, "theta_6.0")
    return SEED + stable_offset(pair, channel, "theta_6.0")


def _dn_seed(pair: str, channel: str, match_stack: bool) -> int:
    if match_stack:
        return SEED + stable_offset(pair, "D_vs_neither")
    return SEED + stable_offset(pair, channel, "D_vs_neither")


def _wp_seed(pair: str, channel: str, match_stack: bool) -> int:
    if match_stack:
        return SEED + stable_offset(pair, "wrong_pocket")
    return SEED + stable_offset(pair, channel, "wrong_pocket")


def summarize_channel(channel: str, packs: dict, label_rule: str, extra_note: str,
                      match_stack: bool = False, full_si: bool = True) -> tuple[list, list, list]:
    t2, wp, cal = [], [], []
    nonempty = {}
    for spec in PAIRS:
        pair = spec["pair"]
        recs = packs.get(pair) or []
        if not recs:
            t2.append(
                {
                    "channel": channel,
                    "pair": pair,
                    "label_rule": label_rule,
                    "bootstrap": "ligand_non_stratified_dualA_B_pool",
                    "n_dual": 0,
                    "n_A_only": 0,
                    "n_B_only": 0,
                    "n_neither": 0,
                    "n_scored_both_ends": 0,
                    "auroc_D_vs_A_pocketB": "",
                    "auroc_D_vs_B_pocketA": "",
                    "summary_min": "",
                    "ci_lo": "",
                    "ci_hi": "",
                    "n_boot_ok": 0,
                    "D_vs_neither_mean": "",
                    "D_vs_neither_ci_lo": "",
                    "D_vs_neither_ci_hi": "",
                    "note": extra_note + "; no both-end scores",
                }
            )
            continue
        nonempty[pair] = recs
        da, db, sm, nD, nA, nB = directional(recs)
        lo, hi, n_ok = boot_pm_ci(recs, seed=_t2_seed(pair, channel, match_stack))
        D = [r for r in recs if r["cls"] == "dual"]
        N = [r for r in recs if r["cls"] == "neither"]
        if D and N:
            dn, dn_lo, dn_hi = boot_single(
                [r["vina_mean"] for r in D],
                [r["vina_mean"] for r in N],
                seed=_dn_seed(pair, channel, match_stack),
            )
        else:
            dn = dn_lo = dn_hi = float("nan")
        t2.append(
            {
                "channel": channel,
                "pair": pair,
                "label_rule": label_rule,
                "bootstrap": "ligand_non_stratified_dualA_B_pool",
                "n_dual": nD,
                "n_A_only": nA,
                "n_B_only": nB,
                "n_neither": len(N),
                "n_scored_both_ends": len(recs),
                "auroc_D_vs_A_pocketB": r4(da),
                "auroc_D_vs_B_pocketA": r4(db),
                "summary_min": r4(sm),
                "ci_lo": r4(lo),
                "ci_hi": r4(hi),
                "n_boot_ok": n_ok,
                "D_vs_neither_mean": r4(dn),
                "D_vs_neither_ci_lo": r4(dn_lo),
                "D_vs_neither_ci_hi": r4(dn_hi),
                "note": extra_note,
            }
        )
        if full_si:
            wp_row = boot_paired(
                recs,
                "vina_B",
                "vina_A",
                "vina_A",
                "vina_B",
                seed=_wp_seed(pair, channel, match_stack),
            )
            wp.append({"channel": channel, "pair": pair, **wp_row})
    if full_si and nonempty:
        cal_rows = caliper_table(nonempty)
        for row in cal_rows:
            row["channel"] = channel
        cal.extend(cal_rows)
    return t2, wp, cal


def cell(rec) -> str:
    if not rec or rec.get("summary_min") in (None, ""):
        return "—"
    return f"{rec['summary_min']} [{rec.get('ci_lo', '')}, {rec.get('ci_hi', '')}]"


def _side(x) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "na"
    if v > 0.5:
        return "above"
    if v < 0.5:
        return "below"
    return "at"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)

    completeness = []
    all_t2, all_wp, all_cal = [], [], []

    def add(channel, packs, label, note, missing_by_pair=None, status_counts=None,
            match_stack=False, full_si=True):
        t2, wp, cal = summarize_channel(
            channel, packs, label, note, match_stack=match_stack, full_si=full_si
        )
        all_t2.extend(t2)
        all_wp.extend(wp)
        all_cal.extend(cal)
        for spec in PAIRS:
            pair = spec["pair"]
            recs = packs.get(pair) or []
            completeness.append(
                {
                    "channel": channel,
                    "pair": pair,
                    "n_scored_both_ends": len(recs),
                    "n_missing_both_ends": len((missing_by_pair or {}).get(pair, [])),
                    "missing_ligands": ";".join((missing_by_pair or {}).get(pair, [])[:20]),
                    "status_counts": json.dumps(status_counts or {}, sort_keys=True) if spec == PAIRS[0] else "",
                }
            )

    # production Vina (reference; identical to seed 20260727)
    vina, vina_st = load_status_scores(LOCAL / "tables" / "scores_vina_mode1_v1.csv", "pair", "score_S")
    packs, miss = {}, {}
    for spec in PAIRS:
        recs, missing = recs_from_panel(spec, vina.get(spec["pair"], {}))
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "vina_20260727",
        packs,
        "theta_6.0",
        "primary production Vina; canonical Table 2 row for this pair",
        miss,
        dict(Counter(r["status"] for r in vina_st)),
        match_stack=True,
        full_si=True,
    )
    vina_min = {r["pair"]: r["summary_min"] for r in all_t2 if r["channel"] == "vina_20260727"}

    # five-seed: per-seed summary_min only (S54-style). Full SI is reserved
    # for the per-ligand median channel below.
    seed_mins = defaultdict(list)
    seed_maps = {}
    for seed in SEEDS:
        path = LOCAL / "tables" / "multiseed" / f"scores_vina_mode1_seed{seed}.csv"
        smap, st = load_status_scores(path, "pair", "score_S")
        seed_maps[seed] = smap
        packs, miss = {}, {}
        for spec in PAIRS:
            recs, missing = recs_from_panel(spec, smap.get(spec["pair"], {}))
            packs[spec["pair"]] = recs
            miss[spec["pair"]] = missing
        ch = f"vina_seed{seed}"
        add(
            ch,
            packs,
            "theta_6.0",
            "five-seed sensitivity; primary remains 20260727; do not average into Table 2",
            miss,
            dict(Counter(r["status"] for r in st)),
            match_stack=False,
            full_si=False,
        )
        for spec in PAIRS:
            recs = packs[spec["pair"]]
            if recs:
                _da, _db, sm, _nD, _nA, _nB = directional(recs)
                seed_mins[spec["pair"]].append(sm)

    # Per-ligand median S across the five seeds, then Table-2 stats.
    med_packs, med_miss = {}, {}
    for spec in PAIRS:
        pair = spec["pair"]
        ta, tb = spec["target_a"], spec["target_b"]
        ligs = set()
        for smap in seed_maps.values():
            ligs.update(smap.get(pair, {}))
        med_map = {}
        for lig in ligs:
            va, vb = [], []
            for smap in seed_maps.values():
                rec = smap.get(pair, {}).get(lig, {})
                if ta in rec:
                    va.append(rec[ta])
                if tb in rec:
                    vb.append(rec[tb])
            out = {}
            if len(va) == 5:
                out[ta] = float(np.median(va))
            if len(vb) == 5:
                out[tb] = float(np.median(vb))
            if out:
                med_map[lig] = out
        recs, missing = recs_from_panel(spec, med_map)
        med_packs[pair] = recs
        med_miss[pair] = missing
    add(
        "vina_fiveseed_per_ligand_median",
        med_packs,
        "theta_6.0",
        "per-ligand median of five Vina seeds; not a Table 2 replacement; primary remains 20260727",
        med_miss,
        {"n_seeds": 5, "seeds": ",".join(map(str, SEEDS))},
        match_stack=False,
        full_si=True,
    )

    seed_agg = []
    med_by_pair = {r["pair"]: r for r in all_t2 if r["channel"] == "vina_fiveseed_per_ligand_median"}
    for spec in PAIRS:
        vals = [v for v in seed_mins[spec["pair"]] if v == v]
        med_row = med_by_pair.get(spec["pair"], {})
        seed_agg.append(
            {
                "pair": spec["pair"],
                "n_seeds": len(vals),
                "primary_20260727": r4(vals[0]) if vals else "",
                "median_of_per_seed_summary_min": r4(float(sorted(vals)[len(vals) // 2])) if vals else "",
                "min_five_seed": r4(min(vals)) if vals else "",
                "max_five_seed": r4(max(vals)) if vals else "",
                "range_of_per_seed_summary_min": r4(max(vals) - min(vals)) if vals else "",
                "per_ligand_median_summary_min": med_row.get("summary_min", ""),
                "per_ligand_median_ci_lo": med_row.get("ci_lo", ""),
                "per_ligand_median_ci_hi": med_row.get("ci_hi", ""),
                "crosses_0.5": int(bool(vals) and min(vals) < 0.5 < max(vals)),
                "note": "median_of_per_seed_summary_min is S54-style; per_ligand_median is a single score channel",
            }
        )

    # RTM
    rtm, rtm_st = load_status_scores(LOCAL / "tables" / "scores_rtm_best9_v1.csv", "pair", "rtm_best")
    packs, miss = {}, {}
    for spec in PAIRS:
        recs, missing = recs_from_panel(spec, rtm.get(spec["pair"], {}))
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "rtm_best9",
        packs,
        "theta_6.0",
        "RTMScore model1 best-of-9 on production Vina poses; higher is better; not an engine bake-off",
        miss,
        dict(Counter(r["status"] for r in rtm_st)),
        match_stack=False,
        full_si=True,
    )

    # GNINA CNN — ligand ID join
    cnn_path = LOCAL / "tables" / "scores_gnina_cnn_best9_v1.csv"
    aff_map, n_aff, aff_dups = load_cnn_by_ligand(cnn_path, "cnn_affinity")
    sc_map, n_sc, sc_dups = load_cnn_by_ligand(cnn_path, "cnn_score")
    packs, miss = {}, {}
    for spec in PAIRS:
        recs, missing = recs_from_panel(spec, aff_map)
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "gnina_cnn_affinity",
        packs,
        "theta_6.0",
        "PRIMARY CNN readout (pK-like). Joined by ligand ID, not pair_hint.",
        miss,
        {"n_ligand_target_rows": n_aff, "duplicate_ligand_target": aff_dups},
        match_stack=False,
        full_si=True,
    )
    packs, miss = {}, {}
    for spec in PAIRS:
        recs, missing = recs_from_panel(spec, sc_map)
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "gnina_cnn_score",
        packs,
        "theta_6.0",
        "SENSITIVITY only. Do not replace cnn_affinity because one pair looks better.",
        miss,
        {"n_ligand_target_rows": n_sc, "duplicate_ligand_target": sc_dups},
        match_stack=False,
        full_si=True,
    )

    # independent GNINA JAK1/TYK2
    gni, gni_st = load_status_scores(
        LOCAL / "tables" / "scores_gnina_independent_jak1_tyk2_v1.csv",
        "pair",
        "score_S",
        ok=("success",),
    )
    packs, miss = {}, {}
    for spec in PAIRS:
        if spec["pair"] != "JAK1/TYK2":
            packs[spec["pair"]] = []
            miss[spec["pair"]] = []
            continue
        recs, missing = recs_from_panel(spec, gni.get(spec["pair"], {}))
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "gnina_independent_jak1_tyk2",
        packs,
        "theta_6.0",
        "independent GNINA search, JAK1/TYK2 only; formulation-gap pair; not a bake-off",
        miss,
        dict(Counter(r["status"] for r in gni_st)),
        match_stack=False,
        full_si=True,
    )

    # holdout
    ho, ho_st = load_status_scores(HO_SCORES, "pair", "score_S")
    packs, miss = {}, {}
    for spec in PAIRS:
        recs, missing = recs_from_holdout(spec, ho.get(spec["pair"], {}))
        packs[spec["pair"]] = recs
        miss[spec["pair"]] = missing
    add(
        "holdout_vina_20260727",
        packs,
        "strict_6.5_5.5",
        "unused-pool holdout; NOT external validation; JAK1/JAK2 drawn 20/20/18; no neither class",
        miss,
        dict(Counter(r["status"] for r in ho_st)),
        match_stack=False,
        full_si=True,
    )

    write_csv(OUT / "table2_comparable_by_channel_v1.csv", all_t2)
    write_csv(OUT / "wrong_pocket_by_channel_v1.csv", all_wp)
    write_csv(OUT / "property_caliper_by_channel_v1.csv", all_cal)
    write_csv(OUT / "fiveseed_summary_min_aggregate_v1.csv", seed_agg)
    write_csv(OUT / "completeness_v1.csv", completeness)
    write_csv(OUT / "independent_gnina_job_status_v1.csv", gni_st)
    write_csv(OUT / "holdout_job_status_v1.csv", ho_st)

    # Reproduce frozen production numbers
    frozen_t2 = {r["pair"]: r for r in csv.DictReader(FROZEN_T2.open())}
    frozen_wp = {r["pair"]: r for r in csv.DictReader(FROZEN_WP.open())} if FROZEN_WP.exists() else {}
    checks = []
    ok_all = True
    for spec in PAIRS:
        pair = spec["pair"]
        got = next(r for r in all_t2 if r["channel"] == "vina_20260727" and r["pair"] == pair)
        exp = frozen_t2[pair]
        for col in ("summary_min", "ci_lo", "ci_hi", "auroc_D_vs_A_pocketB",
                    "auroc_D_vs_B_pocketA", "D_vs_neither_mean"):
            exp_key = "D_vs_neither_vina_mean" if col == "D_vs_neither_mean" else col
            gv, ev = got.get(col, ""), exp.get(exp_key, "")
            match = (gv == "") == (ev == "")
            if gv != "" and ev != "":
                match = abs(float(gv) - float(ev)) < 1.5e-4
            checks.append(
                {
                    "check": f"t2_{col}",
                    "pair": pair,
                    "got": gv,
                    "expected": ev,
                    "match": int(match),
                }
            )
            ok_all = ok_all and match
        gwp = next((r for r in all_wp if r["channel"] == "vina_20260727" and r["pair"] == pair), {})
        ewp = frozen_wp.get(pair, {})
        if ewp:
            for col in ("delta_matched_minus_wrong", "matched_summary_min", "wrong_summary_min"):
                gv, ev = gwp.get(col, ""), ewp.get(col, "")
                match = gv == ev or (
                    gv != "" and ev != "" and abs(float(gv) - float(ev)) < 1.5e-4
                )
                checks.append(
                    {
                        "check": f"wp_{col}",
                        "pair": pair,
                        "got": gv,
                        "expected": ev,
                        "match": int(match),
                    }
                )
                ok_all = ok_all and match
    checks.append({"check": "all_production_reproduced", "pair": "ALL", "got": int(ok_all),
                   "expected": 1, "match": int(ok_all)})
    write_csv(OUT / "sanity_reproduce_stack_vina_v1.csv", checks)
    if not ok_all:
        raise SystemExit("production Vina did not reproduce frozen five_pair_stack_v1 table")

    by = {(r["channel"], r["pair"]): r for r in all_t2}
    wp_by = {(r["channel"], r["pair"]): r for r in all_wp}
    cal_by = {
        (r["channel"], r["pair"], r["contrast"], float(r["caliper_sd"])): r
        for r in all_cal
    }
    merged = []
    deltas = []
    for spec in PAIRS:
        pair = spec["pair"]
        row = {"pair": pair, "primary_vina_summary_min": vina_min.get(pair, "")}
        for ch, short in (
            ("vina_20260727", "vina"),
            ("vina_fiveseed_per_ligand_median", "fiveseed_med"),
            ("rtm_best9", "rtm"),
            ("gnina_cnn_affinity", "cnn_aff"),
            ("gnina_cnn_score", "cnn_score"),
            ("gnina_independent_jak1_tyk2", "gnina_indep"),
            ("holdout_vina_20260727", "holdout"),
        ):
            r = by.get((ch, pair), {})
            row[f"{short}_min"] = r.get("summary_min", "")
            row[f"{short}_ci"] = (
                f"[{r.get('ci_lo', '')}, {r.get('ci_hi', '')}]" if r.get("summary_min") not in (None, "") else ""
            )
            row[f"{short}_Dn"] = r.get("D_vs_neither_mean", "")
            row[f"{short}_nDAB"] = (
                f"{r.get('n_dual', '')}/{r.get('n_A_only', '')}/{r.get('n_B_only', '')}"
                if r.get("n_scored_both_ends") else ""
            )
            w = wp_by.get((ch, pair), {})
            row[f"{short}_wp_delta"] = w.get("delta_matched_minus_wrong", "")
            ca = cal_by.get((ch, pair, "D_vs_A_pocketB", 1.0), {})
            cb = cal_by.get((ch, pair, "D_vs_B_pocketA", 1.0), {})
            row[f"{short}_cal1_DA"] = ca.get("auroc_matched", "")
            row[f"{short}_cal1_DB"] = cb.get("auroc_matched", "")
            if ch != "vina_20260727" and r.get("summary_min") not in (None, "") and vina_min.get(pair) not in (None, ""):
                deltas.append(
                    {
                        "pair": pair,
                        "channel": ch,
                        "primary_vina_min": vina_min[pair],
                        "channel_min": r["summary_min"],
                        "delta_channel_minus_vina": r4(float(r["summary_min"]) - float(vina_min[pair])),
                        "channel_side_of_0.5": _side(r["summary_min"]),
                        "vina_side_of_0.5": _side(vina_min[pair]),
                        "flips_side_of_0.5": int(_side(r["summary_min"]) != _side(vina_min[pair])),
                        "channel_ci_excludes_0.5": int(
                            r.get("ci_lo") not in (None, "")
                            and r.get("ci_hi") not in (None, "")
                            and (float(r["ci_hi"]) < 0.5 or float(r["ci_lo"]) > 0.5)
                        ),
                    }
                )
        sa = next((x for x in seed_agg if x["pair"] == pair), {})
        row["fiveseed_median_of_aurocs"] = sa.get("median_of_per_seed_summary_min", "")
        row["fiveseed_range"] = sa.get("range_of_per_seed_summary_min", "")
        row["fiveseed_crosses_0.5"] = sa.get("crosses_0.5", "")
        row["note"] = (
            "cnn_affinity is the CNN primary. Holdout is unused-pool, not external. "
            "Does not replace Table 2."
        )
        merged.append(row)
    write_csv(OUT / "merged_channel_summary_v1.csv", merged)
    write_csv(OUT / "channel_minus_primary_vina_v1.csv", deltas)

    # Data-driven bullets
    bullets = []
    cross = [s for s in seed_agg if s["crosses_0.5"]]
    if not cross:
        bullets.append(
            "- Five-seed per-seed `summary_min` ranges do not cross 0.5 on any pair. "
            "Primary seed 20260727 is kept."
        )
    else:
        bullets.append(
            "- Five-seed ranges cross 0.5 on: "
            + ", ".join(s["pair"] for s in cross)
            + ". Primary seed is still 20260727; do not pick a prettier seed."
        )
    g = by.get(("gnina_independent_jak1_tyk2", "JAK1/TYK2"), {})
    v = by[("vina_20260727", "JAK1/TYK2")]
    if g.get("summary_min") not in (None, ""):
        bullets.append(
            f"- Independent GNINA on JAK1/TYK2: `summary_min` {cell(g)}, "
            f"Dual-vs-neither {g.get('D_vs_neither_mean', '')} "
            f"(Vina {cell(v)}, Dual-vs-neither {v.get('D_vs_neither_mean', '')}). "
            "The formulation gap (weak min, stronger Dual-vs-neither) remains."
        )
    rtm_p = by.get(("rtm_best9", "PPARG/PPARA"), {})
    v_p = by[("vina_20260727", "PPARG/PPARA")]
    if rtm_p.get("summary_min") not in (None, ""):
        bullets.append(
            f"- RTM on PPARG/PPARA: `summary_min` {cell(rtm_p)} vs Vina {cell(v_p)}; "
            f"Dual-vs-neither {rtm_p.get('D_vs_neither_mean', '')} vs "
            f"{v_p.get('D_vs_neither_mean', '')}. "
            "A high Dual-vs-neither with a collapsed min is the same negative-class "
            "failure mode, not a PPAR success."
        )
    lifts = []
    for spec in PAIRS:
        rr = by.get(("rtm_best9", spec["pair"]), {})
        vv = by[("vina_20260727", spec["pair"])]
        if rr.get("summary_min") not in (None, "") and vv.get("summary_min") not in (None, ""):
            dlt = float(rr["summary_min"]) - float(vv["summary_min"])
            lifts.append((spec["pair"], dlt, rr["summary_min"], vv["summary_min"]))
    if lifts:
        up = [x for x in lifts if x[1] > 0.05]
        down = [x for x in lifts if x[1] < -0.05]
        bits = []
        if up:
            bits.append("lifts " + ", ".join(f"{p} ({a} vs {b})" for p, _, a, b in up))
        if down:
            bits.append("drops " + ", ".join(f"{p} ({a} vs {b})" for p, _, a, b in down))
        if bits:
            bullets.append("- RTM vs Vina `summary_min`: " + "; ".join(bits) + ".")
    c_j = by.get(("gnina_cnn_affinity", "JAK1/TYK2"), {})
    c_p = by.get(("gnina_cnn_affinity", "PPARG/PPARA"), {})
    bullets.append(
        f"- CNN affinity (primary CNN readout): JAK1/TYK2 `{c_j.get('summary_min', '')}` "
        f"[{c_j.get('ci_lo', '')}, {c_j.get('ci_hi', '')}]; "
        f"PPARG/PPARA `{c_p.get('summary_min', '')}` "
        f"[{c_p.get('ci_lo', '')}, {c_p.get('ci_hi', '')}]. "
        "`cnn_score` is sensitivity only and is not promoted."
    )
    h_p = by.get(("holdout_vina_20260727", "PPARG/PPARA"), {})
    if h_p.get("summary_min") not in (None, ""):
        bullets.append(
            f"- Holdout PPARG/PPARA `summary_min` {cell(h_p)} vs main-panel {cell(v_p)}. "
            "The main-panel PPARG min does **not** replicate on the unused-pool draw "
            "(membership-sensitive). Holdout is not external validation."
        )
    h_j = by.get(("holdout_vina_20260727", "JAK1/JAK2"), {})
    v_j = by[("vina_20260727", "JAK1/JAK2")]
    if h_j.get("summary_min") not in (None, ""):
        bullets.append(
            f"- Holdout JAK1/JAK2 stays same-direction ({cell(h_j)}; main {cell(v_j)}). "
            "Drawn 20/20/18; do not relax the Murcko cap."
        )
    weak_h = []
    for spec in PAIRS:
        hh = by.get(("holdout_vina_20260727", spec["pair"]), {})
        if hh.get("summary_min") not in (None, "") and float(hh["summary_min"]) < 0.5:
            weak_h.append(f"{spec['pair']} {hh['summary_min']}")
    if weak_h:
        bullets.append("- Holdout stays weak (<0.5) on: " + ", ".join(weak_h) + ".")
    wp_ex = [r for r in all_wp if r.get("ci_excludes_zero") in (True, "True", "true", 1, "1")]
    if not wp_ex:
        bullets.append(
            "- Wrong-pocket paired Δ CIs include 0 on every reported channel/pair "
            "(same qualitative result as production Vina)."
        )
    else:
        bullets.append(
            "- Wrong-pocket Δ CI excludes 0 on: "
            + ", ".join(f"{r['channel']} {r['pair']}" for r in wp_ex)
            + ". Those are **same-pose rescores**, not independent docking. "
            "Production Vina Δ CIs all include 0. Do not read rescore Δ as "
            "pocket-geometry proof. `cnn_score` on F2/F10 is sensitivity only."
        )

    lines = [
        "# Five-pair local channels (same article, not Table 2)\n\n",
        "Scores were already computed locally. This file is **statistics only**: ",
        "Table-2-comparable ligand-level non-stratified `summary_min` CIs ",
        "(B=2000, seed 20260729 + SHA offset), wrong-pocket paired Δ, and ",
        "property caliper (1 SD). CNN joined by **ligand ID**, not `pair_hint`. ",
        "Primary CNN readout = `cnn_affinity`. Holdout is unused-pool, **not** ",
        "external validation.\n\n",
        "Production Vina (`vina_20260727`) **reproduced** ",
        "`tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv` ",
        "(same estimand, same bootstrap seeds).\n\n",
        "## Primary Vina (seed 20260727) vs local channels (`summary_min` [95% CI])\n\n",
        "| pair | Vina | five-seed per-seed median (range) | per-ligand median | RTM | CNN affinity | Holdout | indep. GNINA |\n",
        "|---|---|---|---|---|---|---|---|\n",
    ]
    for spec in PAIRS:
        pair = spec["pair"]
        v = by[("vina_20260727", pair)]
        r = by.get(("rtm_best9", pair), {})
        c = by.get(("gnina_cnn_affinity", pair), {})
        h = by.get(("holdout_vina_20260727", pair), {})
        g = by.get(("gnina_independent_jak1_tyk2", pair), {})
        m = by.get(("vina_fiveseed_per_ligand_median", pair), {})
        sa = next((x for x in seed_agg if x["pair"] == pair), {})
        gcell = cell(g) if pair == "JAK1/TYK2" else "not run"
        lines.append(
            f"| {pair} | {cell(v)} | {sa.get('median_of_per_seed_summary_min', '')} "
            f"({sa.get('min_five_seed', '')}–{sa.get('max_five_seed', '')}) | "
            f"{cell(m)} | {cell(r)} | {cell(c)} | {cell(h)} | {gcell} |\n"
        )
    channel_table = "".join(lines[lines.index(
        "## Primary Vina (seed 20260727) vs local channels (`summary_min` [95% CI])\n\n"
    ):])
    lines += [
        "\n## Completeness that matters\n\n",
        "- Production Vina / RTM / CNN / per-ligand five-seed median use the same ",
        "both-end set as the zero-dock stack (F2 107, JAK1/TYK2 109, JAK1/JAK2 110, ",
        "PPARG 109, PPARA/PPARD 110).\n",
        "- Later Vina seeds recovered a few production misses (`F2F10_105/106`, ",
        "`J1TYK2_092`, `PGPA_030`). Per-seed AUROCs are therefore not on an identical ",
        "ligand set. The per-ligand median requires all five seeds and stays on the ",
        "production intersection. Do not prefer a later seed because it scored more ligands.\n",
        "- Independent GNINA JAK1/TYK2: 210/220 jobs ok; four CG0-fail ligands ",
        "(`J1TYK2_033/034/039/065`) plus `J1TYK2_092` timeout are dropped (not imputed).\n",
        "- Holdout: `HOF2F10_045` both-end timeout (dual; 19/20/20); `HOPGPA_001` ",
        "missing 9V8H (A-only; 20/19/20). JAK1/JAK2 remains 20/20/18.\n",
        "- CNN `pair_hint` concatenates shared receptors; join is by ligand ID. ",
        f"Duplicate ligand+target rows: {aff_dups}.\n\n",
        "## What changes relative to production Vina\n\n",
    ]
    lines += [b + "\n" for b in bullets]
    lines += [
        "\n## Still true\n\n",
        "Internal four-state formulation audit. Not database-external. ",
        "Does not restock Table 2 or change the title. ",
        "Does not hard-dock BindingDB. JAK1/JAK2 holdout remains 20/20/18. ",
        "Do not promote `cnn_score`. Do not average channels into one master AUROC.\n",
    ]
    (AN / "FIVE_PAIR_LOCAL_CHANNELS_V1.md").write_text("".join(lines), encoding="utf-8")

    merged_md = [
        "# Five-pair results merged (zero-dock + dump-gated + local channels)\n\n",
        "Same article, 8-row destination after withdrawing PIK3CA/PIK3CB ",
        "(`PROJECT_IDENTITY_LOCK_V1.md`). Five pairs added after the ChEMBL 37 ",
        "census; not a 2026-07-23 freeze. **Does not restock Table 2 or retitle.**\n\n",
        "This note merges three already-run layers. It does not re-dock.\n\n",
        "## 1. Primary numbers (production Vina, θ = 6.0)\n\n",
        "Source: `FIVE_PAIR_STACK_V1.md`. Estimand = ligand-level non-stratified ",
        "`summary_min` of Dual-vs-A (pocket B) and Dual-vs-B (pocket A). ",
        "Reproduced by this channel run.\n\n",
        "| pair | n D/A/B | summary_min [95% CI] | Dual-vs-neither | best descriptor |\n",
        "|---|---:|---|---:|---|\n",
        "| F2/F10 | 31/32/32 | 0.3448 [0.2109, 0.4773] | 0.5188 | clogp 0.5151 |\n",
        "| JAK1/TYK2 | 31/32/32 | 0.3649 [0.2306, 0.503] | 0.7696 | clogp 0.5796 |\n",
        "| JAK1/JAK2 | 32/32/32 | 0.5884 [0.4444, 0.7246] | 0.7299 | heavy 0.5781 |\n",
        "| PPARG/PPARA | 32/31/32 | 0.6492 [0.5045, 0.7508] | 0.6853 | tpsa 0.6274 |\n",
        "| PPARA/PPARD | 32/32/32 | 0.4463 [0.2958, 0.5841] | 0.5647 | clogp 0.5635 |\n\n",
        "Only PPARG/PPARA has a `summary_min` CI entirely above 0.5. ",
        "JAK1/JAK2 point > 0.5 but CI includes 0.5. The other three are ",
        "anti-directional or straddle 0.5. Dual-vs-neither is **not** the ",
        "primary claim (it can look good while a weak arm fails).\n\n",
        "## 2. Dump-gated and leftover holdout IDs\n\n",
        "Source: `FIVE_PAIR_DUMP_GATED_V1.md`. sqlite SHA matches the freeze.\n\n",
        "- Dump max vs panel pChEMBL: **0/545** mismatches; **0** missing dump ends.\n",
        "- Leftover vs frozen panel summary: exact MATCH. Holdout IDs frozen ",
        "(seed 20260731, Murcko cap 3). JAK1/JAK2 **20/20/18**.\n",
        "- 1 max→median class flip: PPARA/PPARD CHEMBL121. Do not switch labels.\n",
        "- 2018 time split (earliest document year): reportable only on ",
        "JAK1/TYK2 (min 0.3681) and JAK1/JAK2 (min 0.6531). Shared JAK1; ",
        "**not** external validation. 2020 fails n≥10 everywhere.\n\n",
        "## 3. BindingDB / PubChem count-only\n\n",
        "Source: `FIVE_PAIR_CROSSDB_V1.md`. `equal_only` does **not** flip the ",
        "≥50 hard-neg gate. Do not hard-dock BindingDB.\n\n",
        "## 4. Local score channels (this run)\n\n",
        "Source: `FIVE_PAIR_LOCAL_CHANNELS_V1.md` and ",
        "`tables/five_pair_local_channels_v1/`.\n\n",
        channel_table,
        "\n",
    ]
    merged_md += [
        "## 5. What the channels do to the story\n\n",
    ]
    merged_md += [b + "\n" for b in bullets]
    merged_md += [
        "\n## 6. Honest combined reading\n\n",
        "- The five new pairs are **not** a second PPAR-success paper. ",
        "PPARG/PPARA is the only main-panel Vina min whose CI clears 0.5, and ",
        "that number is membership-sensitive (holdout) and formulation-sensitive (RTM).\n",
        "- JAK1/JAK2 is the most internally consistent Vina pair (main + holdout ",
        "same direction; five-seed stable), but the CI still includes 0.5 and ",
        "it is not an independent confirmation of JAK1/TYK2.\n",
        "- F2/F10, JAK1/TYK2, and PPARA/PPARD stay weak-armed under Vina. ",
        "RTM/CNN can move a point estimate without turning the stack into a ",
        "pocket-matched dual classifier.\n",
        "- Independent GNINA on the gap pair does not close the EGFR-like gap.\n",
        "- Unused-pool holdout and 2018 year-split are **internal** confirmations. ",
        "The frozen BindingDB external contract still has 0 pairs. Do not relax it.\n",
        "- Stack completeness is now a **precondition** for a later 8-row Table 2 ",
        "draft. That retitle/restock is a later manuscript step, not this commit.\n\n",
        "## 7. Still forbidden\n\n",
        "Restock Table 2 today; retitle; average eight pairs; hard-dock BindingDB; ",
        "promote `cnn_score`; re-draw holdout IDs; relax Murcko cap; run independent ",
        "GNINA on the other four new pairs; LigPrep / seed 42; CTSK ordinary Vina; ",
        "PIK3CB 2Y3A swap.\n",
    ]
    (AN / "FIVE_PAIR_MERGED_V1.md").write_text("".join(merged_md), encoding="utf-8")

    stack_path = AN / "FIVE_PAIR_STACK_V1.md"
    stack_txt = stack_path.read_text(encoding="utf-8")
    old = (
        "## Still local (user machine)\n\n"
        "- **five_seed_vina** — `local_recompute`: user will submit locally; seeds 20260727 + 20260811–20260814; see LOCAL_RECOMPUTE_PACK_V1.md\n"
        "- **rtm_best_of_9** — `local_recompute`: poses gitignored; regenerate 9 modes then rtmscore_model1\n"
        "- **gnina_cnn_rescore** — `local_recompute`: same poses; --cnn_scoring rescore --minimize\n"
        "- **independent_gnina_search** — `local_rule_subset`: JAK1/TYK2 only (EGFR-like formulation gap)\n"
    )
    new = (
        "## Local channels (scores uploaded; stats in this repo)\n\n"
        "Local A–E scores are in `local_track_b_v0/tables/` "
        "(`ALLPAIRS_STACK_RUN_COMPLETE_V1.md`). "
        "Table-2-comparable statistics are in `FIVE_PAIR_LOCAL_CHANNELS_V1.md` "
        "and the merged reading in `FIVE_PAIR_MERGED_V1.md`. "
        "Primary manuscript numbers remain production Vina 20260727. "
        "Does not restock Table 2.\n"
    )
    if old in stack_txt:
        stack_path.write_text(stack_txt.replace(old, new), encoding="utf-8")
    elif "## Local channels (scores uploaded; stats in this repo)" not in stack_txt:
        raise SystemExit("FIVE_PAIR_STACK_V1.md still-local block not found")

    print("wrote", OUT)
    print("".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
