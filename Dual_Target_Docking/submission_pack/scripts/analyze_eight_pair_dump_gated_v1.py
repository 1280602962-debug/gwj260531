#!/usr/bin/env python3
"""Join all eight Table 2 pairs to a local ChEMBL 37 sqlite with one rule set.

Requires --sqlite. Does not dock. Does not restock Table 2. Does not draw
holdout IDs.

Same operations as analyze_five_pair_dump_gated_v1.py on the production
ligands:

  - STANDARD_OK harvest (IC50, Ki, Kd, EC50, Potency, IC50app, Ki app)
  - dump max versus median pChEMBL at θ = 6.0 on frozen Vina
  - dump-max versus production-panel match (tol 0.015)
  - unused-pool leftover counts at the strict 6.5/5.5 rule
  - dump-level both-end / θ=6.0 / 6.5/5.5 / strict-small-molecule pool counts
  - earliest-document year split (archived; not formal SI)

EGFR/HER2, AChE/BChE, and PIK3CA/mTOR production tables remain the
2026-07-23 REST harvest. This script only asks whether those frozen
pChEMBL values agree with ChEMBL 37. Panel-construction history is not
rewritten.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
UNIV = SCRIPTS.parent
DATA = UNIV.parent
DT = DATA.parent
LOCAL = UNIV / "local_track_b_v0"
TAB = UNIV / "tables"
OUT = LOCAL / "tables" / "eight_pair_dump_gated_v1"
AN = LOCAL / "analysis"
sys.path.insert(0, str(SCRIPTS))
from analyze_five_pair_dump_gated_v1 import (  # noqa: E402
    PRIMARY_YEAR_CUT,
    SEED,
    THETA,
    YEAR_CUTS,
    assign_theta,
    boot_pm_ci,
    directional,
    harvest_activities,
    r4,
    stable_offset,
    write_csv,
)
from pair_ligand_identity_qc_v1 import (  # noqa: E402
    HI,
    LO,
    STANDARD_OK,
    classify,
    connect,
    harvest,
    mol_properties,
    resolve_targets,
)

MATCH_TOL = 0.015
LOCKED_MIN = {
    "EGFR/HER2": 0.430,
    "AChE/BChE": 0.606,
    "PIK3CA/mTOR": 0.692,
    "F2/F10": 0.345,
    "JAK1/TYK2": 0.365,
    "JAK1/JAK2": 0.588,
    "PPARG/PPARA": 0.649,
    "PPARA/PPARD": 0.446,
}
EXPECTED_LEFTOVER = {
    "F2/F10": {"dual": 312, "A_only": 76, "B_only": 245},
    "JAK1/TYK2": {"dual": 1874, "A_only": 59, "B_only": 80},
    "JAK1/JAK2": {"dual": 5953, "A_only": 76, "B_only": 21},
    "PPARG/PPARA": {"dual": 408, "A_only": 50, "B_only": 59},
    "PPARA/PPARD": {"dual": 187, "A_only": 50, "B_only": 68},
}
REST_PAIRS = {"EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"}

PAIRS = [
    {
        "pair": "EGFR/HER2",
        "genes": ("EGFR", "ERBB2"),
        "panel": DT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        "scores": DT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "id_panel": "panel_id",
        "id_score": "ligand",
        "pA": "pchembl_EGFR",
        "pB": "pchembl_HER2",
        "vina_A": "vina_3POZ_hb",
        "vina_B": "vina_3RCD_hb",
        "smiles": "smiles",
        "source": "2026-07-23 ChEMBL REST",
    },
    {
        "pair": "AChE/BChE",
        "genes": ("ACHE", "BCHE"),
        "panel": DT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "scores": DT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "id_panel": "ligand",
        "id_score": "ligand",
        "pA": "pchembl_ACHE",
        "pB": "pchembl_BCHE",
        "vina_A": "vina_ACHE_hb",
        "vina_B": "vina_BCHE_hb",
        "smiles": "smiles",
        "source": "2026-07-23 ChEMBL REST",
    },
    {
        "pair": "PIK3CA/mTOR",
        "genes": ("PIK3CA", "MTOR"),
        "panel": DT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "scores": DT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        "id_panel": "panel_id",
        "id_score": "ligand",
        "pA": "pchembl_PIK3CA",
        "pB": "pchembl_MTOR",
        "vina_A": "vina_4L23_hb",
        "vina_B": "vina_4JT6_hb",
        "smiles": "smiles",
        "source": "2026-07-23 ChEMBL REST",
    },
    {
        "pair": "F2/F10",
        "genes": ("F2", "F10"),
        "panel": TAB / "track_b_panels" / "panel_F2_F10_v1.csv",
        "id_panel": "panel_id",
        "vina_targets": ("4UDW", "2JKH"),
        "pA": "pchembl_A",
        "pB": "pchembl_B",
        "smiles": "canonical_smiles",
        "source": "ChEMBL 37 dump extract",
    },
    {
        "pair": "JAK1/TYK2",
        "genes": ("JAK1", "TYK2"),
        "panel": TAB / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "id_panel": "panel_id",
        "vina_targets": ("6N7A", "3LXP"),
        "pA": "pchembl_A",
        "pB": "pchembl_B",
        "smiles": "canonical_smiles",
        "source": "ChEMBL 37 dump extract",
    },
    {
        "pair": "JAK1/JAK2",
        "genes": ("JAK1", "JAK2"),
        "panel": TAB / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "id_panel": "panel_id",
        "vina_targets": ("6N7A", "8BXH"),
        "pA": "pchembl_A",
        "pB": "pchembl_B",
        "smiles": "canonical_smiles",
        "source": "ChEMBL 37 dump extract",
    },
    {
        "pair": "PPARG/PPARA",
        "genes": ("PPARG", "PPARA"),
        "panel": TAB / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "id_panel": "panel_id",
        "vina_targets": ("9V8H", "6LXA"),
        "pA": "pchembl_A",
        "pB": "pchembl_B",
        "smiles": "canonical_smiles",
        "source": "ChEMBL 37 dump extract",
    },
    {
        "pair": "PPARA/PPARD",
        "genes": ("PPARA", "PPARD"),
        "panel": TAB / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "id_panel": "panel_id",
        "vina_targets": ("6LXA", "5U3Q"),
        "pA": "pchembl_A",
        "pB": "pchembl_B",
        "smiles": "canonical_smiles",
        "source": "ChEMBL 37 dump extract",
    },
]


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def r3_half_up(value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(Decimal(str(float(value))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def load_five_pair_vina():
    path = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
    out = {}
    for row in csv.DictReader(path.open()):
        out.setdefault(row["pair"], {}).setdefault(row["ligand"], {})[row["target"]] = float(row["score_S"])
    return out


def load_score_map(spec, five_vina):
    if "vina_targets" in spec:
        ta, tb = spec["vina_targets"]
        pair_map = five_vina.get(spec["pair"], {})
        out = {}
        for lig, scores in pair_map.items():
            sa, sb = scores.get(ta), scores.get(tb)
            if sa is None or sb is None:
                continue
            out[lig] = (float(sa), float(sb))
        return out
    out = {}
    for row in csv.DictReader(spec["scores"].open()):
        lig = row[spec["id_score"]]
        sa = parse_float(row.get(spec["vina_A"]))
        sb = parse_float(row.get(spec["vina_B"]))
        if sa is None or sb is None:
            continue
        out[lig] = (sa, sb)
    return out


def leftover_buckets(ma, mb, used_chembl, props_by_mol, mol_to_chembl):
    buckets = {"dual": [], "A_only": [], "B_only": []}
    n_strict_sm = {k: 0 for k in buckets}
    leftover = {k: 0 for k in buckets}
    both = set(ma) & set(mb)
    for mol in both:
        pa, pb = ma[mol], mb[mol]
        if pa >= HI and pb >= HI:
            cls = "dual"
        elif pa >= HI and pb <= LO:
            cls = "A_only"
        elif pb >= HI and pa <= LO:
            cls = "B_only"
        else:
            continue
        props = props_by_mol.get(mol) or {}
        good, _ = classify(props)
        if not good:
            continue
        n_strict_sm[cls] += 1
        cid = props.get("chembl_id") or mol_to_chembl.get(mol)
        if cid in used_chembl:
            continue
        leftover[cls] += 1
        buckets[cls].append(cid)
    return n_strict_sm, leftover


def dump_pool_counts(ma, mb, props_by_mol):
    both = set(ma) & set(mb)
    theta6 = {"dual": 0, "A_only": 0, "B_only": 0, "neither": 0}
    strict = {"dual": 0, "A_only": 0, "B_only": 0}
    strict_sm = {"dual": 0, "A_only": 0, "B_only": 0}
    for mol in both:
        pa, pb = ma[mol], mb[mol]
        cls6 = assign_theta(pa, pb)
        if cls6:
            theta6[cls6] += 1
        if pa >= HI and pb >= HI:
            cls = "dual"
        elif pa >= HI and pb <= LO:
            cls = "A_only"
        elif pb >= HI and pa <= LO:
            cls = "B_only"
        else:
            continue
        strict[cls] += 1
        props = props_by_mol.get(mol) or {}
        good, _ = classify(props)
        if good:
            strict_sm[cls] += 1
    return {
        "n_both": len(both),
        "n_theta6_dual": theta6["dual"],
        "n_theta6_A_only": theta6["A_only"],
        "n_theta6_B_only": theta6["B_only"],
        "n_theta6_neither": theta6["neither"],
        "n_strict_dual": strict["dual"],
        "n_strict_A_only": strict["A_only"],
        "n_strict_B_only": strict["B_only"],
        "n_strict_smallmol_dual": strict_sm["dual"],
        "n_strict_smallmol_A_only": strict_sm["A_only"],
        "n_strict_smallmol_B_only": strict_sm["B_only"],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    args = ap.parse_args()
    if not args.sqlite.exists():
        print(f"missing sqlite: {args.sqlite}", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)

    con = connect(args.sqlite)
    genes = {g for spec in PAIRS for g in spec["genes"]}
    meta = resolve_targets(con, genes)
    missing = sorted(genes - set(meta))
    if missing:
        print(f"unresolved genes: {missing}", file=sys.stderr)
        return 2
    for spec in PAIRS:
        ga, gb = spec["genes"]
        print(
            f"resolved {spec['pair']}: {ga}={meta[ga]['chembl']} "
            f"{gb}={meta[gb]['chembl']}",
            flush=True,
        )

    five_vina = load_five_pair_vina()
    maps = harvest(con, {m["tid"] for m in meta.values()})

    maxmed_rows = []
    maxmed_auroc = []
    year_rows = []
    leftover_rows = []
    pool_rows = []
    parity_rows = []

    mol_to_chembl = {}
    for spec in PAIRS:
        pair = spec["pair"]
        ga, gb = spec["genes"]
        tid_a, tid_b = meta[ga]["tid"], meta[gb]["tid"]
        ma, mb = maps[tid_a], maps[tid_b]
        both = set(ma) & set(mb)
        props = mol_properties(con, both)
        for mol, rec in props.items():
            if rec.get("chembl_id"):
                mol_to_chembl[mol] = rec["chembl_id"]
        pool = dump_pool_counts(ma, mb, props)
        pool_rows.append({"pair": pair, "gene_A": ga, "gene_B": gb, **pool})

        panel = list(csv.DictReader(spec["panel"].open()))
        used = {row["molecule_chembl_id"] for row in panel}
        n_strict_sm, leftover = leftover_buckets(ma, mb, used, props, mol_to_chembl)
        exp = EXPECTED_LEFTOVER.get(pair)
        leftover_ok = "" if exp is None else int(all(leftover[c] == exp[c] for c in leftover))
        leftover_rows.append(
            {
                "pair": pair,
                "n_panel_ids": len(used),
                "n_strict_smallmol_dual": n_strict_sm["dual"],
                "n_strict_smallmol_A_only": n_strict_sm["A_only"],
                "n_strict_smallmol_B_only": n_strict_sm["B_only"],
                "leftover_dual": leftover["dual"],
                "leftover_A_only": leftover["A_only"],
                "leftover_B_only": leftover["B_only"],
                "expected_leftover_dual": "" if exp is None else exp["dual"],
                "expected_leftover_A_only": "" if exp is None else exp["A_only"],
                "expected_leftover_B_only": "" if exp is None else exp["B_only"],
                "matches_frozen_five_pair_summary": leftover_ok,
                "note": (
                    "unused-pool leftover after excluding production-panel ChEMBL IDs; "
                    "strict 6.5/5.5 + drug-like small-molecule filter; no holdout drawn"
                ),
            }
        )

        vina_map = load_score_map(spec, five_vina)
        ids = [row["molecule_chembl_id"] for row in panel]
        acts = harvest_activities(con, {tid_a, tid_b}, ids)
        by = defaultdict(lambda: {"A": [], "B": [], "docs": set(), "years": []})
        for act in acts:
            end = "A" if int(act["tid"]) == tid_a else ("B" if int(act["tid"]) == tid_b else None)
            if end is None:
                continue
            by[act["molecule_chembl_id"]][end].append(float(act["pchembl_value"]))
            if act.get("document_chembl_id"):
                by[act["molecule_chembl_id"]]["docs"].add(act["document_chembl_id"])
            if act.get("year") is not None:
                by[act["molecule_chembl_id"]]["years"].append(int(act["year"]))

        recs_max, recs_med = [], []
        pack = []
        n_missing = n_mismatch = n_flip = n_scored = 0
        for row in panel:
            lig = row[spec["id_panel"]]
            cid = row["molecule_chembl_id"]
            sc = vina_map.get(lig)
            if sc is None:
                continue
            sa, sb = sc
            n_scored += 1
            vals = by[cid]
            has_a = bool(vals["A"])
            has_b = bool(vals["B"])
            panel_a = parse_float(row[spec["pA"]])
            panel_b = parse_float(row[spec["pB"]])
            max_a = max(vals["A"]) if has_a else panel_a
            max_b = max(vals["B"]) if has_b else panel_b
            med_a = statistics.median(vals["A"]) if has_a else None
            med_b = statistics.median(vals["B"]) if has_b else None
            match_a = int(has_a and panel_a is not None and abs(max_a - panel_a) < MATCH_TOL)
            match_b = int(has_b and panel_b is not None and abs(max_b - panel_b) < MATCH_TOL)
            if not has_a or not has_b:
                n_missing += 1
            if not (match_a and match_b):
                n_mismatch += 1
            cls_max = assign_theta(max_a, max_b)
            cls_med = assign_theta(med_a, med_b) if med_a is not None and med_b is not None else None
            if cls_med and cls_max != cls_med:
                n_flip += 1
            recs_max.append({"cls": cls_max, "vina_A": sa, "vina_B": sb})
            if cls_med:
                recs_med.append({"cls": cls_med, "vina_A": sa, "vina_B": sb})
            years = [y for y in vals["years"] if y is not None]
            first_year = min(years) if years else None
            pack.append(
                {
                    "cls": cls_max,
                    "vina_A": sa,
                    "vina_B": sb,
                    "first_year": first_year,
                }
            )
            maxmed_rows.append(
                {
                    "pair": pair,
                    "ligand": lig,
                    "molecule_chembl_id": cid,
                    "n_act_A": len(vals["A"]),
                    "n_act_B": len(vals["B"]),
                    "max_A": r4(max_a),
                    "median_A": r4(med_a),
                    "max_B": r4(max_b),
                    "median_B": r4(med_b),
                    "class_max": cls_max or "",
                    "class_median": cls_med or "",
                    "flip_max_to_median": int(cls_max != cls_med) if cls_med else "",
                    "n_documents": len(vals["docs"]),
                    "first_year": first_year if first_year is not None else "",
                    "panel_pA": row[spec["pA"]],
                    "panel_pB": row[spec["pB"]],
                    "dump_missing_A": int(not has_a),
                    "dump_missing_B": int(not has_b),
                    "dump_max_matches_panel_A": match_a,
                    "dump_max_matches_panel_B": match_b,
                    "production_source": spec["source"],
                }
            )

        auroc_by_label = {}
        for label, recs in (("max_pchembl", recs_max), ("median_pchembl", recs_med)):
            da, db, sm, nD, nA, nB = directional(recs)
            lo, hi = boot_pm_ci(recs, SEED + stable_offset(pair, label))
            auroc_by_label[label] = {
                "n_dual": nD,
                "n_A_only": nA,
                "n_B_only": nB,
                "summary_min": sm,
                "ci_lo": lo,
                "ci_hi": hi,
            }
            maxmed_auroc.append(
                {
                    "pair": pair,
                    "aggregation": label,
                    "n_dual": nD,
                    "n_A_only": nA,
                    "n_B_only": nB,
                    "auroc_D_vs_A": r4(da),
                    "auroc_D_vs_B": r4(db),
                    "summary_min": r4(sm),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "production_source": spec["source"],
                    "note": "same frozen Vina; labels from dump max vs median at θ=6.0",
                }
            )

        for cut in YEAR_CUTS:
            for split, pred in (
                ("train_first_year_lt", lambda y, c=cut: y is not None and y < c),
                ("test_first_year_ge", lambda y, c=cut: y is not None and y >= c),
            ):
                sub = [
                    {"cls": r["cls"], "vina_A": r["vina_A"], "vina_B": r["vina_B"]}
                    for r in pack
                    if pred(r.get("first_year"))
                ]
                da, db, sm, nD, nA, nB = directional(sub)
                nN = sum(r["cls"] == "neither" for r in pack if pred(r.get("first_year")))
                n_undated = sum(r.get("first_year") is None for r in pack)
                powered = split.startswith("test") and nD >= 10 and nA >= 10 and nB >= 10
                if powered:
                    lo, hi = boot_pm_ci(sub, SEED + stable_offset(pair, "year", cut, split))
                else:
                    lo = hi = float("nan")
                year_rows.append(
                    {
                        "pair": pair,
                        "cutoff_year": cut,
                        "split": split,
                        "year_definition": "min_document_year_STANDARD_OK_pchembl",
                        "n_dual": nD,
                        "n_A_only": nA,
                        "n_B_only": nB,
                        "n_neither": nN,
                        "n_undated_panel_scored": n_undated,
                        "auroc_reportable": int(powered),
                        "auroc_D_vs_A": r4(da) if powered else "",
                        "auroc_D_vs_B": r4(db) if powered else "",
                        "summary_min": r4(sm) if powered else "",
                        "ci_lo": r4(lo) if powered else "",
                        "ci_hi": r4(hi) if powered else "",
                        "note": (
                            "archived year split; not formal SI; "
                            "test AUROC only if dual/A/B each n≥10"
                        ),
                    }
                )

        dump_min = auroc_by_label["max_pchembl"]["summary_min"]
        locked = LOCKED_MIN[pair]
        parity_rows.append(
            {
                "pair": pair,
                "production_label_source": spec["source"],
                "n_scored_joined": n_scored,
                "n_dump_missing_end": n_missing,
                "n_dump_mismatch": n_mismatch,
                "n_flip_max_to_median": n_flip,
                "dump_max_n_dual": auroc_by_label["max_pchembl"]["n_dual"],
                "dump_max_n_A_only": auroc_by_label["max_pchembl"]["n_A_only"],
                "dump_max_n_B_only": auroc_by_label["max_pchembl"]["n_B_only"],
                "dump_max_summary_min": r4(dump_min),
                "dump_max_summary_min_3dp": r3_half_up(dump_min),
                "locked_table2_summary_min": r3_half_up(locked),
                "dump_max_matches_locked_3dp": int(r3_half_up(dump_min) == r3_half_up(locked)),
                "leftover_strict_unused_D_A_B": (
                    f"{leftover['dual']}/{leftover['A_only']}/{leftover['B_only']}"
                ),
                "matches_frozen_five_pair_leftover": leftover_ok,
                "dump_pool_n_both": pool["n_both"],
                "dump_pool_theta6_D_A_B": (
                    f"{pool['n_theta6_dual']}/{pool['n_theta6_A_only']}/{pool['n_theta6_B_only']}"
                ),
                "dump_pool_strict_D_A_B": (
                    f"{pool['n_strict_dual']}/{pool['n_strict_A_only']}/{pool['n_strict_B_only']}"
                ),
                "dump_pool_strict_smallmol_D_A_B": (
                    f"{pool['n_strict_smallmol_dual']}/"
                    f"{pool['n_strict_smallmol_A_only']}/"
                    f"{pool['n_strict_smallmol_B_only']}"
                ),
            }
        )
        print(
            f"{pair}: scored={n_scored} missing={n_missing} mismatch={n_mismatch} "
            f"flips={n_flip} dump_min={r3_half_up(dump_min)} "
            f"locked={locked} leftover={leftover['dual']}/{leftover['A_only']}/{leftover['B_only']}",
            flush=True,
        )

    write_csv(OUT / "max_vs_median_ligand_v1.csv", maxmed_rows)
    write_csv(OUT / "max_vs_median_auroc_v1.csv", maxmed_auroc)
    write_csv(OUT / "time_split_v1.csv", year_rows)
    write_csv(OUT / "leftover_counts_v1.csv", leftover_rows)
    write_csv(OUT / "dump_pool_counts_v1.csv", pool_rows)
    write_csv(OUT / "parity_v1.csv", parity_rows)

    rest = [r for r in parity_rows if r["pair"] in REST_PAIRS]
    dump5 = [r for r in parity_rows if r["pair"] not in REST_PAIRS]
    n_rest_mismatch = sum(int(r["n_dump_mismatch"]) for r in rest)
    n_rest_missing = sum(int(r["n_dump_missing_end"]) for r in rest)
    n_rest_scored = sum(int(r["n_scored_joined"]) for r in rest)
    n_rest_flip = sum(int(r["n_flip_max_to_median"]) for r in rest)
    n_five_mismatch = sum(int(r["n_dump_mismatch"]) for r in dump5)
    n_five_missing = sum(int(r["n_dump_missing_end"]) for r in dump5)
    leftover5_ok = all(str(r["matches_frozen_five_pair_leftover"]) == "1" for r in dump5)
    all_locked = all(int(r["dump_max_matches_locked_3dp"]) for r in parity_rows)
    n_2018_ok = sum(
        1
        for r in year_rows
        if r["cutoff_year"] == PRIMARY_YEAR_CUT
        and r["split"].startswith("test")
        and int(r["auroc_reportable"])
    )

    checks = {
        "sqlite": str(args.sqlite),
        "standard_ok": list(STANDARD_OK),
        "match_tol": MATCH_TOL,
        "theta": THETA,
        "n_pairs": len(parity_rows),
        "rest_n_scored": n_rest_scored,
        "rest_n_missing": n_rest_missing,
        "rest_n_mismatch": n_rest_mismatch,
        "rest_n_flip_max_to_median": n_rest_flip,
        "five_n_missing": n_five_missing,
        "five_n_mismatch": n_five_mismatch,
        "five_leftover_matches_frozen": leftover5_ok,
        "all_eight_dump_max_matches_locked_3dp": all_locked,
        "n_2018_test_powered": n_2018_ok,
        "same_label_rule": True,
        "same_panel_construction": False,
        "holdout_drawn": 0,
        "replaces_table2": False,
    }
    (OUT / "run_checks_v1.json").write_text(
        json.dumps(checks, indent=2) + "\n",
        encoding="utf-8",
    )

    def dash(rows, key):
        return " / ".join(f"{r['pair']} {r[key]}" for r in rows)

    lines = [
        "# Eight-pair ChEMBL 37 dump-gated check\n\n",
        "Same STANDARD_OK harvest, θ=6.0 max-versus-median, match tolerance ",
        f"{MATCH_TOL}, leftover unused-pool counts, and dump-level pool counts ",
        "on all eight Table 2 pairs. Frozen Vina scores are unchanged. ",
        "This run does **not** restock Table 2 and does **not** draw holdouts.\n\n",
        f"sqlite: `{args.sqlite}` (ChEMBL 37; tarball SHA-256 ",
        "`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).\n\n",
        "## Label-source sameness\n\n",
        "All eight pairs used the same dump join. Production harvest history ",
        "is **not** identical: EGFR/HER2, AChE/BChE, and PIK3CA/mTOR remain ",
        "the 2026-07-23 REST panels; the other five remain the ChEMBL 37 ",
        "extract. Dump agreement does not rewrite that history.\n\n",
        f"- REST-three scored ligands joined: **{n_rest_scored}**\n",
        f"- REST-three dump-missing ends: **{n_rest_missing}**\n",
        f"- REST-three dump-max vs panel mismatches (tol {MATCH_TOL}): **{n_rest_mismatch}**\n",
        f"- REST-three θ=6.0 max→median class flips: **{n_rest_flip}**\n",
        f"- dump-five missing / mismatch: **{n_five_missing}** / **{n_five_mismatch}**\n",
        f"- dump-five leftover vs frozen `track_b_panel_summary_v1.csv`: ",
        f"**{'MATCH' if leftover5_ok else 'MISMATCH'}**\n",
        f"- all eight dump-max `summary_min` vs locked Table 2 (3 d.p. half-up): ",
        f"**{'MATCH' if all_locked else 'MISMATCH'}**\n",
        f"- 2018 test pairs with dual/A/B each n≥10: **{n_2018_ok}** / 8 ",
        "(archived; not formal SI)\n\n",
        "## Pair parity\n\n",
        "| pair | source | scored | missing | mismatch | flips | dump max min | Table 2 | leftover D/A/B |\n",
        "|---|---|---:|---:|---:|---:|---:|---:|---|\n",
    ]
    for r in parity_rows:
        lines.append(
            f"| {r['pair']} | {r['production_label_source']} | "
            f"{r['n_scored_joined']} | {r['n_dump_missing_end']} | "
            f"{r['n_dump_mismatch']} | {r['n_flip_max_to_median']} | "
            f"{r['dump_max_summary_min_3dp']} | {r['locked_table2_summary_min']} | "
            f"{r['leftover_strict_unused_D_A_B']} |\n"
        )
    lines.extend(
        [
            "\nDump-level pool (same harvest, not panel membership):\n\n",
            f"- θ=6.0 D/A/B: {dash(parity_rows, 'dump_pool_theta6_D_A_B')}\n",
            f"- strict 6.5/5.5 D/A/B: {dash(parity_rows, 'dump_pool_strict_D_A_B')}\n\n",
            "The 2026-08-26 API snapshot on the REST three remains a separate ",
            "sensitivity (Table S3). It is not this dump join.\n",
        ]
    )
    note = "".join(lines)
    (AN / "EIGHT_PAIR_DUMP_GATED_V1.md").write_text(note, encoding="utf-8")
    (AN / "THREE_PAIR_DUMP_GATED_V1.md").write_text(
        "# Three-pair ChEMBL 37 dump-gated check\n\n"
        "EGFR/HER2, AChE/BChE, and PIK3CA/mTOR now use the same ChEMBL 37 "
        "dump operations as the other five pairs. Full eight-pair parity is in "
        "`EIGHT_PAIR_DUMP_GATED_V1.md` and "
        "`eight_pair_dump_gated_v1/parity_v1.csv`.\n\n" + note,
        encoding="utf-8",
    )
    print("wrote", OUT)
    print("".join(lines))
    rest_ok = n_rest_missing == 0 and n_rest_mismatch == 0
    five_ok = n_five_missing == 0 and n_five_mismatch == 0 and leftover5_ok
    return 0 if rest_ok and five_ok and all_locked else 1


if __name__ == "__main__":
    raise SystemExit(main())
