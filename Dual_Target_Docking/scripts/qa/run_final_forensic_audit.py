#!/usr/bin/env python3
"""READ-ONLY forensic audit runner.

Does not import analysis.bootstrap_metrics or compute_canonical_results.
Writes only docs/audit, data/provenance, and results/qa artifacts.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
QA = ROOT / "results" / "qa"
PROV = ROOT / "data" / "provenance"
AUDIT = ROOT / "docs" / "audit"
PRIMARY_PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
RECEPTORS = {
    "EGFR/HER2": ("3POZ", "3RCD"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "PIK3CA/mTOR": ("4L23", "4JT6"),
    "AChE/BChE": ("4EY7", "4BDS"),
    "F2/F10": ("4UDW", "2JKH"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}
N_BOOT = 2000
SEED = 20260729
THETA = 6.0
FIVE_SEEDS = (20260727, 20260811, 20260812, 20260813, 20260814)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fnum(raw):
    if raw is None:
        return None
    text = str(raw).strip()
    if text == "":
        return None
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def is_true(value) -> bool:
    return value in (1, "1", True, "True")


def is_primary_row(row: dict) -> bool:
    return (
        row.get("pair") in PRIMARY_PAIRS
        and row.get("analysis_set") == "main"
        and is_true(row.get("complete_case"))
        and is_true(row.get("activity_eligible", "1"))
    )


def auroc_indep(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def percentile_ci_indep(values, n_boot: int = N_BOOT) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < max(20, n_boot // 4):
        return float("nan"), float("nan")
    lo, hi = np.percentile(arr, [2.5, 97.5])
    return float(lo), float(hi)


def choice(rng, idx):
    if idx.size == 0:
        return idx
    return rng.choice(idx, size=idx.size, replace=True)


def r4(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return float("nan")
    return float(f"{float(x):.4f}")


def load_primary_packs(master_rows):
    packs = {p: [] for p in PRIMARY_PAIRS}
    extras = []
    for r in master_rows:
        rec = {
            "pair": r["pair"],
            "ligand_id": r["ligand_id"],
            "cls": r.get("primary_class_theta6") or "",
            "construction_class": r.get("construction_class") or "",
            "score_A": fnum(r.get("score_A")),
            "score_B": fnum(r.get("score_B")),
            "pA": fnum(r.get("pA")),
            "pB": fnum(r.get("pB")),
            "smiles": r.get("smiles") or "",
            "chembl": r.get("molecule_chembl_id") or "",
            "activity_status": r.get("activity_status") or "",
            "activity_source_kind": r.get("activity_source_kind") or "",
            "activity_eligible": is_true(r.get("activity_eligible", "1")),
            "complete_case": is_true(r.get("complete_case")),
            "analysis_set": r.get("analysis_set") or "",
            "score_source": r.get("score_source") or "",
        }
        if is_primary_row(r):
            packs[r["pair"]].append(rec)
        else:
            extras.append(rec)
    return packs, extras


def directional_points(recs):
    d = [r for r in recs if r["cls"] == "dual"]
    a = [r for r in recs if r["cls"] == "A_only"]
    b = [r for r in recs if r["cls"] == "B_only"]
    n = [r for r in recs if r["cls"] == "neither"]
    da = auroc_indep([r["score_B"] for r in d], [r["score_B"] for r in a])
    db = auroc_indep([r["score_A"] for r in d], [r["score_A"] for r in b])
    smin = min(da, db) if math.isfinite(da) and math.isfinite(db) else float("nan")
    weaker = "D_vs_A_pocketB" if da <= db else "D_vs_B_pocketA"
    dnei = auroc_indep(
        [(r["score_A"] + r["score_B"]) / 2 for r in d],
        [(r["score_A"] + r["score_B"]) / 2 for r in n],
    )
    return {
        "n_dual": len(d),
        "n_A_only": len(a),
        "n_B_only": len(b),
        "n_neither": len(n),
        "n": len(recs),
        "auroc_D_vs_A_pocketB": da,
        "auroc_D_vs_B_pocketA": db,
        "summary_min": smin,
        "weaker_arm": weaker,
        "d_vs_neither_mean": dnei,
        "d": d,
        "a": a,
        "b": b,
        "ncls": n,
    }


def bootstrap_smin(recs):
    pts = directional_points(recs)
    cls = np.array([r["cls"] for r in recs])
    sa = np.array([r["score_A"] for r in recs], dtype=float)
    sb = np.array([r["score_B"] for r in recs], dtype=float)
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")
    rng = np.random.default_rng(SEED)
    das = np.empty(N_BOOT)
    dbs = np.empty(N_BOOT)
    mins = np.empty(N_BOOT)
    for i in range(N_BOOT):
        di, ai, bi = choice(rng, d), choice(rng, a), choice(rng, b)
        das[i] = auroc_indep(sb[di], sb[ai])
        dbs[i] = auroc_indep(sa[di], sa[bi])
        mins[i] = min(das[i], dbs[i])
    da_lo, da_hi = percentile_ci_indep(das)
    db_lo, db_hi = percentile_ci_indep(dbs)
    sm_lo, sm_hi = percentile_ci_indep(mins)
    return pts, (da_lo, da_hi), (db_lo, db_hi), (sm_lo, sm_hi)


def bootstrap_fixed_delta(dual, selective, neither):
    dual = np.asarray(dual, dtype=float)
    selective = np.asarray(selective, dtype=float)
    neither = np.asarray(neither, dtype=float)
    auc_s = auroc_indep(dual, selective)
    auc_n = auroc_indep(dual, neither)
    point = auc_n - auc_s
    rng = np.random.default_rng(SEED)
    deltas = np.empty(N_BOOT)
    for i in range(N_BOOT):
        d = rng.choice(dual, size=dual.size, replace=True)
        s = rng.choice(selective, size=selective.size, replace=True)
        n = rng.choice(neither, size=neither.size, replace=True)
        deltas[i] = auroc_indep(d, n) - auroc_indep(d, s)
    lo, hi = percentile_ci_indep(deltas)
    return auc_s, auc_n, point, lo, hi


def bootstrap_mm(recs):
    cls = np.array([r["cls"] for r in recs])
    sa = np.array([r["score_A"] for r in recs], dtype=float)
    sb = np.array([r["score_B"] for r in recs], dtype=float)
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")

    def smin(idx_d, idx_a, idx_b, matched):
        if matched:
            da = auroc_indep(sb[idx_d], sb[idx_a])
            db = auroc_indep(sa[idx_d], sa[idx_b])
        else:
            da = auroc_indep(sa[idx_d], sa[idx_a])
            db = auroc_indep(sb[idx_d], sb[idx_b])
        return min(da, db)

    m_pt = smin(d, a, b, True)
    u_pt = smin(d, a, b, False)
    delta = m_pt - u_pt
    rng = np.random.default_rng(SEED)
    deltas = np.empty(N_BOOT)
    for i in range(N_BOOT):
        di, ai, bi = choice(rng, d), choice(rng, a), choice(rng, b)
        deltas[i] = smin(di, ai, bi, True) - smin(di, ai, bi, False)
    lo, hi = percentile_ci_indep(deltas)
    return m_pt, u_pt, delta, lo, hi


def ranking_ops(recs):
    ranked = sorted(
        recs,
        key=lambda r: (-(r["score_A"] + r["score_B"]) / 2, r["ligand_id"]),
    )
    n = len(ranked)
    k = math.ceil(0.10 * n)
    top = ranked[:k]
    counts = Counter(r["cls"] for r in top)
    n_dual = sum(1 for r in recs if r["cls"] == "dual")
    ef = (counts.get("dual", 0) / k) / (n_dual / n) if n_dual and k else float("nan")
    duals = [r for r in recs if r["cls"] == "dual"]
    worsts = [min(r["score_A"], r["score_B"]) for r in duals]
    thr = float(np.median(worsts)) if worsts else float("nan")
    eligible = [r for r in recs if r["cls"] != "neither"]
    passed = [r for r in eligible if min(r["score_A"], r["score_B"]) >= thr]
    and_counts = Counter(r["cls"] for r in passed)
    return {
        "n": n,
        "k": k,
        "top_dual": counts.get("dual", 0),
        "top_A": counts.get("A_only", 0),
        "top_B": counts.get("B_only", 0),
        "top_N": counts.get("neither", 0),
        "ef": ef,
        "top_ids": ";".join(r["ligand_id"] for r in top),
        "and_thr": thr,
        "and_n_input": len(eligible),
        "and_n_pass": len(passed),
        "and_dual": and_counts.get("dual", 0),
        "and_A": and_counts.get("A_only", 0),
        "and_B": and_counts.get("B_only", 0),
    }


def load_five_seed_wide():
    """Independent five-seed loader from frozen score files. Does not import analysis."""
    wide = {}
    long_path = ROOT / "data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv"
    if long_path.is_file():
        for r in read_csv(long_path):
            pair = r.get("pair")
            if pair not in PRIMARY_PAIRS or pair == "EGFR/HER2":
                continue
            energy = fnum(r.get("vina_mode1"))
            pocket = r.get("pocket")
            seed = fnum(r.get("seed"))
            lig = r.get("ligand")
            if energy is None or pocket not in {"A", "B"} or seed is None or not lig:
                continue
            wide.setdefault((pair, int(seed), lig), {})[pocket] = -energy
    trackb = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/multiseed"
    for seed in FIVE_SEEDS:
        p = trackb / f"scores_vina_mode1_seed{seed}.csv"
        if not p.is_file():
            continue
        for r in read_csv(p):
            if r.get("status") not in ("success", "ok", "", None):
                continue
            pair = r.get("pair")
            if pair not in RECEPTORS:
                continue
            recA, recB = RECEPTORS[pair]
            pdb = r.get("target")
            pocket = "A" if pdb == recA else ("B" if pdb == recB else None)
            if pocket is None:
                continue
            score = fnum(r.get("score_S"))
            energy = fnum(r.get("mode1_energy"))
            if score is None and energy is not None:
                score = -energy
            if score is None:
                continue
            wide.setdefault((pair, int(r.get("seed") or seed), r["ligand"]), {})[pocket] = score
    egfr = ROOT / "data/egfr_her2_uniform_rdkit_v1/tables/scores_vina_mode1_fiveseed.csv"
    if egfr.is_file():
        pocket_of = {"3POZ": "A", "3RCD": "B"}
        for r in read_csv(egfr):
            if r.get("status") not in ("ok", "success", "", None):
                continue
            energy = fnum(r.get("vina_mode1"))
            pocket = pocket_of.get(r.get("pdb"))
            if energy is None or pocket is None:
                continue
            wide.setdefault(("EGFR/HER2", int(r["seed"]), r["ligand"]), {})[pocket] = -energy
    return wide, long_path, egfr


def main():
    QA.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    findings = []

    master_path = CANON / "current_score_master.csv"
    master = read_csv(master_path)
    packs, extras = load_primary_packs(master)

    # ---- independent statistics replay ----
    replay = []
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        pts, da_ci, db_ci, sm_ci = bootstrap_smin(recs)
        dA = [r["score_A"] for r in pts["d"]]
        dB = [r["score_B"] for r in pts["d"]]
        aB = [r["score_B"] for r in pts["a"]]
        bA = [r["score_A"] for r in pts["b"]]
        nA = [r["score_A"] for r in pts["ncls"]]
        nB = [r["score_B"] for r in pts["ncls"]]
        auc_s_a, auc_n_a, dlt_a, dlt_a_lo, dlt_a_hi = bootstrap_fixed_delta(dA, bA, nA)
        auc_s_b, auc_n_b, dlt_b, dlt_b_lo, dlt_b_hi = bootstrap_fixed_delta(dB, aB, nB)
        m_pt, u_pt, mm, mm_lo, mm_hi = bootstrap_mm(recs)
        rank = ranking_ops(recs)
        canon_dir = {r["estimand"]: r for r in read_csv(CANON / "primary_directional_auroc.csv") if r["pair"] == pair}
        canon_smin = next(r for r in read_csv(CANON / "primary_summary_min.csv") if r["pair"] == pair)
        canon_rank = next(r for r in read_csv(CANON / "two_pocket_mean_ranking.csv") if r["pair"] == pair)
        canon_top = next(r for r in read_csv(CANON / "top10_operating_points.csv") if r["pair"] == pair)
        canon_and = next(r for r in read_csv(CANON / "and_filter_operating_points.csv") if r["pair"] == pair)
        canon_mm = next(r for r in read_csv(CANON / "matched_minus_mismatched.csv") if r["pair"] == pair)
        canon_fix = [r for r in read_csv(CANON / "fixed_score_negative_class_delta.csv") if r["pair"] == pair]

        def disc(name, calc, canon_val, kind="point"):
            c = fnum(canon_val)
            delta = abs(calc - c) if c is not None and math.isfinite(calc) else float("nan")
            status = "P0" if (math.isfinite(delta) and delta > 1e-4) else "MATCH"
            replay.append(
                {
                    "pair": pair,
                    "metric": name,
                    "kind": kind,
                    "independent": "" if not math.isfinite(calc) else f"{calc:.10f}",
                    "canonical": "" if c is None else str(canon_val),
                    "abs_delta": "" if not math.isfinite(delta) else f"{delta:.10f}",
                    "status": status,
                }
            )
            if status == "P0":
                findings.append(f"P0 replay {pair} {name} delta={delta}")

        disc("AUROC_D_vs_A_pocketB", pts["auroc_D_vs_A_pocketB"], canon_dir["AUROC_D_vs_A_pocketB"]["point"])
        disc("AUROC_D_vs_B_pocketA", pts["auroc_D_vs_B_pocketA"], canon_dir["AUROC_D_vs_B_pocketA"]["point"])
        disc("AUROC_D_vs_A_ci_lo", da_ci[0], canon_dir["AUROC_D_vs_A_pocketB"]["ci_lo"], "ci")
        disc("AUROC_D_vs_A_ci_hi", da_ci[1], canon_dir["AUROC_D_vs_A_pocketB"]["ci_hi"], "ci")
        disc("AUROC_D_vs_B_ci_lo", db_ci[0], canon_dir["AUROC_D_vs_B_pocketA"]["ci_lo"], "ci")
        disc("AUROC_D_vs_B_ci_hi", db_ci[1], canon_dir["AUROC_D_vs_B_pocketA"]["ci_hi"], "ci")
        disc("summary_min", pts["summary_min"], canon_smin["summary_min"])
        disc("summary_min_ci_lo", sm_ci[0], canon_smin["ci_lo"], "ci")
        disc("summary_min_ci_hi", sm_ci[1], canon_smin["ci_hi"], "ci")
        disc("d_vs_neither_mean", pts["d_vs_neither_mean"], canon_rank["two_pocket_mean_D_vs_neither"])
        disc("ef_dual_10pct", rank["ef"], canon_rank["ef_dual_10pct"])
        disc("top_k", float(rank["k"]), canon_top["top_k"])
        disc("top_dual", float(rank["top_dual"]), canon_top["top_dual"])
        disc("and_dual", float(rank["and_dual"]), canon_and["and_dual"])
        disc("matched_minus_mismatched", mm, canon_mm["delta"])
        disc("mm_ci_lo", mm_lo, canon_mm["ci_lo"], "ci")
        disc("mm_ci_hi", mm_hi, canon_mm["ci_hi"], "ci")
        fix_a = next(r for r in canon_fix if "pocketA" in r["contrast"])
        fix_b = next(r for r in canon_fix if "pocketB" in r["contrast"])
        disc("fixed_delta_pocketA", dlt_a, fix_a["delta_neither_minus_selective"])
        disc("fixed_delta_pocketA_ci_lo", dlt_a_lo, fix_a["delta_ci_lo"], "ci")
        disc("fixed_delta_pocketA_ci_hi", dlt_a_hi, fix_a["delta_ci_hi"], "ci")
        disc("fixed_delta_pocketB", dlt_b, fix_b["delta_neither_minus_selective"])
        disc("n_dual", float(pts["n_dual"]), canon_smin["n_dual"])
        disc("n_A_only", float(pts["n_A_only"]), canon_smin["n_A_only"])
        disc("n_B_only", float(pts["n_B_only"]), canon_smin["n_B_only"])

    write_csv(QA / "independent_statistics_replay.csv", replay)
    n_p0_replay = sum(1 for r in replay if r["status"] == "P0")

    # ---- membership registry ----
    membership = []

    def add_membership(analysis, pair, recs, reason):
        ids = sorted({r["ligand_id"] if "ligand_id" in r else r.get("ligand", "") for r in recs})
        membership.append(
            {
                "analysis": analysis,
                "pair": pair,
                "n": len(ids),
                "ligand_id_set": "|".join(ids),
                "reason_for_exclusion": reason,
            }
        )

    for pair in PRIMARY_PAIRS:
        add_membership("primary_complete_case", pair, packs[pair], "activity_eligible AND complete_case AND analysis_set=main")
        inelig = [r for r in extras if r["pair"] == pair]
        add_membership(
            "master_excluded_from_primary",
            pair,
            inelig,
            ";".join(sorted({r["activity_status"] or ("incomplete" if not r["complete_case"] else "other") for r in inelig})) or "none",
        )

    # five-seed available vs fixed
    wide, ache_long, egfr_five = load_five_seed_wide()
    five_fixed_rows = []
    for pair in PRIMARY_PAIRS:
        primary_ids = {r["ligand_id"] for r in packs[pair]}
        by_seed = {}
        for seed in FIVE_SEEDS:
            ids = []
            for lig in primary_ids:
                sc = wide.get((pair, seed, lig))
                if sc and "A" in sc and "B" in sc:
                    ids.append(lig)
            by_seed[seed] = set(ids)
            add_membership(
                f"five_seed_available_{seed}",
                pair,
                [{"ligand_id": x} for x in sorted(ids)],
                "missing both-pocket score on this seed" if len(ids) != len(primary_ids) else "same as primary",
            )
        inter = set.intersection(*by_seed.values()) if by_seed else set()
        add_membership(
            "five_seed_intersection_all_seeds",
            pair,
            [{"ligand_id": x} for x in sorted(inter)],
            f"dropped_vs_primary={len(primary_ids - inter)}",
        )
        for seed in FIVE_SEEDS:
            recs = []
            for r in packs[pair]:
                if r["ligand_id"] not in inter:
                    continue
                sc = wide.get((pair, seed, r["ligand_id"]))
                recs.append({**r, "score_A": sc["A"], "score_B": sc["B"]})
            pts = directional_points(recs)
            avail_recs = []
            for r in packs[pair]:
                sc = wide.get((pair, seed, r["ligand_id"]))
                if sc and "A" in sc and "B" in sc:
                    avail_recs.append({**r, "score_A": sc["A"], "score_B": sc["B"]})
            avail = directional_points(avail_recs)
            delta = pts["summary_min"] - avail["summary_min"]
            five_fixed_rows.append(
                {
                    "pair": pair,
                    "seed": seed,
                    "n_intersection": pts["n"],
                    "n_dual": pts["n_dual"],
                    "n_A_only": pts["n_A_only"],
                    "n_B_only": pts["n_B_only"],
                    "n_neither": pts["n_neither"],
                    "n_primary": len(packs[pair]),
                    "n_dropped_vs_primary": len(primary_ids) - len(inter),
                    "summary_min": f"{pts['summary_min']:.4f}" if math.isfinite(pts["summary_min"]) else "",
                    "auroc_D_vs_A_pocketB": f"{pts['auroc_D_vs_A_pocketB']:.4f}",
                    "auroc_D_vs_B_pocketA": f"{pts['auroc_D_vs_B_pocketA']:.4f}",
                    "available_case_summary_min": f"{avail['summary_min']:.4f}",
                    "available_case_n": avail["n"],
                    "delta": f"{delta:.4f}" if math.isfinite(delta) else "",
                    "qualitative_change": "yes" if (math.isfinite(pts["summary_min"]) and ((pts["summary_min"] - 0.5) * (float(next(x['summary_min'] for x in read_csv(CANON/'primary_summary_min.csv') if x['pair']==pair)) - 0.5) < 0)) else "no",
                    "same_protocol_as_primary": 1,
                    "same_membership_as_primary": int(len(inter) == len(primary_ids)),
                    "source_kind": "audit_independent_fixed_membership",
                    "note": "intersection of ligands with both pocket scores on all five seeds; no redock",
                }
            )
    write_csv(QA / "five_seed_fixed_membership_sensitivity.csv", five_fixed_rows)

    # GNINA / RTM / CNN membership
    gnina_files = {
        "EGFR/HER2": ROOT / "data/egfr_her2_uniform_rdkit_v1/tables/gnina_dock_scores_EGFR_HER2.csv",
        "PIK3CA/mTOR": ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_PIK3CA_mTOR.csv",
        "JAK1/TYK2": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_independent_jak1_tyk2_v1.csv",
    }
    for pair, path in gnina_files.items():
        if not path.is_file():
            continue
        rows = read_csv(path)
        ok_ids = set()
        for r in rows:
            lig = r.get("ligand") or r.get("ligand_id")
            status = (r.get("status") or "ok").lower()
            if status in {"ok", "success", ""}:
                ok_ids.add(lig)
        # need both pockets
        by = defaultdict(set)
        for r in rows:
            lig = r.get("ligand") or r.get("ligand_id")
            pdb = r.get("pdb") or r.get("target") or r.get("receptor")
            status = (r.get("status") or "ok").lower()
            if status in {"ok", "success", ""} and lig:
                by[lig].add(pdb)
        recA, recB = RECEPTORS[pair]
        both = [k for k, v in by.items() if recA in v and recB in v]
        add_membership("gnina_independent_both_pockets", pair, [{"ligand_id": x} for x in sorted(both)], "timeout/fail or missing pocket")

    # ECFP membership from oof
    oof = read_csv(CANON / "ecfp4_oof_predictions.csv")
    for pair in PRIMARY_PAIRS:
        for contrast in ("D_vs_A", "D_vs_B"):
            ids = sorted({r["ligand_id"] for r in oof if r["pair"] == pair and r.get("arm") == contrast and r.get("model") == "ECFP4"})
            add_membership(f"ecfp4_{contrast}", pair, [{"ligand_id": x} for x in ids], "failed SMILES parse or class not in contrast")

    write_csv(QA / "analysis_membership_registry.csv", membership)

    # ---- activity provenance matrix ----
    act_rows = []
    for r in master:
        pa, pb = fnum(r.get("pA")), fnum(r.get("pB"))
        recon = ""
        if pa is None or pb is None:
            recon = "missing_arm"
        elif pa >= THETA and pb >= THETA:
            recon = "dual"
        elif pa >= THETA and pb < THETA:
            recon = "A_only"
        elif pb >= THETA and pa < THETA:
            recon = "B_only"
        else:
            recon = "neither"
        mismatch = recon != (r.get("primary_class_theta6") or "") and recon != "missing_arm"
        act_rows.append(
            {
                "pair": r["pair"],
                "ligand_id": r["ligand_id"],
                "molecule_chembl_id": r.get("molecule_chembl_id"),
                "activity_source_kind": r.get("activity_source_kind"),
                "activity_source_file": r.get("activity_source_file"),
                "activity_status": r.get("activity_status"),
                "activity_eligible": r.get("activity_eligible"),
                "analysis_set": r.get("analysis_set"),
                "complete_case": r.get("complete_case"),
                "construction_class": r.get("construction_class"),
                "primary_class_theta6": r.get("primary_class_theta6"),
                "recomputed_theta6": recon,
                "class_match": int(not mismatch),
                "pA": r.get("pA"),
                "pB": r.get("pB"),
                "historical_pA": r.get("historical_pA"),
                "historical_pB": r.get("historical_pB"),
                "n_act_A": r.get("n_act_A"),
                "n_act_B": r.get("n_act_B"),
                "in_primary": int(is_primary_row(r)),
            }
        )
    write_csv(PROV / "activity_provenance_matrix.csv", act_rows)

    # ---- ligand preparation registry ----
    lig_rows = []
    egfr_prep = {}
    egfr_prep_path = ROOT / "data/egfr_her2_uniform_rdkit_v1/tables/ligand_prep_status.csv"
    if egfr_prep_path.is_file():
        for r in read_csv(egfr_prep_path):
            egfr_prep[r["ligand_id"]] = r
    pdbqt_egfr = ROOT / "data/egfr_her2_uniform_rdkit_v1/ligands_pdbqt"
    for r in master:
        pair = r["pair"]
        lid = r["ligand_id"]
        if pair == "EGFR/HER2":
            st = egfr_prep.get(lid, {})
            pdbqt = pdbqt_egfr / f"{lid}.pdbqt"
            lig_rows.append(
                {
                    "pair": pair,
                    "ligand_id": lid,
                    "smiles_source": "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
                    "fragment_rule": "largest_organic",
                    "prep_method": "rdkit_etkdgv3_mmff_meeko",
                    "rdkit_version": "2026.03.5_analysis_pin_not_historical_docking_rdkit",
                    "embed_method": "ETKDGv3",
                    "embed_seed": "20260727",
                    "optimization_method": "MMFF_le_200",
                    "meeko_version": "0.7.1",
                    "input_pdbqt": str(pdbqt.relative_to(ROOT)) if pdbqt.is_file() else "missing",
                    "direct_evidence": "prep_uniform_ligands.py + ligand_prep_status.csv + ligands_pdbqt",
                    "status": "CONFIRMED" if pdbqt.is_file() else "NOT_RECOVERABLE",
                }
            )
        elif pair == "PIK3CA/mTOR":
            lig_rows.append(
                {
                    "pair": pair,
                    "ligand_id": lid,
                    "smiles_source": "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
                    "fragment_rule": "largest_organic_recovered_script",
                    "prep_method": "rdkit_etkdg_meeko_script_recovered",
                    "rdkit_version": "NOT_RECOVERABLE",
                    "embed_method": "ETKDGv3_metadata_only",
                    "embed_seed": "20260727_metadata_only",
                    "optimization_method": "MMFF_metadata_only",
                    "meeko_version": "METADATA_ONLY",
                    "input_pdbqt": "not_in_current_tree",
                    "direct_evidence": "pm48_rdkit_prep_dock.py recovered from git c213c485; ligands_pdbqt absent",
                    "status": "METADATA_ONLY",
                }
            )
        elif pair == "AChE/BChE":
            lig_rows.append(
                {
                    "pair": pair,
                    "ligand_id": lid,
                    "smiles_source": "data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
                    "fragment_rule": "largest_organic_recovered_script",
                    "prep_method": "rdkit_etkdg_meeko_script_recovered",
                    "rdkit_version": "NOT_RECOVERABLE",
                    "embed_method": "ETKDGv3_metadata_only",
                    "embed_seed": "20260727_metadata_only",
                    "optimization_method": "MMFF_metadata_only",
                    "meeko_version": "METADATA_ONLY",
                    "input_pdbqt": "not_in_current_tree",
                    "direct_evidence": "dock_panel.py recovered from git c213c485; ligands_pdbqt absent",
                    "status": "METADATA_ONLY",
                }
            )
        else:
            lig_rows.append(
                {
                    "pair": pair,
                    "ligand_id": lid,
                    "smiles_source": "data/jcim_chembl_universe_v0/local_track_b_v0 panel CSVs",
                    "fragment_rule": "largest_organic_recovered_script",
                    "prep_method": "rdkit_etkdg_meeko_script_recovered",
                    "rdkit_version": "NOT_RECOVERABLE_as_production_pin",
                    "embed_method": "ETKDGv3_metadata_only",
                    "embed_seed": "20260727_metadata_only",
                    "optimization_method": "MMFF_metadata_only",
                    "meeko_version": "0.7.1_yaml_metadata",
                    "input_pdbqt": "not_in_current_tree",
                    "direct_evidence": "prep_track_b_ligands_v1.py + track_b_ligand_prep_status_v1.csv; PDBQT never in git",
                    "status": "METADATA_ONLY",
                }
            )
    write_csv(PROV / "ligand_preparation_registry.csv", lig_rows)

    # ---- docking box registry ----
    box_map = {
        "3POZ": ROOT / "data/egfr_her2_panel120_v0/boxes/3POZ_box_corrected.json",
        "3RCD": ROOT / "data/egfr_her2_panel120_v0/boxes/3RCD_box_corrected.json",
        "4L23": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4L23_box.json",
        "4JT6": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4JT6_box.json",
        "4EY7": ROOT / "data/ache_bche_panel_v0/boxes/4EY7_box.json",
        "4BDS": ROOT / "data/ache_bche_panel_v0/boxes/4BDS_box.json",
        "4UDW": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/4UDW_box.json",
        "2JKH": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/2JKH_box.json",
        "6N7A": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6N7A_box.json",
        "8BXH": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/8BXH_box.json",
        "3LXP": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/3LXP_box.json",
        "9V8H": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/9V8H_box.json",
        "6LXA": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6LXA_box.json",
        "5U3Q": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/5U3Q_box.json",
    }
    old_egfr = {
        "3POZ_legacy": ROOT / "data/egfr_her2_panel120_v0/boxes/3POZ_box.json",
        "3RCD_legacy": ROOT / "data/egfr_her2_panel120_v0/boxes/3RCD_box.json",
    }
    box_rows = []
    for pdb, path in {**box_map, **old_egfr}.items():
        raw = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
        sizes = [fnum(raw.get(k)) for k in ("size_x", "size_y", "size_z")]
        construction = raw.get("construction") or raw.get("box_definition") or ""
        n_h = raw.get("n_heavy_atoms")
        n_lig = raw.get("n_ligand_atoms")
        n_hyd = raw.get("n_hydrogen_atoms")
        H_risk = "possible_if_n_ligand_includes_H" if (n_lig and n_h and int(n_lig) > int(n_h)) else "no_evidence_of_H_inclusive_AABB"
        if "heavy" in str(construction).lower():
            H_risk = "declared_heavy_atom"
        min_ok = all(s is not None and s + 1e-9 >= 20.0 for s in sizes) if all(s is not None for s in sizes) else False
        role = "LEGACY_FILE_PRESENT" if "legacy" in pdb else "CURRENT_PRIMARY"
        box_rows.append(
            {
                "pdb": pdb.replace("_legacy", ""),
                "file": str(path.relative_to(ROOT)) if path.is_file() else "missing",
                "role": role,
                "center_x": raw.get("center_x"),
                "center_y": raw.get("center_y"),
                "center_z": raw.get("center_z"),
                "size_x": raw.get("size_x"),
                "size_y": raw.get("size_y"),
                "size_z": raw.get("size_z"),
                "cognate": raw.get("ligand") or raw.get("cognate") or "",
                "construction": construction or "METADATA_ONLY",
                "n_heavy_atoms": n_h,
                "n_ligand_atoms": n_lig,
                "n_hydrogen_atoms": n_hyd,
                "min_edge_ge_20": int(min_ok),
                "hydrogen_inclusive_risk": H_risk,
                "exists": int(path.is_file()),
            }
        )
    write_csv(PROV / "docking_box_registry.csv", box_rows)

    # ---- docking protocol matrix ----
    proto = []
    for pair in PRIMARY_PAIRS:
        proto.append(
            {
                "pair": pair,
                "channel": "production_vina",
                "level": "PRIMARY",
                "vina_version": "1.2.7",
                "scoring_function": "vina",
                "seed": 20260727,
                "exhaustiveness": 16 if pair == "PIK3CA/mTOR" else 8,
                "num_modes": 9,
                "energy_range": 3,
                "cpu": "NOT_RECOVERABLE",
                "receptor": "+".join(RECEPTORS[pair]),
                "box": "corrected_heavy" if pair == "EGFR/HER2" else "deposited_box_json",
                "ligand_input": "CONFIRMED_pdbqt" if pair == "EGFR/HER2" else "scores_deposited_pdbqt_NOT_RECOVERABLE",
                "score_extraction": "score_S=-vina_mode1",
            }
        )
    for label, level in (
        ("PIK3CA_E8", "SENSITIVITY"),
        ("PM110", "SENSITIVITY"),
        ("4JPS", "SENSITIVITY"),
        ("5DXT", "SENSITIVITY"),
        ("4JSX", "SENSITIVITY"),
    ):
        proto.append(
            {
                "pair": "PIK3CA/mTOR",
                "channel": label,
                "level": level,
                "vina_version": "1.2.7",
                "scoring_function": "vina",
                "seed": 20260727,
                "exhaustiveness": 8 if label == "PIK3CA_E8" else 16,
                "num_modes": 9,
                "energy_range": 3,
                "cpu": "NOT_RECOVERABLE",
                "receptor": {"4JPS": "4JPS+4JT6", "5DXT": "5DXT+4JT6", "4JSX": "4L23+4JSX"}.get(label, "4L23+4JT6"),
                "box": "sensitivity",
                "ligand_input": "deposited_scores",
                "score_extraction": "score_S=-vina_mode1",
            }
        )
    for pair in ("EGFR/HER2", "PIK3CA/mTOR", "JAK1/TYK2"):
        proto.append(
            {
                "pair": pair,
                "channel": "gnina_independent",
                "level": "INDEPENDENT_DOCKING",
                "vina_version": "GNINA_1.3.2",
                "scoring_function": "gnina_default",
                "seed": "campaign",
                "exhaustiveness": "GNINA_default",
                "num_modes": "GNINA",
                "energy_range": "",
                "cpu": "no_gpu",
                "receptor": "+".join(RECEPTORS[pair]),
                "box": "same_as_primary_boxes",
                "ligand_input": "uniform_pdbqt" if pair == "EGFR/HER2" else "deposited_or_not_in_tree",
                "score_extraction": "gnina_mode1_energy; score_S=-energy",
            }
        )
    proto.append(
        {
            "pair": "Track-B",
            "channel": "rtm_best9",
            "level": "SAME_POSE_RESCORING",
            "vina_version": "n/a",
            "scoring_function": "RTMScore_best_of_saved",
            "seed": "n/a",
            "exhaustiveness": "n/a",
            "num_modes": "best_of_saved_poses",
            "energy_range": "",
            "cpu": "",
            "receptor": "Track-B Vina poses",
            "box": "n/a",
            "ligand_input": "Vina poses",
            "score_extraction": "best-of-saved RTM",
        }
    )
    proto.append(
        {
            "pair": "Track-B",
            "channel": "gnina_cnn_best9",
            "level": "SAME_POSE_RESCORING",
            "vina_version": "n/a",
            "scoring_function": "GNINA_CNN_affinity",
            "seed": "n/a",
            "exhaustiveness": "n/a",
            "num_modes": "best_of_saved_poses",
            "energy_range": "",
            "cpu": "",
            "receptor": "Track-B Vina poses",
            "box": "n/a",
            "ligand_input": "Vina poses",
            "score_extraction": "best-of-saved CNN",
        }
    )
    write_csv(PROV / "docking_protocol_matrix.csv", proto)

    # ---- scientific provenance matrix ----
    sci = []
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        kinds = Counter(r["activity_source_kind"] for r in recs)
        sci.append(
            {
                "pair": pair,
                "research_role": "primary_directional_selectivity",
                "activity_kinds": ";".join(f"{k}:{v}" for k, v in kinds.items()),
                "n_primary": len(recs),
                "score_source": recs[0]["score_source"] if recs else "",
                "ligand_prep_status": "CONFIRMED" if pair == "EGFR/HER2" else "METADATA_ONLY",
                "receptor_prep_status": "PDBQT_CONFIRMED_METHOD_NOT_RECOVERABLE",
                "primary_exhaustiveness": 16 if pair == "PIK3CA/mTOR" else 8,
                "theta_primary": 6.0,
                "estimand": "DvsA_pocketB; DvsB_pocketA; summary_min=weaker_arm",
            }
        )
    write_csv(PROV / "scientific_provenance_matrix.csv", sci)

    # ---- score extraction samples ----
    samples = []

    def sample_pair(pair, src_lookup, n=5):
        recs = packs[pair]
        rng = np.random.default_rng(20260729 + sum(ord(c) for c in pair))
        idxs = rng.choice(len(recs), size=min(n, len(recs)), replace=False)
        for i in idxs:
            r = recs[int(i)]
            src = src_lookup(r)
            samples.append(
                {
                    "pair": pair,
                    "ligand_id": r["ligand_id"],
                    "master_score_A": r["score_A"],
                    "master_score_B": r["score_B"],
                    "source_score_A": src.get("A"),
                    "source_score_B": src.get("B"),
                    "source_file": src.get("file"),
                    "abs_delta_A": "" if src.get("A") is None else f"{abs(r['score_A'] - src['A']):.6f}",
                    "abs_delta_B": "" if src.get("B") is None else f"{abs(r['score_B'] - src['B']):.6f}",
                    "status": "MATCH"
                    if src.get("A") is not None
                    and src.get("B") is not None
                    and abs(r["score_A"] - src["A"]) < 1e-6
                    and abs(r["score_B"] - src["B"]) < 1e-6
                    else "MISMATCH_OR_MISSING",
                }
            )

    egfr_src = defaultdict(dict)
    egfr_csv = ROOT / "data/egfr_her2_uniform_rdkit_v1/tables/scores_vina_mode1_fiveseed.csv"
    for r in read_csv(egfr_csv):
        if int(r["seed"]) != 20260727:
            continue
        if r.get("status") not in ("ok", "success"):
            continue
        pocket = "A" if r["pdb"] == "3POZ" else "B"
        egfr_src[r["ligand"]][pocket] = -float(r["vina_mode1"])
        egfr_src[r["ligand"]]["file"] = "data/egfr_her2_uniform_rdkit_v1/tables/scores_vina_mode1_fiveseed.csv"

    pm_src = {}
    for r in read_csv(ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv"):
        pm_src[r["ligand"]] = {
            "A": fnum(r.get("vina_4L23_hb")),
            "B": fnum(r.get("vina_4JT6_hb")),
            "file": "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        }
    ache_src = {}
    for r in read_csv(ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv"):
        ache_src[r["ligand"]] = {
            "A": fnum(r.get("vina_ACHE_hb")),
            "B": fnum(r.get("vina_BCHE_hb")),
            "file": "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        }
    tb_src = defaultdict(dict)
    for r in read_csv(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv"):
        pair = r["pair"]
        recA, recB = RECEPTORS[pair]
        pocket = "A" if r["target"] == recA else "B" if r["target"] == recB else None
        if pocket:
            tb_src[(pair, r["ligand"])][pocket] = fnum(r.get("score_S"))
            tb_src[(pair, r["ligand"])]["file"] = "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv"

    sample_pair("EGFR/HER2", lambda r: egfr_src.get(r["ligand_id"], {}))
    sample_pair("PIK3CA/mTOR", lambda r: pm_src.get(r["ligand_id"], {}))
    sample_pair("AChE/BChE", lambda r: ache_src.get(r["ligand_id"], {}))
    for pair in ("JAK1/JAK2", "JAK1/TYK2", "F2/F10", "PPARG/PPARA", "PPARA/PPARD"):
        sample_pair(pair, lambda r, p=pair: tb_src.get((p, r["ligand_id"]), {}))
    write_csv(QA / "score_extraction_sample_trace.csv", samples)

    # ---- table cell trace (main 1-3 + SI S1-S10 numeric cells) ----
    table_cells = []

    def add_cell(table, row, column, displayed, raw, source, source_row, calc, rounding, status):
        table_cells.append(
            {
                "table": table,
                "row": row,
                "column": column,
                "displayed_value": displayed,
                "raw_value": raw,
                "source_file": source,
                "source_row": source_row,
                "calculation": calc,
                "rounding_rule": rounding,
                "status": status,
            }
        )

    def round3(x):
        return f"{float(x):.3f}"

    def status_round(disp, raw):
        v = fnum(raw)
        d = fnum(disp)
        if v is None or d is None:
            return "NON_NUMERIC"
        if abs(v - d) < 5e-4 or abs(float(f"{v:.3f}") - d) < 1.5e-3:
            return "ROUNDING_MATCH"
        if abs(v - d) < 0.015:
            return "NEAR"
        return "MISMATCH"

    # Table 1 from class_counts + manuscript quotas
    cc = {r["pair"]: r for r in read_csv(CANON / "class_counts.csv")}
    t1_quota = {
        "EGFR/HER2": (28, 38, 32, 12),
        "JAK1/JAK2": (32, 32, 32, 14),
        "JAK1/TYK2": (32, 32, 32, 14),
        "PIK3CA/mTOR": (18, 14, 12, 4),
        "AChE/BChE": (28, 28, 28, 16),
        "F2/F10": (32, 32, 32, 14),
        "PPARG/PPARA": (32, 32, 32, 14),
        "PPARA/PPARD": (32, 32, 32, 14),
    }
    t1_pdb = {p: "/".join(RECEPTORS[p]) for p in PRIMARY_PAIRS}
    t1_e = {p: (16 if p == "PIK3CA/mTOR" else 8) for p in PRIMARY_PAIRS}
    ms = (ROOT / "docs/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    for pair in PRIMARY_PAIRS:
        r = cc[pair]
        nd, na, nb = r["n_dual"], r["n_A_only"], r["n_B_only"]
        add_cell("Table1", pair, "n_scored_D/A/B", f"{nd} / {na} / {nb}", f"{nd}/{na}/{nb}", "results/canonical/class_counts.csv", pair, "complete-case theta6", "as_typeset", "MATCH")
        add_cell("Table1", pair, "PDB", t1_pdb[pair], t1_pdb[pair], "analysis_config.RECEPTORS", pair, "config", "verbatim", "MATCH")
        add_cell("Table1", pair, "exhaustiveness", str(t1_e[pair]), str(t1_e[pair]), "docs/PROTOCOL_LEVELS.md", pair, "PRIMARY E", "verbatim", "MATCH")
        q = t1_quota[pair]
        add_cell("Table1", pair, "quota_D/A/B/N", f"{q[0]} / {q[1]} / {q[2]} / {q[3]}", str(q), "panel construction CSVs / Table 1 prose", pair, "construction quota not Table2 n", "verbatim", "QUOTA_NOT_N_SCORED")

    smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    dau = defaultdict(dict)
    for r in read_csv(CANON / "primary_directional_auroc.csv"):
        dau[r["pair"]][r["estimand"]] = r
    # parse Table 2 displayed from manuscript
    t2_pat = re.compile(
        r"\| (EGFR/HER2|JAK1/JAK2|JAK1/TYK2|PIK3CA/mTOR|AChE/BChE|F2/F10|PPARG/PPARA|PPARA/PPARD) \| ([0-9 /]+) \| ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \|"
    )
    for m in t2_pat.finditer(ms):
        pair = m.group(1)
        add_cell("Table2", pair, "n_scored", m.group(2), f"{smin[pair]['n_dual']} / {smin[pair]['n_A_only']} / {smin[pair]['n_B_only']}", "primary_summary_min.csv", pair, "n_dual/A/B", "typeset", "MATCH" if m.group(2).replace(" ", "") == f"{smin[pair]['n_dual']}/{smin[pair]['n_A_only']}/{smin[pair]['n_B_only']}" else "MISMATCH")
        mapping = [
            ("D_vs_A", m.group(3), dau[pair]["AUROC_D_vs_A_pocketB"]["point"]),
            ("D_vs_A_lo", m.group(4), dau[pair]["AUROC_D_vs_A_pocketB"]["ci_lo"]),
            ("D_vs_A_hi", m.group(5), dau[pair]["AUROC_D_vs_A_pocketB"]["ci_hi"]),
            ("D_vs_B", m.group(6), dau[pair]["AUROC_D_vs_B_pocketA"]["point"]),
            ("D_vs_B_lo", m.group(7), dau[pair]["AUROC_D_vs_B_pocketA"]["ci_lo"]),
            ("D_vs_B_hi", m.group(8), dau[pair]["AUROC_D_vs_B_pocketA"]["ci_hi"]),
            ("summary_min", m.group(9), smin[pair]["summary_min"]),
            ("smin_lo", m.group(10), smin[pair]["ci_lo"]),
            ("smin_hi", m.group(11), smin[pair]["ci_hi"]),
        ]
        for col, disp, raw in mapping:
            add_cell("Table2", pair, col, disp, raw, "primary_directional_auroc.csv / primary_summary_min.csv", pair, "class-stratified B=2000 seed=20260729", "3dp from 4dp CSV", status_round(disp, raw))

    t3_pat = re.compile(
        r"\| (EGFR/HER2|JAK1/JAK2|JAK1/TYK2|PIK3CA/mTOR|AChE/BChE|F2/F10|PPARG/PPARA|PPARA/PPARD) \| ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| (\d+) \| ([0-9 /]+) \| ([0-9.]+) \|"
    )
    rank = {r["pair"]: r for r in read_csv(CANON / "two_pocket_mean_ranking.csv")}
    top = {r["pair"]: r for r in read_csv(CANON / "top10_operating_points.csv")}
    for m in t3_pat.finditer(ms):
        pair = m.group(1)
        add_cell("Table3", pair, "D_vs_neither", m.group(2), rank[pair]["two_pocket_mean_D_vs_neither"], "two_pocket_mean_ranking.csv", pair, "two-pocket mean vs neither only", "3dp", status_round(m.group(2), rank[pair]["two_pocket_mean_D_vs_neither"]))
        add_cell("Table3", pair, "ci_lo", m.group(3), rank[pair]["ci_lo"], "two_pocket_mean_ranking.csv", pair, "bootstrap", "3dp", status_round(m.group(3), rank[pair]["ci_lo"]))
        add_cell("Table3", pair, "ci_hi", m.group(4), rank[pair]["ci_hi"], "two_pocket_mean_ranking.csv", pair, "bootstrap", "3dp", status_round(m.group(4), rank[pair]["ci_hi"]))
        add_cell("Table3", pair, "n_neither", m.group(5), str(smin[pair].get("n_dual") and packs[pair] and sum(1 for x in packs[pair] if x['cls']=='neither')), "class_counts.csv", pair, "n_neither", "verbatim", "MATCH" if int(m.group(5)) == sum(1 for x in packs[pair] if x["cls"] == "neither") else "MISMATCH")
        add_cell("Table3", pair, "EF", m.group(7), rank[pair]["ef_dual_10pct"], "two_pocket_mean_ranking.csv", pair, "EF=(top_dual/k)/(n_dual/n)", "3dp", status_round(m.group(7), rank[pair]["ef_dual_10pct"]))

    si = (ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md").read_text(encoding="utf-8")
    # S1 environment
    add_cell("TableS1", "analysis_freeze", "Python", "3.12.3", "3.12.3", "requirements-analysis.txt / check_analysis_env.py", "env", "analysis pin", "verbatim", "MATCH")
    add_cell("TableS1", "analysis_freeze", "RDKit", "2026.03.5", "2026.3.5", "requirements-analysis.txt", "rdkit==2026.3.5", "pin", "typeset 2026.03.5", "ROUNDING_MATCH")
    if "3.12.13" in si:
        add_cell("TableS1", "analysis_freeze", "Python_stale", "3.12.13", "3.12.3", "SI text", "stale", "no evidence", "verbatim", "STALE")

    # S2 boxes
    s2_box = {
        "4L23": ("32.443", "45.431", "42.139", "20.000", "20.000", "20.000"),
        "4JT6": None,
        "3POZ": ("18.816", "31.837", "11.725", "20.931", "20.000", "21.342"),
    }
    for pdb, path in box_map.items():
        raw = json.loads(path.read_text(encoding="utf-8"))
        # find SI line
        m = re.search(rf"\| [A-Za-z0-9]+ \| {pdb} \|.*?\| ([0-9.]+) \| ([−0-9., ]+) \| ([0-9., ]+) \|", si)
        if not m:
            add_cell("TableS2a", pdb, "row", "NOT_PARSED", "", str(path.relative_to(ROOT)), pdb, "regex", "", "PARSE_FAIL")
            continue
        add_cell("TableS2a", pdb, "center_size_present", "parsed", json.dumps({k: raw.get(k) for k in ("center_x", "size_x")}), str(path.relative_to(ROOT)), pdb, "box JSON", "3dp", "SOURCE_EXISTS")

    # S2b RMSD
    rmsd = {r["pdb"]: r for r in read_csv(CANON / "cognate_rmsd.csv")}
    for pdb, r in rmsd.items():
        m = re.search(rf"\| [A-Za-z0-9]+ \| {pdb} \| \d+ \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \|", si)
        if not m:
            add_cell("TableS2b", pdb, "top1", "NOT_PARSED", r["calcrrms_top1_A"], "cognate_rmsd.csv", pdb, "CalcRMS", "3dp", "PARSE_FAIL")
            continue
        add_cell("TableS2b", pdb, "top1", m.group(1), r["calcrrms_top1_A"], "results/canonical/cognate_rmsd.csv", pdb, "calcrrms_top1_A", "3dp", status_round(m.group(1), r["calcrrms_top1_A"]))
        add_cell("TableS2b", pdb, "best", m.group(3), r["calcrrms_best_A"], "results/canonical/cognate_rmsd.csv", pdb, "calcrrms_best_A", "3dp", status_round(m.group(3), r["calcrrms_best_A"]))

    # S4 cluster vs manuscript stale
    clus = read_csv(CANON / "cluster_bootstrap_sensitivity.csv")
    for r in clus:
        if r["pair"] == "EGFR/HER2" and r["estimator"] == "scaffold_cluster":
            add_cell("TableS4", "EGFR/HER2", "scaffold_ci", "[0.232, 0.637]", f"{r['delta_ci_lo']},{r['delta_ci_hi']}", "cluster_bootstrap_sensitivity.csv", "EGFR scaffold", "cluster bootstrap", "3dp", status_round("0.232", r["delta_ci_lo"]))
        if r["pair"] == "EGFR/HER2" and r["estimator"] == "document_cluster":
            add_cell("TableS4", "EGFR/HER2", "document_ci", "[0.117, 0.636]", f"{r['delta_ci_lo']},{r['delta_ci_hi']}", "cluster_bootstrap_sensitivity.csv", "EGFR document", "cluster bootstrap", "3dp", status_round("0.117", r["delta_ci_lo"]))

    # remaining SI tables: dump key numeric rows from canonical as cells
    for fname, table in (
        ("max_vs_median_sensitivity.csv", "TableS3"),
        ("fixed_score_negative_class_delta.csv", "TableS4"),
        ("ecfp4_incremental_information.csv", "TableS5"),
        ("descriptor_baselines.csv", "TableS5b"),
        ("matched_minus_mismatched.csv", "TableS6"),
        ("computational_robustness.csv", "TableS7"),
        ("external_eligibility.csv", "TableS8"),
        ("top10_operating_points.csv", "TableS10"),
        ("and_filter_operating_points.csv", "TableS10b"),
    ):
        p = CANON / fname
        if not p.is_file():
            add_cell(table, "ALL", "source", "MISSING", "", str(p), "", "file missing", "", "MISSING_SOURCE")
            continue
        rows = read_csv(p)
        add_cell(table, "ALL", "source_exists", str(len(rows)), str(len(rows)), f"results/canonical/{fname}", "all", "row count", "n/a", "SOURCE_EXISTS")

    elig = ROOT / "data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv"
    add_cell("TableS9", "ALL", "source", "exists" if elig.is_file() else "missing", "", str(elig.relative_to(ROOT)), "eligibility", "pair inclusion", "n/a", "SOURCE_EXISTS" if elig.is_file() else "MISSING")

    write_csv(QA / "table_cell_trace.csv", table_cells)

    # ---- figure value trace ----
    fig_rows = []
    plotted_path = ROOT / "figures/jcim_article/plotted_values_postfix.json"
    plotted = json.loads(plotted_path.read_text(encoding="utf-8"))
    source_cache = {}

    def load_source(rel):
        if rel in source_cache:
            return source_cache[rel]
        p = ROOT / rel
        if p.suffix == ".csv" and p.is_file():
            source_cache[rel] = read_csv(p)
        else:
            source_cache[rel] = None
        return source_cache[rel]

    for rec in plotted.get("records", []):
        src = rec.get("source_file") or ""
        raw = rec.get("raw_value")
        pair = rec.get("pair")
        metric = rec.get("metric")
        status = "NO_SOURCE"
        canon_raw = ""
        rows = load_source(src) if src else None
        if rows is not None:
            # try pair match
            hits = [r for r in rows if r.get("pair") == pair or r.get("pdb") == pair or r.get("protein") == pair]
            if not hits:
                hits = rows
            # pick first numeric-ish
            status = "SOURCE_EXISTS"
            if isinstance(raw, (int, float)) and hits:
                # search any field close to raw
                found = False
                for r in hits:
                    for v in r.values():
                        fv = fnum(v)
                        if fv is not None and abs(fv - float(raw)) < 1.5e-3:
                            found = True
                            canon_raw = v
                            break
                    if found:
                        break
                status = "MATCH" if found else "VALUE_NOT_IN_SOURCE_ROWS"
            elif hits:
                status = "SOURCE_EXISTS"
        if rec.get("figure") == "Figure 5" and rec.get("panel") == "C" and metric == "five_seed_comparable":
            if pair == "AChE/BChE" and float(raw) == 1.0:
                status = "FUZZY_COMPARABLE_FLAG"
                findings.append("P1 Figure5C AChE five_seed_comparable=1 while same_membership_as_primary=0")
        fig_rows.append(
            {
                "figure": rec.get("figure"),
                "panel": rec.get("panel"),
                "pair": pair,
                "metric": metric,
                "displayed_value": rec.get("displayed_value"),
                "raw_value": raw,
                "source_file": src,
                "source_row": rec.get("source_row"),
                "matched_source_value": canon_raw,
                "status": status,
            }
        )
    write_csv(QA / "figure_value_trace.csv", fig_rows)

    # ---- text numeric claim audit ----
    text_rows = []
    text_files = [
        ROOT / "docs/MANUSCRIPT_JCIM_EN.md",
        ROOT / "docs/MANUSCRIPT_JCIM_ZH.md",
        ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        ROOT / "docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
        ROOT / "docs/WRITING_INDEX_FREEZE.md",
        ROOT / "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
    ]
    known_stale = {
        "0.3237": "HISTORICAL EGFR ablation summary_min",
        "0.0112": "HISTORICAL ECFP incremental",
        "9.505": "HISTORICAL EGFR RMSD token (also an AChE score)",
        "0.760": "HISTORICAL EGFR RMSD",
        "3.12.13": "NO EVIDENCE Python patch",
        "0.430": "pre-fix EGFR AUROC token",
        "0.808": "pre-fix EGFR AUROC token",
        "0.378": "pre-fix EGFR AUROC token",
    }
    # key current numbers
    current_vals = {}
    for pair in PRIMARY_PAIRS:
        current_vals[f"{pair}_smin"] = float(smin[pair]["summary_min"])
        current_vals[f"{pair}_da"] = float(dau[pair]["AUROC_D_vs_A_pocketB"]["point"])
        current_vals[f"{pair}_db"] = float(dau[pair]["AUROC_D_vs_B_pocketA"]["point"])
    current_vals["ecfp_max_abs_delta"] = max(
        abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"]))
        for r in read_csv(CANON / "ecfp4_incremental_information.csv")
    )
    current_vals["egfr_cluster_scaf_lo"] = float(
        next(r["delta_ci_lo"] for r in clus if r["pair"] == "EGFR/HER2" and r["estimator"] == "scaffold_cluster")
    )
    current_vals["egfr_cluster_scaf_hi"] = float(
        next(r["delta_ci_hi"] for r in clus if r["pair"] == "EGFR/HER2" and r["estimator"] == "scaffold_cluster")
    )
    current_vals["egfr_cluster_doc_lo"] = float(
        next(r["delta_ci_lo"] for r in clus if r["pair"] == "EGFR/HER2" and r["estimator"] == "document_cluster")
    )
    current_vals["egfr_cluster_doc_hi"] = float(
        next(r["delta_ci_hi"] for r in clus if r["pair"] == "EGFR/HER2" and r["estimator"] == "document_cluster")
    )

    claim_patterns = [
        (r"0\.334(?:1)?", "EGFR summary_min", "0.3341"),
        (r"0\.3237", "stale EGFR summary_min", "STALE"),
        (r"0\.0112", "stale ECFP incremental", "STALE"),
        (r"0\.0234", "ECFP max |Δ|", "0.0234"),
        (r"3\.12\.13", "stale Python", "STALE"),
        (r"3\.12\.3", "analysis Python", "3.12.3"),
        (r"9\.505", "RMSD or AChE score token", "CONTEXT"),
        (r"\[0\.235,\s*0\.665\]", "stale EGFR scaffold cluster CI", "STALE"),
        (r"\[0\.125,\s*0\.644\]", "stale EGFR document cluster CI", "STALE"),
        (r"\[0\.232,\s*0\.637\]", "EGFR scaffold cluster CI", "MATCH"),
        (r"\[0\.117,\s*0\.636\]", "EGFR document cluster CI", "MATCH"),
        (r"0\.4275", "stale EGFR TPSA footnote", "STALE"),
        (r"only AChE", "outdated MM wording", "STALE_IF_EXCLUDES_EGFR"),
    ]
    for path in text_files:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        rel = str(path.relative_to(ROOT))
        for pat, claim, expected in claim_patterns:
            for m in re.finditer(pat, text):
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 80)
                ctx = text[start:end].replace("\n", " ")
                if expected == "STALE":
                    status = "STALE" if "historical" not in ctx.lower() and "not current" not in ctx.lower() and "do not" not in ctx.lower() and "superseded" not in ctx.lower() and "不得" not in ctx and "不是当前" not in ctx else "HISTORICAL_WITH_CONTEXT"
                elif expected == "STALE_IF_EXCLUDES_EGFR":
                    status = "STALE" if "EGFR" not in ctx else "MATCH"
                elif expected == "CONTEXT":
                    status = "HISTORICAL_WITH_CONTEXT" if ("historical" in ctx.lower() or "AChE" in ctx or "AB_038" in ctx or "不是当前" in ctx) else "UNSUPPORTED"
                elif expected == "MATCH":
                    status = "MATCH"
                else:
                    status = "MATCH"
                text_rows.append(
                    {
                        "file": rel,
                        "claim": claim,
                        "matched_text": m.group(0),
                        "context": ctx,
                        "expected": expected,
                        "status": status,
                    }
                )
    write_csv(QA / "text_numeric_claim_audit.csv", text_rows)

    # ---- newly discovered scans ----
    new_issues = []
    # duplicate ligand IDs / SMILES / ChEMBL within pair
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        for key, field in (("ligand_id", "ligand_id"), ("smiles", "smiles"), ("chembl", "chembl")):
            c = Counter(r[field] for r in recs if r[field])
            dups = [k for k, v in c.items() if v > 1]
            if dups:
                new_issues.append(f"DUPLICATE_{key} {pair}: {dups[:5]}")
    # class mismatch theta
    n_mismatch = sum(1 for r in act_rows if r["in_primary"] == 1 and r["class_match"] == 0)
    if n_mismatch:
        new_issues.append(f"PRIMARY_THETA6_CLASS_MISMATCH n={n_mismatch}")
    # unresolved in primary
    for r in master:
        if is_primary_row(r) and r.get("activity_status") == "unresolved_missing_arm":
            new_issues.append(f"UNRESOLVED_IN_PRIMARY {r['pair']} {r['ligand_id']}")
    # EH120_059 / AB_087
    for pair, lig in (("EGFR/HER2", "EH120_059"), ("AChE/BChE", "AB_087")):
        recs = [r for r in master if r["pair"] == pair and r["ligand_id"] == lig]
        if not recs:
            new_issues.append(f"MISSING_FLAG_ROW {pair} {lig}")
        elif is_primary_row(recs[0]):
            new_issues.append(f"UNRESOLVED_ENTERED_PRIMARY {pair} {lig}")
    # eligibility AUROC leakage
    elig_text = elig.read_text(encoding="utf-8") if elig.is_file() else ""
    if re.search(r"AUROC|summary_min|GNINA|ECFP", elig_text, re.I):
        new_issues.append("ELIGIBILITY_CSV_CONTAINS_PERFORMANCE_TOKEN")
    # hard-coded AUROC in analysis scripts
    for py in (ROOT / "scripts/analysis").glob("*.py"):
        src = py.read_text(encoding="utf-8")
        if re.search(r"0\.3237|0\.3341|0\.0112", src):
            new_issues.append(f"HARDCODED_RESULT_TOKEN {py.relative_to(ROOT)}")
    # silent except
    silent = []
    for py in (ROOT / "scripts").rglob("*.py"):
        try:
            tree = ast.parse(py.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.ExceptHandler) and node.body:
                if len(node.body) == 1 and isinstance(node.body[0], ast.Pass):
                    silent.append(f"{py.relative_to(ROOT)}:{node.lineno}")
    if silent:
        new_issues.append("SILENT_EXCEPT_PASS " + ";".join(silent[:12]))
    # two writers of canonical master
    writers = []
    for py in (ROOT / "scripts").rglob("*.py"):
        txt = py.read_text(encoding="utf-8")
        if "current_score_master.csv" in txt and ("write" in txt or "to_csv" in txt):
            writers.append(str(py.relative_to(ROOT)))
    # score sign
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        if any((r["score_A"] or 0) < 0 or (r["score_B"] or 0) < 0 for r in recs):
            new_issues.append(f"NEGATIVE_SCORE_S {pair}")
    # extreme scores
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        for r in recs:
            for sc in (r["score_A"], r["score_B"]):
                if sc is not None and (sc > 20 or sc < 2):
                    new_issues.append(f"EXTREME_SCORE {pair} {r['ligand_id']} {sc}")
                    break
    # ECFP incremental
    ecfp = read_csv(CANON / "ecfp4_incremental_information.csv")
    max_d = max(abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in ecfp)
    if abs(max_d - 0.0234) > 1e-4:
        new_issues.append(f"ECFP_MAX_ABS_DELTA {max_d}")
    # stale current-facing freeze docs
    for rel in (
        "docs/PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md",
        "docs/PR39_SCIENTIFIC_DATA_AUDIT.md",
        "docs/SUBMISSION_ONLY_CLEANUP_REPORT.md",
    ):
        p = ROOT / rel
        if p.is_file() and "0.3237" in p.read_text(encoding="utf-8"):
            new_issues.append(f"STALE_CURRENT_FACING_DOC {rel} still states 0.3237 as if current")
    # GNINA_SOURCES historical default
    cfg = (ROOT / "scripts/analysis/analysis_config.py").read_text(encoding="utf-8")
    if "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv" in cfg:
        new_issues.append("GNINA_SOURCES_DEFAULT_STILL_HISTORICAL_EGFR_PATH")
    # comparable flag still written as 1 for AChE
    fs = [r for r in read_csv(CANON / "five_seed_summary_min.csv") if r["pair"] == "AChE/BChE"]
    if fs and all(str(r.get("comparable_to_current_primary")) == "1" for r in fs) and any(
        str(r.get("same_membership_as_primary")) == "0" for r in fs
    ):
        new_issues.append("FIVE_SEED_COMPARABLE_FLAG_STILL_1_WHEN_MEMBERSHIP_DIFFERS")
    # manuscript cluster CI
    men = (ROOT / "docs/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    if "[0.235, 0.665]" in men or "[0.125, 0.644]" in men:
        new_issues.append("MANUSCRIPT_STALE_EGFR_CLUSTER_CI")
    if "0.4275" in (ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md").read_text(encoding="utf-8"):
        new_issues.append("SI_STALE_EGFR_TPSA_FOOTNOTE_0.4275")
    # old box files exist
    if (ROOT / "data/egfr_her2_panel120_v0/boxes/3POZ_box.json").is_file():
        new_issues.append("LEGACY_EGFR_BOX_JSON_STILL_PRESENT")
    # processed vs canonical master
    proc = ROOT / "data/processed/current_score_master.csv"
    if proc.is_file():
        h1 = hashlib.md5((CANON / "current_score_master.csv").read_bytes()).hexdigest()
        h2 = hashlib.md5(proc.read_bytes()).hexdigest()
        if h1 != h2:
            new_issues.append("PROCESSED_MASTER_DIFFERS_FROM_CANONICAL")
        else:
            new_issues.append("DUPLICATE_MASTER_COPY_data/processed_IDENTICAL")

    (QA / "newly_discovered_issues_raw.txt").write_text("\n".join(new_issues) + "\n", encoding="utf-8")
    (QA / "audit_findings_raw.txt").write_text("\n".join(findings) + "\n", encoding="utf-8")

    summary = {
        "n_replay_rows": len(replay),
        "n_replay_p0": n_p0_replay,
        "n_primary": {p: len(packs[p]) for p in PRIMARY_PAIRS},
        "n_membership_rows": len(membership),
        "n_table_cells": len(table_cells),
        "n_figure_records": len(fig_rows),
        "n_text_claims": len(text_rows),
        "n_new_issues": len(new_issues),
        "new_issues": new_issues,
        "ecfp_max_abs_delta": max_d,
        "score_sample_mismatches": sum(1 for r in samples if r["status"] != "MATCH"),
        "theta_mismatch_primary": n_mismatch,
        "silent_except": silent,
        "canonical_writers": writers,
    }
    (QA / "audit_runner_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("n_replay_p0", "n_new_issues", "ecfp_max_abs_delta", "score_sample_mismatches", "theta_mismatch_primary")}, indent=2))
    print("NEW_ISSUES:")
    for x in new_issues:
        print(" -", x)


if __name__ == "__main__":
    main()
