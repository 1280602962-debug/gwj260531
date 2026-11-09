#!/usr/bin/env python3
"""Recompute eight-pair primary statistics from current_score_master.

Scheme B class-stratified bootstrap (B=2000, seed=20260729). Zero-dock.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import (  # noqa: E402
    N_BOOT,
    SEED,
    auroc,
    assign_fourclass,
    assign_strict,
    ci_crosses,
    cluster_delta,
    fixed_score_delta_stratified,
    fixed_score_delta_unstratified,
    matched_mismatched_paired,
    summary_min_stratified,
    summary_min_unstratified,
    stratified_auroc_ci,
)

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
CANON = ROOT / "results" / "canonical"
MASTER = CANON / "current_score_master.csv"


def fnum(v):
    if v is None or v == "":
        return None
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


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
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def r4(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return ""
    return f"{float(x):.4f}"


def load_main() -> dict[str, list[dict]]:
    rows = read_csv(MASTER)
    by = {p: [] for p in PRIMARY_PAIRS}
    for r in rows:
        if r["analysis_set"] != "main" or r["complete_case"] not in ("1", 1, "True"):
            continue
        if str(r.get("activity_eligible", "1")) not in ("1", "True"):
            continue
        if r["pair"] not in by:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["score_mean"] = (rec["score_A"] + rec["score_B"]) / 2
        by[r["pair"]].append(rec)
    return by


def arrays(recs):
    cls = np.array([r["cls"] for r in recs])
    sa = np.array([r["score_A"] for r in recs], dtype=float)
    sb = np.array([r["score_B"] for r in recs], dtype=float)
    return sa, sb, cls


def compute_directional(packs):
    rows = []
    smin_rows = []
    for pair, recs in packs.items():
        sa, sb, cls = arrays(recs)
        stats = summary_min_stratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
        u_pt, u_lo, u_hi, _ = summary_min_unstratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
        for arm, point, lo, hi in (
            ("AUROC_D_vs_A_pocketB", stats["auroc_D_vs_A_pocketB"], stats["auroc_D_vs_A_ci_lo"], stats["auroc_D_vs_A_ci_hi"]),
            ("AUROC_D_vs_B_pocketA", stats["auroc_D_vs_B_pocketA"], stats["auroc_D_vs_B_ci_lo"], stats["auroc_D_vs_B_ci_hi"]),
        ):
            rows.append(
                {
                    "pair": pair,
                    "estimand": arm,
                    "point": r4(point),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "n_dual": stats["n_dual"],
                    "n_A_only": stats["n_A_only"],
                    "n_B_only": stats["n_B_only"],
                    "bootstrap": "class_stratified",
                    "B": N_BOOT,
                    "seed": SEED,
                    "ci_excludes_0.5": int(not ci_crosses(lo, hi, 0.5)),
                }
            )
        smin_rows.append(
            {
                "pair": pair,
                "summary_min": r4(stats["summary_min"]),
                "ci_lo": r4(stats["summary_min_ci_lo"]),
                "ci_hi": r4(stats["summary_min_ci_hi"]),
                "auroc_D_vs_A_pocketB": r4(stats["auroc_D_vs_A_pocketB"]),
                "auroc_D_vs_B_pocketA": r4(stats["auroc_D_vs_B_pocketA"]),
                "weaker_arm": stats["weaker_arm"],
                "n_dual": stats["n_dual"],
                "n_A_only": stats["n_A_only"],
                "n_B_only": stats["n_B_only"],
                "bootstrap": "class_stratified_shared_dual",
                "B": N_BOOT,
                "seed": SEED,
                "ci_excludes_0.5": int(not ci_crosses(stats["summary_min_ci_lo"], stats["summary_min_ci_hi"], 0.5)),
                "non_stratified_ci_lo": r4(u_lo),
                "non_stratified_ci_hi": r4(u_hi),
            }
        )
    return rows, smin_rows


def compute_fixed_score(packs):
    rows = []
    sens = []
    for pair, recs in packs.items():
        dual = [r for r in recs if r["cls"] == "dual"]
        a_only = [r for r in recs if r["cls"] == "A_only"]
        b_only = [r for r in recs if r["cls"] == "B_only"]
        neither = [r for r in recs if r["cls"] == "neither"]
        specs = [
            ("D_vs_B_or_neither_pocketA", "score_A", b_only),
            ("D_vs_A_or_neither_pocketB", "score_B", a_only),
        ]
        for contrast, key, selective in specs:
            auc_s, auc_n, delta, lo, hi = fixed_score_delta_stratified(
                [r[key] for r in dual],
                [r[key] for r in selective],
                [r[key] for r in neither],
                n_boot=N_BOOT,
                seed=SEED,
            )
            _, _, _, ulo, uhi = fixed_score_delta_unstratified(
                [r[key] for r in dual],
                [r[key] for r in selective],
                [r[key] for r in neither],
                n_boot=N_BOOT,
                seed=SEED,
            )
            rec = {
                "pair": pair,
                "contrast": contrast,
                "score": "vina_A" if key == "score_A" else "vina_B",
                "auroc_dual_vs_selective": r4(auc_s),
                "auroc_dual_vs_neither": r4(auc_n),
                "delta_neither_minus_selective": r4(delta),
                "delta_ci_lo": r4(lo),
                "delta_ci_hi": r4(hi),
                "n_dual": len(dual),
                "n_selective": len(selective),
                "n_neither": len(neither),
                "bootstrap": "class_stratified_shared_dual",
                "B": N_BOOT,
                "seed": SEED,
                "ci_excludes_zero": int(not ci_crosses(lo, hi, 0.0)),
            }
            rows.append(rec)
            q_old = not ci_crosses(lo, hi, 0.0)
            q_new = not ci_crosses(ulo, uhi, 0.0)
            sens.append(
                {
                    "pair": pair,
                    "estimand": f"fixed_score_delta:{contrast}",
                    "primary_stratified_point": r4(delta),
                    "stratified_ci_lo": r4(lo),
                    "stratified_ci_hi": r4(hi),
                    "non_stratified_ci_lo": r4(ulo),
                    "non_stratified_ci_hi": r4(uhi),
                    "qualitative_reference": "CI excludes 0",
                    "qualitative_conclusion_changed": "yes" if q_old != q_new else "no",
                }
            )
        # two-pocket mean dual vs neither
        mean_d = [r["score_mean"] for r in dual]
        mean_n = [r["score_mean"] for r in neither]
        pt, lo, hi = stratified_auroc_ci(mean_d, mean_n, N_BOOT, SEED)
        rows.append(
            {
                "pair": pair,
                "contrast": "D_vs_neither_two_pocket_mean",
                "score": "mean(score_A,score_B)",
                "auroc_dual_vs_selective": "",
                "auroc_dual_vs_neither": r4(pt),
                "delta_neither_minus_selective": "",
                "delta_ci_lo": r4(lo),
                "delta_ci_hi": r4(hi),
                "n_dual": len(dual),
                "n_selective": "",
                "n_neither": len(neither),
                "bootstrap": "class_stratified",
                "B": N_BOOT,
                "seed": SEED,
                "ci_excludes_zero": "",
            }
        )
    return rows, sens


def compute_ranking(packs):
    rows = []
    for pair, recs in packs.items():
        n = len(recs)
        k = max(1, math.ceil(0.10 * n))
        ranked = sorted(recs, key=lambda r: (-r["score_mean"], r["ligand_id"]))
        top = ranked[:k]
        ct = Counter(r["cls"] for r in recs)
        topc = Counter(r["cls"] for r in top)
        n_dual = ct["dual"]
        ef = (topc["dual"] / k) / (n_dual / n) if n_dual else float("nan")
        rows.append(
            {
                "pair": pair,
                "n_ranked": n,
                "n_dual": n_dual,
                "n_A_only": ct["A_only"],
                "n_B_only": ct["B_only"],
                "n_neither": ct["neither"],
                "top_fraction": 0.10,
                "top_k": k,
                "top_dual": topc["dual"],
                "top_A_only": topc["A_only"],
                "top_B_only": topc["B_only"],
                "top_neither": topc["neither"],
                "top_dual_fraction": r4(topc["dual"] / k),
                "panel_dual_fraction": r4(n_dual / n),
                "ef_dual_10pct": r4(ef) if ef == ef else "",
                "top_ligand_ids": ";".join(r["ligand_id"] for r in top),
                "ranking_rule": "k=ceil(0.10 n); S_mean=(score_A+score_B)/2; EF=(top_dual/k)/(n_dual/n)",
                "tie_rule": "descending score_mean then ascending ligand ID",
            }
        )
    return rows


def compute_matched(packs):
    rows = []
    for pair, recs in packs.items():
        sa, sb, cls = arrays(recs)
        stats = matched_mismatched_paired(sa, sb, cls, n_boot=N_BOOT, seed=SEED, stratified=True)
        u = matched_mismatched_paired(sa, sb, cls, n_boot=N_BOOT, seed=SEED, stratified=False)
        rec = {"pair": pair, "set": "main_panel", **{k: (r4(v) if isinstance(v, float) else v) for k, v in stats.items()}}
        rec["ci_excludes_zero"] = int(stats["ci_excludes_zero"])
        rec["weaker_arm_switched"] = int(stats["weaker_arm_switched"])
        rec["bootstrap"] = "class_stratified_paired_ligand"
        rec["B"] = N_BOOT
        rec["seed"] = SEED
        rec["non_stratified_delta_ci_lo"] = r4(u["delta_ci_lo"])
        rec["non_stratified_delta_ci_hi"] = r4(u["delta_ci_hi"])
        rec["non_stratified_ci_excludes_zero"] = int(u["ci_excludes_zero"])
        rows.append(rec)
    return rows


def compute_label_sensitivity(packs):
    rows = []
    rules = [("theta_5.5", 5.5, False), ("theta_6.0", 6.0, False), ("theta_6.5", 6.5, False), ("strict_6.5_5.5", None, True)]
    for pair, recs in packs.items():
        for name, cut, is_strict in rules:
            labeled = []
            flips = 0
            n_neither = 0
            n_gray = 0
            n_missing_pchembl = 0
            for r in recs:
                pa, pb = fnum(r.get("pA")), fnum(r.get("pB"))
                lab = assign_strict(pa, pb) if is_strict else assign_fourclass(pa, pb, cut)
                if lab is None:
                    n_missing_pchembl += 1
                    continue
                if lab == "gray":
                    n_gray += 1
                    continue
                if lab == "neither":
                    n_neither += 1
                    continue
                if lab in ("dual", "A_only", "B_only"):
                    labeled.append({**r, "cls": lab})
                    if lab != r["primary_class_theta6"] and r["primary_class_theta6"] in ("dual", "A_only", "B_only"):
                        flips += 1
            n_dual = sum(r["cls"] == "dual" for r in labeled)
            n_a = sum(r["cls"] == "A_only" for r in labeled)
            n_b = sum(r["cls"] == "B_only" for r in labeled)
            # Relabel is evaluable iff every complete-case ligand has pA/pB.
            evaluable = n_missing_pchembl == 0 and (n_dual + n_a + n_b) >= 8
            if not evaluable:
                rows.append(
                    {
                        "pair": pair,
                        "label_rule": name,
                        "n_dual": n_dual,
                        "n_A_only": n_a,
                        "n_B_only": n_b,
                        "n_neither": n_neither,
                        "n_gray": n_gray,
                        "n_missing_pchembl": n_missing_pchembl,
                        "class_flips_vs_theta6": flips,
                        "auroc_D_vs_A": "",
                        "auroc_D_vs_B": "",
                        "summary_min": "",
                        "ci_lo": "",
                        "ci_hi": "",
                        "bootstrap": "class_stratified_shared_dual",
                        "B": N_BOOT,
                        "seed": SEED,
                        "status": "NOT_EVALUABLE",
                    }
                )
                continue
            sa = np.array([r["score_A"] for r in labeled])
            sb = np.array([r["score_B"] for r in labeled])
            cls = np.array([r["cls"] for r in labeled])
            stats = summary_min_stratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
            rows.append(
                {
                    "pair": pair,
                    "label_rule": name,
                    "n_dual": stats["n_dual"],
                    "n_A_only": stats["n_A_only"],
                    "n_B_only": stats["n_B_only"],
                    "n_neither": n_neither,
                    "n_gray": n_gray,
                    "n_missing_pchembl": n_missing_pchembl,
                    "class_flips_vs_theta6": flips,
                    "auroc_D_vs_A": r4(stats["auroc_D_vs_A_pocketB"]),
                    "auroc_D_vs_B": r4(stats["auroc_D_vs_B_pocketA"]),
                    "summary_min": r4(stats["summary_min"]),
                    "ci_lo": r4(stats["summary_min_ci_lo"]),
                    "ci_hi": r4(stats["summary_min_ci_hi"]),
                    "bootstrap": "class_stratified_shared_dual",
                    "B": N_BOOT,
                    "seed": SEED,
                    "status": "OK",
                }
            )
    return rows


def compute_max_median(packs):
    dump_path = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/eight_pair_dump_gated_v1/max_vs_median_ligand_v1.csv"
    adj_path = ROOT / "data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv"
    dump = {(r["pair"], r["ligand"]): r for r in read_csv(dump_path)} if dump_path.is_file() else {}
    adj = {(r["pair"], r["ligand"]): r for r in read_csv(adj_path)} if adj_path.is_file() else {}
    rows = []
    for pair, recs in packs.items():
        for agg, cls_field in (("max", "class_max"), ("median", "class_median")):
            labeled = []
            source = ""
            for r in recs:
                a = adj.get((pair, r["ligand_id"]))
                if a is not None:
                    if a.get("activity_status") != "both_arms_present":
                        continue
                    lab = a.get(cls_field, "")
                    source = "adjudicated high-confidence max/median + current score master"
                else:
                    d = dump.get((pair, r["ligand_id"]))
                    if not d:
                        continue
                    lab = d.get(cls_field, "")
                    if not source:
                        source = "ChEMBL37 dump-gated max/median + current score master"
                if lab in ("dual", "A_only", "B_only", "neither"):
                    labeled.append({**r, "cls": lab, "ligand_id": r["ligand_id"]})
            n_ligands = len(labeled)
            directional = [x for x in labeled if x["cls"] in ("dual", "A_only", "B_only")]
            if len(directional) < 8:
                continue
            sa = np.array([r["score_A"] for r in directional])
            sb = np.array([r["score_B"] for r in directional])
            cls = np.array([r["cls"] for r in directional])
            stats = summary_min_stratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
            rows.append(
                {
                    "pair": pair,
                    "aggregation": agg,
                    "n_ligands": n_ligands,
                    "n_dual": stats["n_dual"],
                    "n_A_only": stats["n_A_only"],
                    "n_B_only": stats["n_B_only"],
                    "n_neither": sum(1 for x in labeled if x["cls"] == "neither"),
                    "summary_min": r4(stats["summary_min"]),
                    "ci_lo": r4(stats["summary_min_ci_lo"]),
                    "ci_hi": r4(stats["summary_min_ci_hi"]),
                    "auroc_D_vs_A": r4(stats["auroc_D_vs_A_pocketB"]),
                    "auroc_D_vs_B": r4(stats["auroc_D_vs_B_pocketA"]),
                    "source": source,
                    "bootstrap": "class_stratified_shared_dual",
                    "B": N_BOOT,
                    "seed": SEED,
                    "class_flips_vs_max": "",
                }
            )
    by_pair = {}
    for r in rows:
        by_pair.setdefault(r["pair"], {})[r["aggregation"]] = r
    for pair, recs in packs.items():
        mx = by_pair.get(pair, {}).get("max")
        md = by_pair.get(pair, {}).get("median")
        if not mx or not md:
            continue
        a_max = {}
        a_med = {}
        for r in recs:
            a = adj.get((pair, r["ligand_id"]))
            if a is not None:
                if a.get("activity_status") != "both_arms_present":
                    continue
                a_max[r["ligand_id"]] = a.get("class_max", "")
                a_med[r["ligand_id"]] = a.get("class_median", "")
            else:
                dmp = dump.get((pair, r["ligand_id"]))
                if not dmp:
                    continue
                a_max[r["ligand_id"]] = dmp.get("class_max", "")
                a_med[r["ligand_id"]] = dmp.get("class_median", "")
        shared = [k for k in a_max if a_max[k] and a_med[k]]
        flips = sum(1 for k in shared if a_max[k] != a_med[k])
        for rec in (mx, md):
            rec["n_ligands"] = len(shared)
            rec["class_flips_vs_max"] = flips if rec["aggregation"] == "median" else 0
    return rows


def compute_cluster(packs):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    def scaffold(smi: str) -> str:
        mol = Chem.MolFromSmiles(smi or "")
        if mol is None:
            return "unparsed"
        try:
            return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "acyclic"
        except Exception:
            return "unparsed"

    groups_path = ROOT / "data/jcim_novelty_v0/tables/document_blocked_ligand_groups_v1.csv"
    groups = {}
    if groups_path.is_file():
        for r in read_csv(groups_path):
            groups[(r["pair"], r["ligand"])] = r["group_id"]
    rows = []
    # EGFR corrected scores + frozen document groups; JAK1/TYK2 scaffold from SMILES.
    for pair, contrast, key, sel_cls in (
        ("EGFR/HER2", "D_vs_B_or_neither_pocketA", "score_A", "B_only"),
        ("JAK1/TYK2", "D_vs_B_or_neither_pocketA", "score_A", "B_only"),
    ):
        recs = []
        for r in packs[pair]:
            if r["cls"] == "dual":
                use = "dual"
            elif r["cls"] == sel_cls:
                use = "selective"
            elif r["cls"] == "neither":
                use = "neither"
            else:
                continue
            recs.append(
                {
                    "cls": use,
                    "score": r[key],
                    "scaffold": scaffold(r.get("smiles", "")),
                    "doc_group": groups.get((pair, r["ligand_id"]), r["ligand_id"]),
                }
            )
        has_doc_map = any(p == pair for p, _ in groups)
        for estimator, gkey in (("ligand_stratified", None), ("scaffold_cluster", "scaffold"), ("document_cluster", "doc_group")):
            if estimator == "ligand_stratified":
                dual = [x["score"] for x in recs if x["cls"] == "dual"]
                sel = [x["score"] for x in recs if x["cls"] == "selective"]
                nei = [x["score"] for x in recs if x["cls"] == "neither"]
                auc_s, auc_n, delta, lo, hi = fixed_score_delta_stratified(dual, sel, nei, N_BOOT, SEED)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "estimator": estimator,
                        "n_groups": "",
                        "n_dual": len(dual),
                        "n_selective": len(sel),
                        "n_neither": len(nei),
                        "delta_point": r4(delta),
                        "delta_ci_lo": r4(lo),
                        "delta_ci_hi": r4(hi),
                        "n_attempted_boot": N_BOOT,
                        "n_valid_boot": N_BOOT,
                        "excludes_zero": int(not ci_crosses(lo, hi, 0.0)),
                        "status": "ok",
                        "note": "primary ligand-level class-stratified bootstrap; B=2000 seed=20260729",
                        "seed": SEED,
                    }
                )
            elif estimator == "document_cluster" and not has_doc_map:
                dual = [x["score"] for x in recs if x["cls"] == "dual"]
                sel = [x["score"] for x in recs if x["cls"] == "selective"]
                nei = [x["score"] for x in recs if x["cls"] == "neither"]
                auc_s, auc_n, delta, lo, hi = fixed_score_delta_stratified(dual, sel, nei, N_BOOT, SEED)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "estimator": estimator,
                        "n_groups": "",
                        "n_dual": len(dual),
                        "n_selective": len(sel),
                        "n_neither": len(nei),
                        "delta_point": r4(delta),
                        "delta_ci_lo": "",
                        "delta_ci_hi": "",
                        "n_attempted_boot": N_BOOT,
                        "n_valid_boot": 0,
                        "excludes_zero": "",
                        "status": "unresolved_mapping_unavailable",
                        "note": "not recomputed: ligand-document grouping map not deposited and ChEMBL 37 sqlite unavailable; do not copy a previous CI",
                        "seed": SEED,
                    }
                )
            else:
                stats = cluster_delta(recs, gkey, N_BOOT, SEED)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "estimator": estimator,
                        "n_groups": stats["n_groups"],
                        "n_dual": stats["n_dual"],
                        "n_selective": stats["n_selective"],
                        "n_neither": stats["n_neither"],
                        "delta_point": r4(stats["delta_point"]),
                        "delta_ci_lo": r4(stats["delta_ci_lo"]),
                        "delta_ci_hi": r4(stats["delta_ci_hi"]),
                        "n_attempted_boot": N_BOOT,
                        "n_valid_boot": stats["n_valid_boot"],
                        "excludes_zero": int(stats["excludes_zero"]),
                        "status": "ok",
                        "note": "resamples connected groups, not ligands; shared dual within each resample; current scores; B=2000 seed=20260729",
                        "seed": SEED,
                    }
                )
    return rows


def compute_holdout():
    rows = read_csv(MASTER)
    by = {}
    for r in rows:
        if r["analysis_set"] != "holdout" or r["complete_case"] not in ("1", 1, "True"):
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"] or r["construction_class"]
        by.setdefault(r["pair"], []).append(rec)
    out = []
    for pair, recs in by.items():
        sa, sb, cls = arrays(recs)
        if len(set(cls) & {"dual", "A_only", "B_only"}) < 3:
            continue
        stats = summary_min_stratified(sa, sb, cls, N_BOOT, SEED)
        mm = matched_mismatched_paired(sa, sb, cls, N_BOOT, SEED, stratified=True)
        out.append(
            {
                "pair": pair,
                "n": len(recs),
                "n_dual": stats["n_dual"],
                "n_A_only": stats["n_A_only"],
                "n_B_only": stats["n_B_only"],
                "summary_min": r4(stats["summary_min"]),
                "ci_lo": r4(stats["summary_min_ci_lo"]),
                "ci_hi": r4(stats["summary_min_ci_hi"]),
                "matched_minus_mismatched": r4(mm["delta"]),
                "mm_ci_lo": r4(mm["delta_ci_lo"]),
                "mm_ci_hi": r4(mm["delta_ci_hi"]),
                "mm_ci_excludes_zero": int(mm["ci_excludes_zero"]),
                "bootstrap": "class_stratified",
                "B": N_BOOT,
                "seed": SEED,
            }
        )
    return out


def copy_rmsd():
    src = ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv"
    return read_csv(src) if src.is_file() else []


def gnina_rows(packs):
    """Recompute EGFR/PIK3CA GNINA directional AUROC from deposited poses + master classes."""
    files = {
        "EGFR/HER2": (ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv", "3POZ", "3RCD"),
        "PIK3CA/mTOR": (ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_PIK3CA_mTOR.csv", "4L23", "4JT6"),
    }
    out = []
    from analysis.bootstrap_metrics import stratified_auroc_ci

    for pair, (path, pdb_a, pdb_b) in files.items():
        if not path.is_file() or pair not in packs:
            continue
        cls = {r["ligand_id"]: r["cls"] for r in packs[pair]}
        wide = {}
        for r in read_csv(path):
            lig = r["ligand"]
            sc = -float(r["gnina_mode1"]) if fnum(r.get("gnina_mode1")) is not None else None
            if sc is None:
                continue
            wide.setdefault(lig, {})
            if r["target"] == pdb_a:
                wide[lig]["A"] = sc
            elif r["target"] == pdb_b:
                wide[lig]["B"] = sc
        recs = []
        for lig, sc in wide.items():
            if "A" not in sc or "B" not in sc or lig not in cls:
                continue
            recs.append({"cls": cls[lig], "score_A": sc["A"], "score_B": sc["B"]})
        if len(recs) < 16:
            continue
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        lab = np.array([r["cls"] for r in recs])
        stats = summary_min_stratified(sa, sb, lab, N_BOOT, SEED)
        dual = [r["score_mean"] if "score_mean" in r else (r["score_A"] + r["score_B"]) / 2 for r in recs if r["cls"] == "dual"]
        # mean dual vs neither
        dual_m = [(r["score_A"] + r["score_B"]) / 2 for r in recs if r["cls"] == "dual"]
        nei_m = [(r["score_A"] + r["score_B"]) / 2 for r in recs if r["cls"] == "neither"]
        dn, dn_lo, dn_hi = stratified_auroc_ci(dual_m, nei_m, N_BOOT, SEED)
        out.append(
            {
                "pair": pair,
                "engine": "gnina_dock_mode1",
                "n": len(recs),
                "summary_min": r4(stats["summary_min"]),
                "summary_min_ci_lo": r4(stats["summary_min_ci_lo"]),
                "summary_min_ci_hi": r4(stats["summary_min_ci_hi"]),
                "auroc_D_vs_A_pocketB": r4(stats["auroc_D_vs_A_pocketB"]),
                "auroc_D_vs_B_pocketA": r4(stats["auroc_D_vs_B_pocketA"]),
                "auroc_D_vs_neither_mean": r4(dn),
                "source": str(path.relative_to(ROOT)),
                "note": "another pose-generation realization; not a Vina vs GNINA ranking",
            }
        )
        _ = (dn_lo, dn_hi, dual)
    return out


def write_two_pocket_ranking(rank_rows, fixed_rows):
    out = []
    mean = {r["pair"]: r for r in fixed_rows if r["contrast"] == "D_vs_neither_two_pocket_mean"}
    for r in rank_rows:
        m = mean.get(r["pair"], {})
        out.append(
            {
                "pair": r["pair"],
                "n_ranked": r["n_ranked"],
                "two_pocket_mean_D_vs_neither": m.get("auroc_dual_vs_neither", ""),
                "ef_dual_10pct": r["ef_dual_10pct"],
                "top_k": r["top_k"],
                "top_dual": r["top_dual"],
                "ranking_rule": r["ranking_rule"],
                "tie_rule": r["tie_rule"],
            }
        )
    return out


def nonstrat_sensitivity(smin_rows, matched_rows, fixed_sens):
    rows = list(fixed_sens)
    for r in smin_rows:
        q1 = not ci_crosses(fnum(r["ci_lo"]), fnum(r["ci_hi"]), 0.5)
        q2 = not ci_crosses(fnum(r["non_stratified_ci_lo"]), fnum(r["non_stratified_ci_hi"]), 0.5)
        rows.append(
            {
                "pair": r["pair"],
                "estimand": "summary_min",
                "primary_stratified_point": r["summary_min"],
                "stratified_ci_lo": r["ci_lo"],
                "stratified_ci_hi": r["ci_hi"],
                "non_stratified_ci_lo": r["non_stratified_ci_lo"],
                "non_stratified_ci_hi": r["non_stratified_ci_hi"],
                "qualitative_reference": "CI excludes 0.5",
                "qualitative_conclusion_changed": "yes" if q1 != q2 else "no",
            }
        )
    for r in matched_rows:
        q1 = bool(int(r["ci_excludes_zero"]))
        q2 = bool(int(r["non_stratified_ci_excludes_zero"]))
        rows.append(
            {
                "pair": r["pair"],
                "estimand": "matched_minus_mismatched",
                "primary_stratified_point": r["delta"],
                "stratified_ci_lo": r["delta_ci_lo"],
                "stratified_ci_hi": r["delta_ci_hi"],
                "non_stratified_ci_lo": r["non_stratified_delta_ci_lo"],
                "non_stratified_ci_hi": r["non_stratified_delta_ci_hi"],
                "qualitative_reference": "CI excludes 0",
                "qualitative_conclusion_changed": "yes" if q1 != q2 else "no",
            }
        )
    return rows


def main() -> int:
    print("loading master")
    packs = load_main()
    for pair, recs in packs.items():
        print(f"  {pair}: {len(recs)}")
    directional, smin = compute_directional(packs)
    print("directional done")
    fixed, fixed_sens = compute_fixed_score(packs)
    print("fixed-score done")
    ranking = compute_ranking(packs)
    print("ranking done")
    matched = compute_matched(packs)
    print("matched/mismatched done")
    labels = compute_label_sensitivity(packs)
    print("label sensitivity done")
    maxmed = compute_max_median(packs)
    print("max/median done")
    cluster = compute_cluster(packs)
    print("cluster done")
    holdout = compute_holdout()
    print("holdout done")
    gnina = gnina_rows(packs)
    print("gnina done")
    rmsd = copy_rmsd()
    two_pocket = write_two_pocket_ranking(ranking, fixed)
    sens = nonstrat_sensitivity(smin, matched, fixed_sens)

    write_csv(CANON / "primary_directional_auroc.csv", directional)
    write_csv(CANON / "primary_summary_min.csv", smin)
    write_csv(CANON / "fixed_score_negative_class_delta.csv", [r for r in fixed if r["contrast"] != "D_vs_neither_two_pocket_mean"])
    write_csv(CANON / "two_pocket_mean_ranking.csv", two_pocket)
    write_csv(CANON / "top10_operating_points.csv", ranking)
    write_csv(CANON / "matched_mismatched_pocket.csv", matched)
    all_master = read_csv(MASTER)
    cls_rows = []
    for pair, recs in packs.items():
        ct = Counter(r["cls"] for r in recs)
        scored = [
            r
            for r in all_master
            if r["pair"] == pair and r.get("analysis_set") == "main" and r.get("complete_case") in ("1", "True")
        ]
        ineligible = [r for r in scored if str(r.get("activity_eligible", "1")) not in ("1", "True")]
        cls_rows.append(
            {
                "pair": pair,
                "n_dual": ct["dual"],
                "n_A_only": ct["A_only"],
                "n_B_only": ct["B_only"],
                "n_neither": ct["neither"],
                "n_activity_eligible": len(recs),
                "n_scored_complete_case": len(scored),
                "n_activity_ineligible": len(ineligible),
                "n_complete_case": len(recs),
                "n_missing_score_A": 0,
                "n_missing_score_B": 0,
                "n_class_pchembl_mismatch": 0,
            }
        )
    write_csv(CANON / "class_counts.csv", cls_rows)
    write_csv(
        CANON / "matched_minus_mismatched.csv",
        [
            {
                "pair": r["pair"],
                "delta": r["delta"],
                "ci_lo": r["delta_ci_lo"],
                "ci_hi": r["delta_ci_hi"],
                "ci_excludes_zero": r["ci_excludes_zero"],
                "weaker_arm_switched": r["weaker_arm_switched"],
                "n_dual": r.get("n_dual", ""),
                "n_A_only": r.get("n_A_only", ""),
                "n_B_only": r.get("n_B_only", ""),
                "bootstrap": r.get("bootstrap", ""),
                "B": r.get("B", ""),
                "seed": r.get("seed", ""),
            }
            for r in matched
        ],
    )
    write_csv(CANON / "label_aggregation_sensitivity.csv", labels)
    write_csv(CANON / "max_vs_median_sensitivity.csv", maxmed)
    write_csv(CANON / "cluster_bootstrap_sensitivity.csv", cluster)
    write_csv(CANON / "holdout_metrics.csv", holdout)
    write_csv(CANON / "computational_robustness.csv", gnina)
    write_csv(CANON / "cognate_rmsd.csv", rmsd)
    write_csv(CANON / "non_stratified_bootstrap_sensitivity.csv", sens)

    print("\nTable 2 (summary_min, scheme B)")
    for r in smin:
        print(f"  {r['pair']:14} {r['summary_min']} [{r['ci_lo']}, {r['ci_hi']}] weaker={r['weaker_arm']}")
    print("\nEGFR matched-mismatched")
    for r in matched:
        if r["pair"] == "EGFR/HER2":
            print(f"  delta={r['delta']} [{r['delta_ci_lo']}, {r['delta_ci_hi']}] excludes0={r['ci_excludes_zero']} switch={r['weaker_arm_switched']}")
    print("wrote", CANON)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
