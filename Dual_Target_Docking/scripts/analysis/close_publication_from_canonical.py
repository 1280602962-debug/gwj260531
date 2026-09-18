#!/usr/bin/env python3
"""Derive figure/manuscript CSVs from results/canonical and finish remaining estimands.

Zero-dock. No new pairs. Scheme B: class-stratified B=2000 seed=20260729.

Retired. Current Figure 3 / Table S5 numbers come from
`results/canonical/ecfp4_incremental_information.csv` and
`results/canonical/descriptor_baselines.csv`. Historical publication CSVs
are not rewritten.
"""
from __future__ import annotations

import csv
import math
import shutil
import sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import (  # noqa: E402
    N_BOOT,
    SEED,
    auroc,
    matched_mismatched_paired,
    stratified_auroc_ci,
    summary_min_stratified,
)

CANON = ROOT / "results" / "canonical"
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = fields or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fnum(v):
    if v is None or v == "":
        return None
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def r4(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return ""
    return f"{float(x):.4f}"


def r3(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(float(x))):
        return ""
    d = Decimal(str(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return f"{d:.3f}"


def signed3(x) -> str:
    v = float(x)
    s = r3(abs(v))
    if v > 0:
        return "+" + s if s != "0.000" else s
    if v < 0:
        return "−" + s
    return s


def ci3(lo, hi) -> str:
    lo_s = r3(lo)
    hi_s = r3(hi)
    if lo_s.startswith("-"):
        lo_s = "−" + lo_s[1:]
    if hi_s.startswith("-"):
        hi_s = "−" + hi_s[1:]
    return f"[{lo_s}, {hi_s}]"


def load_main():
    packs = {p: [] for p in PAIRS}
    for r in read_csv(CANON / "current_score_master.csv"):
        if r["analysis_set"] != "main" or r["complete_case"] not in {"1", "True", 1}:
            continue
        if str(r.get("activity_eligible", "1")) not in {"1", "True"}:
            continue
        pair = r["pair"]
        if pair not in packs:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["score_mean"] = (rec["score_A"] + rec["score_B"]) / 2
        rec["score_worst"] = min(rec["score_A"], rec["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["ligand_id"] = r["ligand_id"]
        packs[pair].append(rec)
    return packs


def load_holdout():
    packs = {}
    for r in read_csv(CANON / "current_score_master.csv"):
        if r["analysis_set"] != "holdout" or r["complete_case"] not in {"1", "True", 1}:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"] or r["construction_class"]
        rec["ligand_id"] = r["ligand_id"]
        packs.setdefault(r["pair"], []).append(rec)
    return packs


def gnina_formulation(packs):
    files = {
        "EGFR/HER2": (
            ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv",
            "3POZ",
            "3RCD",
        ),
        "PIK3CA/mTOR": (
            ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_PIK3CA_mTOR.csv",
            "4L23",
            "4JT6",
        ),
        "JAK1/TYK2": (
            ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_independent_jak1_tyk2_v1.csv",
            "6N7A",
            "3LXP",
        ),
    }
    egfr_prod = ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv"
    rows = []
    robust = []
    for pair, (path, pdb_a, pdb_b) in files.items():
        if not path.is_file():
            continue
        cls = {r["ligand_id"]: r["cls"] for r in packs[pair]}
        wide = {}
        for r in read_csv(path):
            lig = r.get("ligand") or r.get("ligand_id")
            raw = fnum(r.get("gnina_mode1"))
            if raw is None:
                continue
            sc = -raw if raw < 0 or "score_S" not in r else float(r.get("score_S") or -raw)
            if fnum(r.get("score_S")) is not None:
                sc = float(r["score_S"])
            elif raw is not None:
                sc = -raw
            wide.setdefault(lig, {})
            tgt = r.get("target")
            if tgt == pdb_a:
                wide[lig]["A"] = sc
            elif tgt == pdb_b:
                wide[lig]["B"] = sc
        recs = []
        for lig, sc in wide.items():
            if "A" not in sc or "B" not in sc or lig not in cls:
                continue
            recs.append({"cls": cls[lig], "score_A": sc["A"], "score_B": sc["B"], "ligand": lig})
        if len(recs) < 16:
            continue
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        lab = np.array([r["cls"] for r in recs])
        stats = summary_min_stratified(sa, sb, lab, N_BOOT, SEED)
        dual = [r for r in recs if r["cls"] == "dual"]
        aonly = [r for r in recs if r["cls"] == "A_only"]
        bonly = [r for r in recs if r["cls"] == "B_only"]
        neither = [r for r in recs if r["cls"] == "neither"]
        da, da_lo, da_hi = stratified_auroc_ci(
            [r["score_B"] for r in dual], [r["score_B"] for r in aonly], N_BOOT, SEED
        )
        db, db_lo, db_hi = stratified_auroc_ci(
            [r["score_A"] for r in dual], [r["score_A"] for r in bonly], N_BOOT, SEED
        )
        dn, dn_lo, dn_hi = stratified_auroc_ci(
            [(r["score_A"] + r["score_B"]) / 2 for r in dual],
            [(r["score_A"] + r["score_B"]) / 2 for r in neither],
            N_BOOT,
            SEED,
        )
        src = str(path.relative_to(ROOT))
        rows.extend(
            [
                {
                    "pair": pair,
                    "engine": "gnina_dock_mode1",
                    "formulation": "dualfourclass_directional",
                    "contrast": "D_vs_A_pocketB",
                    "score": "score_B",
                    "n_pos": len(dual),
                    "n_neg": len(aonly),
                    "auroc": r4(da),
                    "ci_lo": r4(da_lo),
                    "ci_hi": r4(da_hi),
                    "note": "another pose-generation realization; not a Vina vs GNINA ranking",
                },
                {
                    "pair": pair,
                    "engine": "gnina_dock_mode1",
                    "formulation": "dualfourclass_directional",
                    "contrast": "D_vs_B_pocketA",
                    "score": "score_A",
                    "n_pos": len(dual),
                    "n_neg": len(bonly),
                    "auroc": r4(db),
                    "ci_lo": r4(db_lo),
                    "ci_hi": r4(db_hi),
                    "note": "another pose-generation realization; not a Vina vs GNINA ranking",
                },
                {
                    "pair": pair,
                    "engine": "gnina_dock_mode1",
                    "formulation": "dualfourclass_directional",
                    "contrast": "summary_min",
                    "score": "min(D/A,D/B)",
                    "n_pos": len(dual),
                    "n_neg": min(len(aonly), len(bonly)),
                    "auroc": r4(stats["summary_min"]),
                    "ci_lo": r4(stats["summary_min_ci_lo"]),
                    "ci_hi": r4(stats["summary_min_ci_hi"]),
                    "note": "another pose-generation realization; not a Vina vs GNINA ranking",
                },
                {
                    "pair": pair,
                    "engine": "gnina_dock_mode1",
                    "formulation": "conventional_dual_vs_neither",
                    "contrast": "D_vs_neither_mean",
                    "score": "score_mean",
                    "n_pos": len(dual),
                    "n_neg": len(neither),
                    "auroc": r4(dn),
                    "ci_lo": r4(dn_lo),
                    "ci_hi": r4(dn_hi),
                    "note": "another pose-generation realization; not a Vina vs GNINA ranking",
                },
            ]
        )
        robust.append(
            {
                "pair": pair,
                "engine": "gnina_dock_mode1",
                "n": len(recs),
                "summary_min": r4(stats["summary_min"]),
                "summary_min_ci_lo": r4(stats["summary_min_ci_lo"]),
                "summary_min_ci_hi": r4(stats["summary_min_ci_hi"]),
                "auroc_D_vs_A_pocketB": r4(da),
                "auroc_D_vs_B_pocketA": r4(db),
                "auroc_D_vs_neither_mean": r4(dn),
                "d_vs_neither_ci_lo": r4(dn_lo),
                "d_vs_neither_ci_hi": r4(dn_hi),
                "n_neither": len(neither),
                "source": src,
                "note": "another pose-generation realization; not a Vina vs GNINA ranking",
            }
        )
        if pair == "EGFR/HER2" and path != egfr_prod:
            egfr_prod.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, egfr_prod)
    return rows, robust


def receptor_sub(packs):
    recs = packs["PIK3CA/mTOR"]
    by_id = {r["ligand_id"]: r for r in recs}
    specs = [
        ("4JPS", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt4JPS.csv", "A"),
        ("5DXT", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt5DXT.csv", "A"),
        ("4JSX", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt4JSX.csv", "B"),
    ]
    out_rows = []
    for pdb, path, pocket in specs:
        alt = {r["ligand"]: -float(r["vina_mode1"]) for r in read_csv(path) if r.get("status") == "success"}
        labeled = []
        for lig, base in by_id.items():
            if lig not in alt:
                continue
            sa = alt[lig] if pocket == "A" else base["score_A"]
            sb = alt[lig] if pocket == "B" else base["score_B"]
            labeled.append({"cls": base["cls"], "score_A": sa, "score_B": sb})
        sa = np.array([r["score_A"] for r in labeled])
        sb = np.array([r["score_B"] for r in labeled])
        lab = np.array([r["cls"] for r in labeled])
        stats = summary_min_stratified(sa, sb, lab, N_BOOT, SEED)
        ct = Counter(lab)
        row = {
            "replacement": pdb,
            "pocket_replaced": pocket,
            "n": len(labeled),
            "n_dual": ct["dual"],
            "n_A_only": ct["A_only"],
            "n_B_only": ct["B_only"],
            "auroc_D_vs_A": r4(stats["auroc_D_vs_A_pocketB"]),
            "auroc_D_vs_B": r4(stats["auroc_D_vs_B_pocketA"]),
            "summary_min": r4(stats["summary_min"]),
            "summary_min_ci_lo": r4(stats["summary_min_ci_lo"]),
            "summary_min_ci_hi": r4(stats["summary_min_ci_hi"]),
            "primary_summary_min": "0.6921",
            "note": f"pocket {pocket} replaced by alt receptor; other pocket from current score master; scheme B",
        }
        out_rows.append(row)
        dest = ROOT / f"data/jcim_structure_robust_v0/tables/pocket_matched_PM48_alt{pdb}_v1.csv"
        if pocket == "A":
            write_csv(
                dest,
                [
                    {
                        "receptor": pdb,
                        "n": row["n"],
                        "n_dual": row["n_dual"],
                        "n_A_only": row["n_A_only"],
                        "n_B_only": row["n_B_only"],
                        "auroc_D_vs_A": row["auroc_D_vs_A"],
                        "auroc_D_vs_B": row["auroc_D_vs_B"],
                        "summary_min": row["summary_min"],
                        "summary_min_ci_lo": row["summary_min_ci_lo"],
                        "summary_min_ci_hi": row["summary_min_ci_hi"],
                        "primary_summary_min": row["primary_summary_min"],
                        "note": row["note"],
                    }
                ],
            )
        else:
            write_csv(
                dest,
                [
                    {
                        "receptor": pdb,
                        "pocket": "B",
                        "n": row["n"],
                        "n_dual": row["n_dual"],
                        "n_A_only": row["n_A_only"],
                        "n_B_only": row["n_B_only"],
                        "summary_min": row["summary_min"],
                        "auroc_D_vs_B": row["auroc_D_vs_B"],
                        "summary_min_dup": row["summary_min"],
                        "summary_min_ci_lo": row["summary_min_ci_lo"],
                        "summary_min_ci_hi": row["summary_min_ci_hi"],
                        "primary_summary_min": row["primary_summary_min"],
                        "note": row["note"],
                    }
                ],
            )
    write_csv(CANON / "receptor_substitution.csv", out_rows)
    return out_rows


def two_pocket(packs, smin, ranking):
    smin_d = {r["pair"]: r for r in smin}
    rank_d = {r["pair"]: r for r in ranking}
    form = []
    two = []
    for pair, recs in packs.items():
        dual = [r for r in recs if r["cls"] == "dual"]
        neither = [r for r in recs if r["cls"] == "neither"]
        nondual = [r for r in recs if r["cls"] != "dual"]
        pt, lo, hi = stratified_auroc_ci(
            [r["score_mean"] for r in dual], [r["score_mean"] for r in neither], N_BOOT, SEED
        )
        all_pt, all_lo, all_hi = stratified_auroc_ci(
            [r["score_mean"] for r in dual], [r["score_mean"] for r in nondual], N_BOOT, SEED
        )
        worst_pt, wlo, whi = stratified_auroc_ci(
            [r["score_worst"] for r in dual], [r["score_worst"] for r in neither], N_BOOT, SEED
        )
        s = smin_d[pair]
        rk = rank_d[pair]
        two.append(
            {
                "pair": pair,
                "n_ranked": rk["n_ranked"],
                "two_pocket_mean_D_vs_neither": r4(pt),
                "ci_lo": r4(lo),
                "ci_hi": r4(hi),
                "ef_dual_10pct": rk["ef_dual_10pct"],
                "top_k": rk["top_k"],
                "top_dual": rk["top_dual"],
                "ranking_rule": rk["ranking_rule"],
                "tie_rule": rk["tie_rule"],
            }
        )
        form.extend(
            [
                {
                    "pair": pair,
                    "formulation": "dualfourclass_directional",
                    "contrast": "D_vs_A_pocketB",
                    "score": "vina_B",
                    "n_pos": s["n_dual"],
                    "n_neg": s["n_A_only"],
                    "auroc": s["auroc_D_vs_A_pocketB"],
                    "ci_lo": "",
                    "ci_hi": "",
                    "underpowered": 0,
                    "note": "primary: dual vs A-only scored in pocket B; CI in primary_directional_auroc.csv",
                },
                {
                    "pair": pair,
                    "formulation": "dualfourclass_directional",
                    "contrast": "D_vs_B_pocketA",
                    "score": "vina_A",
                    "n_pos": s["n_dual"],
                    "n_neg": s["n_B_only"],
                    "auroc": s["auroc_D_vs_B_pocketA"],
                    "ci_lo": "",
                    "ci_hi": "",
                    "underpowered": 0,
                    "note": "primary: dual vs B-only scored in pocket A; CI in primary_directional_auroc.csv",
                },
                {
                    "pair": pair,
                    "formulation": "conventional_dual_vs_neither",
                    "contrast": "D_vs_neither_mean",
                    "score": "vina_mean",
                    "n_pos": len(dual),
                    "n_neg": len(neither),
                    "auroc": r4(pt),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "underpowered": int(len(neither) < 10),
                    "note": "two-pocket mean dual vs neither; scheme-B class-stratified",
                },
                {
                    "pair": pair,
                    "formulation": "conventional_dual_vs_neither",
                    "contrast": "D_vs_neither_worst",
                    "score": "vina_worst",
                    "n_pos": len(dual),
                    "n_neg": len(neither),
                    "auroc": r4(worst_pt),
                    "ci_lo": r4(wlo),
                    "ci_hi": r4(whi),
                    "underpowered": int(len(neither) < 10),
                    "note": "AND-like dual vs neither using min(pocket scores)",
                },
                {
                    "pair": pair,
                    "formulation": "dual_vs_all_nondual",
                    "contrast": "D_vs_A+B+neither_mean",
                    "score": "vina_mean",
                    "n_pos": len(dual),
                    "n_neg": len(nondual),
                    "auroc": r4(all_pt),
                    "ci_lo": r4(all_lo),
                    "ci_hi": r4(all_hi),
                    "underpowered": 0,
                    "note": "selectives counted as negatives; still not directional",
                },
            ]
        )
    write_csv(CANON / "two_pocket_mean_ranking.csv", two)
    write_csv(ROOT / "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv", form)
    return two, form


def and_filter(packs, ranking):
    rows = []
    for pair, recs in packs.items():
        dual_w = [r["score_worst"] for r in recs if r["cls"] == "dual"]
        thresh = float(np.median(dual_w))
        usable = [r for r in recs if r["cls"] in {"dual", "A_only", "B_only"}]
        kept = [r for r in usable if r["score_worst"] >= thresh]
        ct = Counter(r["cls"] for r in kept)
        n_dual = sum(1 for r in recs if r["cls"] == "dual")
        rk = next(r for r in ranking if r["pair"] == pair)
        rows.append(
            {
                **rk,
                "and_threshold": r4(thresh),
                "and_n_input": len(usable),
                "and_n_pass": len(kept),
                "and_dual": ct["dual"],
                "and_A_only": ct["A_only"],
                "and_B_only": ct["B_only"],
                "and_dual_recall": r4(ct["dual"] / n_dual) if n_dual else "",
                "and_dual_precision": r4(ct["dual"] / len(kept)) if kept else "",
                "filter_rule": "score_worst >= median_dual_score_worst; neither excluded",
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv", rows)
    write_csv(
        ROOT / "data/jcim_novelty_v0/tables/operating_point_examples_review_v1.csv",
        [
            {
                "pair": r["pair"],
                "n_ranked": r["n_ranked"],
                "top_k": r["top_k"],
                "top_dual": r["top_dual"],
                "top_A_only": r["top_A_only"],
                "top_B_only": r["top_B_only"],
                "top_neither": r["top_neither"],
                "top_ligand_ids": r["top_ligand_ids"],
                "threshold": r["and_threshold"],
                "n_filter_input": r["and_n_input"],
                "retained_dual": r["and_dual"],
                "retained_A_only": r["and_A_only"],
                "retained_B_only": r["and_B_only"],
                "dual_recall": r["and_dual_recall"],
                "dual_precision": r["and_dual_precision"],
                "filter_rule": r["filter_rule"],
                "tie_rule": r["tie_rule"],
            }
            for r in rows
        ],
    )
    return rows


def descriptors(packs, desc):
    out = []
    desc_d = {r["pair"]: r for r in desc}
    for pair, recs in packs.items():
        d = desc_d[pair]
        mw_da = mw_db = None
        dual = [r for r in recs if r["cls"] == "dual"]
        aonly = [r for r in recs if r["cls"] == "A_only"]
        bonly = [r for r in recs if r["cls"] == "B_only"]

        def mw_of(r):
            mol = Chem.MolFromSmiles(r.get("smiles") or "")
            return Descriptors.MolWt(mol) if mol is not None else float("nan")

        if dual:
            mw_da = auroc([mw_of(r) for r in dual], [mw_of(r) for r in aonly])
            mw_db = auroc([mw_of(r) for r in dual], [mw_of(r) for r in bonly])
        mw_smin = min(mw_da, mw_db) if mw_da is not None else float("nan")
        cands = {
            "heavy": float(d["heavy_summary_min"]),
            "mw": mw_smin,
            "clogp": float(d["clogp_summary_min"]),
            "tpsa": float(d["tpsa_summary_min"]),
        }
        best_name, best_val = max(cands.items(), key=lambda kv: kv[1])
        out.append(
            {
                "pair": pair,
                "n_dual": d["n_dual"],
                "n_A_only": d["n_A_only"],
                "n_B_only": d["n_B_only"],
                "heavy_D_vs_A": d["heavy_D_vs_A"],
                "heavy_D_vs_B": d["heavy_D_vs_B"],
                "heavy_summary_min": d["heavy_summary_min"],
                "mw_D_vs_A": r4(mw_da),
                "mw_D_vs_B": r4(mw_db),
                "mw_summary_min": r4(mw_smin),
                "clogp_D_vs_A": d["clogp_D_vs_A"],
                "clogp_D_vs_B": d["clogp_D_vs_B"],
                "clogp_summary_min": d["clogp_summary_min"],
                "tpsa_D_vs_A": d["tpsa_D_vs_A"],
                "tpsa_D_vs_B": d["tpsa_D_vs_B"],
                "tpsa_summary_min": d["tpsa_summary_min"],
                "best_single_descriptor": best_name,
                "best_single_descriptor_summary_min": r4(best_val),
                "note": "best_single_descriptor from current score master; scheme-B Vina comparison in descriptor_baselines.csv",
            }
        )
    return out


def export_core(packs, hold):
    smin = read_csv(CANON / "primary_summary_min.csv")
    direc = read_csv(CANON / "primary_directional_auroc.csv")
    labels = read_csv(CANON / "label_aggregation_sensitivity.csv")
    ranking = read_csv(CANON / "top10_operating_points.csv")
    matched = read_csv(CANON / "matched_mismatched_pocket.csv")
    fixed = read_csv(CANON / "fixed_score_negative_class_delta.csv")
    inc = read_csv(CANON / "ecfp4_incremental_information.csv")
    desc = read_csv(CANON / "descriptor_baselines.csv")
    cluster = read_csv(CANON / "cluster_bootstrap_sensitivity.csv")
    rmsd = read_csv(CANON / "cognate_rmsd.csv")
    da = {(r["pair"], r["estimand"]): r for r in direc}

    unified = []
    for r in labels:
        unified.append(
            {
                "pair": r["pair"],
                "label_rule": r["label_rule"],
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "n_neither_excluded": "",
                "underpowered": int(min(int(r["n_dual"]), int(r["n_A_only"]), int(r["n_B_only"])) < 8),
                "pocket_matched_summary_min": r["summary_min"],
                "auroc_D_vs_A": r["auroc_D_vs_A"],
                "auroc_D_vs_B": r["auroc_D_vs_B"],
                "ci_lo": r["ci_lo"],
                "ci_hi": r["ci_hi"],
                "bootstrap": r["bootstrap"],
                "seed": r["seed"],
            }
        )
    write_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv", unified)

    two, form = two_pocket(packs, smin, ranking)
    form_nei = {(r["pair"], r["contrast"]): r for r in form}
    rank_and = and_filter(packs, ranking)
    descriptors(packs, desc)

    five = []
    for pair in ("F2/F10", "JAK1/TYK2", "JAK1/JAK2", "PPARG/PPARA", "PPARA/PPARD"):
        s = next(r for r in smin if r["pair"] == pair)
        d = next(r for r in desc if r["pair"] == pair)
        n = form_nei[(pair, "D_vs_neither_mean")]
        five.append(
            {
                "pair": pair,
                "label_rule": "theta_6.0",
                "bootstrap": "class_stratified_shared_dual",
                "n_dual": s["n_dual"],
                "n_A_only": s["n_A_only"],
                "n_B_only": s["n_B_only"],
                "n_neither": n["n_neg"],
                "n_scored_both_ends": next(r for r in ranking if r["pair"] == pair)["n_ranked"],
                "auroc_D_vs_A_pocketB": s["auroc_D_vs_A_pocketB"],
                "auroc_D_vs_B_pocketA": s["auroc_D_vs_B_pocketA"],
                "summary_min": s["summary_min"],
                "ci_lo": s["ci_lo"],
                "ci_hi": s["ci_hi"],
                "n_boot_ok": 2000,
                "D_vs_neither_vina_mean": n["auroc"],
                "D_vs_neither_ci_lo": n["ci_lo"],
                "D_vs_neither_ci_hi": n["ci_hi"],
                "best_single_descriptor": d["best_descriptor"],
                "best_single_descriptor_summary_min": d["best_descriptor_summary_min"],
                "note": "Derived from results/canonical; scheme-B primary",
            }
        )
    write_csv(
        ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv",
        five,
    )

    equal = []
    for r in fixed:
        equal.append(
            {
                **r,
                "n_dual": r["n_dual"],
                "n_selective": r["n_selective"],
                "n_neither": r["n_neither"],
                "underpowered_neither": int(int(r["n_neither"] or 0) < 10),
                "note": "scheme-B class-stratified shared dual; results/canonical",
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv", equal)
    write_csv(
        ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv",
        [r for r in equal if r["pair"] in {"F2/F10", "JAK1/TYK2", "JAK1/JAK2", "PPARG/PPARA", "PPARA/PPARD"}],
    )

    wp = []
    for r in matched:
        wp.append(
            {
                "pair": r["pair"],
                "set": r["set"],
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "n_boot_ok": 2000,
                "matched_D_vs_A": r["matched_auroc_D_vs_A"],
                "matched_D_vs_B": r["matched_auroc_D_vs_B"],
                "matched_summary_min": r["matched_summary_min"],
                "matched_ci_lo": r.get("matched_summary_min", ""),
                "matched_ci_hi": "",
                "wrong_D_vs_A": r["mismatched_auroc_D_vs_A"],
                "wrong_D_vs_B": r["mismatched_auroc_D_vs_B"],
                "wrong_summary_min": r["mismatched_summary_min"],
                "delta_matched_minus_wrong": r["delta"],
                "delta_ci_lo": r["delta_ci_lo"],
                "delta_ci_hi": r["delta_ci_hi"],
                "ci_excludes_zero": "True" if str(r["ci_excludes_zero"]) in {"1", "True"} else "False",
                "weaker_arm_switched": r["weaker_arm_switched"],
            }
        )
    five_wp = []
    for pair, recs in hold.items():
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        cls = np.array([r["cls"] for r in recs])
        stats = summary_min_stratified(sa, sb, cls, N_BOOT, SEED)
        mm = matched_mismatched_paired(sa, sb, cls, N_BOOT, SEED, stratified=True)
        wp.append(
            {
                "pair": pair,
                "set": "unused_pool_holdout",
                "n_dual": stats["n_dual"],
                "n_A_only": stats["n_A_only"],
                "n_B_only": stats["n_B_only"],
                "n_boot_ok": 2000,
                "matched_D_vs_A": r4(stats["auroc_D_vs_A_pocketB"]),
                "matched_D_vs_B": r4(stats["auroc_D_vs_B_pocketA"]),
                "matched_summary_min": r4(stats["summary_min"]),
                "matched_ci_lo": r4(stats["summary_min_ci_lo"]),
                "matched_ci_hi": r4(stats["summary_min_ci_hi"]),
                "wrong_D_vs_A": r4(mm["mismatched_auroc_D_vs_A"]),
                "wrong_D_vs_B": r4(mm["mismatched_auroc_D_vs_B"]),
                "wrong_summary_min": r4(mm["mismatched_summary_min"]),
                "delta_matched_minus_wrong": r4(mm["delta"]),
                "delta_ci_lo": r4(mm["delta_ci_lo"]),
                "delta_ci_hi": r4(mm["delta_ci_hi"]),
                "ci_excludes_zero": "True" if mm["ci_excludes_zero"] else "False",
                "weaker_arm_switched": int(mm["weaker_arm_switched"]),
            }
        )
        five_wp.append(
            {
                "channel": "holdout_vina_20260727",
                "pair": pair,
                "n_dual": stats["n_dual"],
                "n_A_only": stats["n_A_only"],
                "n_B_only": stats["n_B_only"],
                "n_boot_ok": 2000,
                "matched_D_vs_A": r4(stats["auroc_D_vs_A_pocketB"]),
                "matched_D_vs_B": r4(stats["auroc_D_vs_B_pocketA"]),
                "matched_summary_min": r4(stats["summary_min"]),
                "matched_ci_lo": r4(stats["summary_min_ci_lo"]),
                "matched_ci_hi": r4(stats["summary_min_ci_hi"]),
                "wrong_D_vs_A": r4(mm["mismatched_auroc_D_vs_A"]),
                "wrong_D_vs_B": r4(mm["mismatched_auroc_D_vs_B"]),
                "wrong_summary_min": r4(mm["mismatched_summary_min"]),
                "wrong_ci_lo": "",
                "wrong_ci_hi": "",
                "delta_matched_minus_wrong": r4(mm["delta"]),
                "delta_ci_lo": r4(mm["delta_ci_lo"]),
                "delta_ci_hi": r4(mm["delta_ci_hi"]),
                "ci_excludes_zero": "True" if mm["ci_excludes_zero"] else "False",
            }
        )
    for r in matched:
        five_wp.append(
            {
                "channel": "vina_20260727",
                "pair": r["pair"],
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "n_boot_ok": 2000,
                "matched_D_vs_A": r["matched_auroc_D_vs_A"],
                "matched_D_vs_B": r["matched_auroc_D_vs_B"],
                "matched_summary_min": r["matched_summary_min"],
                "matched_ci_lo": "",
                "matched_ci_hi": "",
                "wrong_D_vs_A": r["mismatched_auroc_D_vs_A"],
                "wrong_D_vs_B": r["mismatched_auroc_D_vs_B"],
                "wrong_summary_min": r["mismatched_summary_min"],
                "wrong_ci_lo": "",
                "wrong_ci_hi": "",
                "delta_matched_minus_wrong": r["delta"],
                "delta_ci_lo": r["delta_ci_lo"],
                "delta_ci_hi": r["delta_ci_hi"],
                "ci_excludes_zero": "True" if str(r["ci_excludes_zero"]) in {"1", "True"} else "False",
            }
        )
    write_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv", wp)
    write_csv(
        ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv",
        five_wp,
    )

    write_csv(ROOT / "data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv", cluster)
    write_csv(ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv", rmsd)

    ml = []
    inc_out = []
    five_ecfp = []
    for r in inc:
        for model, col in (
            ("ECFP4", "cv_auroc_ECFP4"),
            ("docking", "cv_auroc_docking"),
            ("ECFP4+docking", "cv_auroc_ECFP4_docking"),
        ):
            rec = {
                "pair": r["pair"],
                "contrast": r["contrast"],
                "model": model,
                "n": r["n"],
                "n_pos": r["n_pos"],
                "n_neg": r["n_neg"],
                "n_scaffolds": r["n_scaffolds"],
                "n_splits": 5,
                "cv_auroc": r[col],
                "rank_auroc_docking": r["rank_auroc_docking"],
                "delta_ECFP4_plus_docking_minus_ECFP4": r["delta_ECFP4_plus_docking_minus_ECFP4"],
                "note": r["note"],
            }
            inc_out.append(rec)
            five_ecfp.append(rec)
        ml.append(
            {
                "pair": r["pair"],
                "contrast": r["contrast"],
                "method": "ECFP4_logistic_cv",
                "cv_scheme": "scaffold_GroupKFold",
                "n_splits": 5,
                "n_scaffolds": r["n_scaffolds"],
                "n": r["n"],
                "auroc_ml": r["cv_auroc_ECFP4"],
                "auroc_dock_pocket_matched": r["rank_auroc_docking"],
                "delta_ml_minus_dock": r4(float(r["cv_auroc_ECFP4"]) - float(r["rank_auroc_docking"])),
                "note": "primary_no_scaffold_overlap_across_folds; from results/canonical",
            }
        )

    master = read_csv(CANON / "current_score_master.csv")
    memb = []
    for r in master:
        if r["analysis_set"] != "main" or r["complete_case"] not in {"1", "True"}:
            continue
        if str(r.get("activity_eligible", "1")) not in {"1", "True"}:
            continue
        memb.append(
            {
                "pair": r["pair"],
                "ligand": r["ligand_id"],
                "cls": r["primary_class_theta6"],
                "pA": r["pA"],
                "pB": r["pB"],
                "score_A": r["score_A"],
                "score_B": r["score_B"],
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv", memb)
    dest = ROOT / "data/processed/current_score_master.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(CANON / "current_score_master.csv", dest)

    hold_out = []
    hold_metrics = {r["pair"]: r for r in read_csv(CANON / "holdout_metrics.csv")}
    for pair, r in hold_metrics.items():
        hold_out.append(
            {
                "prefix": "",
                "pair": pair,
                "variant": "pocket_matched_vina",
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "auroc_D_vs_A": "",
                "auroc_D_vs_B": "",
                "summary_min": r["summary_min"],
                "summary_min_ci_lo": r["ci_lo"],
                "summary_min_ci_hi": r["ci_hi"],
                "main_panel_pocket_matched_vina": next(x["summary_min"] for x in smin if x["pair"] == pair),
                "delta_vs_main_panel": "",
                "note": "scheme-B class-stratified from current score master holdout rows",
            }
        )
    write_csv(ROOT / "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv", hold_out)
    return smin, direc, two, rank_and, matched, desc, da


def overlay_fiveseed():
    src = ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv"
    rows = read_csv(src)
    corr = ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_EGFR_corrected.csv"
    if not corr.is_file():
        corr = ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_corrected.csv"
    if corr.is_file():
        by = {(r["pair"], r["seed"]): r for r in read_csv(corr)}
        out = []
        for r in rows:
            k = (r["pair"], r["seed"])
            out.append(by[k] if k in by else r)
        write_csv(src, out)
        rows = out
    ache = ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_ACHE_corrected.csv"
    if ache.is_file():
        by = {(r["pair"], r["seed"]): r for r in read_csv(ache)}
        out = []
        for r in rows:
            k = (r["pair"], r["seed"])
            out.append(by[k] if k in by else r)
        write_csv(src, out)


def main() -> int:
    print("close_publication_from_canonical.py is retired.")
    print("Figures and tables read results/canonical. Historical publication CSVs are not rewritten.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
