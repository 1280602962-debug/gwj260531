#!/usr/bin/env python3
"""Zero-dock analysis stack for the five post-census pairs.

Copies the all-pairs items already reported for every then-primary pair in
Methods 2.4–2.6. Destination identity: PROJECT_IDENTITY_LOCK_V1.md
(8-row main table after withdrawing PIK3CA/PIK3CB). Does not restock Table 2
or retitle the article.

Table-2-comparable summary_min CIs resample the dual+A-only+B-only ligand
pool without class stratification (same estimand as
unified_threshold_sensitivity_v2.csv). Existing track_b_directional_auroc_v1
CIs are class-preserving and are not reused here as Table 2 intervals.

Dump-gated items (max→median, document year, document-blocked CV,
document-cluster bootstrap, leftover holdout IDs, BindingDB recount) are
recorded as blocked when no ChEMBL sqlite / BindingDB cache is present.
Independent GNINA, five-seed Vina, RTM, and GNINA CNN are local recompute
and are not run here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DT = ROOT.parent.parent
LOCAL = ROOT / "local_track_b_v0"
TAB = ROOT / "tables"
OUT = LOCAL / "tables" / "five_pair_stack_v1"
AN = LOCAL / "analysis"
N_BOOT = 2000
SEED = 20260729
N_MC = 1000
TRUE_AUCS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)
CALIPERS = (0.5, 1.0)

PAIRS = [
    {
        "pair": "F2/F10",
        "system": "coagulation",
        "panel": TAB / "track_b_panels" / "panel_F2_F10_v1.csv",
        "target_a": "4UDW",
        "target_b": "2JKH",
    },
    {
        "pair": "JAK1/TYK2",
        "system": "JAK",
        "panel": TAB / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "3LXP",
    },
    {
        "pair": "JAK1/JAK2",
        "system": "JAK",
        "panel": TAB / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "8BXH",
    },
    {
        "pair": "PPARG/PPARA",
        "system": "PPAR",
        "panel": TAB / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "target_a": "9V8H",
        "target_b": "6LXA",
    },
    {
        "pair": "PPARA/PPARD",
        "system": "PPAR",
        "panel": TAB / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "target_a": "6LXA",
        "target_b": "5U3Q",
    },
]


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus


def r4(x):
    if x is None or x != x:
        return ""
    return round(float(x), 4)


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def assign_fourclass(pA, pB, cut: float):
    if pA is None or pB is None:
        return None
    a, b = pA >= cut, pB >= cut
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def assign_strict(pA, pB):
    if pA is None or pB is None:
        return None
    if pA >= 6.5 and pB >= 6.5:
        return "dual"
    if pA >= 6.5 and pB <= 5.5:
        return "A_only"
    if pB >= 6.5 and pA <= 5.5:
        return "B_only"
    if pA <= 5.5 and pB <= 5.5:
        return "neither"
    return "gray"


def largest_fragment(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    return max(frags, key=lambda m: m.GetNumHeavyAtoms()) if len(frags) > 1 else frags[0]


def directional(recs, key_da="vina_B", key_db="vina_A"):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r[key_da] for r in D], [r[key_da] for r in A])
    db = auroc([r[key_db] for r in D], [r[key_db] for r in B])
    return da, db, min(da, db), len(D), len(A), len(B)


def boot_pm_ci(recs, key_da="vina_B", key_db="vina_A", n_boot=N_BOOT, seed=SEED):
    """Ligand-level non-stratified bootstrap of summary_min (Table 2 estimand)."""
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    if len(usable) < 8:
        return float("nan"), float("nan"), 0
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        sub = [usable[int(i)] for i in ii]
        _, _, mn, *_ = directional(sub, key_da, key_db)
        if mn == mn:
            mins.append(mn)
    if len(mins) < n_boot // 2:
        return float("nan"), float("nan"), len(mins)
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi), len(mins)


def boot_single(pos, neg, n_boot=N_BOOT, seed=SEED):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = [auroc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)) for _ in range(n_boot)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc(pos, neg)), float(lo), float(hi)


def boot_paired(recs, key_m_da, key_m_db, key_w_da, key_w_db, seed=SEED):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    da_m, db_m, sm_m, nD, nA, nB = directional(usable, key_m_da, key_m_db)
    da_w, db_w, sm_w, *_ = directional(usable, key_w_da, key_w_db)
    point = float(r4(sm_m) - r4(sm_w)) if sm_m == sm_m and sm_w == sm_w else float("nan")
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    deltas, sms, sws = [], [], []
    for _ in range(N_BOOT):
        sub = [usable[int(i)] for i in rng.choice(idx, size=len(idx), replace=True)]
        _, _, m, *_ = directional(sub, key_m_da, key_m_db)
        _, _, w, *_ = directional(sub, key_w_da, key_w_db)
        if m != m or w != w:
            continue
        sms.append(m)
        sws.append(w)
        deltas.append(m - w)
    if len(deltas) < N_BOOT // 2:
        dlo = dhi = mlo = mhi = wlo = whi = float("nan")
    else:
        dlo, dhi = np.percentile(deltas, [2.5, 97.5])
        mlo, mhi = np.percentile(sms, [2.5, 97.5])
        wlo, whi = np.percentile(sws, [2.5, 97.5])
    return {
        "n_dual": nD,
        "n_A_only": nA,
        "n_B_only": nB,
        "n_boot_ok": len(deltas),
        "matched_D_vs_A": r4(da_m),
        "matched_D_vs_B": r4(db_m),
        "matched_summary_min": r4(sm_m),
        "matched_ci_lo": r4(mlo),
        "matched_ci_hi": r4(mhi),
        "wrong_D_vs_A": r4(da_w),
        "wrong_D_vs_B": r4(db_w),
        "wrong_summary_min": r4(sm_w),
        "wrong_ci_lo": r4(wlo),
        "wrong_ci_hi": r4(whi),
        "delta_matched_minus_wrong": r4(point),
        "delta_ci_lo": r4(dlo),
        "delta_ci_hi": r4(dhi),
        "ci_excludes_zero": bool(dlo == dlo and (dlo > 0 or dhi < 0)),
    }


def load_scores():
    rows = list(csv.DictReader((LOCAL / "tables" / "scores_vina_mode1_v1.csv").open()))
    out = {}
    for r in rows:
        out.setdefault(r["pair"], {}).setdefault(r["ligand"], {})[r["target"]] = float(r["score_S"])
    return out


def load_job_status():
    rows = list(csv.DictReader((LOCAL / "tables" / "job_status.csv").open()))
    by = {}
    for r in rows:
        by.setdefault(r["pair"], {}).setdefault(r["ligand"], []).append(r)
    return by


def assemble(spec, scores, jobs, label_col="theta6_class"):
    panel = list(csv.DictReader(spec["panel"].open()))
    sc = scores.get(spec["pair"], {})
    recs, incomplete = [], []
    for row in panel:
        lig = row["panel_id"]
        mol = largest_fragment(row["canonical_smiles"])
        if mol is None:
            continue
        ha = mol.GetNumHeavyAtoms()
        pA = float(row["pchembl_A"])
        pB = float(row["pchembl_B"])
        m = sc.get(lig, {})
        sa = m.get(spec["target_a"])
        sb = m.get(spec["target_b"])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        except Exception:
            scaf = ""
        rec = {
            "pair": spec["pair"],
            "ligand": lig,
            "chembl_id": row["molecule_chembl_id"],
            "smiles": row["canonical_smiles"],
            "cls": row[label_col],
            "panel_class": row["class"],
            "pA": pA,
            "pB": pB,
            "vina_A": sa,
            "vina_B": sb,
            "heavy": float(ha),
            "mw": float(Descriptors.MolWt(mol)),
            "clogp": float(Crippen.MolLogP(mol)),
            "tpsa": float(Descriptors.TPSA(mol)),
            "charge": float(sum(a.GetFormalCharge() for a in mol.GetAtoms())),
            "rotatable": float(Lipinski.NumRotatableBonds(mol)),
            "scaffold": scaf or f"NONE:{lig}",
            "fp": fp,
        }
        if sa is not None and sb is not None:
            rec["vina_mean"] = (sa + sb) / 2.0
            rec["vina_worst"] = min(sa, sb)
            rec["le_A"] = sa / ha
            rec["le_B"] = sb / ha
            recs.append(rec)
        else:
            reasons = []
            for j in jobs.get(spec["pair"], {}).get(lig, []):
                if j["status"] != "success":
                    reasons.append(f"{j['target']}:{j['status']}:{j.get('reason', '')}")
            incomplete.append({**rec, "status_reason": ";".join(reasons)})
    return recs, incomplete


def table2_comparable(packs):
    rows = []
    for pair, recs in packs.items():
        da, db, sm, nD, nA, nB = directional(recs)
        lo, hi, n_ok = boot_pm_ci(recs, seed=SEED + stable_offset(pair, "theta_6.0"))
        N = [r for r in recs if r["cls"] == "neither"]
        D = [r for r in recs if r["cls"] == "dual"]
        dn, dn_lo, dn_hi = boot_single(
            [r["vina_mean"] for r in D],
            [r["vina_mean"] for r in N],
            seed=SEED + stable_offset(pair, "D_vs_neither"),
        )
        nondual = [r for r in recs if r["cls"] != "dual"]
        dall, dall_lo, dall_hi = boot_single(
            [r["vina_mean"] for r in D],
            [r["vina_mean"] for r in nondual],
            seed=SEED + stable_offset(pair, "D_vs_all_nondual"),
        )
        desc_mins = {}
        for d in ("heavy", "mw", "clogp", "tpsa"):
            dda = auroc([r[d] for r in D], [r[d] for r in recs if r["cls"] == "A_only"])
            ddb = auroc([r[d] for r in D], [r[d] for r in recs if r["cls"] == "B_only"])
            desc_mins[d] = min(dda, ddb)
        best = max(desc_mins, key=desc_mins.get)
        rows.append(
            {
                "pair": pair,
                "label_rule": "theta_6.0",
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
                "D_vs_neither_vina_mean": r4(dn),
                "D_vs_neither_ci_lo": r4(dn_lo),
                "D_vs_neither_ci_hi": r4(dn_hi),
                "D_vs_all_nonduals": r4(dall),
                "D_vs_all_nonduals_ci_lo": r4(dall_lo),
                "D_vs_all_nonduals_ci_hi": r4(dall_hi),
                "heavy_summary_min": r4(desc_mins["heavy"]),
                "mw_summary_min": r4(desc_mins["mw"]),
                "clogp_summary_min": r4(desc_mins["clogp"]),
                "tpsa_summary_min": r4(desc_mins["tpsa"]),
                "best_single_descriptor": best,
                "best_single_descriptor_summary_min": r4(desc_mins[best]),
                "note": "Table-2-comparable; does not replace frozen K=4 Table 2",
            }
        )
    return rows


def aggregation_s26(t2):
    rows = []
    for r in t2:
        a, b = float(r["auroc_D_vs_A_pocketB"]), float(r["auroc_D_vs_B_pocketA"])
        g = math.sqrt(a * b)
        h = 2 * a * b / (a + b) if (a + b) else float("nan")
        rows.append(
            {
                "pair": r["pair"],
                "auroc_D_vs_A": r4(a),
                "auroc_D_vs_B": r4(b),
                "summary_min": r4(min(a, b)),
                "summary_mean": r4((a + b) / 2),
                "summary_geometric": r4(g),
                "summary_harmonic": r4(h),
                "dual_vs_neither_vina_mean": r["D_vs_neither_vina_mean"],
                "note": "primary remains summary_min; other aggregates are sensitivity only",
            }
        )
    return rows


def equal_score_s34(packs):
    rows = []
    for pair, recs in packs.items():
        D = [r for r in recs if r["cls"] == "dual"]
        A = [r for r in recs if r["cls"] == "A_only"]
        B = [r for r in recs if r["cls"] == "B_only"]
        N = [r for r in recs if r["cls"] == "neither"]
        for contrast, sel, key, note in (
            ("D_vs_A_or_neither_pocketB", A, "vina_B", "pocket B fixed; A-only versus neither"),
            ("D_vs_B_or_neither_pocketA", B, "vina_A", "pocket A fixed; B-only versus neither"),
        ):
            dual = [r[key] for r in D]
            S = [r[key] for r in sel]
            Nn = [r[key] for r in N]
            rng = np.random.default_rng(SEED + stable_offset(pair, contrast))
            deltas = []
            for _ in range(N_BOOT):
                d = rng.choice(dual, len(dual), True)
                s = rng.choice(S, len(S), True)
                n = rng.choice(Nn, len(Nn), True)
                deltas.append(auroc(d, n) - auroc(d, s))
            lo, hi = np.percentile(deltas, [2.5, 97.5])
            a_s, a_n = auroc(dual, S), auroc(dual, Nn)
            rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "score": key,
                    "n_dual": len(dual),
                    "n_selective": len(S),
                    "n_neither": len(Nn),
                    "auroc_dual_vs_selective": r4(a_s),
                    "auroc_dual_vs_neither": r4(a_n),
                    "delta_neither_minus_selective": r4(a_n - a_s),
                    "delta_ci_lo": r4(lo),
                    "delta_ci_hi": r4(hi),
                    "underpowered_neither": int(len(Nn) < 8),
                    "note": note + "; shared dual resample, independent negatives; exploratory, not multiplicity-adjusted",
                }
            )
    return rows


def and_filter(packs):
    rows = []
    for pair, recs in packs.items():
        duals = [r for r in recs if r["cls"] == "dual"]
        lib = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
        for score in ("vina_worst", "vina_mean"):
            ds = np.array([r[score] for r in duals], float)
            for pct in (10, 25, 50, 75, 90):
                cut = float(np.percentile(ds, pct))
                passed = [r for r in lib if r[score] >= cut]
                n_pass = len(passed)
                n_d = sum(r["cls"] == "dual" for r in passed)
                n_a = sum(r["cls"] == "A_only" for r in passed)
                n_b = sum(r["cls"] == "B_only" for r in passed)
                rows.append(
                    {
                        "pair": pair,
                        "score": score,
                        "dual_percentile": pct,
                        "threshold": r4(cut),
                        "n_library": len(lib),
                        "n_dual_library": len(duals),
                        "n_pass": n_pass,
                        "n_dual_pass": n_d,
                        "n_A_only_pass": n_a,
                        "n_B_only_pass": n_b,
                        "precision_dual": r4(n_d / n_pass if n_pass else float("nan")),
                        "recall_dual": r4(n_d / len(duals) if duals else float("nan")),
                        "hardneg_fraction_pass": r4((n_a + n_b) / n_pass if n_pass else float("nan")),
                        "note": "AND-like filter on Dual+A-only+B-only; neither excluded",
                    }
                )
    return rows


def descriptor_paired_delta(packs):
    rows = []
    for pair, recs in packs.items():
        best = None
        best_min = -1
        for d in ("heavy", "mw", "clogp", "tpsa"):
            da, db, sm, *_ = directional(recs, d, d)
            if sm > best_min:
                best_min, best = sm, d
        stats = boot_paired(recs, "vina_B", "vina_A", best, best, seed=SEED + stable_offset(pair, "desc_delta", best))
        rows.append(
            {
                "pair": pair,
                "best_descriptor": best,
                **{k: stats[k] for k in stats if k.startswith("n_") or k.startswith("matched") or k.startswith("wrong") or k.startswith("delta") or k == "ci_excludes_zero"},
                "vina_summary_min": stats["matched_summary_min"],
                "descriptor_summary_min": stats["wrong_summary_min"],
                "delta_vina_minus_descriptor": stats["delta_matched_minus_wrong"],
                "note": "paired ligand-pool bootstrap of summary_min(Vina) minus summary_min(best single descriptor); exploratory",
            }
        )
    return rows


def incremental_ecfp4(packs):
    rows = []
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, dock_key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
            groups = np.array([r["scaffold"] for r in kept])
            dock = np.array([[r[dock_key]] for r in kept], dtype=float)
            fp = np.vstack([np.asarray(r["fp"]) for r in kept])
            rank_dock = auroc(
                [r[dock_key] for r in kept if r["cls"] == pos_cls],
                [r[dock_key] for r in kept if r["cls"] == neg_cls],
            )
            n_pos, n_neg = int(y.sum()), int((1 - y).sum())
            n_scaf = len(set(groups))
            n_splits = min(5, n_scaf, n_pos, n_neg)
            for name, X in (
                ("docking", dock),
                ("ECFP4", fp),
                ("ECFP4+docking", np.hstack([fp, dock])),
            ):
                auc = float("nan")
                if n_splits >= 2 and n_pos >= 6 and n_neg >= 6:
                    try:
                        cv = GroupKFold(n_splits=n_splits)
                        lr = LogisticRegression(max_iter=4000, C=1.0)
                        prob = cross_val_predict(lr, X, y, cv=cv, groups=groups, method="predict_proba")[:, 1]
                        auc = float(roc_auc_score(y, prob))
                    except Exception:
                        auc = float("nan")
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "model": name,
                        "n": len(kept),
                        "n_pos": n_pos,
                        "n_neg": n_neg,
                        "n_scaffolds": n_scaf,
                        "n_splits": n_splits if n_splits >= 2 else 0,
                        "cv_auroc": r4(auc),
                        "rank_auroc_docking": r4(rank_dock),
                        "note": "scaffold GroupKFold logistic; docking rank AUROC is not the logistic AUROC",
                    }
                )
    return rows


def caliper_table(packs):
    rows = []
    feats = ("mw", "clogp", "tpsa", "heavy")
    for pair, recs in packs.items():
        duals = [r for r in recs if r["cls"] == "dual"]
        for contrast, negs, score_key in (
            ("D_vs_A_pocketB", [r for r in recs if r["cls"] == "A_only"], "vina_B"),
            ("D_vs_B_pocketA", [r for r in recs if r["cls"] == "B_only"], "vina_A"),
        ):
            pool = duals + negs
            mu = np.array([np.mean([r[f] for r in pool]) for f in feats])
            sd = np.array([np.std([r[f] for r in pool], ddof=0) for f in feats])
            sd = np.where(sd < 1e-12, 1.0, sd)

            def z(row):
                return (np.array([row[f] for f in feats]) - mu) / sd

            zd = [z(r) for r in duals]
            zn = [z(r) for r in negs]
            order = np.argsort([r["ligand"] for r in duals])
            for caliper in CALIPERS:
                used = set()
                md, mn = [], []
                for i in order:
                    best = best_j = None
                    for j, vec in enumerate(zn):
                        if j in used:
                            continue
                        dist = float(np.linalg.norm(zd[i] - vec))
                        if dist <= caliper and (best is None or dist < best):
                            best, best_j = dist, j
                    if best_j is not None:
                        used.add(best_j)
                        md.append(duals[i])
                        mn.append(negs[best_j])
                n = min(len(md), len(mn))
                if n == 0:
                    pt = lo = hi = float("nan")
                else:
                    pt, lo, hi = boot_single(
                        [r[score_key] for r in md],
                        [r[score_key] for r in mn],
                        seed=SEED + stable_offset(pair, contrast, caliper),
                    )
                full_pt = auroc([r[score_key] for r in duals], [r[score_key] for r in negs])
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "caliper_sd": caliper,
                        "n_dual_matched": len(md),
                        "n_neg_matched": len(mn),
                        "n_dual_full": len(duals),
                        "n_neg_full": len(negs),
                        "auroc_matched": r4(pt),
                        "ci_lo": r4(lo),
                        "ci_hi": r4(hi),
                        "auroc_full": r4(full_pt),
                        "delta_matched_minus_full": r4(pt - full_pt) if pt == pt else "",
                        "underpowered": int(n < 8),
                        "note": "1:1 greedy match on z-scored MW/cLogP/TPSA/heavy; Euclidean caliper in SD units",
                    }
                )
    return rows


def le_and_subsets(packs):
    rows = []
    for pair, recs in packs.items():
        for name, kda, kdb in (
            ("pocket_matched", "vina_B", "vina_A"),
            ("wrong_pocket", "vina_A", "vina_B"),
            ("le_pocket_matched", "le_B", "le_A"),
            ("worst_pocket", "vina_worst", "vina_worst"),
            ("pooled_mean", "vina_mean", "vina_mean"),
        ):
            da, db, sm, nD, nA, nB = directional(recs, kda, kdb)
            lo, hi, n_ok = boot_pm_ci(recs, kda, kdb, seed=SEED + stable_offset(pair, name))
            rows.append(
                {
                    "pair": pair,
                    "aggregation": name,
                    "subset": "all_scored",
                    "auroc_D_vs_A": r4(da),
                    "auroc_D_vs_B": r4(db),
                    "summary_min": r4(sm),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "n_dual": nD,
                    "n_A_only": nA,
                    "n_B_only": nB,
                    "n_boot_ok": n_ok,
                }
            )
        # potency / size constrained greedy 1:1
        D = [r for r in recs if r["cls"] == "dual"]
        for contrast, others, pot_d, pot_o, kda, kdb in (
            ("D_vs_A", [r for r in recs if r["cls"] == "A_only"], "pB", "pB", "vina_B", "vina_B"),
            ("D_vs_B", [r for r in recs if r["cls"] == "B_only"], "pA", "pA", "vina_A", "vina_A"),
        ):
            used = set()
            kd, ko = [], []
            for d in D:
                best = best_i = None
                for i, o in enumerate(others):
                    if i in used:
                        continue
                    dpot = abs(d[pot_d] - o[pot_o])
                    dsz = abs(d["heavy"] - o["heavy"])
                    if dpot <= 0.5 and dsz <= 2.0:
                        dist = dpot + 0.1 * dsz
                        if best is None or dist < best:
                            best, best_i = dist, i
                if best_i is not None:
                    used.add(best_i)
                    kd.append(d)
                    ko.append(others[best_i])
            n = min(len(kd), len(ko))
            if n == 0:
                pt = lo = hi = float("nan")
            else:
                pt, lo, hi = boot_single(
                    [r[kda] for r in kd],
                    [r[kdb] for r in ko],
                    seed=SEED + stable_offset(pair, contrast, "pot_size"),
                )
            rows.append(
                {
                    "pair": pair,
                    "aggregation": "pocket_matched",
                    "subset": f"{contrast}_potency0.5_size2",
                    "auroc_D_vs_A": r4(pt) if contrast == "D_vs_A" else "",
                    "auroc_D_vs_B": r4(pt) if contrast == "D_vs_B" else "",
                    "summary_min": r4(pt),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "n_dual": n,
                    "n_A_only": n if contrast == "D_vs_A" else "",
                    "n_B_only": n if contrast == "D_vs_B" else "",
                    "n_boot_ok": "",
                    "underpowered": int(n < 8),
                    "note": "|ΔpChEMBL|≤0.5 and |Δheavy|≤2 greedy 1:1",
                }
            )
    return rows


def scaffold_cluster_bootstrap(packs):
    rows = []
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            groups = defaultdict(list)
            for r in kept:
                groups[r["scaffold"]].append(r)
            names = list(groups)
            point = auroc(
                [r[key] for r in kept if r["cls"] == pos_cls],
                [r[key] for r in kept if r["cls"] == neg_cls],
            )
            rng = np.random.default_rng(SEED + stable_offset(pair, contrast, "scaffold"))
            vals = []
            for _ in range(N_BOOT):
                chosen = rng.choice(names, size=len(names), replace=True)
                sub = [x for g in chosen for x in groups[g]]
                pos = [r[key] for r in sub if r["cls"] == pos_cls]
                neg = [r[key] for r in sub if r["cls"] == neg_cls]
                if len(pos) < 2 or len(neg) < 2:
                    continue
                vals.append(auroc(pos, neg))
            if len(vals) < N_BOOT // 2:
                lo = hi = float("nan")
            else:
                lo, hi = np.percentile(vals, [2.5, 97.5])
            rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "n_ligands": len(kept),
                    "n_scaffolds": len(names),
                    "n_boot_ok": len(vals),
                    "auroc": r4(point),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "note": "resample Bemis–Murcko scaffold groups; document-cluster blocked (no document_id on panels)",
                }
            )
    return rows


def threshold_grid(packs):
    rows = []
    rules = [("theta_5.5", 5.5, False), ("theta_6.0", 6.0, False), ("theta_6.5", 6.5, False), ("strict_6.5_5.5", None, True)]
    for pair, recs in packs.items():
        for rule_name, cut, is_strict in rules:
            labeled = []
            n_gray = 0
            for r in recs:
                lab = assign_strict(r["pA"], r["pB"]) if is_strict else assign_fourclass(r["pA"], r["pB"], cut)
                if lab == "gray":
                    n_gray += 1
                    continue
                if lab in ("dual", "A_only", "B_only", "neither"):
                    labeled.append({**r, "cls": lab})
            dir_recs = [r for r in labeled if r["cls"] in ("dual", "A_only", "B_only")]
            da, db, sm, nD, nA, nB = directional(dir_recs) if dir_recs else (float("nan"),) * 3 + (0, 0, 0)
            lo, hi, n_ok = boot_pm_ci(dir_recs, seed=SEED + stable_offset(pair, rule_name)) if dir_recs else (float("nan"), float("nan"), 0)
            nN = sum(r["cls"] == "neither" for r in labeled)
            rows.append(
                {
                    "pair": pair,
                    "label_rule": rule_name,
                    "n_scored_both_ends": len(recs),
                    "n_gray_excluded": n_gray,
                    "n_dual": nD,
                    "n_A_only": nA,
                    "n_B_only": nB,
                    "n_neither": nN,
                    "auroc_D_vs_A": r4(da),
                    "auroc_D_vs_B": r4(db),
                    "summary_min": r4(sm),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "n_boot_ok": n_ok,
                    "note": "same frozen Vina scores; labels recomputed from panel pChEMBL",
                }
            )
    return rows


def rank_extreme(packs, incomplete):
    fail_rows, sens = [], []
    for pair, recs in packs.items():
        inc = incomplete.get(pair, [])
        for r in inc:
            fail_rows.append(
                {
                    "pair": pair,
                    "ligand": r["ligand"],
                    "chembl_id": r["chembl_id"],
                    "class": r["cls"],
                    "score_A_available": int(r["vina_A"] is not None),
                    "score_B_available": int(r["vina_B"] is not None),
                    "heavy_atoms": int(r["heavy"]),
                    "mw": r4(r["mw"]),
                    "clogp": r4(r["clogp"]),
                    "tpsa": r4(r["tpsa"]),
                    "formal_charge": int(r["charge"]),
                    "rotatable_bonds": int(r["rotatable"]),
                    "status_reason": r.get("status_reason", ""),
                    "note": "largest organic fragment; descriptive coverage audit",
                }
            )
        records = []
        for r in recs:
            records.append({"class": r["cls"], "A": r["vina_A"], "B": r["vina_B"]})
        for r in inc:
            records.append({"class": r["cls"], "A": r["vina_A"], "B": r["vina_B"]})
        complete = [r for r in records if r["A"] is not None and r["B"] is not None]
        for contrast, neg_cls, end in (("D_vs_A_pocketB", "A_only", "B"), ("D_vs_B_pocketA", "B_only", "A")):
            cpos = [r[end] for r in complete if r["class"] == "dual"]
            cneg = [r[end] for r in complete if r["class"] == neg_cls]
            apos = [r[end] for r in records if r["class"] == "dual" and r[end] is not None]
            aneg = [r[end] for r in records if r["class"] == neg_cls and r[end] is not None]
            tpos = sum(r["class"] == "dual" for r in records)
            tneg = sum(r["class"] == neg_cls for r in records)
            observed_pairs = len(apos) * len(aneg)
            total_pairs = tpos * tneg
            observed_wins = auroc(apos, aneg) * observed_pairs if observed_pairs else 0.0
            sens.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "complete_case_n_pos": len(cpos),
                    "complete_case_n_neg": len(cneg),
                    "complete_case_auroc": r4(auroc(cpos, cneg)),
                    "arm_available_n_pos": len(apos),
                    "arm_available_n_neg": len(aneg),
                    "arm_available_auroc": r4(auroc(apos, aneg)),
                    "full_panel_n_pos": tpos,
                    "full_panel_n_neg": tneg,
                    "rank_extreme_lower_bound": r4(observed_wins / total_pairs if total_pairs else float("nan")),
                    "rank_extreme_upper_bound": r4(
                        (observed_wins + total_pairs - observed_pairs) / total_pairs if total_pairs else float("nan")
                    ),
                    "note": "bounds assign every missing comparison against/for the claimed direction; not an imputation model",
                }
            )
    return fail_rows, sens


def class_chemistry(packs):
    rows = []
    for pair, recs in packs.items():
        duals = [r for r in recs if r["cls"] == "dual"]
        for cls in ("dual", "A_only", "B_only", "neither"):
            sub = [r for r in recs if r["cls"] == cls]
            if not sub:
                continue
            nn = []
            if cls != "dual" and duals:
                for r in sub:
                    nn.append(max(DataStructs.TanimotoSimilarity(r["fp"], d["fp"]) for d in duals))
            scafs = [r["scaffold"] for r in sub]
            rows.append(
                {
                    "pair": pair,
                    "class": cls,
                    "n": len(sub),
                    "mw_median": r4(statistics.median(r["mw"] for r in sub)),
                    "mw_q1": r4(float(np.percentile([r["mw"] for r in sub], 25))),
                    "mw_q3": r4(float(np.percentile([r["mw"] for r in sub], 75))),
                    "heavy_median": r4(statistics.median(r["heavy"] for r in sub)),
                    "clogp_median": r4(statistics.median(r["clogp"] for r in sub)),
                    "tpsa_median": r4(statistics.median(r["tpsa"] for r in sub)),
                    "charge_median": r4(statistics.median(r["charge"] for r in sub)),
                    "rotatable_median": r4(statistics.median(r["rotatable"] for r in sub)),
                    "n_scaffolds": len(set(scafs)),
                    "singleton_scaffold_fraction": r4(
                        sum(1 for v in Counter(scafs).values() if v == 1) / len(sub)
                    ),
                    "nearest_dual_ecfp4_median": r4(statistics.median(nn)) if nn else "",
                    "note": "post-hoc descriptive; not a battery of unadjusted tests",
                }
            )
    return rows


def batch_auroc(pos_b, neg_b):
    diff = pos_b[:, :, None] - neg_b[:, None, :]
    return ((diff > 0).sum(axis=(1, 2)) + 0.5 * (diff == 0).sum(axis=(1, 2))) / (
        pos_b.shape[1] * neg_b.shape[1]
    )


def mu_from_auc(auc: float) -> float:
    from statistics import NormalDist

    if auc <= 0:
        return -np.inf
    if auc >= 1:
        return np.inf
    return float(np.sqrt(2.0) * NormalDist().inv_cdf(float(auc)))


def detectable_effect(t2):
    from statistics import NormalDist  # noqa: F401

    rng = np.random.default_rng(SEED)
    rows = []
    for r in t2:
        n_d, n_a, n_b, n_n = r["n_dual"], r["n_A_only"], r["n_B_only"], r["n_neither"]
        for true_auc in TRUE_AUCS:
            mu = mu_from_auc(true_auc)
            print(f"detectable {r['pair']} true={true_auc:.2f}", flush=True)
            for contrast, n_pos, n_neg in (
                ("dual_vs_A_only", n_d, n_a),
                ("dual_vs_B_only", n_d, n_b),
                ("dual_vs_neither", n_d, n_n),
            ):
                excl = 0
                for _ in range(N_MC):
                    pos = rng.normal(mu, 1.0, size=n_pos)
                    neg = rng.normal(0.0, 1.0, size=n_neg)
                    pb = pos[rng.integers(0, n_pos, size=(N_BOOT, n_pos))]
                    nb = neg[rng.integers(0, n_neg, size=(N_BOOT, n_neg))]
                    aucs = batch_auroc(pb, nb)
                    lo, hi = np.percentile(aucs, [2.5, 97.5])
                    excl += int(lo > 0.5 or hi < 0.5)
                p = excl / N_MC
                rows.append(
                    {
                        "pair": r["pair"],
                        "contrast": contrast,
                        "n_pos": n_pos,
                        "n_neg": n_neg,
                        "true_auroc": f"{true_auc:.2f}",
                        "n_mc": N_MC,
                        "n_boot": N_BOOT,
                        "p_ci_excludes_0p5": f"{p:.6g}",
                        "se_binomial": f"{math.sqrt(p * (1 - p) / N_MC):.6g}",
                    }
                )
            excl = 0
            for _ in range(N_MC):
                dual_sa = rng.normal(mu, 1.0, size=n_d)
                dual_sb = rng.normal(mu, 1.0, size=n_d)
                a_sb = rng.normal(0.0, 1.0, size=n_a)
                b_sa = rng.normal(0.0, 1.0, size=n_b)
                idd = rng.integers(0, n_d, size=(N_BOOT, n_d))
                ida = rng.integers(0, n_a, size=(N_BOOT, n_a))
                idb = rng.integers(0, n_b, size=(N_BOOT, n_b))
                smin = np.minimum(batch_auroc(dual_sb[idd], a_sb[ida]), batch_auroc(dual_sa[idd], b_sa[idb]))
                lo, hi = np.percentile(smin, [2.5, 97.5])
                excl += int(lo > 0.5 or hi < 0.5)
            p = excl / N_MC
            rows.append(
                {
                    "pair": r["pair"],
                    "contrast": "summary_min",
                    "n_pos": n_d,
                    "n_neg": f"{n_a}/{n_b}",
                    "true_auroc": f"{true_auc:.2f}",
                    "n_mc": N_MC,
                    "n_boot": N_BOOT,
                    "p_ci_excludes_0p5": f"{p:.6g}",
                    "se_binomial": f"{math.sqrt(p * (1 - p) / N_MC):.6g}",
                }
            )
    return rows


def leftover_and_blocked():
    summary = list(csv.DictReader((TAB / "track_b_panel_summary_v1.csv").open()))
    rows = []
    for r in summary:
        leftover_d = int(r["n_strict_dual_smallmol"]) - int(r["n_panel_dual"])
        leftover_a = int(r["n_strict_A_only_smallmol"]) - int(r["n_panel_A_only"])
        leftover_b = int(r["n_strict_B_only_smallmol"]) - int(r["n_panel_B_only"])
        eligible = leftover_d >= 20 and leftover_a >= 20 and leftover_b >= 20
        thin = eligible and min(leftover_d, leftover_a, leftover_b) < 25
        rows.append(
            {
                "pair": r["pair"],
                "n_strict_dual_smallmol": r["n_strict_dual_smallmol"],
                "n_strict_A_only_smallmol": r["n_strict_A_only_smallmol"],
                "n_strict_B_only_smallmol": r["n_strict_B_only_smallmol"],
                "n_panel_dual": r["n_panel_dual"],
                "n_panel_A_only": r["n_panel_A_only"],
                "n_panel_B_only": r["n_panel_B_only"],
                "leftover_dual": leftover_d,
                "leftover_A_only": leftover_a,
                "leftover_B_only": leftover_b,
                "holdout_20_20_20_eligible": int(eligible),
                "holdout_thin_margin": int(thin),
                "note": (
                    "eligible under the original 20/20/20 leftover gate; JAK1/JAK2 leftover B-only=21 is thin (margin=1)"
                    if r["pair"] == "JAK1/JAK2"
                    else "IDs not extracted here: ChEMBL sqlite required"
                ),
            }
        )
    blocked = [
        {
            "item": "max_vs_median_pchembl",
            "status": "blocked_no_sqlite",
            "why": "panel CSV stores one pChEMBL per end; repeat-record graph needs ChEMBL 37 sqlite",
        },
        {
            "item": "document_year_split",
            "status": "blocked_no_sqlite",
            "why": "no document.year on panels; report counts only after dump join; AUROC only if dual/A/B each n≥10 after 2018",
        },
        {
            "item": "document_cluster_bootstrap",
            "status": "blocked_no_sqlite",
            "why": "no document_id; scaffold-cluster bootstrap was computed instead",
        },
        {
            "item": "document_blocked_cv",
            "status": "blocked_no_sqlite",
            "why": "same missing document_id; ECFP4 used Bemis–Murcko GroupKFold",
        },
        {
            "item": "bindingdb_pubchem_count_only",
            "status": "blocked_no_cache",
            "why": "jcim_supply_crossdb_v0 caches only the original K=4 UniProts; five new pairs need a new count-only fetch (no Docker)",
        },
        {
            "item": "holdout_panel_ids",
            "status": "blocked_no_sqlite",
            "why": "leftover counts are known; 20/20/20 member lists need the dump + HOLDOUT_SEED=20260731",
        },
        {
            "item": "five_seed_vina",
            "status": "local_recompute",
            "why": "user will submit locally; seeds 20260727 + 20260811–20260814; see LOCAL_RECOMPUTE_PACK_V1.md",
        },
        {
            "item": "rtm_best_of_9",
            "status": "local_recompute",
            "why": "poses gitignored; regenerate 9 modes then rtmscore_model1",
        },
        {
            "item": "gnina_cnn_rescore",
            "status": "local_recompute",
            "why": "same poses; --cnn_scoring rescore --minimize",
        },
        {
            "item": "independent_gnina_search",
            "status": "local_rule_subset",
            "why": "JAK1/TYK2 only (EGFR-like formulation gap)",
        },
    ]
    return rows, blocked


def write_verdict(t2, leftover, blocked, incr):
    lines = [
        "# Five-pair zero-dock stack (same article, not Track B)\n\n",
        "Destination: 8-row main table after withdrawing PIK3CA/PIK3CB ",
        "(`PROJECT_IDENTITY_LOCK_V1.md`). These five pairs were added after the ",
        "ChEMBL 37 census; they were not frozen on 2026-07-23. This run does ",
        "**not** restock Table 2 or change the title.\n\n",
        "Bootstrap for `summary_min` is the Table 2 estimand: ligand-level ",
        "non-stratified resample of the dual+A-only+B-only pool, B=2000, seed ",
        "20260729 + SHA-256 pair offset. Existing `track_b_directional_auroc_v1.csv` ",
        "CIs remain class-preserving and are not Table 2 intervals.\n\n",
        "## Table-2-comparable Vina (θ = 6.0, both-end scores)\n\n",
        "| pair | n D/A/B | D/A (B) | D/B (A) | summary_min [95% CI] | Dual-vs-neither | best descriptor |\n",
        "|---|---:|---:|---:|---|---:|---|\n",
    ]
    for r in t2:
        lines.append(
            f"| {r['pair']} | {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']} | "
            f"{r['auroc_D_vs_A_pocketB']} | {r['auroc_D_vs_B_pocketA']} | "
            f"{r['summary_min']} [{r['ci_lo']}, {r['ci_hi']}] | "
            f"{r['D_vs_neither_vina_mean']} | {r['best_single_descriptor']} "
            f"{r['best_single_descriptor_summary_min']} |\n"
        )
    lines.append("\n## Holdout leftover (counts only; IDs need sqlite)\n\n")
    lines.append("| pair | leftover D/A/B | 20/20/20 | thin |\n|---|---:|---:|---:|\n")
    for r in leftover:
        lines.append(
            f"| {r['pair']} | {r['leftover_dual']}/{r['leftover_A_only']}/{r['leftover_B_only']} | "
            f"{r['holdout_20_20_20_eligible']} | {r['holdout_thin_margin']} |\n"
        )
    lines.append("\n## Still blocked or local\n\n")
    for b in blocked:
        lines.append(f"- **{b['item']}** — `{b['status']}`: {b['why']}\n")
    deltas = []
    by = defaultdict(dict)
    for r in incr:
        by[(r["pair"], r["contrast"])][r["model"]] = r
    lines.append("\n## ECFP4 scaffold GroupKFold (docking increment)\n\n")
    lines.append("| pair | contrast | ECFP4 | ECFP4+docking | Δ | docking rank |\n|---|---|---:|---:|---:|---:|\n")
    for (pair, contrast), m in by.items():
        e = m.get("ECFP4", {})
        ed = m.get("ECFP4+docking", {})
        try:
            dlt = float(ed.get("cv_auroc") or "nan") - float(e.get("cv_auroc") or "nan")
        except (TypeError, ValueError):
            dlt = float("nan")
        lines.append(
            f"| {pair} | {contrast} | {e.get('cv_auroc', '')} | {ed.get('cv_auroc', '')} | "
            f"{r4(dlt)} | {e.get('rank_auroc_docking', '')} |\n"
        )
        deltas.append(dlt)
    lines.append(
        "\nIndependent GNINA search is **not** in this run. JAK1/TYK2 is the only "
        "new pair that qualifies under the original formulation-gap rule.\n"
    )
    (AN / "FIVE_PAIR_STACK_V1.md").write_text("".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--skip-s31",
        action="store_true",
        help="Skip the long S31 Monte Carlo; run later with --s31-only.",
    )
    ap.add_argument("--s31-only", action="store_true", help="Only run S31 detectable-effect.")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    scores = load_scores()
    jobs = load_job_status()
    packs = {}
    incomplete = {}
    for spec in PAIRS:
        recs, inc = assemble(spec, scores, jobs, "theta6_class")
        packs[spec["pair"]] = recs
        incomplete[spec["pair"]] = inc
        print(
            f"{spec['pair']}: scored={len(recs)} class={dict(Counter(r['cls'] for r in recs))} incomplete={len(inc)}",
            flush=True,
        )

    if args.s31_only:
        t2_path = OUT / "table2_comparable_theta6_v1.csv"
        if t2_path.exists() and t2_path.stat().st_size > 0:
            t2 = list(csv.DictReader(t2_path.open()))
            for r in t2:
                for k in ("n_dual", "n_A_only", "n_B_only", "n_neither"):
                    r[k] = int(r[k])
        else:
            t2 = table2_comparable(packs)
            write_csv(OUT / "table2_comparable_theta6_v1.csv", t2)
        print("running detectable-effect simulation (S31 settings)...", flush=True)
        det = detectable_effect(t2)
        write_csv(OUT / "detectable_effect_s31_v1.csv", det)
        print("wrote", OUT / "detectable_effect_s31_v1.csv")
        return 0

    t2 = table2_comparable(packs)
    write_csv(OUT / "table2_comparable_theta6_v1.csv", t2)
    write_csv(OUT / "aggregation_s26_v1.csv", aggregation_s26(t2))
    write_csv(OUT / "wrong_pocket_paired_delta_v1.csv", [
        {"pair": pair, **boot_paired(recs, "vina_B", "vina_A", "vina_A", "vina_B", seed=SEED + stable_offset(pair, "wrong_pocket"))}
        for pair, recs in packs.items()
    ])
    write_csv(OUT / "equal_score_negative_s34_v1.csv", equal_score_s34(packs))
    write_csv(OUT / "and_filter_s46_v1.csv", and_filter(packs))
    write_csv(OUT / "descriptor_paired_delta_s19_v1.csv", descriptor_paired_delta(packs))
    incr = incremental_ecfp4(packs)
    write_csv(OUT / "ecfp4_incremental_s20s24_v1.csv", incr)
    write_csv(OUT / "property_caliper_s45_v1.csv", caliper_table(packs))
    write_csv(OUT / "le_and_subsets_v1.csv", le_and_subsets(packs))
    write_csv(OUT / "scaffold_cluster_bootstrap_v1.csv", scaffold_cluster_bootstrap(packs))
    write_csv(OUT / "threshold_grid_v1.csv", threshold_grid(packs))
    fails, sens = rank_extreme(packs, incomplete)
    write_csv(OUT / "docking_failed_ligand_properties_v1.csv", fails)
    write_csv(OUT / "docking_failure_rank_extreme_v1.csv", sens)
    write_csv(OUT / "class_chemistry_s38_v1.csv", class_chemistry(packs))
    leftover, blocked = leftover_and_blocked()
    write_csv(OUT / "holdout_leftover_counts_v1.csv", leftover)
    write_csv(OUT / "blocked_or_local_v1.csv", blocked)

    if args.skip_s31:
        print("Skipping S31 (use --s31-only later).", flush=True)
    else:
        print("running detectable-effect simulation (S31 settings)...", flush=True)
        det = detectable_effect(t2)
        write_csv(OUT / "detectable_effect_s31_v1.csv", det)

    write_verdict(t2, leftover, blocked, incr)
    (OUT / "run_meta_v1.json").write_text(
        json.dumps(
            {
                "seed": SEED,
                "n_boot": N_BOOT,
                "n_mc": N_MC,
                "s31_ran": (not args.skip_s31),
                "identity": "PROJECT_IDENTITY_LOCK_V1.md",
                "table2_bootstrap": "ligand_non_stratified_dualA_B_pool",
                "does_not_replace_table2": True,
                "pairs": [p["pair"] for p in PAIRS],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print("wrote", OUT)
    print("wrote", AN / "FIVE_PAIR_STACK_V1.md")
    for r in t2:
        print(
            f"T2-comparable {r['pair']}: {r['summary_min']} [{r['ci_lo']}, {r['ci_hi']}] "
            f"neither={r['D_vs_neither_vina_mean']} best={r['best_single_descriptor']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
