#!/usr/bin/env python3
"""JCIM strengthen supplement — Wave 1 zero-dock analyses (T0.3–T0.9).

Outputs under data/jcim_strengthen_t0t1_v0/tables/
"""
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
except ImportError:
    LogisticRegression = None

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_strengthen_t0t1_v0"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729
CUTOFFS = (5.5, 6.0, 6.5)
STRICT_HI, STRICT_LO = 6.5, 5.5

PAIR_SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        panel="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        pA="pchembl_EGFR",
        pB="pchembl_HER2",
        lig="ligand",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        panel=None,
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        pA="pchembl_ACHE",
        pB="pchembl_BCHE",
        lig="ligand",
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        panel="data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        pA="pchembl_PIK3CA",
        pB="pchembl_PIK3CB",
        lig="ligand",
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        pA="pchembl_PIK3CA",
        pB="pchembl_MTOR",
        lig="ligand",
    ),
}


def fnum(x):
    try:
        if x in ("", None):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load_csv(p: Path):
    with p.open() as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if not pos or not neg:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def assign_fourclass(pA, pB, cut: float) -> str | None:
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


def assign_strict(pA, pB) -> str | None:
    if pA is None or pB is None:
        return None
    if pA >= STRICT_HI and pB >= STRICT_HI:
        return "dual"
    if pA >= STRICT_HI and pB <= STRICT_LO:
        return "A_only"
    if pB >= STRICT_HI and pA <= STRICT_LO:
        return "B_only"
    if pA <= STRICT_LO and pB <= STRICT_LO:
        return "neither"
    return "gray"


def murcko(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def assemble_records() -> dict[str, list[dict]]:
    packs = {}
    for pair, cfg in PAIR_SPEC.items():
        rows = load_csv(ROOT / cfg["scores"])
        meta = {}
        if cfg["panel"]:
            for r in load_csv(ROOT / cfg["panel"]):
                meta[r.get("panel_id") or r.get("ligand")] = r
        out = []
        for r in rows:
            lig = r.get(cfg["lig"]) or r.get("panel_id")
            m = meta.get(lig, {})
            smi = r.get("smiles") or m.get("smiles")
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            va, vb = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
            if va is None or vb is None:
                continue
            pA = fnum(r.get(cfg["pA"]) or m.get(cfg["pA"]))
            pB = fnum(r.get(cfg["pB"]) or m.get(cfg["pB"]))
            ha = mol.GetNumHeavyAtoms()
            out.append(
                {
                    "pair": pair,
                    "ligand": lig,
                    "cls": r.get("class") or m.get("class"),
                    "smiles": smi,
                    "pA": pA,
                    "pB": pB,
                    "heavy": float(ha),
                    "tpsa": float(Descriptors.TPSA(mol)),
                    "clogp": float(Descriptors.MolLogP(mol)),
                    "vina_A": -va,
                    "vina_B": -vb,
                    "vina_mean": -(va + vb) / 2,
                    "vina_worst": min(-va, -vb),
                    "le_A": -va / ha,
                    "le_B": -vb / ha,
                    "le_B_pm": -vb / ha,
                    "le_A_pm": -va / ha,
                }
            )
        packs[pair] = out
    return packs


def directional_pm(recs, key_da="vina_B", key_db="vina_A"):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r[key_da] for r in D], [r[key_da] for r in A])
    db = auroc([r[key_db] for r in D], [r[key_db] for r in B])
    return da, db, min(da, db), len(D), len(A), len(B)


def boot_pm_ci(recs, key_da="vina_B", key_db="vina_A", n_boot=N_BOOT, seed=SEED):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    if len(usable) < 8:
        return None, None
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        sub = [usable[i] for i in ii]
        _, _, mn, _, _, _ = directional_pm(sub, key_da, key_db)
        if mn == mn:
            mins.append(mn)
    if len(mins) < n_boot // 2:
        return None, None
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi)


def nearest_match(duals, others, key_potency, key_size, dp=0.5, ds=2.0):
    """Keep dual+other pairs where |Δpotency|<=dp and |Δheavy|<=ds."""
    kept_d, kept_o = [], []
    used = set()
    for d in duals:
        best, best_dist = None, 1e9
        for i, o in enumerate(others):
            if i in used:
                continue
            if d.get(key_potency) is None or o.get(key_potency) is None:
                continue
            dpot = abs(d[key_potency] - o[key_potency])
            dsz = abs(d["heavy"] - o["heavy"])
            if dpot <= dp and dsz <= ds:
                dist = dpot + 0.1 * dsz
                if dist < best_dist:
                    best_dist, best = dist, i
        if best is not None:
            used.add(best)
            kept_d.append(d)
            kept_o.append(others[best])
    return kept_d, kept_o


def matched_subset_rows(packs):
    rows = []
    for pair, recs in packs.items():
        for match_type, key_pot in (
            ("potency_matched_D_vs_A", "pA"),
            ("potency_matched_D_vs_B", "pB"),
            ("size_matched_D_vs_A", None),
            ("size_matched_D_vs_B", None),
        ):
            duals = [r for r in recs if r["cls"] == "dual"]
            if match_type.startswith("potency"):
                others = [r for r in recs if r["cls"] == ("A_only" if "A" in match_type else "B_only")]
                kd, ko = nearest_match(duals, others, key_pot, "heavy", dp=0.5, ds=999)
            else:
                others = [r for r in recs if r["cls"] == ("A_only" if "A" in match_type else "B_only")]
                kd, ko = nearest_match(duals, others, "pA", "heavy", dp=999, ds=2.0)
            if match_type.endswith("A"):
                sub = kd + ko
                da, _, mn, nd, na, _ = directional_pm(sub, "vina_B", "vina_A")
                n_arm = min(nd, na)
            else:
                sub = kd + ko
                _, db, mn, nd, _, nb = directional_pm(sub, "vina_A", "vina_B")
                da, mn = db, db
                n_arm = min(nd, nb)
            ci_lo, ci_hi = boot_pm_ci(
                sub,
                "vina_B" if "A" in match_type else "vina_A",
                "vina_A" if "A" in match_type else "vina_B",
                seed=SEED + hash((pair, match_type)) % 99991,
            )
            rows.append(
                {
                    "pair": pair,
                    "match_type": match_type,
                    "n_dual_matched": len(kd),
                    "n_other_matched": len(ko),
                    "n_contrast_min": n_arm,
                    "underpowered": int(n_arm < 8),
                    "auroc_contrast": round(mn, 4) if mn == mn else "",
                    "ci_lo": round(ci_lo, 4) if ci_lo is not None else "",
                    "ci_hi": round(ci_hi, 4) if ci_hi is not None else "",
                    "note": f"|Δp|≤0.5 potency or |Δheavy|≤2 size; pocket-matched vina",
                }
            )
        # combined min for pair
        da, db, mn, nd, na, nb = directional_pm(recs)
        ci_lo, ci_hi = boot_pm_ci(recs)
        rows.append(
            {
                "pair": pair,
                "match_type": "full_panel_pocket_matched",
                "n_dual_matched": nd,
                "n_other_matched": na + nb,
                "n_contrast_min": min(nd, na, nb),
                "underpowered": int(min(nd, na, nb) < 8),
                "auroc_contrast": round(mn, 4),
                "ci_lo": round(ci_lo, 4) if ci_lo else "",
                "ci_hi": round(ci_hi, 4) if ci_hi else "",
                "note": "unmatched full panel reference",
            }
        )
    return rows


def covariate_adjusted(packs):
    if LogisticRegression is None:
        return [{"error": "sklearn not available"}]
    rows = []
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, score_key in (
            ("D_vs_A_only", "dual", "A_only", "vina_B"),
            ("D_vs_B_only", "dual", "B_only", "vina_A"),
        ):
            sub = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            if len(sub) < 10:
                continue
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in sub])
            score = np.array([r[score_key] for r in sub]).reshape(-1, 1)
            heavy = np.array([r["heavy"] for r in sub]).reshape(-1, 1)
            tpsa = np.array([r["tpsa"] for r in sub]).reshape(-1, 1)
            X_full = np.hstack([score, heavy, tpsa])
            try:
                m_score = LogisticRegression(max_iter=2000, random_state=SEED)
                m_full = LogisticRegression(max_iter=2000, random_state=SEED)
                m_score.fit(score, y)
                m_full.fit(X_full, y)
                prob_s = m_score.predict_proba(score)[:, 1]
                prob_f = m_full.predict_proba(X_full)[:, 1]
                auc_s = roc_auc_score(y, prob_s)
                auc_f = roc_auc_score(y, prob_f)
                coef = m_full.coef_[0]
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "score_feature": score_key,
                        "n": len(sub),
                        "n_pos": int(y.sum()),
                        "n_neg": int(len(y) - y.sum()),
                        "auroc_score_only": round(auc_s, 4),
                        "auroc_score_plus_covariates": round(auc_f, 4),
                        "delta_auroc": round(auc_f - auc_s, 4),
                        "coef_score": round(float(coef[0]), 4),
                        "coef_heavy": round(float(coef[1]), 4),
                        "coef_tpsa": round(float(coef[2]), 4),
                        "or_score": round(math.exp(coef[0]), 4),
                    }
                )
            except Exception as e:
                rows.append({"pair": pair, "contrast": contrast, "error": str(e)})
    return rows


def aggregation_sensitivity(packs):
    variants = [
        ("pocket_matched", "vina_B", "vina_A"),
        ("wrong_pocket", "vina_A", "vina_B"),
        ("worst_pocket", "vina_worst", "vina_worst"),
        ("pooled_mean", "vina_mean", "vina_mean"),
        ("le_pocket_matched", "le_B", "le_A"),
    ]
    rows = []
    for pair, recs in packs.items():
        for name, kda, kdb in variants:
            da, db, mn, nd, na, nb = directional_pm(recs, kda, kdb)
            ci_lo, ci_hi = boot_pm_ci(recs, kda, kdb, seed=SEED + hash((pair, name)) % 99991)
            rows.append(
                {
                    "pair": pair,
                    "aggregation": name,
                    "auroc_D_vs_A": round(da, 4),
                    "auroc_D_vs_B": round(db, 4),
                    "summary_min": round(mn, 4),
                    "ci_lo": round(ci_lo, 4) if ci_lo else "",
                    "ci_hi": round(ci_hi, 4) if ci_hi else "",
                    "n_dual": nd,
                    "n_A_only": na,
                    "n_B_only": nb,
                }
            )
    return rows


def unified_threshold_v2(packs):
    rows = []
    rules = [(f"theta_{c}", c, False) for c in CUTOFFS] + [("strict_6.5_5.5", None, True)]
    for pair, recs in packs.items():
        base = [r for r in recs if r.get("pA") is not None and r.get("pB") is not None]
        if len(base) < 15:
            continue
        for rule_name, cut, is_strict in rules:
            labeled = []
            for r in base:
                if is_strict:
                    lab = assign_strict(r["pA"], r["pB"])
                else:
                    lab = assign_fourclass(r["pA"], r["pB"], cut)
                if lab in ("dual", "A_only", "B_only"):
                    labeled.append({**r, "cls": lab})
            if len(labeled) < 8:
                continue
            da, db, mn, nd, na, nb = directional_pm(labeled)
            ci_lo, ci_hi = boot_pm_ci(labeled, seed=SEED + hash((pair, rule_name)) % 99991)
            counts = {c: sum(1 for r in labeled if r["cls"] == c) for c in ("dual", "A_only", "B_only", "neither")}
            rows.append(
                {
                    "pair": pair,
                    "label_rule": rule_name,
                    "n_dual": nd,
                    "n_A_only": na,
                    "n_B_only": nb,
                    "n_neither_excluded": counts.get("neither", 0),
                    "pocket_matched_summary_min": round(mn, 4) if mn == mn else "",
                    "auroc_D_vs_A": round(da, 4) if da == da else "",
                    "auroc_D_vs_B": round(db, 4) if db == db else "",
                    "ci_lo": round(ci_lo, 4) if ci_lo else "",
                    "ci_hi": round(ci_hi, 4) if ci_hi else "",
                }
            )
    return rows


def chembl_aggregation_sensitivity(packs):
    """Compare max (mols_*.json) vs min-end aggregation for paired ChEMBL IDs."""
    target_map = {
        "EGFR/HER2": ("EGFR", "HER2"),
        "AChE/BChE": ("ACHE", "BCHE"),
        "PIK3CA/PIK3CB": ("PIK3CA", "PIK3CB"),
        "PIK3CA/mTOR": ("PIK3CA", "MTOR"),
    }
    mols_cache = {}
    rows = []
    for pair, (ta, tb) in target_map.items():
        for t in (ta, tb):
            if t not in mols_cache:
                p = ROOT / f"data/public_pair_selection/mols_{t}.json"
                mols_cache[t] = json.loads(p.read_text()) if p.exists() else {}
        ma, mb = mols_cache[ta], mols_cache[tb]
        common = set(ma.keys()) & set(mb.keys())
        if len(common) < 50:
            continue
        for agg in ("max_per_target", "min_of_two_max"):
            labeled = []
            for cid in common:
                pa, pb = ma[cid], mb[cid]
                if agg == "min_of_two_max":
                    # conservative: use weaker end as both (stress test)
                    v = min(pa, pb)
                    pa, pb = v, v
                lab = assign_strict(pa, pb)
                if lab not in ("dual", "A_only", "B_only"):
                    continue
                labeled.append({"cls": lab, "pA": pa, "pB": pb, "chembl_id": cid})
            if len(labeled) < 20:
                continue
            rows.append(
                {
                    "pair": pair,
                    "aggregation": agg,
                    "n_common_chembl": len(common),
                    "n_labeled_3class": len(labeled),
                    "n_dual": sum(1 for r in labeled if r["cls"] == "dual"),
                    "n_A_only": sum(1 for r in labeled if r["cls"] == "A_only"),
                    "n_B_only": sum(1 for r in labeled if r["cls"] == "B_only"),
                    "note": "ChEMBL supply-level label counts only; no docking scores in mols_*.json",
                }
            )
    return rows


def scaffold_inventory(packs):
    rows = []
    for pair, recs in packs.items():
        for r in recs:
            sc = murcko(r["smiles"])
            rows.append(
                {
                    "pair": pair,
                    "ligand": r["ligand"],
                    "cls": r["cls"],
                    "murcko_scaffold": sc,
                    "heavy": r["heavy"],
                    "pA": r.get("pA"),
                    "pB": r.get("pB"),
                }
            )
    return rows


def scaffold_bootstrap(packs, inv):
    rows = []
    by_pair = defaultdict(lambda: defaultdict(list))
    for r in inv:
        by_pair[r["pair"]][r["murcko_scaffold"]].append(r["ligand"])

    for pair, recs in packs.items():
        lig2rec = {r["ligand"]: r for r in recs}
        clusters = {sc: [lig2rec[l] for l in ligs if l in lig2rec] for sc, ligs in by_pair[pair].items()}
        clusters = {sc: rs for sc, rs in clusters.items() if rs}
        if len(clusters) < 3:
            continue
        cluster_ids = list(clusters.keys())
        rng = np.random.default_rng(SEED + hash(pair) % 99991)
        mins = []
        for _ in range(N_BOOT):
            picked = rng.choice(cluster_ids, size=len(cluster_ids), replace=True)
            sub = []
            for sc in picked:
                sub.extend(clusters[sc])
            if len(sub) < 8:
                continue
            _, _, mn, _, _, _ = directional_pm(sub)
            if mn == mn:
                mins.append(mn)
        if len(mins) < N_BOOT // 2:
            continue
        lo, hi = np.percentile(mins, [2.5, 97.5])
        _, _, point, nd, na, nb = directional_pm(recs)
        rows.append(
            {
                "pair": pair,
                "n_scaffolds": len(clusters),
                "n_molecules": len(recs),
                "pocket_matched_summary_min": round(point, 4),
                "scaffold_boot_ci_lo": round(float(lo), 4),
                "scaffold_boot_ci_hi": round(float(hi), 4),
                "n_boot_ok": len(mins),
            }
        )
    return rows


def ligand_ml_baseline(packs):
    if LogisticRegression is None:
        return []
    rows = []
    for pair, recs in packs.items():
        sub = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
        if len(sub) < 20:
            continue
        mols, y = [], []
        for r in sub:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
            mols.append(np.array(fp))
            y.append(r["cls"])
        if len(mols) < 20:
            continue
        X = np.vstack(mols)
        y_arr = np.array(y)
        for contrast, pos, neg, dock_key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            mask = np.isin(y_arr, [pos, neg])
            y_bin = (y_arr[mask] == pos).astype(int)
            X_sub = X[mask]
            scores = np.array([r[dock_key] for r, m in zip(sub, range(len(sub))) if y_arr[m] in (pos, neg)])
            if len(set(y_bin)) < 2 or len(y_bin) < 12:
                continue
            try:
                cv = StratifiedKFold(n_splits=min(5, min(np.bincount(y_bin))), shuffle=True, random_state=SEED)
                lr = LogisticRegression(max_iter=2000, C=1.0)
                prob = cross_val_predict(lr, X_sub, y_bin, cv=cv, method="predict_proba")[:, 1]
                auc_ml = roc_auc_score(y_bin, prob)
                auc_dock = roc_auc_score(y_bin, scores)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "method": "ECFP4_logistic_nested_cv",
                        "n": len(y_bin),
                        "auroc_ml": round(auc_ml, 4),
                        "auroc_dock_pocket_matched": round(auc_dock, 4),
                        "delta_ml_minus_dock": round(auc_ml - auc_dock, 4),
                    }
                )
            except Exception:
                pass
    return rows


def main():
    packs = assemble_records()
    write_csv(TAB / "matched_subset_directional_v1.csv", matched_subset_rows(packs))
    write_csv(TAB / "covariate_adjusted_v1.csv", covariate_adjusted(packs))
    write_csv(TAB / "aggregation_sensitivity_v1.csv", aggregation_sensitivity(packs))
    write_csv(TAB / "unified_threshold_sensitivity_v2.csv", unified_threshold_v2(packs))
    chembl = chembl_aggregation_sensitivity(packs)
    write_csv(TAB / "chembl_aggregation_sensitivity_v1.csv", chembl)
    inv = scaffold_inventory(packs)
    write_csv(TAB / "scaffold_inventory_v1.csv", inv)
    write_csv(TAB / "scaffold_bootstrap_ci_v1.csv", scaffold_bootstrap(packs, inv))
    ml = ligand_ml_baseline(packs)
    if ml:
        write_csv(TAB / "ligand_ml_baseline_v1.csv", ml)

    skips = []
    skips.append(
        "T0.7 median/confidence≥8/Homo sapiens: mols_*.json stores single float (max pChEMBL) per ChEMBL ID; "
        "no per-assay median or confidence fields available locally."
    )
    if not chembl:
        skips.append("T0.7 chembl_aggregation_sensitivity_v1.csv: insufficient paired ChEMBL overlap for relabel test.")
    (OUT / "analysis").mkdir(exist_ok=True)
    (OUT / "analysis" / "T0_SKIPS.md").write_text("# T0 skips\n\n" + "\n\n".join(f"- {s}" for s in skips) + "\n")
    print("Wrote T0 strengthen tables to", TAB)


if __name__ == "__main__":
    main()
