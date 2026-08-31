#!/usr/bin/env python3
"""A-level novelty analyses on frozen K=4 scores. No new docking.

Compares conventional dual-target evaluation (dual vs neither / dual vs all
non-duals; single-target inhibitor vs noninhibitor) with DualFourClass
directional hard-negative AUROC; adds chemotype-matched hard negatives,
incremental ligand vs docking information, and mixed-library enrichment.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_novelty_v0"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729

SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        panel="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        panel_key="panel_id",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        panel_key="panel_id",
    ),
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
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
        w = csv.DictWriter(
            fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def boot_auroc(pos, neg, n_boot=N_BOOT, seed=SEED):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        p = rng.choice(pos, size=len(pos), replace=True)
        n = rng.choice(neg, size=len(neg), replace=True)
        vals.append(auroc(p, n))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc(pos, neg)), float(lo), float(hi)


def boot_equal_score_negative_contrast(dual, selective, neither, n_boot=N_BOOT, seed=SEED):
    """Contrast negative classes while holding the score definition fixed.

    Dual observations are shared between estimands and therefore use the same
    bootstrap draw. Selective and neither observations are resampled within
    their own classes. This estimates uncertainty of the descriptive AUROC
    difference; it is not a paired ligand-level test because the negative sets
    contain different compounds.
    """
    dual = np.asarray(dual, dtype=float)
    selective = np.asarray(selective, dtype=float)
    neither = np.asarray(neither, dtype=float)
    if min(len(dual), len(selective), len(neither)) == 0:
        return (float("nan"),) * 5
    rng = np.random.default_rng(seed)
    delta = []
    for _ in range(n_boot):
        d = rng.choice(dual, size=len(dual), replace=True)
        s = rng.choice(selective, size=len(selective), replace=True)
        n = rng.choice(neither, size=len(neither), replace=True)
        delta.append(auroc(d, n) - auroc(d, s))
    lo, hi = np.percentile(delta, [2.5, 97.5])
    auc_selective = auroc(dual, selective)
    auc_neither = auroc(dual, neither)
    return auc_selective, auc_neither, auc_neither - auc_selective, float(lo), float(hi)


def morgan(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    return m, fp


def murcko(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return None


def assemble(pair: str, cfg: dict) -> list[dict]:
    rows = load_csv(ROOT / cfg["scores"])
    smimap = {}
    if cfg["panel"]:
        for r in load_csv(ROOT / cfg["panel"]):
            smimap[r[cfg["panel_key"]]] = r.get("smiles")
    out = []
    for r in rows:
        a, b = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
        if a is None or b is None:
            continue
        lig = r.get("ligand") or r.get("panel_id")
        smi = r.get("smiles") or smimap.get(lig)
        if not smi:
            continue
        mol, fp = morgan(smi)
        if mol is None:
            continue
        rec = {
            "pair": pair,
            "ligand": lig,
            "cls": r.get("class"),
            "smiles": smi,
            "heavy": float(mol.GetNumHeavyAtoms()),
            "mw": float(Descriptors.MolWt(mol)),
            "clogp": float(Descriptors.MolLogP(mol)),
            "tpsa": float(Descriptors.TPSA(mol)),
            "vina_A": -a,
            "vina_B": -b,
            "vina_mean": -(a + b) / 2.0,
            "vina_worst": min(-a, -b),
            "fp": fp,
            "scaffold": murcko(smi) or f"__fail_{lig}",
        }
        out.append(rec)
    return out


def scores_of(recs, cls, key):
    return [r[key] for r in recs if r["cls"] == cls]


def add_contrast(rows, pair, formulation, contrast, score, pos, neg, n_pos, n_neg, note, underpowered=None):
    h = hashlib.md5(f"{pair}|{formulation}|{contrast}".encode()).hexdigest()
    pt, lo, hi = boot_auroc(pos, neg, seed=SEED + (int(h[:8], 16) % 99991))
    if underpowered is None:
        underpowered = int(min(n_pos, n_neg) < 8)
    rows.append(
        {
            "pair": pair,
            "formulation": formulation,
            "contrast": contrast,
            "score": score,
            "n_pos": n_pos,
            "n_neg": n_neg,
            "auroc": "" if pt != pt else round(pt, 4),
            "ci_lo": "" if lo != lo else round(lo, 4),
            "ci_hi": "" if hi != hi else round(hi, 4),
            "underpowered": underpowered,
            "note": note,
        }
    )


def formulation_table(packs):
    rows = []
    for pair, recs in packs.items():
        D = [r for r in recs if r["cls"] == "dual"]
        A = [r for r in recs if r["cls"] == "A_only"]
        B = [r for r in recs if r["cls"] == "B_only"]
        N = [r for r in recs if r["cls"] == "neither"]
        # DualFourClass primary (same definition as Table 2)
        add_contrast(
            rows, pair, "dualfourclass_directional", "D_vs_A_pocketB", "vina_B",
            [r["vina_B"] for r in D], [r["vina_B"] for r in A], len(D), len(A),
            "primary: dual vs A-only scored in pocket B",
        )
        add_contrast(
            rows, pair, "dualfourclass_directional", "D_vs_B_pocketA", "vina_A",
            [r["vina_A"] for r in D], [r["vina_A"] for r in B], len(D), len(B),
            "primary: dual vs B-only scored in pocket A",
        )
        da = auroc([r["vina_B"] for r in D], [r["vina_B"] for r in A])
        db = auroc([r["vina_A"] for r in D], [r["vina_A"] for r in B])
        rows.append(
            {
                "pair": pair,
                "formulation": "dualfourclass_directional",
                "contrast": "summary_min",
                "score": "min(D/A,D/B)",
                "n_pos": len(D),
                "n_neg": min(len(A), len(B)),
                "auroc": round(min(da, db), 4),
                "ci_lo": "",
                "ci_hi": "",
                "underpowered": 0,
                "note": "worst-arm aggregation; CI is on the two arms separately",
            }
        )
        # Conventional: dual vs experimental inactives (neither)
        add_contrast(
            rows, pair, "conventional_dual_vs_neither", "D_vs_neither_mean", "vina_mean",
            [r["vina_mean"] for r in D], [r["vina_mean"] for r in N], len(D), len(N),
            "Zhou-like dual vs noninhibitors using pooled score",
        )
        add_contrast(
            rows, pair, "conventional_dual_vs_neither", "D_vs_neither_worst", "vina_worst",
            [r["vina_worst"] for r in D], [r["vina_worst"] for r in N], len(D), len(N),
            "AND-like dual vs neither using min(pocket scores)",
        )
        # Dual vs all non-duals (selectives treated as negatives)
        nondual = A + B + N
        add_contrast(
            rows, pair, "dual_vs_all_nondual", "D_vs_A+B+neither_mean", "vina_mean",
            [r["vina_mean"] for r in D], [r["vina_mean"] for r in nondual], len(D), len(nondual),
            "selectives counted as negatives; still not directional",
        )
        # Single-target analogue (Zhou first stage)
        add_contrast(
            rows, pair, "single_target_analogue", "A_inhibitors_vs_A_noninhibitors", "vina_A",
            [r["vina_A"] for r in D + A], [r["vina_A"] for r in B + N], len(D) + len(A), len(B) + len(N),
            "pocket A: (dual+A-only) vs (B-only+neither)",
        )
        add_contrast(
            rows, pair, "single_target_analogue", "B_inhibitors_vs_B_noninhibitors", "vina_B",
            [r["vina_B"] for r in D + B], [r["vina_B"] for r in A + N], len(D) + len(B), len(A) + len(N),
            "pocket B: (dual+B-only) vs (A-only+neither)",
        )
    return rows


def equal_score_negative_table(packs):
    """Factorial formulation check: same pocket score, different negatives."""
    rows = []
    for pair, recs in packs.items():
        classes = {c: [r for r in recs if r["cls"] == c] for c in ("dual", "A_only", "B_only", "neither")}
        for contrast, selective_class, score_key, note in (
            ("D_vs_A_or_neither_pocketB", "A_only", "vina_B", "pocket B fixed; A-only versus neither negative class"),
            ("D_vs_B_or_neither_pocketA", "B_only", "vina_A", "pocket A fixed; B-only versus neither negative class"),
        ):
            D = [r[score_key] for r in classes["dual"]]
            S = [r[score_key] for r in classes[selective_class]]
            N = [r[score_key] for r in classes["neither"]]
            h = hashlib.md5(f"equal-score|{pair}|{contrast}".encode()).hexdigest()
            a_s, a_n, delta, lo, hi = boot_equal_score_negative_contrast(
                D, S, N, seed=SEED + (int(h[:8], 16) % 99991)
            )
            rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "score": score_key,
                    "n_dual": len(D),
                    "n_selective": len(S),
                    "n_neither": len(N),
                    "auroc_dual_vs_selective": round(a_s, 4),
                    "auroc_dual_vs_neither": round(a_n, 4),
                    "delta_neither_minus_selective": round(delta, 4),
                    "delta_ci_lo": round(lo, 4),
                    "delta_ci_hi": round(hi, 4),
                    "underpowered_neither": int(len(N) < 8),
                    "note": note + "; joint stratified bootstrap, not paired negative ligands",
                }
            )
    return rows


def chemotype_table(packs):
    rows = []
    for pair, recs in packs.items():
        D = [r for r in recs if r["cls"] == "dual"]
        for arm, neg_cls, score_key, t_grid in (
            ("D_vs_A", "A_only", "vina_B", (0.3, 0.4, 0.5)),
            ("D_vs_B", "B_only", "vina_A", (0.3, 0.4, 0.5)),
        ):
            Neg = [r for r in recs if r["cls"] == neg_cls]
            if not D or not Neg:
                continue
            # nearest-neighbor Tanimoto of each dual to the hard-neg pool
            nn = []
            for d in D:
                sims = [DataStructs.TanimotoSimilarity(d["fp"], n["fp"]) for n in Neg]
                nn.append(max(sims) if sims else 0.0)
            # hard-neg nearest dual
            nn_neg = []
            for n in Neg:
                sims = [DataStructs.TanimotoSimilarity(n["fp"], d["fp"]) for d in D]
                nn_neg.append(max(sims) if sims else 0.0)
            rows.append(
                {
                    "pair": pair,
                    "contrast": arm,
                    "subset": "all",
                    "tanimoto_cut": "",
                    "n_pos": len(D),
                    "n_neg": len(Neg),
                    "median_nn_tanimoto_dual": round(float(np.median(nn)), 4),
                    "median_nn_tanimoto_neg": round(float(np.median(nn_neg)), 4),
                    "auroc": round(auroc([d[score_key] for d in D], [n[score_key] for n in Neg]), 4),
                    "note": "unmatched directional contrast (Table 2 arm)",
                }
            )
            for t in t_grid:
                matched_neg = [n for n, s in zip(Neg, nn_neg) if s >= t]
                if not matched_neg:
                    rows.append(
                        {
                            "pair": pair,
                            "contrast": arm,
                            "subset": "chemotype_matched",
                            "tanimoto_cut": t,
                            "n_pos": len(D),
                            "n_neg": 0,
                            "median_nn_tanimoto_dual": "",
                            "median_nn_tanimoto_neg": "",
                            "auroc": "",
                            "note": "no hard-neg with Tanimoto>=cut to any dual",
                        }
                    )
                    continue
                rows.append(
                    {
                        "pair": pair,
                        "contrast": arm,
                        "subset": "chemotype_matched",
                        "tanimoto_cut": t,
                        "n_pos": len(D),
                        "n_neg": len(matched_neg),
                        "median_nn_tanimoto_dual": round(float(np.median(nn)), 4),
                        "median_nn_tanimoto_neg": round(float(np.median([s for s in nn_neg if s >= t])), 4),
                        "auroc": round(auroc([d[score_key] for d in D], [n[score_key] for n in matched_neg]), 4),
                        "note": "hard-negs with NN Tanimoto to dual >= cut",
                    }
                )
            distant = [n for n, s in zip(Neg, nn_neg) if s < 0.3]
            if distant:
                rows.append(
                    {
                        "pair": pair,
                        "contrast": arm,
                        "subset": "chemotype_distant",
                        "tanimoto_cut": 0.3,
                        "n_pos": len(D),
                        "n_neg": len(distant),
                        "median_nn_tanimoto_dual": round(float(np.median(nn)), 4),
                        "median_nn_tanimoto_neg": round(float(np.median([s for s in nn_neg if s < 0.3])), 4),
                        "auroc": round(auroc([d[score_key] for d in D], [n[score_key] for n in distant]), 4),
                        "note": "hard-negs with NN Tanimoto to dual < 0.3",
                    }
                )
    return rows


def _cv_auroc(X, y, groups):
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    if n_pos < 6 or n_neg < 6 or len(y) < 16:
        return float("nan"), 0
    n_scaf = len(set(groups))
    n_splits = min(5, n_scaf, n_pos, n_neg)
    if n_splits < 2:
        return float("nan"), 0
    cv = GroupKFold(n_splits=n_splits)
    lr = LogisticRegression(max_iter=4000, C=1.0)
    try:
        prob = cross_val_predict(lr, X, y, cv=cv, groups=groups, method="predict_proba")[:, 1]
        return float(roc_auc_score(y, prob)), n_splits
    except Exception:
        return float("nan"), n_splits


def incremental_table(packs):
    rows = []
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, dock_key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            if len(kept) < 16:
                continue
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
            groups = np.array([r["scaffold"] for r in kept])
            phys = np.array([[r["heavy"], r["mw"], r["clogp"], r["tpsa"]] for r in kept], dtype=float)
            dock = np.array([[r[dock_key]] for r in kept], dtype=float)
            fp = np.vstack([np.asarray(r["fp"]) for r in kept])
            models = [
                ("physchem", phys),
                ("ECFP4", fp),
                ("docking", dock),
                ("docking+physchem", np.hstack([dock, phys])),
                ("ECFP4+docking", np.hstack([fp, dock])),
                ("ECFP4+physchem", np.hstack([fp, phys])),
                ("ECFP4+docking+physchem", np.hstack([fp, dock, phys])),
            ]
            rank_dock = auroc(
                [r[dock_key] for r in kept if r["cls"] == pos_cls],
                [r[dock_key] for r in kept if r["cls"] == neg_cls],
            )
            for name, X in models:
                auc, n_splits = _cv_auroc(X, y, groups)
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "model": name,
                        "n": len(kept),
                        "n_pos": int(y.sum()),
                        "n_neg": int((1 - y).sum()),
                        "n_scaffolds": len(set(groups)),
                        "n_splits": n_splits,
                        "cv_auroc": "" if auc != auc else round(auc, 4),
                        "rank_auroc_docking": round(rank_dock, 4),
                        "note": "scaffold GroupKFold logistic; docking rank AUROC is not the logistic AUROC",
                    }
                )
    return rows


def enrichment_table(packs):
    """Mixed four-state library ranked by a single score; dual is the hit class."""
    rows = []
    for pair, recs in packs.items():
        n = len(recs)
        n_dual = sum(r["cls"] == "dual" for r in recs)
        if n_dual == 0:
            continue
        random_hit = n_dual / n
        for score_name in ("vina_mean", "vina_worst", "vina_A", "vina_B"):
            ranked = sorted(recs, key=lambda r: r[score_name], reverse=True)
            for frac, lab in ((0.05, "EF5"), (0.10, "EF10"), (None, "Top10")):
                if lab == "Top10":
                    k = min(10, n)
                else:
                    k = max(1, int(round(frac * n)))
                top = ranked[:k]
                n_hit = sum(r["cls"] == "dual" for r in top)
                counts = Counter(r["cls"] for r in top)
                ef = (n_hit / k) / random_hit if random_hit > 0 else float("nan")
                rows.append(
                    {
                        "pair": pair,
                        "score": score_name,
                        "cutoff": lab,
                        "k": k,
                        "n_library": n,
                        "n_dual_library": n_dual,
                        "n_dual_top": n_hit,
                        "n_A_only_top": counts.get("A_only", 0),
                        "n_B_only_top": counts.get("B_only", 0),
                        "n_neither_top": counts.get("neither", 0),
                        "hit_rate": round(n_hit / k, 4),
                        "EF": round(ef, 3),
                        "hardneg_fraction_top": round((counts.get("A_only", 0) + counts.get("B_only", 0)) / k, 4),
                    }
                )
    return rows


def main():
    packs = {pair: assemble(pair, cfg) for pair, cfg in SPEC.items()}
    for pair, recs in packs.items():
        print(pair, Counter(r["cls"] for r in recs), "n=", len(recs))

    form = formulation_table(packs)
    equal_score = equal_score_negative_table(packs)
    chemo = chemotype_table(packs)
    incr = incremental_table(packs)
    enr = enrichment_table(packs)
    write_csv(TAB / "formulation_conventional_vs_directional_v1.csv", form)
    write_csv(TAB / "formulation_equal_score_negative_v1.csv", equal_score)
    write_csv(TAB / "chemotype_matched_hardneg_v1.csv", chemo)
    write_csv(TAB / "incremental_information_v1.csv", incr)
    write_csv(TAB / "mixed_library_enrichment_v1.csv", enr)

    # compact pairwise summary for the novelty claim
    summary = []
    for pair in SPEC:
        def grab(formu, contrast):
            hit = next((r for r in form if r["pair"] == pair and r["formulation"] == formu and r["contrast"] == contrast), None)
            return hit

        sm = grab("dualfourclass_directional", "summary_min")
        dn = grab("conventional_dual_vs_neither", "D_vs_neither_mean")
        dw = grab("conventional_dual_vs_neither", "D_vs_neither_worst")
        stA = grab("single_target_analogue", "A_inhibitors_vs_A_noninhibitors")
        stB = grab("single_target_analogue", "B_inhibitors_vs_B_noninhibitors")
        summary.append(
            {
                "pair": pair,
                "directional_summary_min": sm["auroc"] if sm else "",
                "conventional_D_vs_neither_mean": dn["auroc"] if dn else "",
                "conventional_D_vs_neither_mean_n_neg": dn["n_neg"] if dn else "",
                "conventional_D_vs_neither_mean_underpowered": dn["underpowered"] if dn else "",
                "conventional_D_vs_neither_worst": dw["auroc"] if dw else "",
                "single_target_pocketA": stA["auroc"] if stA else "",
                "single_target_pocketB": stB["auroc"] if stB else "",
                "delta_neither_mean_minus_summary_min": (
                    round(float(dn["auroc"]) - float(sm["auroc"]), 4)
                    if dn and sm and dn["auroc"] != "" and sm["auroc"] != ""
                    else ""
                ),
            }
        )
    write_csv(TAB / "formulation_summary_v1.csv", summary)

    meta = {
        "n_boot": N_BOOT,
        "seed": SEED,
        "no_new_docking": True,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": importlib.metadata.version("numpy"),
            "rdkit": importlib.metadata.version("rdkit"),
            "scipy": importlib.metadata.version("scipy"),
            "scikit_learn": importlib.metadata.version("scikit-learn"),
            "pandas": importlib.metadata.version("pandas"),
        },
        "zhou_2013": "10.1021/ci400065e",
        "pairs": {p: dict(Counter(r["cls"] for r in recs)) for p, recs in packs.items()},
        "summary": summary,
    }
    (OUT / "analysis" ).mkdir(exist_ok=True)
    (OUT / "analysis" / "formulation_run_meta_v1.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
