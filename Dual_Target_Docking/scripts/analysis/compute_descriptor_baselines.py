#!/usr/bin/env python3
"""Single-descriptor baselines from current score master (RDKit descriptors + scheme-B).

Full-panel best-descriptor summary_min is a descriptive univariate screen, not a
selection-adjusted predictive estimate. Nested scaffold-GroupKFold train-only
descriptor selection supplies the OOF AUROC used as the chemistry-control
predictive number. Does not modify the ECFP4 pipeline.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    DESCRIPTOR_NAMES,
    GROUPKFOLD_MAX_SPLITS,
    PRIMARY_PAIRS,
    add_io_args,
    io_paths,
    is_primary_row,
)
from analysis.bootstrap_metrics import (  # noqa: E402
    N_BOOT,
    SEED,
    auroc,
    percentile_ci,
    summary_min_stratified,
)

CANON = ROOT / "results" / "canonical"
MASTER = CANON / "current_score_master.csv"
GETTERS = {
    "tpsa": lambda mol: float(Descriptors.TPSA(mol)),
    "clogp": lambda mol: float(Descriptors.MolLogP(mol)),
    "heavy": lambda mol: float(Descriptors.HeavyAtomCount(mol)),
    "mw": lambda mol: float(Descriptors.MolWt(mol)),
}


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def scaffold_of(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi or "")
    if mol is None:
        return "invalid"
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "acyclic"
    except Exception:
        return "acyclic"


def load():
    by = {p: [] for p in PRIMARY_PAIRS}
    for r in read_csv(MASTER):
        if not is_primary_row(r):
            continue
        mol = Chem.MolFromSmiles(r.get("smiles") or "")
        if mol is None:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["scaffold"] = scaffold_of(r.get("smiles") or "")
        for name in DESCRIPTOR_NAMES:
            rec[name] = GETTERS[name](mol)
        by[r["pair"]].append(rec)
    return by


def directional_smin(vals, cls):
    d = cls == "dual"
    a = cls == "A_only"
    b = cls == "B_only"
    da = auroc(vals[d], vals[a])
    db = auroc(vals[d], vals[b])
    return da, db, min(da, db)


def nested_oof(recs):
    """Train-only descriptor selection under scaffold GroupKFold.

    Contrasts are selected separately (D vs A on pocket-B analogue = the
    descriptor itself; D vs B likewise). OOF scores are concatenated within
    each contrast. Mixed descriptor scales across folds are accepted as the
    nested-selection readout; fold-wise AUROCs are also stored.
    """
    cls = np.array([r["cls"] for r in recs])
    groups = np.array([r["scaffold"] for r in recs])
    fold_rows = []
    oof = {name: np.full(len(recs), np.nan) for name in ("selected_DvsA", "selected_DvsB")}
    selected_da = []
    selected_db = []
    n_pos = int((cls == "dual").sum())
    # Use all recs; GroupKFold needs >=2 groups and both classes in contrast.
    n_scaf = len(set(groups.tolist()))
    n_splits = min(GROUPKFOLD_MAX_SPLITS, n_scaf)
    if n_splits < 2 or n_pos < 6:
        return None, []
    cv = GroupKFold(n_splits=n_splits)
    dummy_x = np.zeros((len(recs), 1))
    dummy_y = (cls == "dual").astype(int)
    fold_id = np.full(len(recs), -1, dtype=int)
    for fold, (tr, te) in enumerate(cv.split(dummy_x, dummy_y, groups)):
        fold_id[te] = fold
        tr_cls = cls[tr]
        te_cls = cls[te]
        best_da_name, best_da = None, -np.inf
        best_db_name, best_db = None, -np.inf
        for name in DESCRIPTOR_NAMES:
            val = np.array([recs[i][name] for i in range(len(recs))], dtype=float)
            da, db, _ = directional_smin(val[tr], tr_cls)
            if da > best_da:
                best_da, best_da_name = da, name
            if db > best_db:
                best_db, best_db_name = db, name
        if best_da_name is None or best_db_name is None:
            continue
        val_da = np.array([recs[i][best_da_name] for i in range(len(recs))], dtype=float)
        val_db = np.array([recs[i][best_db_name] for i in range(len(recs))], dtype=float)
        oof["selected_DvsA"][te] = val_da[te]
        oof["selected_DvsB"][te] = val_db[te]
        te_da, te_db, te_sm = directional_smin(val_da[te], te_cls) if (te_cls == "dual").any() else (np.nan, np.nan, np.nan)
        selected_da.append(best_da_name)
        selected_db.append(best_db_name)
        fold_rows.append(
            {
                "fold": fold,
                "n_train": int(len(tr)),
                "n_test": int(len(te)),
                "selected_D_vs_A": best_da_name,
                "selected_D_vs_B": best_db_name,
                "train_D_vs_A": f"{best_da:.4f}",
                "train_D_vs_B": f"{best_db:.4f}",
                "test_D_vs_A": "" if te_da != te_da else f"{te_da:.4f}",
                "test_D_vs_B": "" if te_db != te_db else f"{te_db:.4f}",
                "test_summary_min": "" if te_sm != te_sm else f"{te_sm:.4f}",
            }
        )
    mask_a = ~np.isnan(oof["selected_DvsA"])
    mask_b = ~np.isnan(oof["selected_DvsB"])
    da = auroc(oof["selected_DvsA"][mask_a & (cls == "dual")], oof["selected_DvsA"][mask_a & (cls == "A_only")])
    db = auroc(oof["selected_DvsB"][mask_b & (cls == "dual")], oof["selected_DvsB"][mask_b & (cls == "B_only")])
    return {
        "n_splits": n_splits,
        "oof_D_vs_A": da,
        "oof_D_vs_B": db,
        "oof_summary_min": min(da, db),
        "mode_D_vs_A": Counter(selected_da).most_common(1)[0][0] if selected_da else "",
        "mode_D_vs_B": Counter(selected_db).most_common(1)[0][0] if selected_db else "",
    }, fold_rows


def main() -> int:
    import argparse

    global CANON, MASTER
    parser = argparse.ArgumentParser(description="Single-descriptor baselines from current score master.")
    add_io_args(parser)
    args = parser.parse_args()
    CANON, MASTER = io_paths(args.outdir, args.master)
    CANON.mkdir(parents=True, exist_ok=True)

    packs = load()
    rows = []
    fold_out = []
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        cls = np.array([r["cls"] for r in recs])
        vina = summary_min_stratified(sa, sb, cls, n_boot=1, seed=SEED)
        desc_smin = {}
        desc_arms = {}
        for name in DESCRIPTOR_NAMES:
            val = np.array([r[name] for r in recs], dtype=float)
            da, db, sm = directional_smin(val, cls)
            desc_smin[name] = sm
            desc_arms[name] = (da, db)
        best = max(desc_smin, key=desc_smin.get)
        d_idx = np.flatnonzero(cls == "dual")
        a_idx = np.flatnonzero(cls == "A_only")
        b_idx = np.flatnonzero(cls == "B_only")
        rng = np.random.default_rng(SEED)
        deltas = []
        best_val = np.array([r[best] for r in recs], dtype=float)
        for _ in range(N_BOOT):
            di = rng.choice(d_idx, size=d_idx.size, replace=True)
            ai = rng.choice(a_idx, size=a_idx.size, replace=True)
            bi = rng.choice(b_idx, size=b_idx.size, replace=True)
            v = min(auroc(sb[di], sb[ai]), auroc(sa[di], sa[bi]))
            u = min(auroc(best_val[di], best_val[ai]), auroc(best_val[di], best_val[bi]))
            deltas.append(v - u)
        lo, hi = percentile_ci(deltas, N_BOOT)
        point = vina["summary_min"] - desc_smin[best]
        nested, fold_rows = nested_oof(recs)
        rec = {
            "pair": pair,
            "n_dual": vina["n_dual"],
            "n_A_only": vina["n_A_only"],
            "n_B_only": vina["n_B_only"],
            "vina_summary_min": f"{vina['summary_min']:.4f}",
            "vina_D_vs_A": f"{vina['auroc_D_vs_A_pocketB']:.4f}",
            "vina_D_vs_B": f"{vina['auroc_D_vs_B_pocketA']:.4f}",
            "best_descriptor": best,
            "best_descriptor_selection": "full_panel_descriptive_univariate_screen",
            "best_descriptor_summary_min": f"{desc_smin[best]:.4f}",
            "vina_minus_best_descriptor": f"{point:.4f}",
            "delta_ci_lo": f"{lo:.4f}",
            "delta_ci_hi": f"{hi:.4f}",
            "nested_scaffold_cv_oof_summary_min": "" if nested is None else f"{nested['oof_summary_min']:.4f}",
            "nested_scaffold_cv_oof_D_vs_A": "" if nested is None else f"{nested['oof_D_vs_A']:.4f}",
            "nested_scaffold_cv_oof_D_vs_B": "" if nested is None else f"{nested['oof_D_vs_B']:.4f}",
            "nested_selected_D_vs_A_mode": "" if nested is None else nested["mode_D_vs_A"],
            "nested_selected_D_vs_B_mode": "" if nested is None else nested["mode_D_vs_B"],
            "nested_n_splits": 0 if nested is None else nested["n_splits"],
            "predictive_descriptor_estimate": "nested_scaffold_cv_oof_summary_min",
            "bootstrap": "class_stratified_shared_dual",
            "B": N_BOOT,
            "seed": SEED,
        }
        for name in DESCRIPTOR_NAMES:
            rec[f"{name}_D_vs_A"] = f"{desc_arms[name][0]:.4f}"
            rec[f"{name}_D_vs_B"] = f"{desc_arms[name][1]:.4f}"
            rec[f"{name}_summary_min"] = f"{desc_smin[name]:.4f}"
        rows.append(rec)
        for fr in fold_rows:
            fold_out.append({"pair": pair, **fr})
        nested_txt = rec["nested_scaffold_cv_oof_summary_min"]
        print(
            f"{pair:14} vina={rec['vina_summary_min']} best={best} {rec['best_descriptor_summary_min']} "
            f"Δ={rec['vina_minus_best_descriptor']} nested_oof={nested_txt}"
        )
    write_csv(CANON / "descriptor_baselines.csv", rows)
    if fold_out:
        write_csv(CANON / "descriptor_nested_scaffold_cv.csv", fold_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
