#!/usr/bin/env python3
"""Single-descriptor baselines from current score master (RDKit descriptors + scheme-B).

Two analyses are written and must not be conflated:

A. Full-panel descriptive univariate screen (TPSA, cLogP, heavy, MW).
   This is not a selection-adjusted predictive estimate.

B. Train-only nested scaffold-GroupKFold predictive baseline:
   inner CV selects one prespecified descriptor on outer-training data only;
   a single-variable logistic regression (StandardScaler fit on outer-train)
   emits held-out P(dual); pooled AUROC uses those probabilities.

Does not modify the ECFP4 pipeline. Does not modify docking scores.
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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    DESCRIPTOR_NAMES,
    GROUPKFOLD_MAX_SPLITS,
    LOGREG_C,
    LOGREG_MAX_ITER,
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
ARMS = (
    ("D_vs_A", "dual", "A_only"),
    ("D_vs_B", "dual", "B_only"),
)


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


def _lr_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "lr",
                LogisticRegression(
                    max_iter=LOGREG_MAX_ITER,
                    C=LOGREG_C,
                    random_state=SEED,
                    solver="lbfgs",
                ),
            ),
        ]
    )


def _inner_mean_auroc(x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> float:
    """Mean inner-validation AUROC of a single-variable scaled logistic model."""
    n_scaf = len(set(groups.tolist()))
    n_pos = int(y.sum())
    n_neg = int((1 - y).sum())
    n_splits = min(GROUPKFOLD_MAX_SPLITS, n_scaf, n_pos, n_neg)
    if n_splits < 2:
        model = _lr_pipeline()
        model.fit(x.reshape(-1, 1), y)
        prob = model.predict_proba(x.reshape(-1, 1))[:, 1]
        return float(auroc(prob[y == 1], prob[y == 0]))
    cv = GroupKFold(n_splits=n_splits)
    dummy = np.zeros((len(y), 1))
    aucs: list[float] = []
    for tr, va in cv.split(dummy, y, groups):
        if y[tr].sum() == 0 or (1 - y[tr]).sum() == 0:
            continue
        if y[va].sum() == 0 or (1 - y[va]).sum() == 0:
            continue
        model = _lr_pipeline()
        model.fit(x[tr].reshape(-1, 1), y[tr])
        prob = model.predict_proba(x[va].reshape(-1, 1))[:, 1]
        aucs.append(auroc(prob[y[va] == 1], prob[y[va] == 0]))
    if not aucs:
        model = _lr_pipeline()
        model.fit(x.reshape(-1, 1), y)
        prob = model.predict_proba(x.reshape(-1, 1))[:, 1]
        return float(auroc(prob[y == 1], prob[y == 0]))
    return float(np.mean(aucs))


def select_descriptor(train_recs: list[dict], y: np.ndarray, groups: np.ndarray) -> str:
    scores = {}
    for name in DESCRIPTOR_NAMES:
        x = np.array([r[name] for r in train_recs], dtype=float)
        scores[name] = _inner_mean_auroc(x, y, groups)
    return max(scores, key=scores.get)


def nested_arm(recs: list[dict], pair: str, arm: str, pos_cls: str, neg_cls: str):
    """Independent nested CV for one pair × arm.

    Outer GroupKFold is on contrast ligands only (dual vs the arm negative).
    Inner GroupKFold on outer-training scaffolds selects one descriptor.
    Outer-train scaler + logistic regression emit held-out P(dual).
    """
    kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
    y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
    groups = np.array([r["scaffold"] for r in kept])
    n_pos = int(y.sum())
    n_neg = int((1 - y).sum())
    n_scaf = len(set(groups.tolist()))
    n_splits = min(GROUPKFOLD_MAX_SPLITS, n_scaf, n_pos, n_neg)
    if n_splits < 2 or n_pos < 6 or n_neg < 6:
        return None
    cv = GroupKFold(n_splits=n_splits)
    dummy = np.zeros((len(kept), 1))
    oof = np.full(len(kept), np.nan)
    fold_id = np.full(len(kept), -1, dtype=int)
    fold_rows = []
    oof_rows = []
    selected_by_fold: list[str] = []
    seen_keys: set[tuple[str, str, int]] = set()
    for fold, (tr, te) in enumerate(cv.split(dummy, y, groups)):
        key = (pair, arm, fold)
        if key in seen_keys:
            raise SystemExit(f"FAIL: duplicate selection key {key}")
        seen_keys.add(key)
        tr_ids = {kept[i]["ligand_id"] for i in tr}
        te_ids = {kept[i]["ligand_id"] for i in te}
        if tr_ids & te_ids:
            raise SystemExit(f"FAIL: train/test ID leak {pair} {arm} fold={fold}")
        tr_scaf = {groups[i] for i in tr}
        te_scaf = {groups[i] for i in te}
        if tr_scaf & te_scaf:
            raise SystemExit(f"FAIL: scaffold leak {pair} {arm} fold={fold}")
        train_recs = [kept[i] for i in tr]
        y_tr = y[tr]
        g_tr = groups[tr]
        chosen = select_descriptor(train_recs, y_tr, g_tr)
        x_all = np.array([r[chosen] for r in kept], dtype=float)
        model = _lr_pipeline()
        model.fit(x_all[tr].reshape(-1, 1), y_tr)
        scaler = model.named_steps["scaler"]
        lr = model.named_steps["lr"]
        coef = float(lr.coef_.ravel()[0])
        intercept = float(lr.intercept_.ravel()[0])
        scaled_te = scaler.transform(x_all[te].reshape(-1, 1)).ravel()
        prob_te = model.predict_proba(x_all[te].reshape(-1, 1))[:, 1]
        if np.any((prob_te < 0.0) | (prob_te > 1.0)):
            raise SystemExit(f"FAIL: OOF probability outside [0,1] {pair} {arm} fold={fold}")
        oof[te] = prob_te
        fold_id[te] = fold
        selected_by_fold.append(chosen)
        if chosen is None:
            raise SystemExit(f"FAIL: empty descriptor selection {pair} {arm} fold={fold}")
        fold_rows.append(
            {
                "pair": pair,
                "arm": arm,
                "outer_fold": fold,
                "selected_descriptor": chosen,
                "selection_source_pair": pair,
                "selection_source_arm": arm,
                "selection_source_fold": fold,
                "test_descriptor_used": chosen,
                "n_train": int(len(tr)),
                "n_test": int(len(te)),
                "n_pos_train": int(y_tr.sum()),
                "n_neg_train": int((1 - y_tr).sum()),
                "n_pos": n_pos,
                "n_neg": n_neg,
                "training_coefficient": f"{coef:.8f}",
                "training_intercept": f"{intercept:.8f}",
            }
        )
        if fold_rows[-1]["selection_source_arm"] != arm:
            raise SystemExit(f"FAIL: selection_source_arm != arm for {pair} {arm} fold={fold}")
        for j, idx in enumerate(te):
            r = kept[idx]
            oof_rows.append(
                {
                    "pair": pair,
                    "arm": arm,
                    "ligand_id": r["ligand_id"],
                    "class": r["cls"],
                    "outer_fold": fold,
                    "scaffold": r["scaffold"],
                    "selected_descriptor": chosen,
                    "selection_source_pair": pair,
                    "selection_source_arm": arm,
                    "selection_source_fold": fold,
                    "raw_descriptor_value": f"{x_all[idx]:.8f}",
                    "training_coefficient": f"{coef:.8f}",
                    "training_intercept": f"{intercept:.8f}",
                    "scaled_test_value": f"{float(scaled_te[j]):.12f}",
                    "oof_probability": f"{float(prob_te[j]):.12f}",
                    "y": int(y[idx]),
                    "n_train": int(len(tr)),
                    "n_test": int(len(te)),
                }
            )
    if np.isnan(oof).any():
        return None
    # AUROC from the same rounded probabilities that are deposited.
    deposited = np.array([float(r["oof_probability"]) for r in oof_rows], dtype=float)
    y_dep = np.array([int(r["y"]) for r in oof_rows], dtype=int)
    auc = auroc(deposited[y_dep == 1], deposited[y_dep == 0])
    oof = deposited
    y = y_dep
    # kept order must match oof_rows for bootstrap ID alignment
    kept = [next(r for r in kept if r["ligand_id"] == row["ligand_id"]) for row in oof_rows]
    return {
        "pair": pair,
        "arm": arm,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "outer_folds": n_splits,
        "oof_auroc": float(auc),
        "descriptors_selected_by_fold": ",".join(selected_by_fold),
        "number_of_unique_selected_descriptors": len(set(selected_by_fold)),
        "mode_descriptor": Counter(selected_by_fold).most_common(1)[0][0],
        "oof": oof,
        "y": y,
        "kept": kept,
        "fold_id": fold_id,
        "fold_rows": fold_rows,
        "oof_rows": oof_rows,
    }


def shared_dual_oof_bootstrap(arm_a: dict, arm_b: dict) -> dict:
    """Ligand-level shared-dual bootstrap of fixed OOF probabilities.

    This is not a re-run of nested model selection.
    """
    dual_a = {r["ligand_id"]: float(p) for r, p in zip(arm_a["kept"], arm_a["oof"]) if r["cls"] == "dual"}
    dual_b = {r["ligand_id"]: float(p) for r, p in zip(arm_b["kept"], arm_b["oof"]) if r["cls"] == "dual"}
    dual_ids = np.array(sorted(set(dual_a) & set(dual_b)))
    a_only = sorted(
        ((r["ligand_id"], float(p)) for r, p in zip(arm_a["kept"], arm_a["oof"]) if r["cls"] == "A_only"),
        key=lambda t: t[0],
    )
    b_only = sorted(
        ((r["ligand_id"], float(p)) for r, p in zip(arm_b["kept"], arm_b["oof"]) if r["cls"] == "B_only"),
        key=lambda t: t[0],
    )
    da = auroc([dual_a[i] for i in dual_ids], [p for _, p in a_only])
    db = auroc([dual_b[i] for i in dual_ids], [p for _, p in b_only])
    point = min(da, db)
    rng = np.random.default_rng(SEED)
    das = np.empty(N_BOOT, dtype=float)
    dbs = np.empty(N_BOOT, dtype=float)
    mins = np.empty(N_BOOT, dtype=float)
    a_probs = np.array([p for _, p in a_only], dtype=float)
    b_probs = np.array([p for _, p in b_only], dtype=float)
    dual_pa = np.array([dual_a[i] for i in dual_ids], dtype=float)
    dual_pb = np.array([dual_b[i] for i in dual_ids], dtype=float)
    for i in range(N_BOOT):
        di = rng.choice(dual_ids.size, size=dual_ids.size, replace=True)
        ai = rng.choice(a_probs.size, size=a_probs.size, replace=True)
        bi = rng.choice(b_probs.size, size=b_probs.size, replace=True)
        das[i] = auroc(dual_pa[di], a_probs[ai])
        dbs[i] = auroc(dual_pb[di], b_probs[bi])
        mins[i] = min(das[i], dbs[i])
    da_lo, da_hi = percentile_ci(das, N_BOOT)
    db_lo, db_hi = percentile_ci(dbs, N_BOOT)
    sm_lo, sm_hi = percentile_ci(mins, N_BOOT)
    return {
        "oof_D_vs_A": float(da),
        "oof_D_vs_A_ci_lo": da_lo,
        "oof_D_vs_A_ci_hi": da_hi,
        "oof_D_vs_B": float(db),
        "oof_D_vs_B_ci_lo": db_lo,
        "oof_D_vs_B_ci_hi": db_hi,
        "oof_summary_min": float(point),
        "oof_summary_min_ci_lo": sm_lo,
        "oof_summary_min_ci_hi": sm_hi,
        "n_shared_dual": int(dual_ids.size),
        "bootstrap_note": "ligand_level_shared_dual_on_fixed_oof_probabilities; not selection-adjusted model uncertainty",
    }


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
    oof_out = []
    arm_summary = []
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

        arm_models = {}
        for arm, pos_cls, neg_cls in ARMS:
            fitted = nested_arm(recs, pair, arm, pos_cls, neg_cls)
            if fitted is None:
                raise SystemExit(f"FAIL: nested descriptor CV produced no OOF for {pair} {arm}")
            arm_models[arm] = fitted
            fold_out.extend(fitted["fold_rows"])
            oof_out.extend(fitted["oof_rows"])
            arm_summary.append(
                {
                    "pair": pair,
                    "arm": arm,
                    "n_pos": fitted["n_pos"],
                    "n_neg": fitted["n_neg"],
                    "outer_folds": fitted["outer_folds"],
                    "oof_auroc": f"{fitted['oof_auroc']:.4f}",
                    "descriptors_selected_by_fold": fitted["descriptors_selected_by_fold"],
                    "number_of_unique_selected_descriptors": fitted["number_of_unique_selected_descriptors"],
                    "mode_descriptor": fitted["mode_descriptor"],
                    "score_scale": "oof_probability",
                    "selection": "inner_scaffold_groupkfold_on_outer_train",
                    "model": "standardscaler_plus_univariate_logistic_regression",
                }
            )

        nested = shared_dual_oof_bootstrap(arm_models["D_vs_A"], arm_models["D_vs_B"])
        if abs(nested["oof_D_vs_A"] - arm_models["D_vs_A"]["oof_auroc"]) > 1e-10:
            raise SystemExit(f"FAIL: D_vs_A OOF AUROC mismatch {pair}")
        if abs(nested["oof_D_vs_B"] - arm_models["D_vs_B"]["oof_auroc"]) > 1e-10:
            raise SystemExit(f"FAIL: D_vs_B OOF AUROC mismatch {pair}")
        if abs(nested["oof_summary_min"] - min(nested["oof_D_vs_A"], nested["oof_D_vs_B"])) > 1e-12:
            raise SystemExit(f"FAIL: summary_min != min(directional) {pair}")

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
            "nested_scaffold_cv_oof_summary_min": f"{nested['oof_summary_min']:.4f}",
            "nested_scaffold_cv_oof_summary_min_ci_lo": f"{nested['oof_summary_min_ci_lo']:.4f}",
            "nested_scaffold_cv_oof_summary_min_ci_hi": f"{nested['oof_summary_min_ci_hi']:.4f}",
            "nested_scaffold_cv_oof_D_vs_A": f"{nested['oof_D_vs_A']:.4f}",
            "nested_scaffold_cv_oof_D_vs_A_ci_lo": f"{nested['oof_D_vs_A_ci_lo']:.4f}",
            "nested_scaffold_cv_oof_D_vs_A_ci_hi": f"{nested['oof_D_vs_A_ci_hi']:.4f}",
            "nested_scaffold_cv_oof_D_vs_B": f"{nested['oof_D_vs_B']:.4f}",
            "nested_scaffold_cv_oof_D_vs_B_ci_lo": f"{nested['oof_D_vs_B_ci_lo']:.4f}",
            "nested_scaffold_cv_oof_D_vs_B_ci_hi": f"{nested['oof_D_vs_B_ci_hi']:.4f}",
            "nested_selected_D_vs_A_mode": arm_models["D_vs_A"]["mode_descriptor"],
            "nested_selected_D_vs_B_mode": arm_models["D_vs_B"]["mode_descriptor"],
            "nested_n_splits": arm_models["D_vs_A"]["outer_folds"],
            "nested_score_scale": "oof_probability",
            "nested_ci_definition": nested["bootstrap_note"],
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
        print(
            f"{pair:14} vina={rec['vina_summary_min']} best={best} {rec['best_descriptor_summary_min']} "
            f"Δ={rec['vina_minus_best_descriptor']} nested_oof={rec['nested_scaffold_cv_oof_summary_min']} "
            f"[{rec['nested_scaffold_cv_oof_summary_min_ci_lo']}, {rec['nested_scaffold_cv_oof_summary_min_ci_hi']}]"
        )

    # Merge pair×arm summary fields onto the 80 fold-unit table requested by the audit.
    arm_by = {(r["pair"], r["arm"]): r for r in arm_summary}
    for fr in fold_out:
        src = arm_by[(fr["pair"], fr["arm"])]
        fr["oof_auroc"] = src["oof_auroc"]
        fr["descriptors_selected_by_fold"] = src["descriptors_selected_by_fold"]
        fr["number_of_unique_selected_descriptors"] = src["number_of_unique_selected_descriptors"]
        fr["outer_folds"] = src["outer_folds"]

    write_csv(CANON / "descriptor_baselines.csv", rows)
    write_csv(CANON / "descriptor_nested_scaffold_cv.csv", fold_out)
    write_csv(CANON / "descriptor_nested_oof_predictions.csv", oof_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
