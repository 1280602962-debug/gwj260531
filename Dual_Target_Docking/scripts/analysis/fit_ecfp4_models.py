#!/usr/bin/env python3
"""Refit eight-pair × two-arm ECFP4 / docking logistic models from the score master.

Primary: unscaled LogisticRegression, scaffold GroupKFold, same folds for
ECFP4-only, docking-only, and ECFP4+docking. Writes OOF predictions and
fold assignments. StandardScaler arm is sensitivity only.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import SEED, auroc  # noqa: E402

CANON = ROOT / "results" / "canonical"
MASTER = CANON / "current_score_master.csv"
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
ARMS = (
    ("D_vs_A", "dual", "A_only", "score_B"),
    ("D_vs_B", "dual", "B_only", "score_A"),
)


def fnum(v):
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
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def morgan(smi: str):
    mol = Chem.MolFromSmiles(smi or "")
    if mol is None:
        return None, None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "acyclic"
    except Exception:
        scaf = "unparsed"
    return mol, np.asarray(fp, dtype=float), scaf


def load_main():
    by = {p: [] for p in PRIMARY_PAIRS}
    for r in read_csv(MASTER):
        if r["analysis_set"] != "main" or r["complete_case"] not in ("1", "True", 1):
            continue
        if str(r.get("activity_eligible", "1")) not in ("1", "True"):
            continue
        mol, fp, scaf = morgan(r.get("smiles", ""))
        if fp is None:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["fp"] = fp
        rec["scaffold"] = scaf
        rec["tpsa"] = float(Descriptors.TPSA(mol))
        rec["clogp"] = float(Descriptors.MolLogP(mol))
        rec["heavy"] = int(Descriptors.HeavyAtomCount(mol))
        rec["mw"] = float(Descriptors.MolWt(mol))
        by[r["pair"]].append(rec)
    return by


def oof_predict(X, y, groups, scaled: bool):
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    n_scaf = len(set(groups.tolist()))
    n_splits = min(5, n_scaf, n_pos, n_neg)
    if n_splits < 2 or n_pos < 6 or n_neg < 6:
        return None, 0, np.full(len(y), np.nan), np.full(len(y), -1)
    cv = GroupKFold(n_splits=n_splits)
    oof = np.full(len(y), np.nan)
    fold_id = np.full(len(y), -1, dtype=int)
    for fold, (tr, te) in enumerate(cv.split(X, y, groups)):
        fold_id[te] = fold
        if scaled:
            model = Pipeline(
                [
                    ("scaler", StandardScaler()),
                    ("lr", LogisticRegression(max_iter=4000, C=1.0, random_state=SEED)),
                ]
            )
        else:
            model = LogisticRegression(max_iter=4000, C=1.0, random_state=SEED)
        model.fit(X[tr], y[tr])
        oof[te] = model.predict_proba(X[te])[:, 1]
    if np.isnan(oof).any():
        return None, n_splits, oof, fold_id
    return float(roc_auc_score(y, oof)), n_splits, oof, fold_id


def main() -> int:
    packs = load_main()
    fold_rows = []
    oof_rows = []
    inc_rows = []
    scaler_rows = []
    max_abs = 0.0
    max_loc = ""
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        for arm, pos_cls, neg_cls, dock_key in ARMS:
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            if len(kept) < 16:
                continue
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
            groups = np.array([r["scaffold"] for r in kept])
            fp = np.vstack([r["fp"] for r in kept])
            dock = np.array([[r[dock_key]] for r in kept], dtype=float)
            rank_dock = auroc(
                [r[dock_key] for r in kept if r["cls"] == pos_cls],
                [r[dock_key] for r in kept if r["cls"] == neg_cls],
            )
            models = {
                "ECFP4": fp,
                "docking": dock,
                "ECFP4+docking": np.hstack([fp, dock]),
            }
            primary = {}
            fold_ref = None
            for name, X in models.items():
                auc, n_splits, oof, fold_id = oof_predict(X, y, groups, scaled=False)
                primary[name] = auc
                if fold_ref is None:
                    fold_ref = fold_id
                    for i, r in enumerate(kept):
                        fold_rows.append(
                            {
                                "pair": pair,
                                "arm": arm,
                                "ligand_id": r["ligand_id"],
                                "scaffold": r["scaffold"],
                                "fold_id": int(fold_id[i]),
                                "class": r["cls"],
                                "n_splits": n_splits,
                            }
                        )
                for i, r in enumerate(kept):
                    oof_rows.append(
                        {
                            "pair": pair,
                            "arm": arm,
                            "model": name,
                            "scaled": 0,
                            "ligand_id": r["ligand_id"],
                            "class": r["cls"],
                            "y": int(y[i]),
                            "oof_prob": "" if oof is None or not math.isfinite(oof[i]) else f"{oof[i]:.6f}",
                            "dock_score": r[dock_key],
                        }
                    )
                sauc, _, soof, _ = oof_predict(X, y, groups, scaled=True)
                scaler_rows.append(
                    {
                        "pair": pair,
                        "arm": arm,
                        "model": name,
                        "unscaled_cv_auroc": "" if auc is None else f"{auc:.4f}",
                        "scaled_cv_auroc": "" if sauc is None else f"{sauc:.4f}",
                        "rank_auroc_docking": f"{rank_dock:.4f}",
                    }
                )
            inc = None
            if primary.get("ECFP4") is not None and primary.get("ECFP4+docking") is not None:
                inc = primary["ECFP4+docking"] - primary["ECFP4"]
                if abs(inc) > abs(max_abs):
                    max_abs = inc
                    max_loc = f"{pair} {arm}"
            inc_rows.append(
                {
                    "pair": pair,
                    "contrast": arm,
                    "n": len(kept),
                    "n_pos": int(y.sum()),
                    "n_neg": int((1 - y).sum()),
                    "n_scaffolds": len(set(groups.tolist())),
                    "cv_auroc_ECFP4": "" if primary.get("ECFP4") is None else f"{primary['ECFP4']:.4f}",
                    "cv_auroc_docking": "" if primary.get("docking") is None else f"{primary['docking']:.4f}",
                    "cv_auroc_ECFP4_docking": "" if primary.get("ECFP4+docking") is None else f"{primary['ECFP4+docking']:.4f}",
                    "delta_ECFP4_plus_docking_minus_ECFP4": "" if inc is None else f"{inc:.4f}",
                    "rank_auroc_docking": f"{rank_dock:.4f}",
                    "note": "unscaled logistic GroupKFold OOF; docking-only is logistic OOF, not rank AUROC",
                }
            )
    write_csv(CANON / "model_fold_assignments.csv", fold_rows)
    write_csv(CANON / "ecfp4_oof_predictions.csv", oof_rows)
    write_csv(CANON / "ecfp4_incremental_information.csv", inc_rows)
    write_csv(CANON / "ecfp4_scaler_sensitivity.csv", scaler_rows)
    print(f"max |ΔAUROC_incremental| = {abs(max_abs):.4f} at {max_loc}")
    for r in inc_rows:
        print(f"  {r['pair']:14} {r['contrast']:6} ECFP={r['cv_auroc_ECFP4']} +dock={r['cv_auroc_ECFP4_docking']} Δ={r['delta_ECFP4_plus_docking_minus_ECFP4']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
