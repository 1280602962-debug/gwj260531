#!/usr/bin/env python3
"""AChE/BChE ECFP4 OOF + incremental AUROC on the corrected complete-case panel.

Protocol matches independent_metrics_audit_v1 / incremental_information:
  Morgan r=2 2048-bit, LogisticRegression C=1 max_iter=4000
  GroupKFold n_splits=min(5, n_scaffolds, n_pos, n_neg)
  D_vs_A uses vina_B; D_vs_B uses vina_A
  AChE SMILES used as deposited (not largest-fragment Track B rule)
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
PANEL = ROOT / "remediation_outputs/phase2_ache_bche/corrected_main_complete_case.csv"
OUT = ROOT / "remediation_outputs/phase2_ache_bche/ache_ecfp4_incremental_corrected.csv"


def fp_array(mol) -> np.ndarray:
    bv = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    arr = np.zeros((2048,), dtype=np.int8)
    AllChem.DataStructs.ConvertToNumpyArray(bv, arr)
    return arr


def grouped_oof(X, y, groups):
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    n_scaf = len(set(groups.tolist()))
    n_splits = min(5, n_scaf, n_pos, n_neg)
    if n_splits < 2 or n_pos < 6 or n_neg < 6:
        return float("nan"), n_splits, n_scaf, 0
    cv = GroupKFold(n_splits=n_splits)
    leaked = 0
    pred = np.zeros(len(y), dtype=float)
    for train, test in cv.split(X, y, groups):
        if set(groups[train]) & set(groups[test]):
            leaked += 1
        model = LogisticRegression(max_iter=4000, C=1.0)
        model.fit(X[train], y[train])
        pred[test] = model.predict_proba(X[test])[:, 1]
    return float(roc_auc_score(y, pred)), n_splits, n_scaf, leaked


def descriptor_auroc(recs, key, pos_cls, neg_cls):
    pos = [r[key] for r in recs if r["primary_class_theta6"] == pos_cls]
    neg = [r[key] for r in recs if r["primary_class_theta6"] == neg_cls]
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def main() -> int:
    recs = []
    with PANEL.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            mol = Chem.MolFromSmiles(r["canonical_smiles"])
            if mol is None:
                continue
            try:
                scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            except Exception:
                scaf = ""
            recs.append(
                {
                    **r,
                    "_fp": fp_array(mol),
                    "_scaf": scaf or f"__fail_{r['ligand_id']}",
                    "tpsa": Descriptors.TPSA(mol),
                    "mw": Descriptors.MolWt(mol),
                    "clogp": Descriptors.MolLogP(mol),
                    "score_A": float(r["score_A"]),
                    "score_B": float(r["score_B"]),
                }
            )
    rows = []
    pre = {
        ("D_vs_A", "ECFP4"): 0.8948148148148148,
        ("D_vs_A", "ECFP4+dock"): 0.8933333333333333,
        ("D_vs_A", "increment"): -0.001481481481481528,
        ("D_vs_B", "ECFP4"): 0.8214285714285714,
        ("D_vs_B", "ECFP4+dock"): 0.8082010582010581,
        ("D_vs_B", "increment"): -0.013227513227513255,
    }
    for contrast, pos_cls, neg_cls, dock_key in (
        ("D_vs_A", "dual", "A_only", "score_B"),
        ("D_vs_B", "dual", "B_only", "score_A"),
    ):
        kept = [r for r in recs if r["primary_class_theta6"] in (pos_cls, neg_cls)]
        y = np.array([1.0 if r["primary_class_theta6"] == pos_cls else 0.0 for r in kept])
        groups = np.array([r["_scaf"] for r in kept])
        Xfp = np.vstack([r["_fp"] for r in kept])
        dock = np.array([r[dock_key] for r in kept], dtype=float).reshape(-1, 1)
        Xboth = np.hstack([Xfp, dock])
        a_fp, ns, nsc, leak = grouped_oof(Xfp, y, groups)
        a_both, _, _, leak2 = grouped_oof(Xboth, y, groups)
        inc = a_both - a_fp
        tpsa = descriptor_auroc(recs, "tpsa", pos_cls, neg_cls)
        rows.append(
            {
                "pair": "AChE/BChE",
                "contrast": contrast,
                "n_pos": int(y.sum()),
                "n_neg": int((1 - y).sum()),
                "n_splits": ns,
                "n_scaffolds": nsc,
                "scaffold_leakage": leak + leak2,
                "AUROC_ECFP4_OOF": a_fp,
                "AUROC_ECFP4_plus_docking_OOF": a_both,
                "delta_AUROC_incremental": inc,
                "AUROC_TPSA": tpsa,
                "pre_ECFP4": pre[(contrast, "ECFP4")],
                "pre_ECFP4_dock": pre[(contrast, "ECFP4+dock")],
                "pre_increment": pre[(contrast, "increment")],
                "interpretation_ligand_only": "substantial" if a_fp >= 0.70 else "weak",
                "interpretation_increment": "limited" if abs(inc) < 0.05 else "material",
                "ligand_only_changed": "no" if (a_fp >= 0.70) == (pre[(contrast, "ECFP4")] >= 0.70) else "yes",
                "increment_changed": "no" if (abs(inc) < 0.05) == (abs(pre[(contrast, "increment")]) < 0.05) else "yes",
            }
        )
        print(contrast, f"ECFP4={a_fp:.4f} +dock={a_both:.4f} d={inc:.4f} TPSA={tpsa:.4f} leak={leak}")
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
