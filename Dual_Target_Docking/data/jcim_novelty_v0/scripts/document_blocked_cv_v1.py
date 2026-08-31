#!/usr/bin/env python3
"""Document-blocked grouped CV and document-cluster bootstrap.

Same ChEMBL document cannot contribute ligands to both train and test.
Ligands that share any retained high-confidence document are connected into
one group. Grouping is frozen before scores are inspected; a worse result is
kept. If a pair has too few groups or folds lack both classes, the pair is
reported as not stably estimable rather than regrouped.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729
MAX_SPLITS = 5
MIN_VALID_FOLDS = 2

TARGETS = {
    "EGFR/HER2": ("CHEMBL203", "CHEMBL1824"),
    "AChE/BChE": ("CHEMBL220", "CHEMBL1914"),
    "PIK3CA/PIK3CB": ("CHEMBL4005", "CHEMBL3145"),
    "PIK3CA/mTOR": ("CHEMBL4005", "CHEMBL2842"),
}

PANEL_SMILES = {
    "EGFR/HER2": "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
    "AChE/BChE": "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
    "PIK3CA/PIK3CB": "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
    "PIK3CA/mTOR": "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
}

CONTRASTS = (
    ("D_vs_A", "dual", "A_only", "vina_B"),
    ("D_vs_B", "dual", "B_only", "vina_A"),
)


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        if value in ("", None):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    p = np.asarray(pos, float).ravel()
    n = np.asarray(neg, float).ravel()
    if p.size == 0 or n.size == 0:
        return float("nan")
    delta = p[:, None] - n[None, :]
    return float(((delta > 0).sum() + 0.5 * (delta == 0).sum()) / (p.size * n.size))


def union_find(items):
    parent = {item: item for item in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    return find, union


def load_smiles() -> dict[tuple[str, str], str]:
    out = {}
    for pair, relative in PANEL_SMILES.items():
        path = ROOT / relative
        if not path.exists():
            continue
        for row in read_csv(path):
            ligand = row.get("ligand") or row.get("panel_id")
            smiles = row.get("smiles") or ""
            if ligand and smiles:
                out[(pair, ligand)] = smiles
    return out


def load_ligands():
    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    activities = [
        row
        for row in read_csv(TAB / "high_confidence_activity_audit_v1.csv")
        if row.get("keep") in ("1", "True", "true")
    ]
    smiles_map = load_smiles()
    by_mol_target: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in activities:
        by_mol_target[(row["molecule_chembl_id"], row["target_chembl_id"])].append(row)

    packs: dict[str, list[dict]] = defaultdict(list)
    for row in labels:
        pair = row["pair"]
        target_a, target_b = TARGETS[pair]
        mol = row["molecule_chembl_id"]
        docs = {
            rec.get("document_chembl_id")
            for target in (target_a, target_b)
            for rec in by_mol_target[(mol, target)]
            if rec.get("document_chembl_id")
        }
        smiles = smiles_map.get((pair, row["ligand"]), "")
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        rec = {
            "pair": pair,
            "ligand": row["ligand"],
            "molecule_chembl_id": mol,
            "cls": row["frozen_class"],
            "documents": sorted(docs),
            "smiles": smiles,
            "vina_A": fnum(row["vina_A"]),
            "vina_B": fnum(row["vina_B"]),
            "heavy": float(molecule.GetNumHeavyAtoms()) if molecule is not None else None,
            "tpsa": float(Descriptors.TPSA(molecule)) if molecule is not None else None,
            "clogp": float(Descriptors.MolLogP(molecule)) if molecule is not None else None,
            "fp": None,
        }
        if molecule is not None:
            rec["fp"] = np.asarray(AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=2048))
        packs[pair].append(rec)

    for pair, recs in packs.items():
        nodes = []
        for rec in recs:
            nodes.append(("lig", rec["ligand"]))
            nodes.extend(("doc", doc) for doc in rec["documents"])
        find, union = union_find(nodes)
        for rec in recs:
            if rec["documents"]:
                first = rec["documents"][0]
                union(("lig", rec["ligand"]), ("doc", first))
                for doc in rec["documents"][1:]:
                    union(("doc", first), ("doc", doc))
        roots = {}
        next_id = 1
        for rec in recs:
            root = find(("lig", rec["ligand"]))
            if root not in roots:
                roots[root] = f"G{next_id:03d}"
                next_id += 1
            rec["group_id"] = roots[root]
    return packs


def subset(recs, pos_cls, neg_cls):
    return [
        rec
        for rec in recs
        if rec["cls"] in (pos_cls, neg_cls) and rec["vina_A"] is not None and rec["vina_B"] is not None
    ]


def document_cluster_bootstrap(recs, score_key, pos_cls, neg_cls, seed):
    usable = subset(recs, pos_cls, neg_cls)
    groups = sorted({rec["group_id"] for rec in usable})
    by_group = defaultdict(list)
    for rec in usable:
        by_group[rec["group_id"]].append(rec)
    if len(groups) < 2:
        return None, None, 0
    rng = np.random.default_rng(seed)
    values = []
    group_ids = np.asarray(groups)
    for _ in range(N_BOOT):
        chosen = rng.choice(group_ids, size=len(groups), replace=True)
        sample = [item for gid in chosen for item in by_group[gid]]
        pos = [rec[score_key] for rec in sample if rec["cls"] == pos_cls]
        neg = [rec[score_key] for rec in sample if rec["cls"] == neg_cls]
        value = auroc(pos, neg)
        if value == value:
            values.append(value)
    if len(values) < N_BOOT // 2:
        return None, None, len(values)
    lo, hi = np.percentile(values, [2.5, 97.5])
    return float(lo), float(hi), len(values)


def logistic_oof(X, y, groups, train_idx, test_idx):
    if len(set(y[train_idx])) < 2:
        return None
    model = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED)
    model.fit(X[train_idx], y[train_idx])
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X[test_idx])[:, 1]
    return model.decision_function(X[test_idx])


def evaluate_contrast(pair, recs, contrast, pos_cls, neg_cls, score_key):
    usable = subset(recs, pos_cls, neg_cls)
    n_pos = sum(rec["cls"] == pos_cls for rec in usable)
    n_neg = sum(rec["cls"] == neg_cls for rec in usable)
    groups = [rec["group_id"] for rec in usable]
    n_groups = len(set(groups))
    n_docs = len({doc for rec in usable for doc in rec["documents"]})
    class_docs = {
        cls: len({doc for rec in usable if rec["cls"] == cls for doc in rec["documents"]})
        for cls in (pos_cls, neg_cls)
    }
    group_sizes = Counter(groups)
    top_group_n = max(group_sizes.values()) if group_sizes else 0
    rank_point = auroc(
        [rec[score_key] for rec in usable if rec["cls"] == pos_cls],
        [rec[score_key] for rec in usable if rec["cls"] == neg_cls],
    )
    boot_lo, boot_hi, n_valid_boot = document_cluster_bootstrap(
        recs, score_key, pos_cls, neg_cls, SEED + stable_offset(pair, contrast)
    )

    status = "ok"
    notes = []
    if n_pos < 2 or n_neg < 2:
        status = "cannot_stably_estimate"
        notes.append("fewer than two ligands in a class")
    if n_groups < MIN_VALID_FOLDS:
        status = "cannot_stably_estimate"
        notes.append(f"only {n_groups} document-connected groups")

    fold_rows = []
    model_rows = []
    summary = {
        "pair": pair,
        "contrast": contrast,
        "n_ligands": len(usable),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_documents": n_docs,
        "n_pos_documents": class_docs[pos_cls],
        "n_neg_documents": class_docs[neg_cls],
        "n_groups": n_groups,
        "largest_group_n": top_group_n,
        "largest_group_fraction": round(top_group_n / len(usable), 4) if usable else "",
        "n_splits_requested": "",
        "n_valid_folds": 0,
        "rank_auroc_full": round(rank_point, 4) if rank_point == rank_point else "",
        "rank_auroc_doc_blocked_mean_fold": "",
        "rank_auroc_doc_blocked_oof": "",
        "ecfp4_auroc_oof": "",
        "physchem_auroc_oof": "",
        "docking_logistic_auroc_oof": "",
        "doc_cluster_boot_lo": round(boot_lo, 4) if boot_lo is not None else "",
        "doc_cluster_boot_hi": round(boot_hi, 4) if boot_hi is not None else "",
        "n_valid_boot": n_valid_boot,
        "status": status,
        "note": "; ".join(notes) or "document-connected groups frozen; worse AUROC retained",
    }
    if status != "ok":
        return summary, fold_rows, model_rows

    n_splits = min(MAX_SPLITS, n_pos, n_neg, n_groups)
    if n_splits < MIN_VALID_FOLDS:
        summary["status"] = "cannot_stably_estimate"
        summary["n_splits_requested"] = n_splits
        summary["note"] = "requested splits < 2 after group/class constraints"
        return summary, fold_rows, model_rows

    y = np.array([1 if rec["cls"] == pos_cls else 0 for rec in usable])
    group_arr = np.array(groups)
    scores = np.array([rec[score_key] for rec in usable], float)
    fp_ok = all(rec["fp"] is not None for rec in usable)
    phys_ok = all(rec[k] is not None for rec in usable for k in ("heavy", "tpsa", "clogp"))
    X_fp = np.vstack([rec["fp"] for rec in usable]) if fp_ok else None
    X_phys = (
        np.array([[rec["heavy"], rec["tpsa"], rec["clogp"]] for rec in usable], float)
        if phys_ok
        else None
    )
    X_dock = scores.reshape(-1, 1)

    cv = GroupKFold(n_splits=n_splits)
    valid = []
    for fold_id, (train_idx, test_idx) in enumerate(cv.split(scores, y, group_arr), 1):
        if len(set(y[test_idx])) < 2 or len(set(y[train_idx])) < 2:
            fold_rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "fold": fold_id,
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    "n_test_pos": int(y[test_idx].sum()),
                    "n_test_neg": int(len(test_idx) - y[test_idx].sum()),
                    "n_test_groups": int(len(set(group_arr[test_idx]))),
                    "n_test_documents": len(
                        {doc for i in test_idx for doc in usable[i]["documents"]}
                    ),
                    "train_test_group_overlap": int(
                        bool(set(group_arr[train_idx]) & set(group_arr[test_idx]))
                    ),
                    "rank_auroc": "",
                    "kept": 0,
                    "note": "fold dropped: a class is missing in train or test",
                }
            )
            continue
        valid.append((fold_id, train_idx, test_idx))

    summary["n_splits_requested"] = n_splits
    summary["n_valid_folds"] = len(valid)
    if len(valid) < MIN_VALID_FOLDS:
        summary["status"] = "cannot_stably_estimate"
        summary["note"] = (
            f"only {len(valid)} folds contained both classes; grouping rule not changed"
        )
        return summary, fold_rows, model_rows

    oof_rank = np.full(len(usable), np.nan)
    oof_fp = np.full(len(usable), np.nan)
    oof_phys = np.full(len(usable), np.nan)
    oof_dock = np.full(len(usable), np.nan)
    fold_aurocs = []
    for fold_id, train_idx, test_idx in valid:
        fold_auroc = auroc(scores[test_idx][y[test_idx] == 1], scores[test_idx][y[test_idx] == 0])
        fold_aurocs.append(fold_auroc)
        oof_rank[test_idx] = scores[test_idx]
        fold_rows.append(
            {
                "pair": pair,
                "contrast": contrast,
                "fold": fold_id,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "n_test_pos": int(y[test_idx].sum()),
                "n_test_neg": int(len(test_idx) - y[test_idx].sum()),
                "n_test_groups": int(len(set(group_arr[test_idx]))),
                "n_test_documents": len(
                    {doc for i in test_idx for doc in usable[i]["documents"]}
                ),
                "train_test_group_overlap": int(
                    bool(set(group_arr[train_idx]) & set(group_arr[test_idx]))
                ),
                "rank_auroc": round(fold_auroc, 4) if fold_auroc == fold_auroc else "",
                "kept": 1,
                "note": "both classes present; document groups do not overlap",
            }
        )
        if X_fp is not None:
            pred = logistic_oof(X_fp, y, group_arr, train_idx, test_idx)
            if pred is not None:
                oof_fp[test_idx] = pred
        if X_phys is not None:
            pred = logistic_oof(X_phys, y, group_arr, train_idx, test_idx)
            if pred is not None:
                oof_phys[test_idx] = pred
        pred = logistic_oof(X_dock, y, group_arr, train_idx, test_idx)
        if pred is not None:
            oof_dock[test_idx] = pred

    def oof_auc(values):
        mask = np.isfinite(values)
        if mask.sum() < 4 or len(set(y[mask])) < 2:
            return ""
        return round(float(roc_auc_score(y[mask], values[mask])), 4)

    summary["rank_auroc_doc_blocked_mean_fold"] = (
        round(float(np.nanmean(fold_aurocs)), 4) if fold_aurocs else ""
    )
    summary["rank_auroc_doc_blocked_oof"] = oof_auc(oof_rank)
    summary["ecfp4_auroc_oof"] = oof_auc(oof_fp)
    summary["physchem_auroc_oof"] = oof_auc(oof_phys)
    summary["docking_logistic_auroc_oof"] = oof_auc(oof_dock)
    for method, values in (
        ("rank_vina", oof_rank),
        ("ECFP4_logistic", oof_fp),
        ("physchem_logistic", oof_phys),
        ("docking_logistic", oof_dock),
    ):
        auc = oof_auc(values)
        model_rows.append(
            {
                "pair": pair,
                "contrast": contrast,
                "method": method,
                "cv_scheme": "document_connected_GroupKFold",
                "n_valid_folds": len(valid),
                "n_groups": n_groups,
                "n_documents": n_docs,
                "n_pos": n_pos,
                "n_neg": n_neg,
                "auroc_oof": auc,
                "note": "same folds for all methods; no regrouping after seeing AUROC",
            }
        )
    return summary, fold_rows, model_rows


def main() -> None:
    packs = load_ligands()
    group_rows = []
    summary_rows = []
    fold_rows = []
    model_rows = []
    boot_rows = []
    for pair, recs in packs.items():
        for rec in recs:
            group_rows.append(
                {
                    "pair": pair,
                    "ligand": rec["ligand"],
                    "molecule_chembl_id": rec["molecule_chembl_id"],
                    "frozen_class": rec["cls"],
                    "group_id": rec["group_id"],
                    "n_documents": len(rec["documents"]),
                    "documents": ";".join(rec["documents"]),
                }
            )
        for contrast, pos_cls, neg_cls, score_key in CONTRASTS:
            summary, folds, models = evaluate_contrast(
                pair, recs, contrast, pos_cls, neg_cls, score_key
            )
            summary_rows.append(summary)
            fold_rows.extend(folds)
            model_rows.extend(models)
            boot_rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "estimator": "document_cluster_bootstrap",
                    "n_groups": summary["n_groups"],
                    "n_documents": summary["n_documents"],
                    "auroc_point": summary["rank_auroc_full"],
                    "ci_lo": summary["doc_cluster_boot_lo"],
                    "ci_hi": summary["doc_cluster_boot_hi"],
                    "n_valid_boot": summary["n_valid_boot"],
                    "note": "resamples document-connected ligand groups, not individual ligands",
                }
            )

    write_csv(TAB / "document_blocked_ligand_groups_v1.csv", group_rows)
    write_csv(TAB / "document_blocked_cv_summary_v1.csv", summary_rows)
    write_csv(TAB / "document_blocked_cv_folds_v1.csv", fold_rows)
    write_csv(TAB / "document_blocked_cv_methods_v1.csv", model_rows)
    write_csv(TAB / "document_cluster_bootstrap_v1.csv", boot_rows)

    lines = [
        "# Document-blocked CV verdict",
        "",
        "Grouping rule: ligands sharing any retained high-confidence ChEMBL `document_id` are one group.",
        "The same folds are used for ECFP4, physicochemical, and docking logistic models.",
        "Grouping was not changed after seeing AUROC.",
        "",
        "| pair | contrast | n_pos/n_neg | groups | docs | valid folds | rank full | rank blocked mean-fold | ECFP4 | physchem | dock logistic | status |",
        "|------|----------|-------------|-------:|-----:|------------:|----------:|-----------------------:|------:|---------:|--------------:|--------|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['pair']} | {row['contrast']} | {row['n_pos']}/{row['n_neg']} | "
            f"{row['n_groups']} | {row['n_documents']} | {row['n_valid_folds']} | "
            f"{row['rank_auroc_full']} | {row['rank_auroc_doc_blocked_mean_fold']} | "
            f"{row['ecfp4_auroc_oof']} | {row['physchem_auroc_oof']} | "
            f"{row['docking_logistic_auroc_oof']} | {row['status']} |"
        )
    n_ok = sum(row["status"] == "ok" for row in summary_rows)
    lines += [
        "",
        f"Estimable directional arms: {n_ok}/{len(summary_rows)}.",
        "If an arm is not stably estimable, that is a result, not a reason to regroup.",
        "",
    ]
    (ANALYSIS / "DOCUMENT_BLOCKED_CV_VERDICT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary_rows, indent=2))


if __name__ == "__main__":
    main()
