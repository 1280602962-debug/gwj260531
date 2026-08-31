#!/usr/bin/env python3
"""Scaffold-cluster bootstrap for primary directional AUROCs on frozen K=4 pairs.

Resamples Bemis-Murcko scaffold groups (not individual ligands) using the same
frozen high-confidence Vina pocket-matched scores as document_cluster_bootstrap.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729
MIN_GROUPS = 2
MIN_CLASS = 2

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


def murcko_scaffold(smiles: str, ligand: str) -> str:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return f"__parse_failure_{ligand}"
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return f"__parse_failure_{ligand}"
    return scaffold or "__acyclic__"


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


def load_ligands() -> dict[str, list[dict]]:
    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    smiles_map = load_smiles()
    packs: dict[str, list[dict]] = defaultdict(list)
    for row in labels:
        pair = row["pair"]
        smiles = smiles_map.get((pair, row["ligand"]), "")
        rec = {
            "pair": pair,
            "ligand": row["ligand"],
            "molecule_chembl_id": row["molecule_chembl_id"],
            "cls": row["frozen_class"],
            "smiles": smiles,
            "vina_A": fnum(row["vina_A"]),
            "vina_B": fnum(row["vina_B"]),
            "scaffold_id": murcko_scaffold(smiles, row["ligand"]),
        }
        packs[pair].append(rec)
    return packs


def subset(recs, pos_cls, neg_cls):
    return [
        rec
        for rec in recs
        if rec["cls"] in (pos_cls, neg_cls) and rec["vina_A"] is not None and rec["vina_B"] is not None
    ]


def three_class_subset(recs):
    return [
        rec
        for rec in recs
        if rec["cls"] in ("dual", "A_only", "B_only") and rec["vina_A"] is not None and rec["vina_B"] is not None
    ]


def scaffold_groups(recs) -> tuple[list[str], dict[str, list[dict]]]:
    by_group = defaultdict(list)
    for rec in recs:
        by_group[rec["scaffold_id"]].append(rec)
    groups = sorted(by_group)
    return groups, by_group


def scaffolds_per_class(by_group, pos_cls, neg_cls) -> tuple[int, int]:
    pos_groups = sum(
        1 for items in by_group.values() if any(rec["cls"] == pos_cls for rec in items)
    )
    neg_groups = sum(
        1 for items in by_group.values() if any(rec["cls"] == neg_cls for rec in items)
    )
    return pos_groups, neg_groups


def stably_estimable(recs, pos_cls, neg_cls) -> tuple[bool, str]:
    usable = subset(recs, pos_cls, neg_cls)
    n_pos = sum(rec["cls"] == pos_cls for rec in usable)
    n_neg = sum(rec["cls"] == neg_cls for rec in usable)
    groups, by_group = scaffold_groups(usable)
    pos_groups, neg_groups = scaffolds_per_class(by_group, pos_cls, neg_cls)
    if n_pos < MIN_CLASS or n_neg < MIN_CLASS:
        return False, "fewer than two ligands in a class"
    if len(groups) < MIN_GROUPS:
        return False, f"only {len(groups)} scaffold groups"
    if pos_groups < MIN_GROUPS or neg_groups < MIN_GROUPS:
        return False, (
            f"too few scaffolds with both classes "
            f"(pos_scaffolds={pos_groups}, neg_scaffolds={neg_groups})"
        )
    return True, ""


def scaffold_cluster_bootstrap(recs, score_key, pos_cls, neg_cls, seed):
    usable = subset(recs, pos_cls, neg_cls)
    groups, by_group = scaffold_groups(usable)
    if len(groups) < MIN_GROUPS:
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


def summary_min_point(recs):
    da = auroc(
        [rec["vina_B"] for rec in recs if rec["cls"] == "dual"],
        [rec["vina_B"] for rec in recs if rec["cls"] == "A_only"],
    )
    db = auroc(
        [rec["vina_A"] for rec in recs if rec["cls"] == "dual"],
        [rec["vina_A"] for rec in recs if rec["cls"] == "B_only"],
    )
    if da != da or db != db:
        return float("nan")
    return min(da, db)


def summary_min_scaffold_bootstrap(recs, seed):
    usable = three_class_subset(recs)
    groups, by_group = scaffold_groups(usable)
    if len(groups) < MIN_GROUPS:
        return None, None, 0
    rng = np.random.default_rng(seed)
    values = []
    group_ids = np.asarray(groups)
    for _ in range(N_BOOT):
        chosen = rng.choice(group_ids, size=len(groups), replace=True)
        sample = [item for gid in chosen for item in by_group[gid]]
        value = summary_min_point(sample)
        if value == value:
            values.append(value)
    if len(values) < N_BOOT // 2:
        return None, None, len(values)
    lo, hi = np.percentile(values, [2.5, 97.5])
    return float(lo), float(hi), len(values)


def summary_min_estimable(recs) -> tuple[bool, str]:
    ok_a, note_a = stably_estimable(recs, "dual", "A_only")
    ok_b, note_b = stably_estimable(recs, "dual", "B_only")
    if not ok_a and not ok_b:
        return False, note_a or note_b
    if not ok_a:
        return False, f"D_vs_A arm: {note_a}"
    if not ok_b:
        return False, f"D_vs_B arm: {note_b}"
    return True, ""


def evaluate_contrast(pair, recs, contrast, pos_cls, neg_cls, score_key):
    usable = subset(recs, pos_cls, neg_cls)
    groups, _ = scaffold_groups(usable)
    point = auroc(
        [rec[score_key] for rec in usable if rec["cls"] == pos_cls],
        [rec[score_key] for rec in usable if rec["cls"] == neg_cls],
    )
    ok, reason = stably_estimable(recs, pos_cls, neg_cls)
    boot_lo, boot_hi, n_valid_boot = scaffold_cluster_bootstrap(
        recs, score_key, pos_cls, neg_cls, SEED + stable_offset(pair, contrast)
    )
    if not ok or boot_lo is None:
        return {
            "pair": pair,
            "contrast": contrast,
            "estimator": "scaffold_cluster_bootstrap",
            "n_groups": len(groups),
            "n_ligands": len(usable),
            "auroc_point": round(point, 4) if point == point else "",
            "ci_lo": "",
            "ci_hi": "",
            "n_valid_boot": n_valid_boot,
            "note": "not_stably_estimable",
        }
    return {
        "pair": pair,
        "contrast": contrast,
        "estimator": "scaffold_cluster_bootstrap",
        "n_groups": len(groups),
        "n_ligands": len(usable),
        "auroc_point": round(point, 4) if point == point else "",
        "ci_lo": round(boot_lo, 4),
        "ci_hi": round(boot_hi, 4),
        "n_valid_boot": n_valid_boot,
        "note": "resamples Bemis-Murcko scaffold groups, not individual ligands",
    }


def evaluate_summary_min(pair, recs):
    usable = three_class_subset(recs)
    groups, _ = scaffold_groups(usable)
    point = summary_min_point(usable)
    ok, _ = summary_min_estimable(recs)
    boot_lo, boot_hi, n_valid_boot = summary_min_scaffold_bootstrap(
        recs, SEED + stable_offset(pair, "summary_min")
    )
    if not ok or boot_lo is None:
        return {
            "pair": pair,
            "contrast": "summary_min",
            "estimator": "scaffold_cluster_bootstrap",
            "n_groups": len(groups),
            "n_ligands": len(usable),
            "auroc_point": round(point, 4) if point == point else "",
            "ci_lo": "",
            "ci_hi": "",
            "n_valid_boot": n_valid_boot,
            "note": "not_stably_estimable",
        }
    return {
        "pair": pair,
        "contrast": "summary_min",
        "estimator": "scaffold_cluster_bootstrap",
        "n_groups": len(groups),
        "n_ligands": len(usable),
        "auroc_point": round(point, 4) if point == point else "",
        "ci_lo": round(boot_lo, 4),
        "ci_hi": round(boot_hi, 4),
        "n_valid_boot": n_valid_boot,
        "note": "resamples Bemis-Murcko scaffold groups, not individual ligands",
    }


def ligand_bootstrap(recs, score_key, pos_cls, neg_cls, seed):
    pos = [rec for rec in recs if rec["cls"] == pos_cls and rec.get(score_key) is not None]
    neg = [rec for rec in recs if rec["cls"] == neg_cls and rec.get(score_key) is not None]
    if len(pos) < MIN_CLASS or len(neg) < MIN_CLASS:
        return None, None, 0
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(N_BOOT):
        p = [pos[i] for i in rng.integers(0, len(pos), len(pos))]
        n = [neg[i] for i in rng.integers(0, len(neg), len(neg))]
        value = auroc([rec[score_key] for rec in p], [rec[score_key] for rec in n])
        if value == value:
            values.append(value)
    if len(values) < N_BOOT // 2:
        return None, None, len(values)
    lo, hi = np.percentile(values, [2.5, 97.5])
    return float(lo), float(hi), len(values)


def load_reference_tables():
    document = {}
    for row in read_csv(TAB / "document_cluster_bootstrap_v1.csv"):
        document[(row["pair"], row["contrast"])] = row
    threshold = {}
    for row in read_csv(
        ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    ):
        if row["label_rule"] == "theta_6.0":
            threshold[row["pair"]] = row
    return document, threshold


def write_verdict(rows, packs, document_ref, threshold_ref):
    lines = [
        "# Cluster uncertainty verdict (v1)",
        "",
        "Compares ligand-level, document-connected, and Bemis-Murcko scaffold-cluster",
        "bootstrap 95% CIs on frozen high-confidence Vina pocket-matched directional AUROCs.",
        "B = 2000, seed 20260729. Scaffold grouping frozen before score inspection.",
        "",
        "## EGFR/HER2 weak arm (D_vs_B, pocket A)",
        "",
        "| estimator | AUROC | CI lo | CI hi | groups | note |",
        "|-----------|------:|------:|------:|-------:|------|",
    ]
    pair = "EGFR/HER2"
    contrast = "D_vs_B"
    recs = packs[pair]
    lig_lo, lig_hi, _ = ligand_bootstrap(
        recs, "vina_A", "dual", "B_only", SEED + stable_offset(pair, contrast, "ligand")
    )
    doc = document_ref[(pair, contrast)]
    scaf = next(r for r in rows if r["pair"] == pair and r["contrast"] == contrast)
    lines.append(
        f"| ligand bootstrap | {scaf['auroc_point']} | {lig_lo:.4f} | "
        f"{lig_hi:.4f} | — | class-preserving resample |"
    )
    lines.append(
        f"| document cluster | {doc['auroc_point']} | {doc['ci_lo']} | {doc['ci_hi']} | "
        f"{doc['n_groups']} docs | {doc['note']} |"
    )
    lines.append(
        f"| scaffold cluster | {scaf['auroc_point']} | {scaf['ci_lo'] or '—'} | "
        f"{scaf['ci_hi'] or '—'} | {scaf['n_groups']} scaffolds | {scaf['note']} |"
    )
    lines += [
        "",
        "The weak EGFR/HER2 arm (D_vs_B = 0.4297) stays near chance under all three",
        "resampling schemes; scaffold and document CIs are wider than ligand bootstrap",
        "because correlated chemotypes/documents are kept together.",
        "",
        "## PIK3CA/mTOR issues",
        "",
        "| contrast | ligand CI | document CI | scaffold CI | document status | scaffold note |",
        "|----------|-----------|-------------|-------------|-----------------|---------------|",
    ]
    pair = "PIK3CA/mTOR"
    thresh = threshold_ref[pair]
    for contrast in ("D_vs_A", "D_vs_B", "summary_min"):
        doc = document_ref.get((pair, contrast))
        scaf = next(r for r in rows if r["pair"] == pair and r["contrast"] == contrast)
        if contrast == "D_vs_A":
            lig_lo, lig_hi, _ = ligand_bootstrap(
                packs[pair], "vina_B", "dual", "A_only", SEED + stable_offset(pair, contrast, "ligand")
            )
        elif contrast == "D_vs_B":
            lig_lo, lig_hi, _ = ligand_bootstrap(
                packs[pair], "vina_A", "dual", "B_only", SEED + stable_offset(pair, contrast, "ligand")
            )
        else:
            lig_lo, lig_hi = float(thresh["ci_lo"]), float(thresh["ci_hi"])
        lig_ci = f"[{lig_lo:.4f}, {lig_hi:.4f}]"
        doc_ci = (
            f"[{doc['ci_lo']}, {doc['ci_hi']}]"
            if doc and doc["ci_lo"] != ""
            else "—"
        )
        scaf_ci = (
            f"[{scaf['ci_lo']}, {scaf['ci_hi']}]"
            if scaf["ci_lo"] != ""
            else "—"
        )
        doc_status = "cannot_stably_estimate" if contrast == "D_vs_B" else "ok"
        if contrast == "D_vs_B":
            doc_status = "cannot_stably_estimate (1 valid doc-blocked fold)"
        lines.append(
            f"| {contrast} | {lig_ci} | {doc_ci} | {scaf_ci} | {doc_status} | {scaf['note']} |"
        )
    lines += [
        "",
        "## Cross-estimator summary (all K=4 pairs)",
        "",
        "| pair | contrast | ligand CI width | document CI width | scaffold CI width |",
        "|------|----------|----------------:|------------------:|------------------:|",
    ]
    for pair in TARGETS:
        recs = packs[pair]
        for contrast, score_key, pos_cls, neg_cls in (
            ("D_vs_A", "vina_B", "dual", "A_only"),
            ("D_vs_B", "vina_A", "dual", "B_only"),
        ):
            doc = document_ref[(pair, contrast)]
            scaf = next(r for r in rows if r["pair"] == pair and r["contrast"] == contrast)
            lig_lo, lig_hi, _ = ligand_bootstrap(
                recs, score_key, pos_cls, neg_cls, SEED + stable_offset(pair, contrast, "ligand")
            )
            lig_w = float(lig_hi) - float(lig_lo)
            doc_w = float(doc["ci_hi"]) - float(doc["ci_lo"]) if doc["ci_lo"] else float("nan")
            scaf_w = (
                float(scaf["ci_hi"]) - float(scaf["ci_lo"])
                if scaf["ci_lo"] != ""
                else float("nan")
            )
            doc_w_str = f"{doc_w:.3f}" if doc_w == doc_w else "—"
            scaf_w_str = f"{scaf_w:.3f}" if scaf_w == scaf_w else "—"
            lines.append(
                f"| {pair} | {contrast} | {lig_w:.3f} | {doc_w_str} | {scaf_w_str} |"
            )
    lines += [
        "",
        "Scaffold-cluster bootstrap generally widens CIs relative to ligand bootstrap and",
        "is comparable to or slightly narrower than document-cluster bootstrap depending on",
        "whether document co-reporting or scaffold reuse dominates correlation structure.",
        "PIK3CA/mTOR remains the most uncertainty-sensitive pair because of small n and",
        "heavy document/scaffold concentration; document-blocked CV already flags D_vs_B as",
        "not stably estimable under leave-document-out folds.",
        "",
    ]
    (ANALYSIS / "CLUSTER_UNCERTAINTY_VERDICT_V1.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    packs = load_ligands()
    rows = []
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, score_key in CONTRASTS:
            rows.append(evaluate_contrast(pair, recs, contrast, pos_cls, neg_cls, score_key))
        rows.append(evaluate_summary_min(pair, recs))

    write_csv(TAB / "scaffold_cluster_bootstrap_v1.csv", rows)
    document_ref, threshold_ref = load_reference_tables()
    write_verdict(rows, packs, document_ref, threshold_ref)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
