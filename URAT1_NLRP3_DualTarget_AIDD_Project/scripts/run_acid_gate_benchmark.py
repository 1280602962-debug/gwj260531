#!/usr/bin/env python3
"""Retrospective URAT1 acid-gate benchmark: carboxylate actives vs true decoys.

Uses frozen phase-1 gnina 9-mode SDFs (no re-docking).
Compares A1 vs A2 pose selection on the same poses.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import fisher_exact

PROJECT_ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from c1_acid_pose_selection import evaluate_urat1_acid_sdf, load_ref_centroid  # noqa: E402

COOH = Chem.MolFromSmarts("[CX3](=O)[OH,O-]")
POSE_ROOT = PROJECT_ROOT / "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/poses/gnina_sdf"
MAP_CSV = PROJECT_ROOT / "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/meta/mol_index_map.csv"
OUT = PROJECT_ROOT / "data/campaigns/c1/05_metrics/acid_gate_retrospective_benchmark"


def has_acid(smi: str) -> bool:
    m = Chem.MolFromSmiles(smi)
    return bool(m and COOH and m.HasSubstructMatch(COOH))


def bootstrap_or_ci(y_true, y_pred, n_boot: int = 2000, seed: int = 42) -> tuple[float, float]:
    rng = random.Random(seed)
    n = len(y_true)
    ors = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        yt = [y_true[i] for i in idx]
        yp = [y_pred[i] for i in idx]
        a = sum(1 for t, p in zip(yt, yp) if t == 1 and p == 1)
        b = sum(1 for t, p in zip(yt, yp) if t == 1 and p == 0)
        c = sum(1 for t, p in zip(yt, yp) if t == 0 and p == 1)
        d = sum(1 for t, p in zip(yt, yp) if t == 0 and p == 0)
        if min(a, b, c, d) == 0:
            ors.append(float("inf") if a * d > 0 else 0.0)
        else:
            ors.append((a * d) / (b * c))
    lo, hi = np.percentile(ors, [2.5, 97.5])
    return float(lo), float(hi)


def metrics_from_pass(labels: list[int], passes: list[bool]) -> dict:
    y = np.array(labels, dtype=int)
    p = np.array(passes, dtype=bool)
    tp = int(((y == 1) & p).sum())
    fn = int(((y == 1) & ~p).sum())
    fp = int(((y == 0) & p).sum())
    tn = int(((y == 0) & ~p).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    table = [[tp, fn], [fp, tn]]
    oddsr, pval = fisher_exact(table)
    ci_lo, ci_hi = bootstrap_or_ci(labels, passes)
    return {
        "n": int(len(labels)),
        "n_active": int((y == 1).sum()),
        "n_decoy": int((y == 0).sum()),
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "sensitivity": sens,
        "specificity": spec,
        "ppv": ppv,
        "odds_ratio": float(oddsr),
        "fisher_exact_p": float(pval),
        "or_bootstrap_ci95_lo": ci_lo,
        "or_bootstrap_ci95_hi": ci_hi,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="debug subset")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    act = pd.read_csv(PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/actives.csv")
    dec = pd.read_csv(PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/true_decoys.csv")
    act = act[act.canonical_smiles.map(has_acid)].copy()
    dec = dec[dec.canonical_smiles.map(has_acid)].copy()
    act["label"] = 1
    dec["label"] = 0
    pool = pd.concat([act, dec], ignore_index=True)
    pool = pool.drop_duplicates(subset=["canonical_smiles"])

    mmap = pd.read_csv(MAP_CSV)
    smi2mol = dict(zip(mmap.canonical_smiles, mmap.mol_id))
    pool["mol_id"] = pool.canonical_smiles.map(smi2mol)
    pool = pool[pool.mol_id.notna()].copy()
    if args.limit:
        pool = pool.head(args.limit)

    arg_json = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs/arg477_coords.json"
    ref_com = load_ref_centroid(
        PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs/lesinurad_crystal_ref.sdf"
    )

    rows = []
    for _, r in pool.iterrows():
        sdf = POSE_ROOT / f"{r['mol_id']}_out.sdf"
        if not sdf.exists():
            continue
        for rule in ("a1", "a2"):
            ev = evaluate_urat1_acid_sdf(sdf, arg_json, ref_com, r["mol_id"], 42, rule=rule)
            ev["label"] = int(r["label"])
            ev["canonical_smiles"] = r["canonical_smiles"]
            ev["set_role"] = "active" if r["label"] == 1 else "true_decoy"
            rows.append(ev)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "acid_gate_benchmark_per_mol.csv", index=False)

    summary = {"subset": "carboxylate-compatible actives vs true decoys", "pose_source": str(POSE_ROOT)}
    for rule in ("a1", "a2"):
        sub = df[df.pose_selection_rule == rule]
        summary[rule] = metrics_from_pass(sub.label.tolist(), sub.keep_urat1_acid.tolist())

    (OUT / "acid_gate_benchmark_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
