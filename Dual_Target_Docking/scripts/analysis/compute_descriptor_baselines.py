#!/usr/bin/env python3
"""Single-descriptor baselines from current score master (RDKit descriptors + scheme-B)."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import (  # noqa: E402
    N_BOOT,
    SEED,
    auroc,
    percentile_ci,
    summary_min_stratified,
)

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
DESCS = ("tpsa", "clogp", "heavy")


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load():
    by = {p: [] for p in PRIMARY_PAIRS}
    for r in read_csv(MASTER):
        if r["analysis_set"] != "main" or r["complete_case"] not in ("1", 1, "True"):
            continue
        if str(r.get("activity_eligible", "1")) not in ("1", "True"):
            continue
        mol = Chem.MolFromSmiles(r.get("smiles") or "")
        if mol is None:
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["tpsa"] = float(Descriptors.TPSA(mol))
        rec["clogp"] = float(Descriptors.MolLogP(mol))
        rec["heavy"] = float(Descriptors.HeavyAtomCount(mol))
        by[r["pair"]].append(rec)
    return by


def paired_delta(vina_smin, desc_smin, n_boot=N_BOOT, seed=SEED):
    """Bootstrap CI for vina_summary_min - descriptor_summary_min with shared dual draws.

    vina_smin/desc_smin are callables on resampled index lists — we instead
    resample ligands and recompute both summary_mins.
    """
    return vina_smin, desc_smin


def main() -> int:
    packs = load()
    rows = []
    rng_seed = SEED
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        cls = np.array([r["cls"] for r in recs])
        vina = summary_min_stratified(sa, sb, cls, n_boot=1, seed=rng_seed)
        desc_smin = {}
        desc_arms = {}
        for name in DESCS:
            val = np.array([r[name] for r in recs], dtype=float)
            st = summary_min_stratified(val, val, cls, n_boot=1, seed=rng_seed)
            # directional: D vs A uses descriptor as "pocket B" analogue — both arms same feature
            d = cls == "dual"
            a = cls == "A_only"
            b = cls == "B_only"
            da = auroc(val[d], val[a])
            db = auroc(val[d], val[b])
            desc_smin[name] = min(da, db)
            desc_arms[name] = (da, db)
        best = max(desc_smin, key=desc_smin.get)
        # paired bootstrap of vina_smin - best_desc_smin
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
        rec = {
            "pair": pair,
            "n_dual": vina["n_dual"],
            "n_A_only": vina["n_A_only"],
            "n_B_only": vina["n_B_only"],
            "vina_summary_min": f"{vina['summary_min']:.4f}",
            "vina_D_vs_A": f"{vina['auroc_D_vs_A_pocketB']:.4f}",
            "vina_D_vs_B": f"{vina['auroc_D_vs_B_pocketA']:.4f}",
            "best_descriptor": best,
            "best_descriptor_summary_min": f"{desc_smin[best]:.4f}",
            "vina_minus_best_descriptor": f"{point:.4f}",
            "delta_ci_lo": f"{lo:.4f}",
            "delta_ci_hi": f"{hi:.4f}",
            "bootstrap": "class_stratified_shared_dual",
            "B": N_BOOT,
            "seed": SEED,
        }
        for name in DESCS:
            rec[f"{name}_D_vs_A"] = f"{desc_arms[name][0]:.4f}"
            rec[f"{name}_D_vs_B"] = f"{desc_arms[name][1]:.4f}"
            rec[f"{name}_summary_min"] = f"{desc_smin[name]:.4f}"
        rows.append(rec)
        print(f"{pair:14} vina={rec['vina_summary_min']} best={best} {rec['best_descriptor_summary_min']} Δ={rec['vina_minus_best_descriptor']} [{rec['delta_ci_lo']},{rec['delta_ci_hi']}]")
    write_csv(CANON / "descriptor_baselines.csv", rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
