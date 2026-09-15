#!/usr/bin/env python3
"""Independent GNINA metrics for corrected EGFR/HER2 boxes."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SCORES = ROOT / "remediation_outputs/phase_gnina_independent/gnina_dock_scores_EGFR_HER2.csv"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
OUT = ROOT / "remediation_outputs/phase_gnina_independent/independent_dock_formulation_EGFR_HER2.csv"
N_BOOT = 2000
SEED = 20260729


def auc(pos, neg):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def boot(pos, neg, rng):
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    vals = []
    for _ in range(N_BOOT):
        try:
            vals.append(auc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)))
        except ValueError:
            continue
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> int:
    if not SCORES.exists():
        print("waiting", SCORES)
        return 2
    cls = {r["panel_id"]: r["class"] for r in csv.DictReader(PANEL.open())}
    by = {}
    with SCORES.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("status") not in {"ok", "cached"}:
                continue
            try:
                by.setdefault(r["ligand"], {})[r["target"]] = -float(r["gnina_mode1"])
            except (TypeError, ValueError, KeyError):
                continue
    recs = []
    for lig, d in by.items():
        if "3POZ" not in d or "3RCD" not in d:
            continue
        recs.append(
            {
                "cls": cls[lig],
                "A": d["3POZ"],
                "B": d["3RCD"],
                "mean": 0.5 * (d["3POZ"] + d["3RCD"]),
            }
        )
    dual = [r for r in recs if r["cls"] == "dual"]
    aonly = [r for r in recs if r["cls"] == "A_only"]
    bonly = [r for r in recs if r["cls"] == "B_only"]
    nei = [r for r in recs if r["cls"] == "neither"]
    rng = np.random.default_rng(SEED)
    da = auc([r["B"] for r in dual], [r["B"] for r in aonly])
    db = auc([r["A"] for r in dual], [r["A"] for r in bonly])
    da_lo, da_hi = boot([r["B"] for r in dual], [r["B"] for r in aonly], rng)
    db_lo, db_hi = boot([r["A"] for r in dual], [r["A"] for r in bonly], rng)
    smin = min(da, db)
    nei_auc = auc([r["mean"] for r in dual], [r["mean"] for r in nei])
    nei_lo, nei_hi = boot([r["mean"] for r in dual], [r["mean"] for r in nei], rng)
    rows = [
        {
            "pair": "EGFR/HER2",
            "engine": "gnina_dock_mode1",
            "formulation": "dualfourclass_directional",
            "contrast": "D_vs_A_pocketB",
            "n_pos": len(dual),
            "n_neg": len(aonly),
            "auroc": da,
            "ci_lo": da_lo,
            "ci_hi": da_hi,
        },
        {
            "pair": "EGFR/HER2",
            "engine": "gnina_dock_mode1",
            "formulation": "dualfourclass_directional",
            "contrast": "D_vs_B_pocketA",
            "n_pos": len(dual),
            "n_neg": len(bonly),
            "auroc": db,
            "ci_lo": db_lo,
            "ci_hi": db_hi,
        },
        {
            "pair": "EGFR/HER2",
            "engine": "gnina_dock_mode1",
            "formulation": "dualfourclass_directional",
            "contrast": "summary_min",
            "n_pos": len(dual),
            "n_neg": min(len(aonly), len(bonly)),
            "auroc": smin,
            "ci_lo": "",
            "ci_hi": "",
        },
        {
            "pair": "EGFR/HER2",
            "engine": "gnina_dock_mode1",
            "formulation": "conventional_dual_vs_neither",
            "contrast": "D_vs_neither_mean",
            "n_pos": len(dual),
            "n_neg": len(nei),
            "auroc": nei_auc,
            "ci_lo": nei_lo,
            "ci_hi": nei_hi,
        },
    ]
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(
        f"n={len(recs)} D/A/B/N={len(dual)}/{len(aonly)}/{len(bonly)}/{len(nei)} "
        f"DA={da:.4f} DB={db:.4f} smin={smin:.4f} DvsN={nei_auc:.4f}"
    )
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
