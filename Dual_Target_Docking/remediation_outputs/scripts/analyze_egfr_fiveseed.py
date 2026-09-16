#!/usr/bin/env python3
"""Analyze EGFR five-seed Vina after canonical-box redock."""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SCORES = ROOT / "remediation_outputs/phase_fiveseed/scores_vina_mode1_fiveseed.csv"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
OUT = ROOT / "remediation_outputs/phase_fiveseed/multiseed_auroc_by_seed_corrected.csv"


def auc(pos, neg):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def main() -> int:
    if not SCORES.exists():
        print("waiting", SCORES)
        return 2
    cls = {r["panel_id"]: r["class"] for r in csv.DictReader(PANEL.open())}
    by = defaultdict(lambda: defaultdict(dict))
    with SCORES.open() as fh:
        for r in csv.DictReader(fh):
            if r["status"] not in {"ok", "cached", "production_corrected_box"}:
                continue
            try:
                by[int(r["seed"])][r["ligand"]][r["pdb"]] = -float(r["vina_mode1"])
            except (TypeError, ValueError):
                continue
    rows = []
    for seed, ligs in sorted(by.items()):
        recs = []
        for lig, d in ligs.items():
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
        da = auc([r["B"] for r in dual], [r["B"] for r in aonly])
        db = auc([r["A"] for r in dual], [r["A"] for r in bonly])
        smin = min(da, db)
        nei_auc = auc([r["mean"] for r in dual], [r["mean"] for r in nei])
        rows.append(
            {
                "pair": "EGFR/HER2",
                "seed": seed,
                "n_complete": len(recs),
                "n_dual": len(dual),
                "n_A_only": len(aonly),
                "n_B_only": len(bonly),
                "n_neither": len(nei),
                "auroc_dual_vs_A_only": round(da, 4),
                "auroc_dual_vs_B_only": round(db, 4),
                "summary_min": round(smin, 4),
                "auroc_dual_vs_neither_vina_mean": round(nei_auc, 4),
                "formulation_gap_neither_minus_summary_min": round(nei_auc - smin, 4),
            }
        )
        print(rows[-1])
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
