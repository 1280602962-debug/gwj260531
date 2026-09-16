#!/usr/bin/env python3
"""EGFR/HER2 Vina-pose GNINA CNN metrics after canonical-box redock.

Primary readout: best-of-9 cnn_affinity.
Sensitivity: best-of-9 cnn_score (do not promote if one pair looks better).
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
LONG = ROOT / "remediation_outputs/phase_rtm/tables/scores_gnina_long.csv"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
OUT = ROOT / "remediation_outputs/phase_rtm/gnina_cnn_directional_metrics.csv"
BEST_AFF = ROOT / "remediation_outputs/phase_rtm/gnina_cnn_affinity_best.csv"


def auc(pos, neg):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    return float(roc_auc_score(y, s))


def best_of_k(rows, field):
    best = {}
    for r in rows:
        if r.get("status") != "success":
            continue
        try:
            val = float(r[field])
        except (TypeError, ValueError, KeyError):
            continue
        key = (r["ligand"], r["target"])
        prev = best.get(key)
        if prev is None or val > prev[0]:
            best[key] = (val, r.get("mode", ""))
    return best


def directional(cls, by_lig):
    recs = []
    for lig, pockets in by_lig.items():
        if lig not in cls or "3POZ" not in pockets or "3RCD" not in pockets:
            continue
        recs.append({"cls": cls[lig], "A": pockets["3POZ"], "B": pockets["3RCD"]})
    dual = [r for r in recs if r["cls"] == "dual"]
    aonly = [r for r in recs if r["cls"] == "A_only"]
    bonly = [r for r in recs if r["cls"] == "B_only"]
    nei = [r for r in recs if r["cls"] == "neither"]
    da = auc([r["B"] for r in dual], [r["B"] for r in aonly])
    db = auc([r["A"] for r in dual], [r["A"] for r in bonly])
    mean_d = [0.5 * (r["A"] + r["B"]) for r in dual]
    mean_n = [0.5 * (r["A"] + r["B"]) for r in nei]
    return {
        "n": len(recs),
        "n_dual": len(dual),
        "n_A": len(aonly),
        "n_B": len(bonly),
        "n_neither": len(nei),
        "AUROC_D_vs_A_pocketB": da,
        "AUROC_D_vs_B_pocketA": db,
        "summary_min": min(da, db),
        "AUROC_mean_D_vs_neither": auc(mean_d, mean_n),
    }


def main() -> int:
    if not LONG.exists():
        print("waiting", LONG)
        return 2
    cls = {r["panel_id"]: r["class"] for r in csv.DictReader(PANEL.open())}
    rows = list(csv.DictReader(LONG.open()))
    out_rows = []
    best_aff = best_of_k(rows, "cnn_affinity")
    with BEST_AFF.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand", "target", "mode", "cnn_affinity"])
        w.writeheader()
        for (lig, tgt), (val, mode) in sorted(best_aff.items()):
            w.writerow({"ligand": lig, "target": tgt, "mode": mode, "cnn_affinity": val})
    for field, channel in (
        ("cnn_affinity", "gnina_cnn_affinity_best_of_9"),
        ("cnn_score", "gnina_cnn_score_best_of_9_sensitivity"),
    ):
        best = best_of_k(rows, field)
        by_lig = defaultdict(dict)
        for (lig, tgt), (val, _mode) in best.items():
            by_lig[lig][tgt] = val
        m = directional(cls, by_lig)
        out_rows.append({"pair": "EGFR/HER2", "channel": channel, **m})
        print(
            f"{channel} n={m['n']} D/A/B/N={m['n_dual']}/{m['n_A']}/{m['n_B']}/{m['n_neither']} "
            f"DA={m['AUROC_D_vs_A_pocketB']:.4f} DB={m['AUROC_D_vs_B_pocketA']:.4f} "
            f"smin={m['summary_min']:.4f} DvsN={m['AUROC_mean_D_vs_neither']:.4f}"
        )
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print("wrote", OUT)
    print("wrote", BEST_AFF)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
