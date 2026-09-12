#!/usr/bin/env python3
"""Fixed-membership five-seed summary for Track B pairs only.

A ligand enters the intersection only if every seed has a finite mode-1
score on both ends of that pair. n_intersection includes neither.
Directional AUROCs use only dual / A-only / B-only.

This isolates search-seed variation from membership changes caused by
timeouts. It is not a label-source or assay-heterogeneity check.
EGFR/HER2, AChE/BChE, and PIK3CA/mTOR are not in this table.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
MS = LOCAL / "tables" / "multiseed"
PANELS = ROOT / "tables" / "track_b_panels"

SEEDS = [20260727, 20260811, 20260812, 20260813, 20260814]
PAIRS = {
    "F2/F10": ("4UDW", "2JKH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}


def load_labels():
    lab = {}
    for p in PANELS.glob("panel_*_v1.csv"):
        for r in csv.DictReader(p.open(encoding="utf-8", newline="")):
            # Primary analysis uses θ=6.0 class; Track B panels store it as class
            # and, when present, theta6_class (identical on current panels).
            lab[(r["pair"], r["panel_id"])] = r.get("theta6_class") or r["class"]
    return lab


def load_seed(seed: int):
    path = MS / f"scores_vina_mode1_seed{seed}.csv"
    if not path.exists() and seed == 20260727:
        path = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
    rows = defaultdict(dict)
    for r in csv.DictReader(path.open(encoding="utf-8", newline="")):
        if r.get("mode1_energy") in (None, ""):
            continue
        rows[(r["pair"], r["ligand"])][r["target"]] = -float(r["mode1_energy"])
    return rows


def auroc(y, s):
    y = np.asarray(y)
    s = np.asarray(s, float)
    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def weaker(a: float, b: float) -> str:
    if np.isnan(a) or np.isnan(b):
        return ""
    if a < b:
        return "D_vs_A"
    if b < a:
        return "D_vs_B"
    return "tie"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=LOCAL / "tables" / "multiseed_fixed_membership_v1.csv")
    ap.add_argument(
        "--membership-output",
        type=Path,
        default=LOCAL / "tables" / "multiseed_fixed_membership_ids_v1.csv",
    )
    args = ap.parse_args()
    labels = load_labels()
    by_seed = {s: load_seed(s) for s in SEEDS}
    out = []
    members = []
    for pair, (rec_a, rec_b) in PAIRS.items():
        complete = None
        for seed, tab in by_seed.items():
            keys = {k[1] for k, v in tab.items() if k[0] == pair and rec_a in v and rec_b in v}
            complete = keys if complete is None else (complete & keys)
        complete = sorted(complete or [])
        n_by_class = defaultdict(int)
        for lig in complete:
            n_by_class[labels.get((pair, lig), "")] += 1
            members.append(
                {
                    "pair": pair,
                    "ligand_id": lig,
                    "class": labels.get((pair, lig), ""),
                    "in_directional_auroc": int(labels.get((pair, lig), "") in {"dual", "A_only", "B_only"}),
                }
            )
        prod_weak = ""
        for seed, tab in by_seed.items():
            da_y, da_s, db_y, db_s = [], [], [], []
            for lig in complete:
                cls = labels.get((pair, lig), "")
                sc = tab[(pair, lig)]
                if cls in {"dual", "A_only"}:
                    da_y.append(1 if cls == "dual" else 0)
                    da_s.append(sc[rec_b])
                if cls in {"dual", "B_only"}:
                    db_y.append(1 if cls == "dual" else 0)
                    db_s.append(sc[rec_a])
            auc_da = auroc(da_y, da_s)
            auc_db = auroc(db_y, db_s)
            weak = weaker(auc_da, auc_db)
            if seed == SEEDS[0]:
                prod_weak = weak
            out.append(
                {
                    "pair": pair,
                    "seed": seed,
                    "n_intersection": len(complete),
                    "n_intersection_includes_neither": 1,
                    "n_neither": n_by_class.get("neither", 0),
                    "n_dual": n_by_class.get("dual", 0),
                    "n_A_only": n_by_class.get("A_only", 0),
                    "n_B_only": n_by_class.get("B_only", 0),
                    "n_directional": n_by_class.get("dual", 0)
                    + n_by_class.get("A_only", 0)
                    + n_by_class.get("B_only", 0),
                    "auc_D_vs_A": None if np.isnan(auc_da) else round(auc_da, 4),
                    "auc_D_vs_B": None if np.isnan(auc_db) else round(auc_db, 4),
                    "summary_min": None
                    if (np.isnan(auc_da) or np.isnan(auc_db))
                    else round(min(auc_da, auc_db), 4),
                    "weaker_arm": weak,
                    "weaker_switched_from_production": int(bool(prod_weak) and weak != prod_weak),
                    "label_source": "theta6_class_or_class",
                    "membership": "intersection_all_five_seeds",
                }
            )
        print(pair, "intersection", len(complete), dict(n_by_class))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    with args.membership_output.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(members[0]))
        w.writeheader()
        w.writerows(members)
    print("wrote", args.output)
    print("wrote", args.membership_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
