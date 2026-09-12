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
import math
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


def parse_finite_energy(raw: str | None) -> tuple[float | None, str | None]:
    """Return (-energy, None) only when the raw score is a finite number."""
    if raw is None or str(raw).strip() == "":
        return None, "empty"
    text = str(raw).strip()
    try:
        energy = float(text)
    except (TypeError, ValueError):
        return None, "non_numeric"
    if not math.isfinite(energy):
        return None, "not_finite"
    return -energy, None


def load_seed(seed: int):
    path = MS / f"scores_vina_mode1_seed{seed}.csv"
    if not path.exists() and seed == 20260727:
        path = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
    rows = defaultdict(dict)
    exclusions = []
    for r in csv.DictReader(path.open(encoding="utf-8", newline="")):
        score, reason = parse_finite_energy(r.get("mode1_energy"))
        if reason is not None:
            exclusions.append(
                {
                    "seed": seed,
                    "pair": r.get("pair", ""),
                    "ligand": r.get("ligand", ""),
                    "target": r.get("target", ""),
                    "raw_mode1_energy": r.get("mode1_energy", ""),
                    "reason": reason,
                }
            )
            continue
        rows[(r["pair"], r["ligand"])][r["target"]] = score
    return rows, exclusions


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


def _check_finite_parser() -> None:
    assert parse_finite_energy("") == (None, "empty")
    assert parse_finite_energy(None) == (None, "empty")
    assert parse_finite_energy("nan") == (None, "not_finite")
    assert parse_finite_energy("NaN") == (None, "not_finite")
    assert parse_finite_energy("inf") == (None, "not_finite")
    assert parse_finite_energy("-inf") == (None, "not_finite")
    assert parse_finite_energy("not-a-number") == (None, "non_numeric")
    score, reason = parse_finite_energy("-8.2")
    assert reason is None and score == 8.2


def main() -> int:
    _check_finite_parser()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=LOCAL / "tables" / "multiseed_fixed_membership_v1.csv")
    ap.add_argument(
        "--membership-output",
        type=Path,
        default=LOCAL / "tables" / "multiseed_fixed_membership_ids_v1.csv",
    )
    ap.add_argument(
        "--exclusion-output",
        type=Path,
        default=LOCAL / "tables" / "multiseed_fixed_membership_exclusions_v1.csv",
    )
    args = ap.parse_args()
    labels = load_labels()
    by_seed = {}
    exclusions = []
    for seed in SEEDS:
        tab, skipped = load_seed(seed)
        by_seed[seed] = tab
        exclusions.extend(skipped)
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
    excl_fields = ["seed", "pair", "ligand", "target", "raw_mode1_energy", "reason"]
    with args.exclusion_output.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=excl_fields, lineterminator="\n")
        w.writeheader()
        w.writerows(exclusions)
    print("wrote", args.output)
    print("wrote", args.membership_output)
    print("wrote", args.exclusion_output, "n_excluded", len(exclusions))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
