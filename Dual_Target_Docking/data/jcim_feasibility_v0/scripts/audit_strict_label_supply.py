#!/usr/bin/env python3
"""Hard-negative supply audit for the four-class dual-target task.

Uses the already-fetched ChEMBL pChEMBL dictionaries under
data/public_pair_selection/mols_*.json (no network, no docking) and asks, per
candidate target pair, how many molecules fall into each class under

  theta labelling      : active if pChEMBL >= 6.0
  strict/margin labels : dual  = both ends >= 6.5
                         A_only = A >= 6.5 and B <= 5.5
                         B_only = B >= 6.5 and A <= 5.5
                         gray   = measured on both ends but inside the buffer

The point is that panel size for a strict four-class evaluation is capped by
public-data supply of the *rarer* hard-negative class, not by docking budget.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data" / "public_pair_selection"
OUT = Path(__file__).resolve().parents[1] / "tables"

THETA = 6.0
HI = 6.5
LO = 5.5

PAIRS = [
    ("EGFR/HER2", "EGFR", "HER2"),
    ("PIK3CA/MTOR", "PIK3CA", "MTOR"),
    ("ACHE/BCHE", "ACHE", "BCHE"),
    ("PIK3CA/PIK3CB", "PIK3CA", "PIK3CB"),
    ("MCL1/BCL2L1_BclxL", "MCL1", "BCL2L1_BclxL"),
    ("MCL1/BCL2", "MCL1", "BCL2"),
    ("AKT1/RPS6KB1_p70S6K", "AKT1", "RPS6KB1_p70S6K"),
    ("BRD4/HDAC1", "BRD4", "HDAC1"),
    ("BRD4/HDAC6", "BRD4", "HDAC6"),
    ("JAK2/HDAC1", "JAK2", "HDAC1"),
    ("PARP1/MET", "PARP1", "MET"),
    ("CDK6/BRD4", "CDK6", "BRD4"),
]

MIN_HARDNEG_STRICT = 50


def load(target: str) -> dict[str, float]:
    with (SRC / f"mols_{target}.json").open() as fh:
        return json.load(fh)


def audit(a_vals: dict[str, float], b_vals: dict[str, float]) -> dict[str, int]:
    both = set(a_vals) & set(b_vals)
    rec = {
        "n_both_measured": len(both),
        "theta_dual": 0,
        "theta_A_only": 0,
        "theta_B_only": 0,
        "strict_dual": 0,
        "strict_A_only": 0,
        "strict_B_only": 0,
        "strict_neither": 0,
        "gray": 0,
    }
    for mol in both:
        x, y = a_vals[mol], b_vals[mol]
        if x >= THETA and y >= THETA:
            rec["theta_dual"] += 1
        elif x >= THETA:
            rec["theta_A_only"] += 1
        elif y >= THETA:
            rec["theta_B_only"] += 1

        if x >= HI and y >= HI:
            rec["strict_dual"] += 1
        elif x >= HI and y <= LO:
            rec["strict_A_only"] += 1
        elif y >= HI and x <= LO:
            rec["strict_B_only"] += 1
        elif x <= LO and y <= LO:
            rec["strict_neither"] += 1
        else:
            rec["gray"] += 1
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, a, b in PAIRS:
        try:
            rec = audit(load(a), load(b))
        except FileNotFoundError:
            continue
        rec["pair"] = name
        rec["gray_frac"] = (
            round(rec["gray"] / rec["n_both_measured"], 3) if rec["n_both_measured"] else None
        )
        rec["min_strict_hardneg"] = min(rec["strict_A_only"], rec["strict_B_only"])
        rec["supports_strict_panel"] = rec["min_strict_hardneg"] >= MIN_HARDNEG_STRICT
        rows.append(rec)

    rows.sort(key=lambda r: -r["min_strict_hardneg"])
    fields = [
        "pair",
        "n_both_measured",
        "theta_dual",
        "theta_A_only",
        "theta_B_only",
        "strict_dual",
        "strict_A_only",
        "strict_B_only",
        "strict_neither",
        "gray",
        "gray_frac",
        "min_strict_hardneg",
        "supports_strict_panel",
    ]
    path = OUT / "strict_label_supply.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {path}")
    print(
        f"{'pair':22s} {'both':>6s} {'θ d/A/B':>16s} {'strict d/A/B':>16s} "
        f"{'gray':>6s} {'ok':>4s}"
    )
    for r in rows:
        print(
            f"{r['pair']:22s} {r['n_both_measured']:6d} "
            f"{r['theta_dual']:5d}/{r['theta_A_only']:4d}/{r['theta_B_only']:4d} "
            f"{r['strict_dual']:6d}/{r['strict_A_only']:4d}/{r['strict_B_only']:4d} "
            f"{(r['gray_frac'] or 0):6.2f} "
            f"{'Y' if r['supports_strict_panel'] else '-':>4s}"
        )
    n_ok = sum(r["supports_strict_panel"] for r in rows)
    print(
        f"\npairs supporting a strict four-class panel "
        f"(both hard-negative classes >= {MIN_HARDNEG_STRICT}): {n_ok}/{len(rows)}"
    )


if __name__ == "__main__":
    main()
