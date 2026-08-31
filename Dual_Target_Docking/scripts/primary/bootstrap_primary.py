#!/usr/bin/env python3
"""Read the locked Table 2 bootstrap. Do not re-bootstrap.

Canonical source: unified_threshold_sensitivity_v2.csv, label_rule=theta_6.0.
Figure 3, Table S4, the validator, and MASTER must all read this file.
A second hash-offset bootstrap (pocket_matched_directional_v1.csv) is not Table 2.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANONICAL = (
    ROOT
    / "data"
    / "jcim_strengthen_t0t1_v0"
    / "tables"
    / "unified_threshold_sensitivity_v2.csv"
)

LOCKED = {
    "EGFR/HER2": (0.4297, 0.2818, 0.5775),
    "AChE/BChE": (0.6058, 0.4370, 0.7303),
    "PIK3CA/PIK3CB": (0.5000, 0.3502, 0.6495),
    "PIK3CA/mTOR": (0.6921, 0.4702, 0.8133),
}


def load_table2():
    rows = list(csv.DictReader(CANONICAL.open(encoding="utf-8")))
    out = {}
    for r in rows:
        if r.get("label_rule") != "theta_6.0":
            continue
        out[r["pair"]] = (
            float(r["pocket_matched_summary_min"]),
            float(r["ci_lo"]),
            float(r["ci_hi"]),
        )
    return out


def main() -> int:
    got = load_table2()
    for pair, expected in LOCKED.items():
        point, lo, hi = got[pair]
        if any(abs(a - b) > 5e-4 for a, b in zip((point, lo, hi), expected)):
            raise SystemExit(f"{pair}: {got[pair]} != locked {expected}")
        print(
            f"{pair}: summary_min={point:.4f} [{lo:.4f}, {hi:.4f}]  "
            f"(canonical {CANONICAL.relative_to(ROOT)})"
        )
    print("Table 2 lock: PASS (read-only; no new bootstrap)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
