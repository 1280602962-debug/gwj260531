#!/usr/bin/env python3
"""Read the locked Table 2 bootstrap. Do not re-bootstrap.

Canonical sources:
- EGFR/HER2, AChE/BChE, PIK3CA/mTOR:
  unified_threshold_sensitivity_v2.csv, label_rule=theta_6.0
- F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD:
  five_pair_stack_v1/table2_comparable_theta6_v1.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
THREE = (
    ROOT
    / "data"
    / "jcim_strengthen_t0t1_v0"
    / "tables"
    / "unified_threshold_sensitivity_v2.csv"
)
FIVE = (
    ROOT
    / "data"
    / "jcim_chembl_universe_v0"
    / "local_track_b_v0"
    / "tables"
    / "five_pair_stack_v1"
    / "table2_comparable_theta6_v1.csv"
)

CANON = ROOT / "results" / "canonical" / "primary_summary_min.csv"


def _locked_from_canonical() -> dict[str, tuple[float, float, float]]:
    out = {}
    with CANON.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            out[row["pair"]] = (
                float(row["summary_min"]),
                float(row["ci_lo"]),
                float(row["ci_hi"]),
            )
    return out


LOCKED = None  # filled from results/canonical/primary_summary_min.csv in main()

THREE_PAIRS = {"EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"}


def _load_three() -> dict[str, tuple[float, float, float]]:
    out = {}
    with THREE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("label_rule") != "theta_6.0":
                continue
            out[row["pair"]] = (
                float(row["pocket_matched_summary_min"]),
                float(row["ci_lo"]),
                float(row["ci_hi"]),
            )
    return out


def _load_five() -> dict[str, tuple[float, float, float]]:
    out = {}
    with FIVE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            out[row["pair"]] = (
                float(row["summary_min"]),
                float(row["ci_lo"]),
                float(row["ci_hi"]),
            )
    return out


def load_table2() -> dict[str, tuple[float, float, float]]:
    got = _load_three()
    got.update(_load_five())
    return got


def main() -> int:
    global LOCKED
    LOCKED = _locked_from_canonical()
    got = load_table2()
    missing = [pair for pair in LOCKED if pair not in got]
    if missing:
        raise SystemExit(f"missing Table 2 pairs: {missing}")
    extra = sorted(set(got) - set(LOCKED) - {"PIK3CA/PIK3CB"})
    if extra:
        raise SystemExit(f"unexpected Table 2 pairs: {extra}")
    for pair, expected in LOCKED.items():
        point, lo, hi = got[pair]
        if any(abs(a - b) > 5e-4 for a, b in zip((point, lo, hi), expected)):
            raise SystemExit(f"{pair}: {got[pair]} != canonical {expected}")
        src = THREE if pair in THREE_PAIRS else FIVE
        print(
            f"{pair}: summary_min={point:.4f} [{lo:.4f}, {hi:.4f}]  "
            f"(figure source {src.relative_to(ROOT)}; lock=results/canonical/primary_summary_min.csv)"
        )
    print("Table 2 lock: PASS (figure-source CSVs match class-stratified canonical)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
