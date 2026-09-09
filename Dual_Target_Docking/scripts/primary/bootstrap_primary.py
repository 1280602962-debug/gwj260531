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

LOCKED = {
    "EGFR/HER2": (0.4297, 0.2818, 0.5775),
    "AChE/BChE": (0.6058, 0.4370, 0.7303),
    "PIK3CA/mTOR": (0.6921, 0.4702, 0.8133),
    "F2/F10": (0.3448, 0.2109, 0.4773),
    "JAK1/TYK2": (0.3649, 0.2306, 0.5030),
    "JAK1/JAK2": (0.5884, 0.4444, 0.7246),
    "PPARG/PPARA": (0.6492, 0.5045, 0.7508),
    "PPARA/PPARD": (0.4463, 0.2958, 0.5841),
}

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
            raise SystemExit(f"{pair}: {got[pair]} != locked {expected}")
        src = THREE if pair in THREE_PAIRS else FIVE
        print(
            f"{pair}: summary_min={point:.4f} [{lo:.4f}, {hi:.4f}]  "
            f"(canonical {src.relative_to(ROOT)})"
        )
    print("Table 2 lock: PASS (read-only; no new bootstrap; eight pairs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
