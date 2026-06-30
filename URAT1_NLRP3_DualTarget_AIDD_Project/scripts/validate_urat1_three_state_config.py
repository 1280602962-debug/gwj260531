#!/usr/bin/env python3
"""Cross-check URAT1 three-state PDB mapping across yaml and CSV."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
YAML_PATH = PROJECT_ROOT / "config" / "docking_ensemble.yaml"
CSV_PATH = PROJECT_ROOT / "data" / "structures" / "docking_ensemble_pdb.csv"

EXPECTED = {
    "inward_open": "9DKB",
    "occluded": "9B1K",
    "outward_open": "9B1L",
}


def main() -> int:
    with open(YAML_PATH) as f:
        cfg = yaml.safe_load(f)
    three = cfg["urat1_ensemble"]["three_state_primary"]
    errors: list[str] = []

    for state, pdb in EXPECTED.items():
        got = three[state]["pdb_id"]
        if got != pdb:
            errors.append(f"yaml three_state_primary.{state}: expected {pdb}, got {got}")

    df = pd.read_csv(CSV_PATH)
    urat1 = df[df["target"] == "URAT1"]
    for state, pdb in EXPECTED.items():
        role = f"three_state_primary.{state.split('_')[0] if state == 'inward_open' else state}"
        if state == "inward_open":
            role = "three_state_primary.inward"
        rows = urat1[urat1["pdb_id"] == pdb]
        if rows.empty:
            errors.append(f"csv missing PDB {pdb} for {state}")
        elif not str(rows.iloc[0]["ensemble_role"]).startswith("three_state_primary"):
            errors.append(f"csv {pdb} ensemble_role not three_state_primary: {rows.iloc[0]['ensemble_role']}")

    jdz = urat1[urat1["pdb_id"] == "9JDZ"]
    if not jdz.empty and "occluded" in str(jdz.iloc[0]["notes"]).lower():
        errors.append("9JDZ notes must NOT claim occluded/outward")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print(f"  - {e}")
        return 1

    print("URAT1 three-state config OK:")
    for state, pdb in EXPECTED.items():
        print(f"  {state}: {pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
