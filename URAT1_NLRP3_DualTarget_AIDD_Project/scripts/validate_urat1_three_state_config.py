#!/usr/bin/env python3
"""Cross-check URAT1 three-state PDB mapping across yaml, CSV, and targets.yaml."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
YAML_PATH = PROJECT_ROOT / "config" / "docking_ensemble.yaml"
CSV_PATH = PROJECT_ROOT / "data" / "structures" / "docking_ensemble_pdb.csv"
TARGETS_PATH = PROJECT_ROOT / "config" / "targets.yaml"

EXPECTED = {
    "inward_open": "9DKB",
    "occluded": "9B1K",
    "outward_open": "9B1L",
}

EXPECTED_STATES = set(EXPECTED.keys())
PRIMARY_ROLES = {state: f"three_state_primary.{state}" for state in EXPECTED}


def main() -> int:
    with open(YAML_PATH) as f:
        cfg = yaml.safe_load(f)
    with open(TARGETS_PATH) as f:
        targets = yaml.safe_load(f)

    urat1_cfg = cfg["urat1_ensemble"]
    three = urat1_cfg["three_state_primary"]
    errors: list[str] = []

    # 1) three_state_primary PDB IDs
    for state, pdb in EXPECTED.items():
        got = three[state]["pdb_id"]
        if got != pdb:
            errors.append(f"yaml three_state_primary.{state}: expected {pdb}, got {got}")
        if three[state].get("state") != state:
            errors.append(
                f"yaml three_state_primary.{state}.state mismatch: {three[state].get('state')}"
            )

    # 2) pdb_structures must list the same three primaries with matching roles
    pdb_by_id = {s["pdb_id"]: s for s in urat1_cfg["pdb_structures"]}
    for state, pdb in EXPECTED.items():
        if pdb not in pdb_by_id:
            errors.append(f"yaml pdb_structures missing primary {pdb}")
            continue
        entry = pdb_by_id[pdb]
        if entry.get("state") != state:
            errors.append(f"yaml pdb_structures {pdb} state: expected {state}, got {entry.get('state')}")
        if entry.get("ensemble_role") != PRIMARY_ROLES[state]:
            errors.append(
                f"yaml pdb_structures {pdb} role: expected {PRIMARY_ROLES[state]}, "
                f"got {entry.get('ensemble_role')}"
            )

    # 3) teacher_gate redock grid must be inward primary
    gate = urat1_cfg.get("teacher_gate", {})
    inward_pdb = EXPECTED["inward_open"]
    if gate.get("lesinurad_redock_pdb") != inward_pdb:
        errors.append(
            f"teacher_gate.lesinurad_redock_pdb must be {inward_pdb}, got {gate.get('lesinurad_redock_pdb')}"
        )

    # 4) pi_scoring formula consistency
    pi = urat1_cfg.get("pi_scoring", {})
    if "pi_in + pi_occ - pi_out" not in str(pi.get("formula", "")):
        errors.append("pi_scoring.formula must encode S_pi = pi_in + pi_occ - pi_out")

    # 5) CSV mirror of primaries
    df = pd.read_csv(CSV_PATH)
    urat1 = df[df["target"] == "URAT1"]
    for state, pdb in EXPECTED.items():
        rows = urat1[urat1["pdb_id"] == pdb]
        if rows.empty:
            errors.append(f"csv missing PDB {pdb} for {state}")
            continue
        row = rows.iloc[0]
        if row["conformation_state"] != state:
            errors.append(f"csv {pdb} conformation_state: expected {state}, got {row['conformation_state']}")
        if row["ensemble_role"] != PRIMARY_ROLES[state]:
            errors.append(
                f"csv {pdb} ensemble_role: expected {PRIMARY_ROLES[state]}, got {row['ensemble_role']}"
            )

    # 6) 9JDZ is supplementary inward only — never a transport-cycle state grid
    jdz = urat1[urat1["pdb_id"] == "9JDZ"]
    if jdz.empty:
        errors.append("csv missing supplementary 9JDZ entry")
    else:
        row = jdz.iloc[0]
        if row["conformation_state"] != "inward_open":
            errors.append(f"9JDZ must be inward_open, got {row['conformation_state']}")
        notes = str(row["notes"]).lower()
        if "occluded" in notes or "outward" in notes:
            errors.append("9JDZ notes must NOT claim occluded/outward as usable grids")
        if str(row["ensemble_role"]).startswith("three_state_primary"):
            errors.append("9JDZ must not be tagged three_state_primary")

    primary_pdbs = set(EXPECTED.values())
    if primary_pdbs & {"9JDZ"}:
        errors.append("9JDZ must not appear in three_state_primary")

    # 7) targets.yaml benchmark redock PDBs
    urat1_bench = {b["name"]: b.get("pdb_complex") for b in targets["benchmark_compounds"]["urat1"]}
    if urat1_bench.get("lesinurad") != inward_pdb:
        errors.append(f"targets.yaml lesinurad pdb_complex must be {inward_pdb}")
    if urat1_bench.get("benzbromarone") != "9DKA":
        errors.append("targets.yaml benzbromarone pdb_complex must be 9DKA (not 9DKB)")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print(f"  - {e}")
        return 1

    print("URAT1 three-state config OK:")
    for state, pdb in EXPECTED.items():
        print(f"  {state}: {pdb} ({PRIMARY_ROLES[state]})")
    print(f"  teacher_gate redock: {gate.get('lesinurad_redock_pdb')}")
    print(f"  pi_scoring: {pi.get('formula')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
