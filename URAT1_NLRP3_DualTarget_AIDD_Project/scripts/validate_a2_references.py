#!/usr/bin/env python3
"""Validate A2 pose selection on known URAT1 reference ligands (forced-recovery SDFs)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from c1_acid_pose_selection import (  # noqa: E402
    ARG_THRESH_A,
    evaluate_urat1_acid_sdf,
    load_ref_centroid,
)
from parse_c1_sdf_readouts import load_poses, min_acid_arg_dist  # noqa: E402

REFS = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
SDF_DIR = PROJECT_ROOT / "data/campaigns/c1/03_forced_recovery/urat1_9dkb/seed42"
OUT = PROJECT_ROOT / "data/campaigns/c1/03_forced_recovery/a2_reference_validation_seed42"


def best_acid_among_modes(sdf: Path, arg_json: Path) -> float | None:
    arg = json.loads(arg_json.read_text())
    best = None
    for pose in load_poses(sdf):
        d = min_acid_arg_dist(pose, arg["atoms"])
        if d is not None:
            best = d if best is None else min(best, d)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    arg_json = REFS / "arg477_coords.json"
    ref_com = load_ref_centroid(REFS / "lesinurad_crystal_ref.sdf")

    rows = []
    for sdf in sorted(SDF_DIR.glob("*_out.sdf")):
        name = sdf.stem.replace("_out", "")
        for rule in ("a1", "a2"):
            r = evaluate_urat1_acid_sdf(sdf, arg_json, ref_com, name, args.seed, rule=rule)
            r["ligand"] = name
            r["best_carboxylate_arg_among_modes_A"] = best_acid_among_modes(sdf, arg_json)
            rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a2_reference_validation.csv", index=False)

    carboxylate_refs = [
        "lesinurad",
        "verinurad",
        "probenecid",
        "puliginurad",
        "SHR-4640",
        "GSK-3008348",
    ]
    sub = df[df.ligand.isin(carboxylate_refs)].copy()
    summary = {
        "seed": args.seed,
        "arg_threshold_A": ARG_THRESH_A,
        "carboxylate_references": carboxylate_refs,
        "a1_pass": int(sub[sub.pose_selection_rule == "a1"]["keep_urat1_acid"].sum()),
        "a2_pass": int(sub[sub.pose_selection_rule == "a2"]["keep_urat1_acid"].sum()),
        "n_carboxylate_refs": len(carboxylate_refs),
        "trigger_for_A2": "a1 fails to recover known URAT1 carboxylate ligands despite sampling near-crystal poses",
        "outputs": {
            "csv": str((OUT / "a2_reference_validation.csv").relative_to(PROJECT_ROOT)),
        },
    }
    (OUT / "a2_reference_validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    md = [
        "# A2 reference validation (seed 42)",
        "",
        f"- Arg threshold: **{ARG_THRESH_A} Å**",
        f"- A1 carboxylate refs pass: **{summary['a1_pass']}/{summary['n_carboxylate_refs']}**",
        f"- A2 carboxylate refs pass: **{summary['a2_pass']}/{summary['n_carboxylate_refs']}**",
        "",
        "| ligand | rule | keep | Arg477 (Å) | best-acid among modes (Å) | status |",
        "|--------|------|-----:|-----------:|--------------------------:|--------|",
    ]
    for lig in carboxylate_refs:
        for rule in ("a1", "a2"):
            r = sub[(sub.ligand == lig) & (sub.pose_selection_rule == rule)].iloc[0]
            md.append(
                f"| {lig} | {rule} | {r['keep_urat1_acid']} | "
                f"{r['acid_arg477_min_A']:.2f} | {r['best_carboxylate_arg_among_modes_A']:.2f} | "
                f"{r['pose_selection_status']} |"
            )
    (OUT / "a2_reference_validation.md").write_text("\n".join(md) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
