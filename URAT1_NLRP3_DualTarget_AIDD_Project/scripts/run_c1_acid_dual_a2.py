#!/usr/bin/env python3
"""Rescore or re-dock clinical Acid pool under Amendment A2 (geometry-first URAT1)."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from c1_acid_pose_selection import (  # noqa: E402
    ARG_THRESH_A,
    evaluate_nlrp3_pose_sdf,
    evaluate_urat1_acid_sdf,
    load_ref_centroid,
)

REFS = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
CRYSTAL_COM = {
    "urat1_9dkb": REFS / "lesinurad_crystal_ref.sdf",
    "nlrp3_7alv": REFS / "NP3-146_RM5_crystal_ref.sdf",
}


def aggregate_dual(rows: list[dict], seed: int, out_dir: Path) -> None:
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / f"acid_pose_metrics_a2_seed{seed}.csv", index=False)
    u = df[df.target == "urat1_9dkb"][
        ["ligand_id", "keep_urat1_acid", "acid_arg477_min_A", "CNNscore", "CNNaffinity", "pose_selection_status"]
    ].rename(columns={"CNNscore": "u_CNNscore", "CNNaffinity": "u_CNNaffinity"})
    n = df[df.target == "nlrp3_7alv"][
        ["ligand_id", "keep_nlrp3_pose", "CNNscore", "CNNaffinity"]
    ].rename(columns={"CNNscore": "n_CNNscore", "CNNaffinity": "n_CNNaffinity"})
    dual = u.merge(n, on="ligand_id", how="outer")
    dual["keep_dual_acid_geometry"] = dual["keep_urat1_acid"].fillna(False) & dual["keep_nlrp3_pose"].fillna(False)
    dual.to_csv(out_dir / f"acid_dual_keep_a2_seed{seed}.csv", index=False)
    summary = {
        "seed": seed,
        "amendment": "A2",
        "pose_selection_urat1": "geometry_first_then_cnnscore",
        "n_ligands": int(dual.ligand_id.nunique()),
        "n_keep_urat1_arg": int(dual.keep_urat1_acid.fillna(False).sum()),
        "n_keep_nlrp3_pose": int(dual.keep_nlrp3_pose.fillna(False).sum()),
        "n_keep_dual": int(dual.keep_dual_acid_geometry.sum()),
        "arg_threshold_A": ARG_THRESH_A,
        "percentile_used": False,
    }
    (out_dir / f"acid_dual_summary_a2_seed{seed}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--engine-config", type=Path, default=PROJECT_ROOT / "config/docking_c1_cpu.yaml")
    ap.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--seeds", nargs="*", type=int, default=None)
    ap.add_argument("--metrics-only", action="store_true")
    ap.add_argument("--sdf-source-dir", type=Path, default=None, help="read SDFs from another run dir (e.g. A1 for seed42)")
    args = ap.parse_args()

    seeds = args.seeds if args.seeds else [args.seed]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(args.manifest)
    man = man[man.status == "prepared"]
    arg_json = REFS / "arg477_coords.json"
    ref_com_u = load_ref_centroid(CRYSTAL_COM["urat1_9dkb"])
    ref_com_n = load_ref_centroid(CRYSTAL_COM["nlrp3_7alv"])
    eng = yaml.safe_load(args.engine_config.read_text())

    for seed in seeds:
        rows = []
        sdf_root = args.sdf_source_dir if args.sdf_source_dir and seed == 42 else args.output_dir
        for _, r in man.iterrows():
            rid = r["repurposing_id"]
            for tkey in ("urat1_9dkb", "nlrp3_7alv"):
                out_sdf = sdf_root / tkey / f"seed{seed}" / f"{rid}_out.sdf"
                metrics_sdf = args.output_dir / tkey / f"seed{seed}" / f"{rid}_out.sdf"
                if not out_sdf.exists() or out_sdf.stat().st_size == 0:
                    rows.append(
                        {
                            "ligand_id": rid,
                            "target": tkey,
                            "seed": seed,
                            "error": "missing_sdf",
                            "keep_urat1_acid": False,
                            "keep_nlrp3_pose": False,
                        }
                    )
                    continue
                if tkey == "urat1_9dkb":
                    rows.append(
                        evaluate_urat1_acid_sdf(out_sdf, arg_json, ref_com_u, rid, seed, rule="a2")
                    )
                else:
                    rows.append(evaluate_nlrp3_pose_sdf(out_sdf, ref_com_n, rid, seed))
                # record resolved sdf path for traceability
                rows[-1]["sdf"] = str(metrics_sdf if metrics_sdf.exists() else out_sdf)
        aggregate_dual(rows, seed, args.output_dir)


if __name__ == "__main__":
    main()
