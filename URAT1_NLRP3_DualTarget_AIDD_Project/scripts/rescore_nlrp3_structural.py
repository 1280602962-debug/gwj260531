#!/usr/bin/env python3
"""Rescore existing Acid dual NLRP3 SDFs with structural metrics (IFP / overlap / key residues).

Does not call gnina. Reads SDFs from acid_dual (seed42) and/or acid_dual_a2 (seeds 43/44).
Writes nlrp3_structural_metrics and upgraded dual keep tables under acid_dual_a2/.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from c1_nlrp3_pose_metrics import (  # noqa: E402
    crystal_reference_ifp,
    evaluate_nlrp3_structural,
    load_key_map,
    load_receptor_heavy,
)


def find_sdf(rid: str, seed: int, roots: list[Path]) -> Path | None:
    for root in roots:
        p = root / "nlrp3_7alv" / f"seed{seed}" / f"{rid}_out.sdf"
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        type=Path,
        default=PROJECT_ROOT
        / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv",
    )
    ap.add_argument(
        "--a1-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual",
    )
    ap.add_argument(
        "--a2-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2",
    )
    ap.add_argument("--seeds", nargs="*", type=int, default=[42])
    args = ap.parse_args()

    man = pd.read_csv(args.manifest)
    man = man[man.status == "prepared"]
    key_map = load_key_map()
    ref_ifp, ref_heavy, ref_com = crystal_reference_ifp()
    receptor_heavy = load_receptor_heavy()
    roots = [args.a2_dir, args.a1_dir]

    for seed in args.seeds:
        rows = []
        for _, r in man.iterrows():
            rid = r["repurposing_id"]
            sdf = find_sdf(rid, seed, roots)
            if sdf is None:
                rows.append(
                    {
                        "ligand_id": rid,
                        "seed": seed,
                        "error": "missing_sdf",
                        "keep_nlrp3_pose": False,
                        "keep_nlrp3_structural": False,
                    }
                )
                continue
            rows.append(
                evaluate_nlrp3_structural(
                    sdf,
                    rid,
                    seed,
                    key_map=key_map,
                    ref_heavy=ref_heavy,
                    ref_com=ref_com,
                    receptor_heavy=receptor_heavy,
                    ref_ifp=ref_ifp,
                )
            )
        df = pd.DataFrame(rows)
        out = args.a2_dir / f"nlrp3_structural_metrics_seed{seed}.csv"
        args.a2_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)

        # join with URAT1 A2 keep if present
        u_path = args.a2_dir / f"acid_dual_keep_a2_seed{seed}.csv"
        if not u_path.exists() and seed == 42:
            # build from A2 metrics if needed
            u_path = args.a2_dir / "acid_dual_keep_a2_seed42.csv"
        summary = {
            "seed": seed,
            "n": int(len(df)),
            "n_keep_loose": int(df.keep_nlrp3_pose.fillna(False).sum()),
            "n_keep_structural": int(df.keep_nlrp3_structural.fillna(False).sum()),
            "mean_ifp_jaccard": float(df.ifp_jaccard_vs_np3146.dropna().mean())
            if "ifp_jaccard_vs_np3146" in df
            else None,
            "mean_overlap": float(df.pocket_overlap_frac.dropna().mean())
            if "pocket_overlap_frac" in df
            else None,
        }
        if u_path.exists():
            u = pd.read_csv(u_path)
            dual = u.merge(
                df[
                    [
                        "ligand_id",
                        "keep_nlrp3_pose",
                        "keep_nlrp3_structural",
                        "pocket_overlap_frac",
                        "ifp_jaccard_vs_np3146",
                        "n_key_contacts",
                        "key_recovery_frac",
                    ]
                ],
                on="ligand_id",
                how="left",
                suffixes=("_old", ""),
            )
            # prefer structural columns from df
            if "keep_nlrp3_pose_old" in dual.columns:
                dual["keep_nlrp3_pose"] = dual["keep_nlrp3_pose"].fillna(dual["keep_nlrp3_pose_old"])
            dual["keep_dual_acid_geometry"] = dual["keep_urat1_acid"].fillna(False) & dual[
                "keep_nlrp3_pose"
            ].fillna(False)
            dual["keep_dual_acid_structural"] = dual["keep_urat1_acid"].fillna(False) & dual[
                "keep_nlrp3_structural"
            ].fillna(False)
            dual.to_csv(args.a2_dir / f"acid_dual_keep_structural_seed{seed}.csv", index=False)
            summary["n_keep_dual_loose"] = int(dual.keep_dual_acid_geometry.sum())
            summary["n_keep_dual_structural"] = int(dual.keep_dual_acid_structural.sum())

        (args.a2_dir / f"nlrp3_structural_summary_seed{seed}.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
