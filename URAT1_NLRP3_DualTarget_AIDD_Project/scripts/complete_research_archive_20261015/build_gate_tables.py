#!/usr/bin/env python3
"""Build the structural-gate tables for the complete research archive.

Two distinct things are kept apart on purpose:
  * the HISTORICAL decision, read from the deposited per-seed gate tables that
    were written when the campaign ran;
  * the RETROSPECTIVE re-evaluation, read from the stage-3B analysis that
    re-applied the same A2 and NLRP3 structural functions to the same poses.
They are reported side by side, never merged into one verdict column.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
ARCHIVE = PROJECT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
OUT = ARCHIVE / "06_STAGE3"
A2 = PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2"
S3B = PROJECT / "results/nlrp3_model_v2/stage3b_local"
SEEDS = (42, 43, 44)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------- historical per-seed gate decisions ----------------
    hist = []
    for seed in SEEDS:
        urat1 = pd.read_csv(
            A2 / (f"acid_pose_metrics_a2_seed{seed}.csv" if seed == 42
                  else f"acid_pose_metrics_seed{seed}.csv")
        )
        urat1 = urat1[urat1.target == "urat1_9dkb"]
        nlrp3 = pd.read_csv(A2 / f"nlrp3_structural_metrics_seed{seed}.csv")
        dual = pd.read_csv(A2 / f"acid_dual_keep_structural_seed{seed}.csv")
        m = urat1[[
            "ligand_id", "pose_selection_rule", "pose_selection_status", "n_poses",
            "n_geometry_pass_modes", "selected_mode", "CNNscore", "CNNaffinity",
            "acid_arg477_min_A", "has_carboxylate_oxy_in_pose",
            "centroid_to_crystal_lig_A", "keep_urat1_acid",
        ]].rename(columns={
            "pose_selection_status": "urat1_pose_selection_status",
            "n_poses": "urat1_n_poses",
            "n_geometry_pass_modes": "urat1_n_geometry_pass_modes",
            "selected_mode": "urat1_selected_mode",
            "CNNscore": "urat1_CNNscore", "CNNaffinity": "urat1_CNNaffinity",
            "acid_arg477_min_A": "urat1_acid_arg477_min_A",
            "centroid_to_crystal_lig_A": "urat1_centroid_A",
            "keep_urat1_acid": "urat1_A2_pass",
        })
        n = nlrp3[[
            "ligand_id", "n_poses", "selected_mode", "CNNscore", "CNNaffinity",
            "centroid_to_crystal_lig_A", "pocket_overlap_frac", "n_key_contacts",
            "key_recovery_frac", "ifp_jaccard_vs_np3146", "ifp_bits",
            "keep_nlrp3_pose", "keep_nlrp3_structural",
        ]].rename(columns={
            "n_poses": "nlrp3_n_poses", "selected_mode": "nlrp3_selected_mode",
            "CNNscore": "nlrp3_CNNscore", "CNNaffinity": "nlrp3_CNNaffinity",
            "centroid_to_crystal_lig_A": "nlrp3_centroid_A",
            "pocket_overlap_frac": "nlrp3_pocket_overlap_frac",
            "ifp_jaccard_vs_np3146": "nlrp3_ifp_jaccard",
            "keep_nlrp3_pose": "nlrp3_loose_pass",
            "keep_nlrp3_structural": "nlrp3_structural_pass",
        })
        n["nlrp3_clash_any"] = n.ifp_bits.fillna("").str.contains("CLASH:ANY")
        d = dual[["ligand_id", "keep_dual_acid_geometry", "keep_dual_acid_structural"]].rename(
            columns={"keep_dual_acid_geometry": "dual_loose_pass",
                     "keep_dual_acid_structural": "dual_structural_pass"}
        )
        row = m.merge(n, on="ligand_id", how="outer").merge(d, on="ligand_id", how="outer")
        row["seed"] = seed
        row["record_type"] = "historical_decision_as_run"
        hist.append(row)
    histdf = pd.concat(hist, ignore_index=True)

    failures = [json.loads(line) for line in (A2 / "dock_failures.jsonl").read_text().splitlines() if line.strip()]
    fdf = pd.DataFrame(failures)
    fdf["historical_run_status"] = "timeout"
    fdf.to_csv(OUT / "HISTORICAL_DOCK_FAILURES.csv", index=False)
    key = set(zip(fdf.ligand_id, fdf.target, fdf.seed))
    histdf["historical_run_status"] = [
        "timeout_1800s_no_pose_file" if (l, "nlrp3_7alv", s) in key else "dock_completed"
        for l, s in zip(histdf.ligand_id, histdf.seed)
    ]
    histdf.to_csv(OUT / "ALL_STRUCTURAL_GATE_RESULTS.csv", index=False)

    # ---------------- retrospective unified re-evaluation ----------------
    retro = pd.concat([
        pd.read_csv(S3B / "analysis/historical_per_seed_rescored.csv").assign(
            cohort="historical_deposited_156"),
        pd.read_csv(S3B / "analysis/new32_per_seed.csv").assign(cohort="stage3b_new32"),
    ], ignore_index=True)
    retro["record_type"] = "retrospective_unified_rule_20261015"
    retro.to_csv(OUT / "UNIFIED_CANDIDATE_EVALUATION.csv", index=False)

    # ---------------- historical vs retrospective comparison ----------------
    h = histdf[["ligand_id", "seed", "urat1_A2_pass", "nlrp3_structural_pass",
                "dual_structural_pass", "historical_run_status"]].rename(columns={
        "urat1_A2_pass": "hist_urat1_pass",
        "nlrp3_structural_pass": "hist_nlrp3_pass",
        "dual_structural_pass": "hist_dual_pass",
    })
    r = retro[retro.cohort == "historical_deposited_156"][[
        "ligand_id", "seed", "urat1_pass", "nlrp3_structural_pass", "dual_seed_pass",
        "nlrp3_pose_origin",
    ]].rename(columns={
        "urat1_pass": "retro_urat1_pass",
        "nlrp3_structural_pass": "retro_nlrp3_pass",
        "dual_seed_pass": "retro_dual_pass",
    })
    cmp = h.merge(r, on=["ligand_id", "seed"], how="outer")
    for a, b, name in [("hist_urat1_pass", "retro_urat1_pass", "urat1"),
                       ("hist_nlrp3_pass", "retro_nlrp3_pass", "nlrp3"),
                       ("hist_dual_pass", "retro_dual_pass", "dual")]:
        cmp[f"{name}_agrees"] = cmp[a].astype("object") == cmp[b].astype("object")
    cmp.to_csv(OUT / "HISTORICAL_VS_RETROSPECTIVE_GATE.csv", index=False)

    # ---------------- stage 3A scope ----------------
    s3a = pd.read_csv(PROJECT / "results/nlrp3_model_v2/stage3a/new_molecule_docking_decisions.csv")
    s3a.to_csv(OUT / "STAGE3A_DOCKING_DECISIONS.csv", index=False)
    pd.read_csv(PROJECT / "results/nlrp3_model_v2/stage3a/chemistry_pass_157_annotation.csv").to_csv(
        OUT / "STAGE3A_CHEMISTRY_PASS_157.csv", index=False
    )

    print("historical gate rows:", len(histdf), "retrospective rows:", len(retro))
    print("historical timeouts:", len(fdf))
    for name in ("urat1", "nlrp3", "dual"):
        col = f"{name}_agrees"
        print(f"{name} agreement:", int(cmp[col].sum()), "of", len(cmp))
    print("\nhistorical dual_structural_pass per seed:")
    print(histdf.groupby("seed").dual_structural_pass.sum().to_string())


if __name__ == "__main__":
    main()
