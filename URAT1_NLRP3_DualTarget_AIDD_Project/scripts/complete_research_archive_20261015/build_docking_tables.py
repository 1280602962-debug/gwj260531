#!/usr/bin/env python3
"""Build the docking protocol registry, run ledger, pose-level scores and the
ranking-method comparison for the complete research archive.

Pose-level extraction covers the gate-relevant campaigns only: the 156-molecule
acid clinical campaign (both targets, seeds 42/43/44), the stage-3B set (32 new
molecules plus the REP_07837 gap fill) and the seven reference redocks. The
9839-molecule protocol benchmark and the 1583-molecule rank track are reported
at ligand level because their deposited score tables are already one row per
molecule (the rank track was run with num_modes 1).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "scripts"))
from parse_c1_sdf_readouts import load_poses  # noqa: E402

ARCHIVE = PROJECT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
OUT = ARCHIVE / "05_DOCKING"
EXPORT = PROJECT / "docking_export_20260820"
P1 = EXPORT / "01_phase1_benchmark_URAT1_9DKB"
P2 = EXPORT / "02_phase2_dualtarget_9DKB_7ALV"
C1 = PROJECT / "data/campaigns/c1"
S3B = PROJECT / "results/nlrp3_model_v2/stage3b_local"
SEEDS = (42, 43, 44)
PROPS = ("CNNscore", "CNNaffinity", "CNN_VS", "minimizedAffinity")


def fprop(mol, key: str):
    if mol.HasProp(key):
        try:
            return float(mol.GetProp(key))
        except ValueError:
            return None
    return None


def pose_rows(sdf: Path, protocol: str, ligand: str, target: str, seed, origin: str):
    rows = []
    for i, m in enumerate(load_poses(sdf), start=1):
        row = {
            "protocol_id": protocol, "ligand_id": ligand, "target": target,
            "seed": seed, "pose_mode": i, "pose_origin": origin,
            "source_file": str(sdf.relative_to(PROJECT)),
        }
        for p in PROPS:
            row[p] = fprop(m, p)
        if row["CNN_VS"] is None and row["CNNscore"] is not None and row["CNNaffinity"] is not None:
            row["CNN_VS"] = row["CNNscore"] * row["CNNaffinity"]
        rows.append(row)
    return rows


def acid_sdf(ligand: str, target: str, seed: int) -> Path:
    folder = "acid_dual" if seed == 42 else "acid_dual_a2"
    return C1 / "07_clinical_dock" / folder / target / f"seed{seed}" / f"{ligand}_out.sdf"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------- protocol registry ----------------
    registry = [
        {
            "protocol_id": "P_BENCH_9DKB_VINA_GNINA",
            "stage": "S5_protocol_benchmark",
            "purpose": "choose a scoring channel by retrospective enrichment",
            "target": "URAT1", "pdb": "9DKB",
            "receptor_file": "Maestro-prepared 9DKB receptor as deposited in "
                             "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/receptor/",
            "ligand_source": "Schrodinger LigPrep 3D conformers converted by Meeko",
            "engine": "AutoDock Vina and GNINA",
            "engine_version": "see docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/meta/",
            "box_center": "99.980, 102.958, 105.657", "box_size": "20 x 20.01 x 20",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42",
            "scoring": "Vina affinity, CNNscore, CNNaffinity, minimizedAffinity, RTMScore rescoring",
            "hardware": "see logs", "n_molecules": 9839, "n_jobs": "9839 per engine",
            "raw_output": "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/",
            "status": "complete_formal_output",
            "comparability_note": "different receptor preparation and different ligand "
                                  "preparation from the acid campaign; scores are not poolable",
        },
        {
            "protocol_id": "P_RANK_TRACK_9DKB_7ALV",
            "stage": "S7_docking",
            "purpose": "rank the 1588 NLRP3 pool on both targets (Rank track)",
            "target": "URAT1 and NLRP3", "pdb": "9DKB and 7ALV",
            "receptor_file": "docking_export_20260820/02_phase2_dualtarget_9DKB_7ALV/receptors/",
            "ligand_source": "scripts/prepare_ligands_vina.py, no pH 7.4 microspecies enumeration",
            "engine": "GNINA", "engine_version": "see phase2 logs",
            "box_center": "see config/docking_production_p2.yaml",
            "box_size": "see config/docking_production_p2.yaml",
            "exhaustiveness": 32, "num_modes": 1, "seeds": "42",
            "scoring": "CNNaffinity reported as dock_score",
            "hardware": "GPU server", "n_molecules": 1583, "n_jobs": 3166,
            "raw_output": "docking_export_20260820/02_phase2_dualtarget_9DKB_7ALV/",
            "status": "complete_but_superseded",
            "comparability_note": "num_modes 1 means no pose ensemble; the Rank track was "
                                  "closed and replaced by the Acid track",
        },
        {
            "protocol_id": "P_ACID_DUAL_156",
            "stage": "S7_docking",
            "purpose": "formal dual-target docking of the 156 acid clinical molecules",
            "target": "URAT1 and NLRP3", "pdb": "9DKB and 7ALV",
            "receptor_file": "data/structures/prepared/9DKB_receptor.pdbqt and 7ALV_receptor.pdbqt",
            "ligand_source": "scripts/prepare_ligands_c1.py, Dimorphite-DL pH 7.4 then Meeko",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "URAT1 99.966, 102.967, 105.699 / NLRP3 16.756, 35.449, 125.714",
            "box_size": "URAT1 22 / NLRP3 20",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42, 43, 44",
            "scoring": "cnn_scoring rescore", "hardware": "CPU, --no_gpu",
            "n_molecules": 156, "n_jobs": "934 deposited SDFs",
            "raw_output": "data/campaigns/c1/07_clinical_dock/acid_dual/ and acid_dual_a2/",
            "status": "complete_formal_output_with_two_missing_nlrp3_sdfs",
            "comparability_note": "this is the production protocol; stage 3B reproduces it exactly",
        },
        {
            "protocol_id": "P_STAGE3B_NEW32",
            "stage": "S7_docking",
            "purpose": "supplementary docking of 32 new carboxylates plus one gap fill",
            "target": "URAT1 and NLRP3", "pdb": "9DKB and 7ALV",
            "receptor_file": "data/structures/prepared/9DKB_receptor.pdbqt and 7ALV_receptor.pdbqt "
                             "(byte-identical to the 156-molecule campaign)",
            "ligand_source": "scripts/prepare_ligands_c1.py, same code and parameters",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "URAT1 99.966, 102.967, 105.699 / NLRP3 16.756, 35.449, 125.714",
            "box_size": "URAT1 22 / NLRP3 20",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42, 43, 44",
            "scoring": "cnn_scoring rescore", "hardware": "CPU, --no_gpu, 2 jobs x 4 threads",
            "n_molecules": 32, "n_jobs": 194,
            "raw_output": "results/nlrp3_model_v2/stage3b_local/new32/ and gap_REP_07837/",
            "status": "complete_formal_output",
            "comparability_note": "directly comparable to P_ACID_DUAL_156",
        },
        {
            "protocol_id": "P_REFERENCE7_REPRODUCTION",
            "stage": "S6_reference_redock",
            "purpose": "strict numerical reproduction check before any new docking",
            "target": "URAT1 and NLRP3", "pdb": "9DKB and 7ALV",
            "receptor_file": "byte-identical copies of the production receptors",
            "ligand_source": "byte-identical copies of the historical ligand PDBQTs",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "as production", "box_size": "as production",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42, 43",
            "scoring": "cnn_scoring rescore", "hardware": "CPU, --no_gpu",
            "n_molecules": 4, "n_jobs": 7,
            "raw_output": "results/nlrp3_model_v2/stage3b_local/reference7/",
            "status": "complete_formal_output",
            "comparability_note": "matched-pose CNNscore difference 0 and heavy-atom RMSD 0 "
                                  "for all seven jobs on seeds 42 and 43; seed 44 not rerun",
        },
        {
            "protocol_id": "P_SELFDOCK_L2",
            "stage": "S6_reference_redock",
            "purpose": "reference-ligand self-docking and RMSD validation",
            "target": "URAT1 and NLRP3", "pdb": "9DKB, 9DKA, 7ALV",
            "receptor_file": "data/structures/prepared/",
            "ligand_source": "data/campaigns/c1/01_ligand_prep/ and c5/01_ligand_prep/",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "as production", "box_size": "as production",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42, 43, 44",
            "scoring": "cnn_scoring rescore", "hardware": "CPU",
            "n_molecules": "lesinurad, benzbromarone, NP3-146",
            "n_jobs": "see data/campaigns/c1/02_selfdock/ and c5/01_crossdock/",
            "raw_output": "data/campaigns/c1/02_selfdock/l2_selfdock_metrics.csv",
            "status": "complete_formal_output",
            "comparability_note": "defines the crystal anchor used by the A2 geometry rule",
        },
        {
            "protocol_id": "P_CAPABILITY_D1D2",
            "stage": "S7_docking",
            "purpose": "local URAT1 capability probe",
            "target": "URAT1", "pdb": "9DKB",
            "receptor_file": "data/structures/prepared/9DKB_receptor.pdbqt",
            "ligand_source": "data/campaigns/c5/06_capability_dock/",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "as production", "box_size": "as production",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "see config/docking_capability.yaml",
            "scoring": "cnn_scoring rescore", "hardware": "CPU",
            "n_molecules": "20 jobs", "n_jobs": 20,
            "raw_output": "data/campaigns/c5/06_capability_dock/",
            "status": "exploratory_does_not_change_the_frozen_lists",
            "comparability_note": "capability probe only; not part of the nomination path",
        },
        {
            "protocol_id": "P_NLRP3_W4_PANEL",
            "stage": "S8_structural_gate",
            "purpose": "NLRP3 gate validation on a literature panel",
            "target": "NLRP3", "pdb": "7ALV",
            "receptor_file": "data/structures/prepared/7ALV_receptor.pdbqt",
            "ligand_source": "data/campaigns/c5/02_nlrp3_panel/",
            "engine": "GNINA", "engine_version": "v1.3.2 master:f23dd2b",
            "box_center": "16.756, 35.449, 125.714", "box_size": "20",
            "exhaustiveness": 32, "num_modes": 9, "seeds": "42",
            "scoring": "cnn_scoring rescore", "hardware": "CPU",
            "n_molecules": 146, "n_jobs": 146,
            "raw_output": "data/campaigns/c1/05_metrics/nlrp3_structural_panel/ and c5/02_nlrp3_panel/",
            "status": "complete_formal_output",
            "comparability_note": "gate validation panel; scaffold-concentrated positives",
        },
    ]
    pd.DataFrame(registry).to_csv(OUT / "DOCKING_PROTOCOL_REGISTRY.csv", index=False)

    # ---------------- pose-level results ----------------
    acid_ids = pd.read_csv(
        C1 / "07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv"
    ).repurposing_id.tolist()
    rows: list[dict] = []
    missing: list[dict] = []
    for lig in acid_ids:
        for target in ("urat1_9dkb", "nlrp3_7alv"):
            for seed in SEEDS:
                p = acid_sdf(lig, target, seed)
                if p.exists() and p.stat().st_size > 0:
                    rows += pose_rows(p, "P_ACID_DUAL_156", lig, target, seed,
                                      "historical_deposited")
                else:
                    gap = S3B / "gap_REP_07837/rerun" / f"seed{seed}" / "poses_out.sdf"
                    if lig == "REP_07837" and target == "nlrp3_7alv" and gap.exists():
                        rows += pose_rows(gap, "P_STAGE3B_NEW32", lig, target, seed,
                                          "stage3b_gap_fill")
                    else:
                        missing.append({"protocol_id": "P_ACID_DUAL_156", "ligand_id": lig,
                                        "target": target, "seed": seed,
                                        "expected_file": str(p.relative_to(PROJECT)),
                                        "status": "sdf_absent_run_never_completed"})
    new32_ids = pd.read_csv(S3B / "analysis/new32_summary.csv").ligand_id.tolist()
    for lig in new32_ids:
        for target in ("urat1_9dkb", "nlrp3_7alv"):
            for seed in SEEDS:
                p = S3B / "new32/poses" / target / f"seed{seed}" / lig / "poses_out.sdf"
                if p.exists():
                    rows += pose_rows(p, "P_STAGE3B_NEW32", lig, target, seed, "stage3b_new")
                else:
                    missing.append({"protocol_id": "P_STAGE3B_NEW32", "ligand_id": lig,
                                    "target": target, "seed": seed,
                                    "expected_file": str(p.relative_to(PROJECT)),
                                    "status": "sdf_absent"})
    for d in sorted((S3B / "reference7").glob("0*")):
        sdf = d / "poses_out.sdf"
        if not sdf.exists():
            continue
        parts = d.name.split("_")
        lig = "_".join(parts[1:-2])
        pdb = parts[-2]
        seed = int(parts[-1].replace("seed", ""))
        target = "urat1_9dkb" if pdb == "9DKB" else "nlrp3_7alv"
        rows += pose_rows(sdf, "P_REFERENCE7_REPRODUCTION", lig, target, seed,
                          "stage3b_reference_reproduction")
    poses = pd.DataFrame(rows)
    poses.to_csv(OUT / "ALL_POSE_LEVEL_RESULTS.csv", index=False)
    pd.DataFrame(missing).to_csv(OUT / "MISSING_POSE_FILES.csv", index=False)

    # ---------------- long-form per-run summary ----------------
    g = poses.groupby(["protocol_id", "ligand_id", "target", "seed", "pose_origin"], dropna=False)
    long = g.agg(
        n_poses=("pose_mode", "count"),
        max_CNNscore=("CNNscore", "max"),
        max_CNNaffinity=("CNNaffinity", "max"),
        min_minimizedAffinity=("minimizedAffinity", "min"),
        max_CNN_VS=("CNN_VS", "max"),
    ).reset_index()
    sel = poses.loc[poses.groupby(
        ["protocol_id", "ligand_id", "target", "seed"], dropna=False
    )["CNNscore"].idxmax()]
    sel = sel[["protocol_id", "ligand_id", "target", "seed", "pose_mode",
               "CNNaffinity", "minimizedAffinity", "source_file"]].rename(columns={
        "pose_mode": "mode_of_max_CNNscore",
        "CNNaffinity": "CNNaffinity_at_max_CNNscore",
        "minimizedAffinity": "minimizedAffinity_at_max_CNNscore",
    })
    long = long.merge(sel, on=["protocol_id", "ligand_id", "target", "seed"], how="left")
    long["num_modes_requested"] = 9
    long["fewer_than_requested_is_not_a_failure"] = True
    long["run_status"] = "dock_completed"
    long.to_csv(OUT / "ALL_DOCKING_SCORES_LONG.csv", index=False)

    # ---------------- wide, ligand level, all protocols ----------------
    wide = long.pivot_table(
        index="ligand_id",
        columns=["protocol_id", "target", "seed"],
        values=["max_CNNscore", "CNNaffinity_at_max_CNNscore", "max_CNNaffinity"],
        aggfunc="first",
    )
    wide.columns = ["__".join(str(c) for c in col) for col in wide.columns]
    wide = wide.reset_index()
    bench = pd.read_csv(P1 / "scores/mol_protocol_scores.csv").rename(columns={"mol_id": "ligand_id"})
    bench = bench.assign(protocol_id="P_BENCH_9DKB_VINA_GNINA")
    bench.to_csv(OUT / "BENCH_9DKB_LIGAND_SCORES.csv", index=False)
    du = pd.read_csv(P2 / "docking_9dkb_URAT1/docking_9dkb_gnina.csv")[
        ["repurposing_id", "dock_score", "docking_status"]
    ].rename(columns={"repurposing_id": "ligand_id",
                      "dock_score": "P_RANK_TRACK__urat1_9dkb__dock_score",
                      "docking_status": "P_RANK_TRACK__urat1_9dkb__status"})
    dn = pd.read_csv(P2 / "docking_7alv_NLRP3/docking_7alv_gnina.csv")[
        ["repurposing_id", "dock_score", "docking_status"]
    ].rename(columns={"repurposing_id": "ligand_id",
                      "dock_score": "P_RANK_TRACK__nlrp3_7alv__dock_score",
                      "docking_status": "P_RANK_TRACK__nlrp3_7alv__status"})
    wide = wide.merge(du, on="ligand_id", how="outer").merge(dn, on="ligand_id", how="outer")
    wide.to_csv(OUT / "ALL_DOCKING_SCORES_WIDE.csv", index=False)

    # ---------------- ranking method comparison ----------------
    pm = pd.read_csv(P1 / "metrics/protocol_metrics.csv")
    pm["population"] = "URAT1 9DKB benchmark, 9839 molecules, 469 actives"
    pm["used_in_production"] = pm.ID.map({
        "P2": "yes, CNNaffinity was the rank-track production channel",
        "P0": "yes, CNNscore is the pose-selection channel of the acid campaign",
        "P5": "no, RTMScore rescoring was not carried into production",
        "P4": "no, incomplete coverage 8983/9839",
        "P1": "no, Vina affinity rejected on enrichment",
        "P3": "no, minimizedAffinity rejected on enrichment",
    })
    pm["rejection_reason"] = pm.ID.map({
        "P2": "", "P0": "",
        "P5": "best enrichment but adds an external rescoring dependency and zero "
              "random-decoy EF@1%; not adopted",
        "P4": "coverage incomplete and EF@1% 0.647",
        "P1": "EF@1% 0.431, close to random",
        "P3": "EF@1% 0.423, close to random",
    })
    pm.to_csv(OUT / "RANKING_METHOD_COMPARISON.csv", index=False)
    # Confidence intervals live in a separate deposited file with two benchmark
    # rows per protocol; copying it in as-is avoids duplicating the headline table.
    ci = PROJECT / "data/si/protocol_enrichment_ci/protocol_ef_ci.csv"
    if ci.exists():
        pd.read_csv(ci).to_csv(OUT / "RANKING_METHOD_ENRICHMENT_CI.csv", index=False)

    # ---------------- reference redocking ----------------
    ref_rows = []
    u = pd.read_csv(PROJECT / "data/manuscript/tables/urat1_redocking_summary.csv")
    n = pd.read_csv(PROJECT / "data/manuscript/tables/nlrp3_selfdocking_summary.csv")
    for _, r in pd.concat([u, n], ignore_index=True).iterrows():
        ref_rows.append({
            "reference_ligand": r.ligand, "target": r.target, "seed": r.seed,
            "protocol_id": "P_SELFDOCK_L2",
            "top1_pose_rmsd_A": r.top1_pose_rmsd_A,
            "best_of_9_pose_rmsd_A": r.best_of_9_pose_rmsd_A,
            "search_found_near_native_le_2A": bool(r.best_of_9_pose_rmsd_A <= 2.0),
            "scoring_ranked_near_native_first": bool(r.top1_le_2A),
            "rmsd_definition": r.rmsd_definition,
            "raw_output_exists": True,
            "source_file": r.source_file,
        })
    rc = pd.read_csv(S3B / "reference7/task_comparison.csv")
    rc.to_csv(OUT / "REFERENCE7_STRICT_REPRODUCTION.csv", index=False)
    pd.DataFrame(ref_rows).to_csv(OUT / "REFERENCE_REDOCKING.csv", index=False)

    print("pose rows:", len(poses))
    print("missing pose files:", len(missing))
    print("long rows:", len(long), "wide rows:", len(wide))
    print(poses.groupby("protocol_id").size().to_string())


if __name__ == "__main__":
    main()
