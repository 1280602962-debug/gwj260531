#!/usr/bin/env python3
"""Emit the remaining archive tables from deposited files and verified constants.

Nothing here recomputes a docking score, a model probability, or an MD metric.
Derived ranks and joins are labelled as such.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
A = PROJECT / "COMPLETE_RESEARCH_ARCHIVE_20261015"


def write(rel: str, frame: pd.DataFrame) -> None:
    path = A / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    print(rel, len(frame))


def main() -> None:
    # ----- model performance, one row per reported evaluation -----
    write("03_MACHINE_LEARNING/MODEL_PERFORMANCE.csv", pd.DataFrame([
        {
            "model": "URAT1_XGBoost_regressor_v1",
            "artifact": "data/models/urat1_model.joblib",
            "metrics_file": "data/models/training_report.json",
            "evaluation": "pooled_out_of_fold_5fold_murcko_scaffold",
            "n": 822, "auroc": "", "auprc": "",
            "r2": 0.5079205459521259, "spearman": 0.7258124808529701,
            "rmse": 0.6634725892592518, "mae": 0.49847748712321566,
            "ef_10pct_pge6": 1.70991731239274,
            "used_for_clinical_ranking": False,
            "note": "SI and benchmark only. Script header forbids clinical ranking.",
        },
        {
            "model": "NLRP3_XGBoost_classifier_v1",
            "artifact": "data/models/nlrp3_model.joblib",
            "metrics_file": "data/models/training_report.json",
            "evaluation": "mean_of_five_scaffold_fold_metrics_not_pooled_oof",
            "n": 513, "auroc": 0.8934017784193573, "auprc": 0.9143469018214805,
            "r2": "", "spearman": "", "rmse": "", "mae": "",
            "ef_10pct_pge6": 1.5706359287488578,
            "used_for_clinical_ranking": True,
            "note": "This is the 0.8934 number in MANUSCRIPT_FACTS.yaml. It is a mean of fold AUROCs.",
        },
        {
            "model": "NLRP3_v2_structure_only_calibration_none",
            "artifact": "results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib",
            "metrics_file": "results/nlrp3_model_v2/performance_summary.csv",
            "evaluation": "pooled_outer_test_nested_scaffold_cv",
            "n": 513, "auroc": 0.8895916097250914, "auprc": 0.9131075396059838,
            "r2": "", "spearman": "", "rmse": "", "mae": "",
            "ef_10pct_pge6": 1.5575,
            "used_for_clinical_ranking": False,
            "note": "Not the historical 1588-pool model. Recall 0.8742, specificity 0.8128, Brier 0.1261.",
        },
        {
            "model": "NLRP3_v2_n513_production_rescore",
            "artifact": "results/nlrp3_model_v2/nlrp3_v2_structure_only_n513_production.joblib",
            "metrics_file": "results/nlrp3_model_v2/stage2/library_summary.json",
            "evaluation": "no_held_out_metric_full_513_fit",
            "n": 513, "auroc": "", "auprc": "",
            "r2": "", "spearman": "", "rmse": "", "mae": "",
            "ef_10pct_pge6": "",
            "used_for_clinical_ranking": False,
            "note": "Used only to rescore the library in stage 2. Did not replace the frozen 1588 pool.",
        },
    ]))

    write("04_PROTEIN_LIGAND_PREP/RECEPTOR_PREPARATION.csv", pd.DataFrame([
        {"pdb": "9DKB", "target": "URAT1", "uniprot": "Q96S37", "species": "Homo sapiens",
         "state": "inward-open, ligand-bound cryo-EM", "resolution_A": 2.74,
         "crystal_ligand": "A1AIL lesinurad", "chain": "A",
         "preparation": "gemmi protein-only, drop waters and HETATM, Open Babel AddHydrogens pH 7.4, rigid PDBQT",
         "file": "data/structures/prepared/9DKB_receptor.pdbqt", "bytes": 620900, "atom_lines": 7759,
         "partial_charges": "all 0.000 in the deposited PDBQT",
         "his_tautomer": "not recorded", "role": "production",
         "residue_note": "literature Arg477 is prepared ARG A 476"},
        {"pdb": "9DKB", "target": "URAT1", "uniprot": "Q96S37", "species": "Homo sapiens",
         "state": "same structure, Maestro-written receptor", "resolution_A": 2.74,
         "crystal_ligand": "A1AIL lesinurad", "chain": "A",
         "preparation": "Maestro; exact ProtAssign/Epik command cannot be recovered",
         "file": "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/receptor/9DKB_receptor.pdbqt",
         "bytes": 375040, "atom_lines": 4688,
         "partial_charges": "nonzero Maestro charges",
         "his_tautomer": "not recorded", "role": "phase1_benchmark_only",
         "residue_note": "do not pool scores with the production receptor"},
        {"pdb": "7ALV", "target": "NLRP3", "uniprot": "Q96P20", "species": "Homo sapiens",
         "state": "NACHT, X-ray, ADP and inhibitor bound", "resolution_A": 2.84,
         "crystal_ligand": "RM5 NP3-146", "chain": "source chain A; prepared file is resegmented A-I",
         "preparation": "same gemmi plus Open Babel pH 7.4 path as 9DKB",
         "file": "data/structures/prepared/7ALV_receptor.pdbqt", "bytes": 1480113, "atom_lines": 4519,
         "partial_charges": "all 0.000", "his_tautomer": "not recorded", "role": "production",
         "residue_note": "use residue_mapping.csv, not native numbering"},
        {"pdb": "8ETR", "target": "NLRP3", "uniprot": "Q96P20", "species": "Homo sapiens",
         "state": "ligand-bound cryo-EM, GDC-2394", "resolution_A": 3.5,
         "crystal_ligand": "WTN GDC-2394", "chain": "not prepared",
         "preparation": "no prepared PDBQT in the tree",
         "file": "data/structures/pdb/8ETR.cif", "bytes": "", "atom_lines": "",
         "partial_charges": "", "his_tautomer": "", "role": "rigid_transfer_only_not_production",
         "residue_note": ""},
        {"pdb": "9DKA", "target": "URAT1", "uniprot": "Q96S37", "species": "Homo sapiens",
         "state": "benzbromarone-bound", "resolution_A": 3.0,
         "crystal_ligand": "R75 benzbromarone", "chain": "A",
         "preparation": "same Open Babel path", "file": "data/structures/prepared/9DKA_receptor.pdbqt",
         "bytes": 624500, "atom_lines": "", "partial_charges": "0.000",
         "his_tautomer": "not recorded", "role": "crossdock_and_capability_not_production",
         "residue_note": ""},
    ]))

    write("04_PROTEIN_LIGAND_PREP/LIGAND_PREPARATION.csv", pd.DataFrame([
        {"pipeline_id": "L1_dimorphite_c1", "script": "scripts/prepare_ligands_c1.py",
         "protonation": "Dimorphite-DL pH 7.4, precision 1.0, max 8 variants",
         "microspecies_rule": "prefer [CX3](=O)[O-], else [CX3](=O)[OH], else first variant",
         "stereo": "kept as in SMILES; unspecified centres not assigned",
         "embed": "ETKDGv3 seed 0xC0FFEE", "minimize": "MMFF maxIters 200, failure swallowed",
         "writer": "Meeko MoleculePreparation defaults plus PDBQTWriterLegacy",
         "consumed_by": "C1 acid 156, stage 3B new 32, reference ligands",
         "n": "156 and 32", "n_failed": 0,
         "comparability": "same code and same selection rule; package versions at the historical run were not recorded"},
        {"pipeline_id": "L2_vina_no_dimorphite", "script": "scripts/prepare_ligands_vina.py",
         "protonation": "none; RDKit AddHs on the input SMILES ionization",
         "microspecies_rule": "none", "stereo": "as in SMILES",
         "embed": "ETKDGv3 seed 0xC0FFEE", "minimize": "MMFF maxIters 200",
         "writer": "Meeko defaults",
         "consumed_by": "phase-2 rank track, 1588 to 1583",
         "n": 1583, "n_failed": 5,
         "comparability": "not comparable to L1; the rank track was later closed"},
        {"pipeline_id": "L3_ligprep_phase1", "script": "Schrodinger structconvert then split_ligands_meeko.py",
         "protonation": "LigPrep; pH and tautomer settings cannot be recovered",
         "microspecies_rule": "cannot be recovered", "stereo": "cannot be recovered",
         "embed": "LigPrep conformer kept, no re-embed", "minimize": "cannot be recovered",
         "writer": "Meeko",
         "consumed_by": "phase-1 9839-molecule benchmark, paired with the Maestro 9DKB receptor",
         "n": 9839, "n_failed": "10 dropped by LigPrep; mol_09300 Meeko charge fail",
         "comparability": "not comparable to L1 or L2"},
    ]))

    # stereo and protonation problem list from deposited reports
    prep = pd.read_csv(PROJECT / "results/nlrp3_model_v2/stage3b_local/new32/ligand_prep/ligand_prep_report.csv")
    rows = []
    for _, r in prep.iterrows():
        issues = []
        if r.selected_n_unassigned_stereo and r.selected_n_unassigned_stereo > 0:
            issues.append(f"unassigned_stereo_{int(r.selected_n_unassigned_stereo)}")
        if r.n_microspecies and r.n_microspecies > 1:
            issues.append(f"n_microspecies_{int(r.n_microspecies)}")
        if r.formal_charge == 0:
            issues.append("net_charge_0_zwitterion_carboxylate_still_deprotonated")
        if "HEMISUCCINATE" in str(r["name"]):
            issues.append("name_implies_salt_smiles_is_single_component")
        if issues:
            rows.append({
                "repurposing_id": r.repurposing_id, "name": r["name"], "cohort": "stage3b_new32",
                "formal_charge": r.formal_charge, "n_microspecies": r.n_microspecies,
                "unassigned_stereo": r.selected_n_unassigned_stereo,
                "issues": ";".join(issues),
                "source": "results/nlrp3_model_v2/stage3b_local/new32/ligand_prep/ligand_prep_report.csv",
                "action_taken": "recorded_not_reselected",
            })
    man = pd.read_csv(PROJECT / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv")
    # unassigned stereo is not a column there; flag multi-microspecies only
    if "n_microspecies" in man.columns:
        for _, r in man[man.n_microspecies > 1].iterrows():
            rows.append({
                "repurposing_id": r.repurposing_id, "name": r.get("name", ""),
                "cohort": "acid_clinical_156", "formal_charge": "",
                "n_microspecies": r.n_microspecies, "unassigned_stereo": "",
                "issues": f"n_microspecies_{int(r.n_microspecies)}",
                "source": "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv",
                "action_taken": "carboxylate_preference_applied_other_states_not_docked",
            })
    write("04_PROTEIN_LIGAND_PREP/PROTONATION_AND_STEREO_AUDIT.csv", pd.DataFrame(rows))

    write("01_RESEARCH_DESIGN/SCREENING_DECISION_HISTORY.csv", pd.read_csv(
        PROJECT / "data/manuscript/tables/screening_funnel.csv"
    ).assign(record_type="historical_funnel_as_frozen", retrospective_revision="none"))

    # ADMET: population summary plus candidate-name slice. SA was never computed.
    admet = pd.read_csv(PROJECT / "data/repurposing/p2/admet_pool.csv")
    summary = pd.DataFrame([{
        "population": "rank_track_dual_dock_pool",
        "n": len(admet),
        "source": "data/repurposing/p2/admet_pool.csv",
        "tool": "RDKit descriptors in scripts/10_admet_druglikeness.py",
        "sa_score": "not_computed",
        "commercial_admet": "not_run",
        "full_library_8319": "not_run",
        "mw_median": admet.mw.median(), "clogp_median": admet.clogp.median(),
        "tpsa_median": admet.tpsa.median(), "qed_median": admet.qed.median(),
        "lipinski_pass_n": int(admet.lipinski_pass.sum()),
        "veber_pass_n": int(admet.veber_pass.sum()),
        "note": "This is the superseded rank-track pool, not the acid-campaign 156.",
    }])
    write("08_ADMET_AND_MD/ADMET_SA_RESULTS.csv", summary)

    write("08_ADMET_AND_MD/MD_SOURCE_INVENTORY.csv", pd.DataFrame([
        {"item": "git_protocol_lock", "in_github": True,
         "path": "config/md_protocol.yaml",
         "status": "planned_not_authorized",
         "what_it_contains": "3 replicas x 200 ns, ff19SB, TIP3P, results_available false"},
        {"item": "git_manuscript_facts", "in_github": True,
         "path": "MANUSCRIPT_FACTS.yaml",
         "status": "declares_md_not_run",
         "what_it_contains": "md_authorized false, md_results_available false"},
        {"item": "local_c5_w5_analysis", "in_github": False,
         "path": "/mnt/d/CADD paper exercise/NLRP3_URAT1/viz_pack_20260907/viz_pack/",
         "status": "completed_replica_1_only_local",
         "what_it_contains": "4 systems x 200 ns analysis, stripped trajectories, MM-GBSA last 50 ns"},
        {"item": "full_solvated_trajectories", "in_github": False,
         "path": "/home/hww/gwj/NLRP3_URAT1/md_c5_w5 (symlink target absent on this machine)",
         "status": "not_accessible",
         "what_it_contains": "original solvated trajectories and prmtop are missing here"},
        {"item": "replicas_2_and_3", "in_github": False, "path": "",
         "status": "not_run", "what_it_contains": "readme says replicas are not submitted"},
    ]))

    write("10_SCIENTIFIC_AUDIT/CLAIM_EVIDENCE_MATRIX.csv", pd.DataFrame([
        {"claim": "NLRP3 classifier mean-of-fold AUROC is 0.8934",
         "kind": "machine_learning_metric", "supported": "yes_for_v1_mean_of_folds",
         "not_the_same_as": "v2 pooled outer-test AUROC 0.8896, and not a probability of true activity",
         "evidence": "data/models/training_report.json"},
        {"claim": "1588 clinical compounds passed p>=0.5",
         "kind": "deposited_membership", "supported": "yes_the_table_exists",
         "not_the_same_as": "a regenerable score; SI retrain gave 1377",
         "evidence": "data/repurposing/screening/docking_pool_p05.csv"},
        {"claim": "URAT1 regression ranked the clinical library",
         "kind": "process_claim", "supported": "no",
         "not_the_same_as": "",
         "evidence": "scripts/02_train_asymmetric_models.py header; no URAT1 column in clinical scores"},
        {"claim": "Production docking retrieves URAT1 actives",
         "kind": "docking_claim", "supported": "no_rank_track_was_closed",
         "not_the_same_as": "P2 EF@1% 2.54 on a decoy benchmark is not clinical retrieval",
         "evidence": "config/campaign_c1.yaml rank_track_status closed"},
        {"claim": "Lesinurad production self-dock is native-like at top-1",
         "kind": "reference_redock", "supported": "no",
         "not_the_same_as": "seed 42 best-of-9 RMSD 1.18 A shows search can succeed while scoring selects 4.30 A",
         "evidence": "data/manuscript/tables/urat1_redocking_summary.csv"},
        {"claim": "NP3-146 self-dock is native-like at top-1",
         "kind": "reference_redock", "supported": "yes_all_three_seeds",
         "not_the_same_as": "",
         "evidence": "data/manuscript/tables/nlrp3_selfdocking_summary.csv"},
        {"claim": "A2 geometry gate enriches URAT1 actives",
         "kind": "gate_validation", "supported": "no_lr_plus_about_1",
         "not_the_same_as": "A1 LR+ 2.20 is weak-to-moderate and was not the final admissibility rule",
         "evidence": "data/manuscript/tables/urat1_gate_validation.csv"},
        {"claim": "Frozen 12 and 21 still pass the unified structural rule",
         "kind": "retrospective_rescore", "supported": "yes_33_of_33",
         "not_the_same_as": "a new nomination; membership was not changed",
         "evidence": "results/nlrp3_model_v2/stage3b_local/analysis/counts.json"},
        {"claim": "13 new molecules are additional structural-and-chemistry passes",
         "kind": "stage3b_result", "supported": "yes_under_the_unified_rule",
         "not_the_same_as": "membership of the frozen 12/21, and not experimental activity",
         "evidence": "results/nlrp3_model_v2/stage3b_local/analysis/new32_summary.csv"},
        {"claim": "Any nominated molecule inhibits URAT1 and NLRP3",
         "kind": "biological_activity", "supported": "no_direct_evidence",
         "not_the_same_as": "pocket compatibility under a declared geometry rule",
         "evidence": "no wet-lab table exists"},
        {"claim": "PF-03882845 MD proves a stable dual-target complex",
         "kind": "md_observation", "supported": "not_as_stated",
         "not_the_same_as": "local replica-1 analysis exists but is not in git; URAT1 control fails the pre-registered salt-bridge occupancy rule",
         "evidence": "config/md_protocol.yaml and local viz_pack summary.json"},
    ]))

    print("done")


if __name__ == "__main__":
    main()
