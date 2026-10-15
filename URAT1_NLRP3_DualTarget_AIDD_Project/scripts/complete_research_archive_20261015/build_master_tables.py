#!/usr/bin/env python3
"""Build the compound-level master index, model-score table, ranking
transitions and screening-decision ledger for the complete research archive.

Every column is either copied from an existing deposited file or derived by a
rule stated in the `derivation` column of SOURCE_PROVENANCE.csv. No score is
recomputed from structures and no missing value is imputed.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
ARCHIVE = PROJECT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
EXPORT = PROJECT / "docking_export_20260820"
P2 = EXPORT / "02_phase2_dualtarget_9DKB_7ALV"
C1 = PROJECT / "data/campaigns/c1"
C5 = PROJECT / "data/campaigns/c5"
V2 = PROJECT / "results/nlrp3_model_v2"
S3B = V2 / "stage3b_local"


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def rank_desc(series: pd.Series) -> pd.Series:
    """Dense competition rank, 1 = best, higher score is better, NaN stays NaN."""
    return series.rank(ascending=False, method="min", na_option="keep")


def main() -> None:
    (ARCHIVE / "03_MACHINE_LEARNING").mkdir(parents=True, exist_ok=True)
    (ARCHIVE / "07_RANKING_AND_SELECTION").mkdir(parents=True, exist_ok=True)

    ml = read(P2 / "screening_pool/nlrp3_ml_scores_clinical_all.csv")
    cl = read(V2 / "stage2/clinical_library_scores.csv")
    pool = read(P2 / "screening_pool/docking_pool_p05.csv")
    acid303 = read(C1 / "07_clinical_dock/acid_pool/acid_clinical_pool.csv")
    acid156 = read(C1 / "07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv")
    dock_u = read(P2 / "docking_9dkb_URAT1/docking_9dkb_gnina.csv")
    dock_n = read(P2 / "docking_7alv_NLRP3/docking_7alv_gnina.csv")
    pareto = read(P2 / "pareto/pareto_merged_scores.csv")
    nominated = read(P2 / "candidates/nominated_candidates.csv")
    shortlist7 = read(P2 / "candidates/nominated_shortlist_diverse.csv")
    ledger = read(PROJECT / "data/manuscript/tables/screening_decision_ledger.csv")
    primary = read(C5 / "04_shortlist_frozen/primary_candidates.csv")
    backup = read(C5 / "04_shortlist_frozen/backup_candidates.csv")
    reportable = read(C5 / "04_shortlist_frozen/reportable_candidates.csv")
    descr = read(PROJECT / "data/manuscript/tables/acid_descriptors.csv")
    stage3a = read(V2 / "stage3a/new_molecule_docking_decisions.csv")
    anno157 = read(V2 / "stage3a/chemistry_pass_157_annotation.csv")
    new32 = read(S3B / "analysis/new32_summary.csv")
    hist_resc = read(S3B / "analysis/historical_summary_rescored.csv")

    # ---------------- model scores and ranks ----------------
    m = ml[[
        "repurposing_id", "chembl_id", "name", "canonical_smiles", "inchi_key",
        "mw", "max_phase", "scaffold", "p_active_nlrp3", "nlrp3_percentile",
        "nlrp3_pred_active", "max_tanimoto_nlrp3_train", "rank_nlrp3",
    ]].rename(columns={
        "p_active_nlrp3": "nlrp3_p_active_historical_v1",
        "nlrp3_percentile": "nlrp3_percentile_historical_v1",
        "rank_nlrp3": "nlrp3_rank_historical_v1_deposited",
        "max_tanimoto_nlrp3_train": "max_tanimoto_nlrp3_train_historical",
    })
    v2cols = cl[[
        "repurposing_id", "p_historical", "p_n392", "p_n513",
        "in_historical_pool", "in_v2_pool", "in_n392_pool",
        "max_tanimoto_train", "ad_threshold", "in_applicability_domain",
        "in_training_set",
    ]].rename(columns={
        "p_historical": "nlrp3_p_rescored_historical_model",
        "p_n392": "nlrp3_p_v2_n392",
        "p_n513": "nlrp3_p_v2_n513_production",
        "max_tanimoto_train": "max_tanimoto_train_v2",
    })
    model = m.merge(v2cols, on="repurposing_id", how="left", validate="one_to_one")
    model["nlrp3_rank_historical_v1_recomputed"] = rank_desc(model["nlrp3_p_active_historical_v1"]).astype("Int64")
    model["nlrp3_rank_v2_n513"] = rank_desc(model["nlrp3_p_v2_n513_production"]).astype("Int64")
    model["nlrp3_threshold_0p5_pass_historical_v1"] = model["nlrp3_p_active_historical_v1"] >= 0.5
    model["nlrp3_threshold_0p5_pass_v2_n513"] = model["nlrp3_p_v2_n513_production"] >= 0.5
    model["urat1_regression_prediction"] = pd.NA
    model["urat1_regression_rank"] = pd.NA
    model["urat1_model_role"] = "not_applied_to_clinical_library"
    model["model_version_used_for_production_pool"] = "historical_v1_classifier_commit_72b13228"
    model["ranking_population"] = "clinical_library_8319"
    model["ranking_direction"] = "higher_probability_is_better"
    model["tie_handling"] = "min_rank_competition"
    model["missing_value_handling"] = "rank_left_empty_when_score_missing"
    model.to_csv(ARCHIVE / "03_MACHINE_LEARNING/ALL_MODEL_SCORES_AND_RANKS.csv", index=False)

    # ---------------- master compound index ----------------
    idx = model[[
        "repurposing_id", "name", "chembl_id", "canonical_smiles", "inchi_key",
        "max_phase", "mw", "scaffold",
        "nlrp3_p_active_historical_v1", "nlrp3_percentile_historical_v1",
        "nlrp3_rank_historical_v1_recomputed", "nlrp3_p_v2_n513_production",
        "nlrp3_rank_v2_n513", "in_applicability_domain", "in_training_set",
    ]].copy()

    chem = cl[[
        "repurposing_id", "acid_equivalent", "acid_reject_reason", "has_carboxylate",
        "has_tetrazole", "has_acyl_sulfonamide", "pass_veber", "pass_ro5_hb",
        "pass_mw_200_550", "pass_chemistry_soft", "in_deposited_a2_docking_set",
        "docking_complete_a2_three_seed", "historical_name_exclusion_vecabrutinib",
    ]]
    idx = idx.merge(chem, on="repurposing_id", how="left")

    idx["in_clinical_library_8319"] = True
    idx["in_nlrp3_pool_1588"] = idx.repurposing_id.isin(pool.repurposing_id)
    idx["in_acid_pool_303"] = idx.repurposing_id.isin(acid303.repurposing_id)
    idx["in_chemistry_soft_156"] = idx.repurposing_id.isin(acid156.repurposing_id)
    idx["in_rank_track_dock_1583"] = idx.repurposing_id.isin(dock_u.repurposing_id)
    idx["in_rank_track_pareto_1580"] = idx.repurposing_id.isin(pareto.repurposing_id)
    nom_ids = set(nominated.chembl_id.dropna())
    idx["in_rank_track_dual_gate_51"] = idx.chembl_id.isin(nom_ids)
    idx["in_rank_track_diverse_shortlist_7"] = idx.chembl_id.isin(set(shortlist7.chembl_id.dropna()))
    idx["in_frozen_primary_12"] = idx.repurposing_id.isin(primary.ligand_id)
    idx["in_frozen_backup_21"] = idx.repurposing_id.isin(backup.ligand_id)
    idx["in_frozen_reportable_33"] = idx.repurposing_id.isin(reportable.ligand_id)
    idx["in_chemistry_pass_panel_157"] = idx.repurposing_id.isin(anno157.repurposing_id)
    idx["in_stage3a_new_39"] = idx.repurposing_id.isin(stage3a.repurposing_id)
    idx["in_stage3b_new32_docked"] = idx.repurposing_id.isin(new32.ligand_id)

    # descriptors recomputed for the 303 acid pool only, by the manuscript pipeline
    d = descr[[
        "ligand_id", "mw_recomputed", "clogp", "tpsa", "hbd", "hba",
        "rotatable_bonds", "qed_recomputed", "murcko_scaffold", "acid_types",
        "rdkit_version",
    ]].rename(columns={"ligand_id": "repurposing_id"})
    idx = idx.merge(d, on="repurposing_id", how="left")
    idx["descriptor_scope"] = idx.mw_recomputed.notna().map(
        {True: "recomputed_for_acid_pool_303", False: "not_computed"}
    )

    # ligand preparation formal charge, stage 3B only
    prep = read(S3B / "new32/ligand_prep/ligand_prep_report.csv")[[
        "repurposing_id", "formal_charge", "selected_n_unassigned_stereo", "n_microspecies",
    ]].rename(columns={
        "formal_charge": "prep_formal_charge_stage3b",
        "selected_n_unassigned_stereo": "prep_unassigned_stereo_stage3b",
        "n_microspecies": "prep_n_microspecies_stage3b",
    })
    idx = idx.merge(prep, on="repurposing_id", how="left")

    # structural-gate outcome under the unified rules
    gate = pd.concat([
        new32[["ligand_id", "n_dual_seeds", "structural_pass", "nomination_state"]].assign(
            unified_gate_cohort="stage3b_new32"
        ),
        hist_resc[["ligand_id", "n_dual_seeds", "structural_pass", "nomination_state"]].assign(
            unified_gate_cohort="historical_deposited_156"
        ),
    ], ignore_index=True).rename(columns={
        "ligand_id": "repurposing_id",
        "n_dual_seeds": "unified_dual_pass_seeds",
        "structural_pass": "unified_structural_pass",
        "nomination_state": "unified_nomination_state",
    })
    idx = idx.merge(gate, on="repurposing_id", how="left")

    # chemistry alerts, frozen flags take precedence for the historical cohort
    alerts = anno157[[
        "repurposing_id", "pains_any", "brenk", "nih", "beta_lactam",
        "alert_clean_pains_brenk", "aggregation_risk_heuristic",
    ]]
    idx = idx.merge(alerts, on="repurposing_id", how="left")

    def final_status(r: pd.Series) -> str:
        if r.in_frozen_primary_12:
            return "frozen_primary_candidate"
        if r.in_frozen_backup_21:
            return "frozen_backup_candidate"
        if r.in_stage3b_new32_docked:
            return "stage3b_new_" + str(r.unified_nomination_state)
        if r.in_stage3a_new_39:
            return "stage3a_excluded_rule_not_applicable"
        if r.in_chemistry_soft_156:
            return "docked_156_not_nominated"
        if r.in_acid_pool_303:
            return "acid_pool_not_docked_chemistry_soft_fail"
        if r.in_nlrp3_pool_1588:
            return "nlrp3_pool_not_acid_equivalent"
        return "clinical_library_below_nlrp3_threshold"

    idx["final_status"] = idx.apply(final_status, axis=1)

    def evidence_level(r: pd.Series) -> str:
        if r.final_status in {"frozen_primary_candidate", "frozen_backup_candidate"}:
            return "ml_prediction_plus_three_seed_dual_structural_geometry"
        if r.in_stage3b_new32_docked:
            return "ml_prediction_plus_three_seed_dual_structural_geometry"
        if r.in_chemistry_soft_156:
            return "ml_prediction_plus_three_seed_docking"
        if r.in_rank_track_dock_1583:
            return "ml_prediction_plus_single_pose_rank_track_docking_superseded"
        if r.in_nlrp3_pool_1588:
            return "ml_prediction_only"
        return "ml_prediction_only_below_threshold"

    idx["evidence_level"] = idx.apply(evidence_level, axis=1)
    idx["activity_claim"] = "computational_prediction_only_no_experimental_activity_measured"
    idx.to_csv(ARCHIVE / "07_RANKING_AND_SELECTION/ALL_COMPOUNDS_MASTER_INDEX.csv", index=False)

    # ---------------- ranking transitions ----------------
    tr = idx[[
        "repurposing_id", "name", "chembl_id",
        "nlrp3_p_active_historical_v1", "nlrp3_rank_historical_v1_recomputed",
        "nlrp3_p_v2_n513_production", "nlrp3_rank_v2_n513",
    ]].copy()
    du = dock_u[["repurposing_id", "dock_score", "docking_status"]].rename(
        columns={"dock_score": "rank_track_urat1_9dkb_cnnaffinity",
                 "docking_status": "rank_track_urat1_status"})
    dn = dock_n[["repurposing_id", "dock_score", "docking_status"]].rename(
        columns={"dock_score": "rank_track_nlrp3_7alv_cnnaffinity",
                 "docking_status": "rank_track_nlrp3_status"})
    tr = tr.merge(du, on="repurposing_id", how="left").merge(dn, on="repurposing_id", how="left")
    # rank track dock_score is CNNaffinity reported as a negative-signed value in
    # the deposited table; the deposited percentile columns define the direction.
    pp = pareto[["repurposing_id", "s_u_percentile", "s_n_dock_percentile",
                 "s_n_ml_percentile", "s_n_percentile", "pareto_front"]]
    tr = tr.merge(pp, on="repurposing_id", how="left")
    tr["rank_track_urat1_rank"] = rank_desc(tr["s_u_percentile"]).astype("Int64")
    tr["rank_track_nlrp3_dock_rank"] = rank_desc(tr["s_n_dock_percentile"]).astype("Int64")
    tr = tr.merge(
        idx[["repurposing_id", "unified_dual_pass_seeds", "unified_structural_pass",
             "unified_nomination_state", "final_status"]],
        on="repurposing_id", how="left",
    )
    tr["rank_population_ml"] = "clinical_library_8319"
    tr["rank_population_rank_track"] = "pareto_merged_1580"
    tr["rank_population_structural"] = "acid_docked_156_plus_stage3b_32"
    tr["note_not_comparable"] = (
        "ML probability rank, rank-track docking percentile and structural gate "
        "use different populations and different preparations; do not read across"
    )
    tr.to_csv(ARCHIVE / "07_RANKING_AND_SELECTION/RANKING_TRANSITIONS.csv", index=False)

    # ---------------- screening decisions, long form ----------------
    rows = []
    for _, r in idx.iterrows():
        rid = r.repurposing_id
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S3_clinical_library",
                     "decision": "included", "reason": "deposited clinical library",
                     "rule": "Broad repurposing library as deposited", "executed": True})
        passed = bool(r.in_nlrp3_pool_1588)
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S3_nlrp3_activity_pool",
                     "decision": "pass" if passed else "excluded",
                     "reason": "frozen historical classifier p_active >= 0.5" if passed
                     else "frozen historical classifier p_active < 0.5",
                     "rule": "p_active_nlrp3 >= 0.5", "executed": True})
        if not passed:
            continue
        acid_ok = bool(r.in_acid_pool_303)
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S3_acid_equivalent",
                     "decision": "pass" if acid_ok else "excluded",
                     "reason": "acid equivalent SMARTS matched" if acid_ok
                     else str(r.acid_reject_reason),
                     "rule": "deposited acid SMARTS", "executed": True})
        if not acid_ok:
            continue
        soft_ok = bool(r.in_chemistry_soft_156)
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S3_chemistry_soft",
                     "decision": "pass" if soft_ok else "excluded",
                     "reason": "MW 200-550, TPSA<=140, rotatable<=10, HBD<=5, HBA<=10"
                     if soft_ok else "failed one or more soft physicochemical bounds",
                     "rule": "chemistry_soft", "executed": True})
        if not soft_ok:
            continue
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S7_dock_156_three_seed",
                     "decision": "executed", "reason": "docked on 9DKB and 7ALV, seeds 42/43/44",
                     "rule": "GNINA exhaustiveness 32, num_modes 9, cnn_scoring rescore",
                     "executed": True})
        st = r.unified_structural_pass
        rows.append({"repurposing_id": rid, "name": r["name"],
                     "stage": "S8_unified_structural_gate",
                     "decision": "pass" if st is True else ("fail" if st is False else "not_evaluated"),
                     "reason": f"dual-pass seeds = {r.unified_dual_pass_seeds}",
                     "rule": "URAT1 A2 and NLRP3 structural on the same seed, >=2 of 3 seeds",
                     "executed": bool(pd.notna(st))})
        rows.append({"repurposing_id": rid, "name": r["name"], "stage": "S9_frozen_nomination",
                     "decision": "primary" if r.in_frozen_primary_12 else
                     ("backup" if r.in_frozen_backup_21 else "not_nominated"),
                     "reason": r.final_status,
                     "rule": "C5 freeze: A1 seed42 primary branch; backup parallel branch",
                     "executed": True})
    for _, r in stage3a.iterrows():
        rows.append({"repurposing_id": r.repurposing_id, "name": r["name"],
                     "stage": "S9_stage3a_scope",
                     "decision": "dock" if r.docking_decision == "dock" else "rule_not_applicable",
                     "reason": r.reason, "rule": "carboxylate required for the Arg477 geometry rule",
                     "executed": True})
    pd.DataFrame(rows).to_csv(
        ARCHIVE / "07_RANKING_AND_SELECTION/ALL_SCREENING_DECISIONS.csv", index=False
    )

    # ---------------- final candidate evidence ----------------
    cand = []
    for _, r in reportable.iterrows():
        g = idx[idx.repurposing_id == r.ligand_id]
        cand.append({
            "ligand_id": r.ligand_id, "name": r["name"], "chembl_id": r.chembl_id,
            "max_phase": r.max_phase, "cohort": "frozen_2026_09_05",
            "tier": r.tier,
            "group": "primary_12" if r.ligand_id in set(primary.ligand_id) else "backup_21",
            "nlrp3_p_active_historical_v1": r.p_active_nlrp3,
            "nlrp3_percentile": r.nlrp3_percentile,
            "frozen_n_seed_pass": r.n_seed_pass,
            "frozen_acid_arg477_seed42_A": r.acid_arg477_min_A_seed42,
            "frozen_ifp_jaccard_seed42": r.ifp_jaccard_seed42,
            "frozen_n_key_contacts_seed42": r.n_key_contacts_seed42,
            "pains_any": r.pains_any, "brenk": r.brenk,
            "beta_lactam_flag": r.beta_lactam_flag,
            "is_structural_control_not_candidate": r.is_structural_control_not_candidate,
            "frozen_excluded_reason": r.excluded_reason,
            "unified_dual_pass_seeds": (g.unified_dual_pass_seeds.iloc[0] if len(g) else pd.NA),
            "unified_structural_pass": (g.unified_structural_pass.iloc[0] if len(g) else pd.NA),
            "unified_nomination_state": (g.unified_nomination_state.iloc[0] if len(g) else pd.NA),
            "md_executed": False,
            "experimental_activity": "none",
            "evidence_level": "ml_prediction_plus_three_seed_dual_structural_geometry",
        })
    for _, r in new32.iterrows():
        cand.append({
            "ligand_id": r.ligand_id, "name": r["name"], "chembl_id": r.chembl_id,
            "max_phase": r.max_phase, "cohort": "stage3b_new_2026_10_15",
            "tier": "", "group": "stage3b_new32",
            "nlrp3_p_active_historical_v1": pd.NA, "nlrp3_percentile": pd.NA,
            "frozen_n_seed_pass": pd.NA, "frozen_acid_arg477_seed42_A": pd.NA,
            "frozen_ifp_jaccard_seed42": pd.NA, "frozen_n_key_contacts_seed42": pd.NA,
            "pains_any": pd.NA, "brenk": pd.NA,
            "beta_lactam_flag": r.beta_lactam,
            "is_structural_control_not_candidate": False,
            "frozen_excluded_reason": "",
            "unified_dual_pass_seeds": r.n_dual_seeds,
            "unified_structural_pass": r.structural_pass,
            "unified_nomination_state": r.nomination_state,
            "md_executed": False,
            "experimental_activity": "none",
            "evidence_level": "ml_prediction_plus_three_seed_dual_structural_geometry",
        })
    pd.DataFrame(cand).to_csv(
        ARCHIVE / "07_RANKING_AND_SELECTION/FINAL_CANDIDATE_EVIDENCE.csv", index=False
    )

    print("master index rows:", len(idx))
    print("final_status tally:\n", idx.final_status.value_counts().to_string())
    print("decisions rows:", len(rows))


if __name__ == "__main__":
    main()
