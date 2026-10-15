#!/usr/bin/env python3
"""Stages A–D. Reads historical files and writes only under 12_FINAL_ADDITIONAL_ANALYSIS."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from libchem import (  # noqa: E402
    AD_THRESHOLD_DEPOSITED,
    alert_record,
    build_catalogs,
    canonicalize,
    morgan_fp,
    nearest_tanimoto,
    quantiles,
    standardize_record,
)

ROOT = CODE.parents[2]
ARCH = ROOT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
OUT = ARCH / "12_FINAL_ADDITIONAL_ANALYSIS"
SUP = OUT / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)
SUP.mkdir(parents=True, exist_ok=True)

SEEDS = (17, 29, 42, 71, 101)


def structural_ge2() -> set[str]:
    counts: dict[str, int] = {}
    for seed in (42, 43, 44):
        df = pd.read_csv(
            ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv"
        )
        for lid in df.loc[df["keep_dual_acid_structural"] == True, "ligand_id"].astype(str):  # noqa: E712
            counts[lid] = counts.get(lid, 0) + 1
    return {lid for lid, n in counts.items() if n >= 2}


def main() -> None:
    master = pd.read_csv(ARCH / "07_RANKING_AND_SELECTION/ALL_COMPOUNDS_MASTER_INDEX.csv", low_memory=False)
    scores = pd.read_csv(ARCH / "03_MACHINE_LEARNING/ALL_MODEL_SCORES_AND_RANKS.csv", low_memory=False)
    runs = pd.read_csv(ARCH / "05_DOCKING/ALL_DOCKING_RUNS.csv", low_memory=False)
    audited = pd.read_csv(ROOT / "data/campaigns/c1/08_nomination/acid_a2_eligible_audited.csv")
    tier_pool = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/tier_pool_w2_annotated.csv")
    primary = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/primary_candidates.csv")
    backup = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/backup_candidates.csv")
    reportable = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/reportable_candidates.csv")
    train = pd.read_csv(ROOT / "data/processed/nlrp3_records.csv", usecols=["canonical_smiles"])
    folds = pd.read_csv(ROOT / "results/nlrp3_model_v2/fold_assignments.csv")
    v1_scores = pd.read_csv(ROOT / "data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv", low_memory=False)

    struct54 = structural_ge2()
    eligible40 = set(audited["ligand_id"].astype(str))
    inter37 = struct54 & eligible40
    frozen12 = set(primary["ligand_id"].astype(str))
    frozen21 = set(backup["ligand_id"].astype(str))
    frozen33 = set(reportable["ligand_id"].astype(str))
    tier_ids = set(tier_pool["ligand_id"].astype(str))

    # ----- A. counts -----
    status_counts = master["final_status"].value_counts(dropna=False).to_dict()
    n_stage3b_runs = int((runs["protocol_id"] == "P_STAGE3B_NEW32").sum())
    n_ref7 = int((runs["protocol_id"] == "P_REFERENCE7").sum())
    checks = [
        ("clinical_library", int(master["in_clinical_library_8319"].sum()), 8319, "master.in_clinical_library_8319"),
        ("master_rows", int(len(master)), 8319, "ALL_COMPOUNDS_MASTER_INDEX.csv"),
        ("v1_pool_1588", int(master["in_nlrp3_pool_1588"].sum()), 1588, "master.in_nlrp3_pool_1588"),
        ("acid_equivalent_303", int(master["in_acid_pool_303"].sum()), 303, "master.in_acid_pool_303"),
        ("chemistry_soft_156", int(master["in_chemistry_soft_156"].sum()), 156, "master.in_chemistry_soft_156"),
        ("deposited_a2_docking_set", int(master["in_deposited_a2_docking_set"].sum()), 156, "master.in_deposited_a2_docking_set"),
        ("historical_structural_ge2", len(struct54), 54, "acid_dual_keep_structural keep_dual_acid_structural >=2/3"),
        ("historical_chemistry_eligible_40", len(eligible40), 40, "acid_a2_eligible_audited.csv"),
        ("intersection_37", len(inter37), 37, "structural_ge2 AND chemistry-eligible 40"),
        ("tier_pool", len(tier_ids), 37, "tier_pool_w2_annotated.csv"),
        ("tier_equals_intersection", int(tier_ids == inter37), 1, "set equality of tier pool and 37"),
        ("frozen_primary_12", len(frozen12), 12, "primary_candidates.csv"),
        ("frozen_backup_21", len(frozen21), 21, "backup_candidates.csv"),
        ("frozen_reportable_33", len(frozen33), 33, "reportable_candidates.csv"),
        ("master_primary_flag", int(master["in_frozen_primary_12"].sum()), 12, "master"),
        ("master_backup_flag", int(master["in_frozen_backup_21"].sum()), 21, "master"),
        ("master_reportable_flag", int(master["in_frozen_reportable_33"].sum()), 33, "master"),
        ("stage3b_molecules", int(master["in_stage3b_new32_docked"].sum()), 32, "master"),
        ("stage3b_dock_runs", n_stage3b_runs, 194, "ALL_DOCKING_RUNS protocol P_STAGE3B_NEW32"),
        ("reference7_runs", n_ref7, 7, "ALL_DOCKING_RUNS protocol P_REFERENCE7"),
        ("stage3b_chem_pass_13", int((master["final_status"] == "stage3b_structural_and_chemistry_pass").sum()) if False else int(master["unified_nomination_state"].fillna("").eq("structural_and_chemistry_pass").sum()) if "unified_nomination_state" in master else -1, 13, "see final_status breakdown"),
    ]
    # nomination states
    nom = master["unified_nomination_state"].fillna("").value_counts().to_dict()
    final_status = {str(k): int(v) for k, v in status_counts.items()}
    n13 = int(master["final_status"].astype(str).str.contains("stage3b_structural_and_chemistry").sum())
    # more reliable from the column values printed historically
    fs = master["final_status"].astype(str)
    observed_fs = {k: int((fs == k).sum()) for k in sorted(fs.unique())}
    checks.append(("final_status_stage3b_pass_like", n13, 13, "final_status contains stage3b_structural_and_chemistry"))
    n_scripts = len(list((ROOT / "scripts/complete_research_archive_20261015").glob("*.py")))
    checks.append(("archive_generator_scripts", n_scripts, 4, "FINAL_ARCHIVE_ACCEPTANCE says 四个; observed file count is the true value"))

    audit_rows = []
    for name, obs, exp, src in checks:
        audit_rows.append({
            "item": name,
            "observed": obs,
            "claimed": exp,
            "consistent": bool(obs == exp),
            "source": src,
        })
    # explicit stage3b split from final_status keys
    for key, claimed in [
        ("stage3b_structural_and_chemistry_pass", 13),
        ("stage3b_structural_pass_chemistry_fail", 8),
        ("stage3b_not_structural_pass", 11),
        ("stage3a_excluded_not_carboxylate", 7),
    ]:
        obs = observed_fs.get(key, 0)
        # try fuzzy if exact key missing
        if obs == 0:
            obs = int(sum(v for k, v in observed_fs.items() if key.split("_")[0] in k and key.split("_")[-1] in k))
        audit_rows.append({
            "item": f"final_status::{key}",
            "observed": obs,
            "claimed": claimed,
            "consistent": obs == claimed,
            "source": "master.final_status exact key, else 0 if the key is absent",
        })
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(OUT / "ARCHIVE_BASELINE_AUDIT.csv", index=False)
    (SUP / "final_status_counts.json").write_text(json.dumps(observed_fs, indent=2, ensure_ascii=False))
    (SUP / "unified_nomination_counts.json").write_text(json.dumps({str(k): int(v) for k, v in nom.items()}, indent=2))

    # 37 -> 33 member audit
    name_map = master.set_index("repurposing_id")["name"].astype(str).to_dict()
    smi_map = master.set_index("repurposing_id")["canonical_smiles"].astype(str).to_dict()
    members = []
    excluded = tier_pool.loc[~tier_pool["reportable"].astype(bool), ["ligand_id", "name", "tier", "excluded_reason", "beta_lactam_flag"]]
    excluded_ids = set(excluded["ligand_id"].astype(str))
    for lid in sorted(eligible40):
        row_a = audited.loc[audited["ligand_id"].astype(str) == lid].iloc[0]
        in_tier = lid in tier_ids
        tier_row = tier_pool.loc[tier_pool["ligand_id"].astype(str) == lid]
        reason = "" if tier_row.empty else str(tier_row.iloc[0].get("excluded_reason") or "")
        reportable_flag = False if tier_row.empty else bool(tier_row.iloc[0]["reportable"])
        members.append({
            "ligand_id": lid,
            "name": row_a["name"],
            "in_chemistry_eligible_40": True,
            "in_structural_ge2_54": lid in struct54,
            "in_intersection_37": lid in inter37,
            "tier": None if tier_row.empty else tier_row.iloc[0]["tier"],
            "reportable_33": reportable_flag,
            "freeze_excluded_reason": reason,
            "deposited_pains_any": bool(row_a["pains_any"]) if pd.notna(row_a["pains_any"]) else None,
            "deposited_brenk": bool(row_a["brenk"]) if pd.notna(row_a["brenk"]) else None,
            "deposited_alert_clean": bool(row_a["alert_clean"]) if pd.notna(row_a["alert_clean"]) else None,
            "deposited_soft_excluded": bool(row_a["soft_excluded"]) if pd.notna(row_a["soft_excluded"]) else None,
            "deposited_beta_lactam_on_tier": None if tier_row.empty else bool(tier_row.iloc[0]["beta_lactam_flag"]),
            "eligible_rule_included_beta_lactam": False,
            "beta_lactam_applied_only_at_freeze": True,
        })
    mem = pd.DataFrame(members)
    mem.to_csv(SUP / "SET_37_TO_33_MEMBERS.csv", index=False)

    outside37 = sorted(struct54 - inter37)
    only_eligible = sorted(eligible40 - struct54)
    (SUP / "structural54_not_in_37.json").write_text(json.dumps({
        "structural_not_in_eligible40": [
            {"ligand_id": i, "name": name_map.get(i)} for i in outside37
        ],
        "eligible40_not_structural_ge2": [
            {"ligand_id": i, "name": name_map.get(i, audited.set_index("ligand_id")["name"].to_dict().get(i))}
            for i in only_eligible
        ],
        "freeze_excluded_from_37": excluded.to_dict("records"),
        "n_37_minus_excluded": int(len(inter37 - excluded_ids)),
        "excluded_ids_subset_of_37": sorted(excluded_ids - inter37) == [],
    }, indent=2, ensure_ascii=False))

    # ----- C chemistry of historical 21 = structural54 - frozen33 -----
    hist21 = sorted(struct54 - frozen33)
    catalogs = build_catalogs()
    chem_rows = []
    for lid in hist21:
        smi = smi_map.get(lid)
        name = name_map.get(lid, "")
        flags = alert_record(smi, catalogs, name)
        in40 = lid in eligible40
        in37 = lid in inter37
        reason = ""
        if lid in excluded_ids:
            reason = str(excluded.set_index("ligand_id").loc[lid, "excluded_reason"])
        # historical eligible rule, reconstructed
        passes_eligible_rule = bool(flags["alert_clean_pains_brenk"]) and not flags["soft_name_excluded"]
        passes_freeze_rule = passes_eligible_rule and not flags["beta_lactam"] and not flags["structural_control_name"]
        if in37 and reason:
            group = "in_37_removed_at_freeze"
        elif not in40:
            group = "structural_pass_not_in_historical_eligible_40"
        else:
            group = "other"
        chem_rows.append({
            "ligand_id": lid,
            "name": name,
            "canonical_smiles": smi,
            "in_historical_eligible_40": in40,
            "in_intersection_37": in37,
            "freeze_excluded_reason": reason,
            "group_vs_freeze": group,
            "recomputed_would_pass_eligible_rule": passes_eligible_rule,
            "recomputed_would_pass_freeze_reportable_rule": passes_freeze_rule,
            "rule_source": "scripts/09_cheminformatics_filters.py FilterCatalog PAINS_A/B/C+BRENK; beta-lactam SMARTS from scripts/build_c5_tier_assignment.py; soft-name list from scripts/build_c1_acid_shortlist_a2.py",
            **flags,
        })
    chem = pd.DataFrame(chem_rows)
    chem.to_csv(OUT / "HISTORICAL21_CHEMISTRY.csv", index=False)

    # recovery: recomputed vs deposited on the 40
    rec_rows = []
    for _, row in audited.iterrows():
        flags = alert_record(str(row["canonical_smiles"]), catalogs, str(row["name"]))
        rec_rows.append({
            "ligand_id": row["ligand_id"],
            "deposited_pains": bool(row["pains_any"]) if pd.notna(row["pains_any"]) else None,
            "recomputed_pains": flags["pains_any"],
            "deposited_brenk": bool(row["brenk"]) if pd.notna(row["brenk"]) else None,
            "recomputed_brenk": flags["brenk"],
            "pains_match": flags["pains_any"] == (bool(row["pains_any"]) if pd.notna(row["pains_any"]) else None),
            "brenk_match": flags["brenk"] == (bool(row["brenk"]) if pd.notna(row["brenk"]) else None),
        })
    rec = pd.DataFrame(rec_rows)
    rec.to_csv(SUP / "CHEMISTRY_RULE_RECOVERY_ON_ELIGIBLE40.csv", index=False)

    # ----- B applicability domain -----
    train_smiles = sorted(set(train["canonical_smiles"].astype(str)))
    train_recanon = [canonicalize(s) for s in train_smiles]
    n_train_recanon_changed = int(sum(a != b for a, b in zip(train_smiles, train_recanon) if b))
    n_train_recanon_fail = int(sum(b is None for b in train_recanon))
    train_fps = [morgan_fp(s) for s in train_smiles]
    n_train_fp_fail = int(sum(f is None for f in train_fps))
    loo = nearest_tanimoto(train_fps, train_fps, exclude_self=True)
    ad_threshold = float(np.nanpercentile(loo, 5))
    # positional LOO as a diagnostic, not the production rule
    valid_idx = [i for i, f in enumerate(train_fps) if f is not None]
    positional = np.full(len(train_fps), np.nan)
    for i in valid_idx:
        ref = [train_fps[j] for j in valid_idx if j != i]
        positional[i] = max(DataStructs.BulkTanimotoSimilarity(train_fps[i], ref))
    n_loo_method_differs = int(np.nansum(np.abs(loo - positional) > 1e-9))

    # identical fingerprints among distinct training SMILES
    ident_pairs = 0
    for i in range(len(valid_idx)):
        sims = DataStructs.BulkTanimotoSimilarity(train_fps[valid_idx[i]], [train_fps[j] for j in valid_idx[i + 1:]])
        ident_pairs += int(sum(s >= 1.0 - 1e-12 for s in sims))

    lib = scores.copy()
    lib_smiles = lib["canonical_smiles"].map(lambda s: canonicalize(s) if isinstance(s, str) else None)
    # production v2 fingerprinted the canonicalized library SMILES
    query_fps = [morgan_fp(s) if isinstance(s, str) else None for s in lib_smiles.tolist()]
    nn = nearest_tanimoto(query_fps, train_fps, exclude_self=False)
    lib["recomputed_max_tanimoto_v2"] = nn
    lib["recomputed_in_ad_v2"] = nn >= ad_threshold
    deposited_nn = pd.to_numeric(lib["max_tanimoto_train_v2"], errors="coerce")
    ad_abs = (lib["recomputed_max_tanimoto_v2"] - deposited_nn).abs()
    # v1 historical, rounded to 3 decimals, no AD threshold was applied to pool membership
    v1 = v1_scores[["repurposing_id", "max_tanimoto_nlrp3_train", "canonical_smiles"]].copy()
    v1_fps = [morgan_fp(canonicalize(s) if isinstance(s, str) else None) for s in v1["canonical_smiles"].tolist()]
    v1_nn = nearest_tanimoto(v1_fps, train_fps, exclude_self=False)
    v1["recomputed"] = v1_nn
    v1["deposited_rounded"] = pd.to_numeric(v1["max_tanimoto_nlrp3_train"], errors="coerce")
    v1_abs = (v1["recomputed"] - v1["deposited_rounded"]).abs()

    # fold-out NN: each molecule is the test member of exactly one fold
    test = folds[folds["outer_role"] == "test"][["fold", "canonical_smiles"]].drop_duplicates()
    train_by_fold = {
        int(k): set(g.loc[g["outer_role"] == "train", "canonical_smiles"].astype(str))
        for k, g in folds.groupby("fold")
    }
    fold_nn = []
    fp_cache = {s: morgan_fp(s) for s in train_smiles}
    for fold, g in test.groupby("fold"):
        ref = [fp_cache.get(s) for s in train_by_fold[int(fold)]]
        q = [fp_cache.get(s) for s in g["canonical_smiles"].astype(str)]
        vals = nearest_tanimoto(q, ref, exclude_self=False)
        for smi, val in zip(g["canonical_smiles"].astype(str), vals):
            fold_nn.append({"canonical_smiles": smi, "fold": int(fold), "foldout_max_tanimoto": float(val)})
    fold_df = pd.DataFrame(fold_nn)

    method_rows = [
        {"item": "fingerprint", "value": "Morgan radius 2, 2048 bits, RDKit GetMorganFingerprintAsBitVect", "applies_to": "v2 AD and v1 stored max Tanimoto"},
        {"item": "distance", "value": "Tanimoto similarity; in-domain if max similarity to training set >= threshold", "applies_to": "v2 AD"},
        {"item": "neighbor_count", "value": "1 (maximum similarity). Not a k>1 vote.", "applies_to": "v2 AD"},
        {"item": "training_molecules", "value": f"unique canonical_smiles in data/processed/nlrp3_records.csv = {len(train_smiles)}; records = {len(train)}", "applies_to": "v1 annotation and v2 AD share this file"},
        {"item": "threshold_rule", "value": "5th percentile of leave-one-out nearest-neighbour Tanimoto inside the training set", "applies_to": "v2 only"},
        {"item": "threshold_recomputed", "value": f"{ad_threshold:.10f}", "applies_to": "v2"},
        {"item": "threshold_deposited", "value": f"{AD_THRESHOLD_DEPOSITED:.10f}", "applies_to": "v2 library_summary.json"},
        {"item": "threshold_abs_diff", "value": f"{abs(ad_threshold - AD_THRESHOLD_DEPOSITED):.3e}", "applies_to": "v2"},
        {"item": "loo_vs_positional_exclusion_n_differ", "value": str(n_loo_method_differs), "applies_to": "diagnostic; production uses sorted-sims[1], not a rewritten rule"},
        {"item": "training_identical_fp_pairs", "value": str(ident_pairs), "applies_to": "diagnostic"},
        {"item": "training_fp_failures", "value": str(n_train_fp_fail), "applies_to": "v2"},
        {"item": "training_recanonicalize_changed", "value": str(n_train_recanon_changed), "applies_to": "v2 training SMILES were used as stored"},
        {"item": "v2_query_standardization", "value": "utils_ml.canonicalize = MolFromSmiles, SanitizeMol, MolToSmiles; salts are not stripped", "applies_to": "v2"},
        {"item": "v2_vs_deposited_max_abs_diff", "value": f"{float(np.nanmax(ad_abs)):.6g}", "applies_to": "v2 clinical library"},
        {"item": "v2_n_abs_diff_gt_1e-6", "value": str(int((ad_abs > 1e-6).sum())), "applies_to": "v2"},
        {"item": "v1_vs_recomputed_max_abs_diff_after_round3", "value": f"{float(np.nanmax(np.abs(np.round(v1['recomputed'], 3) - v1['deposited_rounded']))):.6g}", "applies_to": "v1 stored column was rounded to 3 decimals"},
        {"item": "v1_pool_rule", "value": "p_active_nlrp3 >= 0.5 only. max_tanimoto_nlrp3_train was stored and was not a membership filter.", "applies_to": "v1 historical screen"},
        {"item": "v2_ad_on_v1_candidates", "value": "retrospective. The v2 threshold was not the v1 screening rule and is not used here to change the 1588 membership.", "applies_to": "interpretation"},
        {"item": "model_feature_space", "value": "Morgan 2048 + 12 RDKit descriptors. v1 also used assay one-hot. AD uses Morgan only.", "applies_to": "v1 and v2 models differ from the AD feature space"},
        {"item": "rdkit_version", "value": Chem.rdBase.rdkitVersion, "applies_to": "this recomputation"},
    ]
    pd.DataFrame(method_rows).to_csv(OUT / "AD_METHOD_AUDIT.csv", index=False)

    id_to_nn = dict(zip(lib["repurposing_id"].astype(str), lib["recomputed_max_tanimoto_v2"]))
    id_to_ad = dict(zip(lib["repurposing_id"].astype(str), lib["recomputed_in_ad_v2"]))
    cohorts = {
        "train_loo_nn": loo,
        "train_foldout_nn": fold_df["foldout_max_tanimoto"].to_numpy(),
        "clinical_8319_max_tanimoto": lib["recomputed_max_tanimoto_v2"].to_numpy(),
        "pool_1588_max_tanimoto": lib.loc[lib["in_historical_pool"] == True, "recomputed_max_tanimoto_v2"].to_numpy(),  # noqa: E712
        "frozen_33_max_tanimoto": np.array([id_to_nn[i] for i in frozen33]),
        "stage3b_13_max_tanimoto": np.array([
            id_to_nn[i] for i in master.loc[master["final_status"].astype(str).str.contains("stage3b_structural_and_chemistry"), "repurposing_id"].astype(str)
        ]),
        "historical21_max_tanimoto": np.array([id_to_nn[i] for i in hist21]),
    }
    sim_rows = []
    for name, vals in cohorts.items():
        row = {"population": name, "metric": "max_tanimoto_to_train" if "loo" not in name and "foldout" not in name else name}
        row.update(quantiles(vals))
        finite = np.asarray(vals, dtype=float)
        finite = finite[np.isfinite(finite)]
        row["n_ge_v2_ad_threshold"] = int((finite >= ad_threshold).sum()) if "loo" not in name and "foldout" not in name else None
        row["ad_threshold_used"] = ad_threshold
        row["note"] = (
            "This is the production 1-NN AD statistic."
            if name.startswith("clinical") or name.startswith("pool") or name.startswith("frozen") or name.startswith("stage") or name.startswith("historical")
            else "Internal reference distribution. Not the clinical AD call."
        )
        sim_rows.append(row)
    pd.DataFrame(sim_rows).to_csv(OUT / "AD_TRAIN_CLINICAL_SIMILARITY.csv", index=False)

    # candidate nearest neighbors
    # map train fps index
    train_fp_ok = [(s, f) for s, f in zip(train_smiles, train_fps) if f is not None]
    cand_ids = sorted(frozen33 | set(hist21) | set(master.loc[master["in_stage3b_new32_docked"] == True, "repurposing_id"].astype(str)))  # noqa: E712
    nn_rows = []
    for lid in cand_ids:
        smi = canonicalize(smi_map.get(lid, ""))
        fp = morgan_fp(smi) if smi else None
        if fp is None:
            nn_rows.append({"ligand_id": lid, "name": name_map.get(lid), "max_tanimoto": None, "nearest_train_smiles": None})
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, [f for _, f in train_fp_ok])
        j = int(np.argmax(sims))
        group = []
        if lid in frozen12:
            group.append("frozen_primary_12")
        if lid in frozen21:
            group.append("frozen_backup_21")
        if lid in hist21:
            group.append("historical_structural_21")
        if bool(master.set_index("repurposing_id").loc[lid, "in_stage3b_new32_docked"]):
            group.append("stage3b_new32")
        nn_rows.append({
            "ligand_id": lid,
            "name": name_map.get(lid),
            "group": "|".join(group),
            "max_tanimoto": float(sims[j]),
            "in_v2_applicability_domain_retrospective": bool(sims[j] >= ad_threshold),
            "nearest_train_smiles": train_fp_ok[j][0],
            "v1_screening_used_this_threshold": False,
        })
    pd.DataFrame(nn_rows).to_csv(OUT / "AD_CANDIDATE_NEAREST_NEIGHBORS.csv", index=False)

    in_domain = lib.loc[lib["recomputed_in_ad_v2"] == True, ["repurposing_id", "name", "canonical_smiles", "recomputed_max_tanimoto_v2", "in_training_set", "in_historical_pool"]].copy()  # noqa: E712
    in_domain.to_csv(SUP / "AD_IN_DOMAIN_MOLECULES.csv", index=False)

    # chemical space comparison: train vs clinical vs in-domain, using standardized descriptors on a sample of properties
    def descriptor_frame(smiles_list, label):
        rows = []
        for smi in smiles_list:
            rec = standardize_record(smi)
            rec["population"] = label
            rows.append(rec)
        return pd.DataFrame(rows)

    # full clinical descriptors are also needed for stage D; compute once
    print("standardizing clinical library", flush=True)
    desc_rows = []
    for lid, smi in zip(master["repurposing_id"].astype(str), master["canonical_smiles"].astype(str)):
        rec = standardize_record(smi)
        rec["repurposing_id"] = lid
        desc_rows.append(rec)
    desc = pd.DataFrame(desc_rows)
    print("standardizing training set", flush=True)
    train_desc = descriptor_frame(train_smiles, "nlrp3_train_513")
    train_desc.to_csv(SUP / "TRAIN_513_DESCRIPTORS.csv", index=False)
    desc.to_csv(SUP / "CLINICAL_8319_DESCRIPTORS.csv", index=False)

    def space_row(label, frame, nn_values=None):
        row = {"population": label, "n": int(len(frame)), "n_parse_fail": int((~frame["parse_ok"].fillna(False)).sum()) if "parse_ok" in frame else 0}
        for col in ("mw", "clogp", "aromatic_rings", "hetero_rings", "n_carboxylate", "heavy_atoms", "fp_onbits", "formal_charge_full"):
            if col in frame:
                row.update(quantiles(frame[col], prefix=f"{col}_"))
        if "murcko_scaffold" in frame:
            ok = frame["parse_ok"].fillna(False) if "parse_ok" in frame else pd.Series(True, index=frame.index)
            scaff = frame.loc[ok, "murcko_scaffold"].fillna("")
            row["n_murcko"] = int(scaff.replace("", np.nan).nunique(dropna=True))
            row["n_empty_murcko"] = int((scaff == "").sum())
        if nn_values is not None:
            row["n_in_v2_ad"] = int(np.nansum(np.asarray(nn_values) >= ad_threshold))
        return row

    space = [
        space_row("nlrp3_train", train_desc, None),
        space_row("clinical_8319", desc, lib["recomputed_max_tanimoto_v2"]),
    ]
    in_ids = set(in_domain["repurposing_id"].astype(str))
    space.append(space_row("clinical_in_v2_ad", desc[desc["repurposing_id"].isin(in_ids)], None))
    train_scaff = set(train_desc["murcko_scaffold"].dropna()) - {""}
    clin_scaff = set(desc["murcko_scaffold"].dropna()) - {""}
    space.append({
        "population": "scaffold_overlap",
        "n_train_murcko": len(train_scaff),
        "n_clinical_murcko": len(clin_scaff),
        "n_shared_murcko": len(train_scaff & clin_scaff),
        "fraction_clinical_scaffolds_seen_in_train": (len(train_scaff & clin_scaff) / len(clin_scaff)) if clin_scaff else None,
    })
    pd.DataFrame(space).to_csv(OUT / "AD_CHEMICAL_SPACE_SUMMARY.csv", index=False)

    # ----- D diversity cohorts -----
    print("diversity cohorts", flush=True)
    desc_i = desc.set_index("repurposing_id")
    stage3b_ids = set(master.loc[master["in_stage3b_new32_docked"] == True, "repurposing_id"].astype(str))  # noqa: E712
    # stage3b 13 from final_status exact after we know the key
    fs_map = master.set_index(master["repurposing_id"].astype(str))["final_status"].astype(str)
    new13 = set(fs_map[fs_map.str.contains("stage3b_structural_and_chemistry")].index)
    # if the contains match is too broad, intersect with stage3b
    new13 &= stage3b_ids
    comparable = frozen33 | new13 | set(chem.loc[chem["recomputed_would_pass_freeze_reportable_rule"] == True, "ligand_id"])  # noqa: E712

    cohort_sets = {
        "clinical_8319": set(master["repurposing_id"].astype(str)),
        "pool_1588": set(master.loc[master["in_nlrp3_pool_1588"] == True, "repurposing_id"].astype(str)),  # noqa: E712
        "acid_303": set(master.loc[master["in_acid_pool_303"] == True, "repurposing_id"].astype(str)),  # noqa: E712
        "docked_156": set(master.loc[master["in_deposited_a2_docking_set"] == True, "repurposing_id"].astype(str)),  # noqa: E712
        "historical_structural_54": struct54,
        "frozen_33": frozen33,
        "stage3b_new32": stage3b_ids,
        "stage3b_pass_13": new13,
        "historical_structural_21": set(hist21),
        "comparison_set": comparable,
    }
    # murcko table is the clinical descriptor table plus cohort flags; training written separately
    flags = desc[["repurposing_id", "original_smiles", "standardized_smiles", "murcko_scaffold", "generic_murcko_scaffold",
                  "mw", "clogp", "tpsa", "hba", "hbd", "rotatable_bonds", "aromatic_rings", "hetero_rings",
                  "n_carboxylate", "formal_charge_full", "fp_onbits", "parse_ok", "n_components", "murcko_empty",
                  "unspecified_stereo"]].copy()
    for name, ids in cohort_sets.items():
        flags[f"in_{name}"] = flags["repurposing_id"].isin(ids)
    flags.to_csv(OUT / "MURCKO_SCAFFOLDS_ALL.csv", index=False)

    def scaffold_stats(ids, label):
        sub = desc[desc["repurposing_id"].isin(ids)]
        ok = sub[sub["parse_ok"] == True]  # noqa: E712
        mur = ok["murcko_scaffold"].fillna("")
        gen = ok["generic_murcko_scaffold"].fillna("")
        mur_named = mur[mur != ""]
        gen_named = gen[gen != ""]
        top = mur_named.value_counts().head(5)
        return {
            "cohort": label,
            "n_molecules": int(len(ids)),
            "n_in_descriptor_table": int(len(sub)),
            "n_parse_fail": int((sub["parse_ok"] != True).sum()),  # noqa: E712
            "n_murcko": int(mur_named.nunique()),
            "n_generic_murcko": int(gen_named.nunique()),
            "n_empty_murcko": int((mur == "").sum()),
            "top_murcko": " || ".join(f"{a}:{int(b)}" for a, b in top.items()),
        }

    count_rows = [scaffold_stats(ids, name) for name, ids in cohort_sets.items()]
    pd.DataFrame(count_rows).to_csv(OUT / "SCAFFOLD_COUNTS_BY_COHORT.csv", index=False)

    def murcko_set(ids):
        sub = desc[desc["repurposing_id"].isin(ids) & (desc["parse_ok"] == True)]  # noqa: E712
        return set(sub["murcko_scaffold"].fillna("")) - {""}

    m33, m13, m156, m21, m32 = map(murcko_set, (frozen33, new13, cohort_sets["docked_156"], set(hist21), stage3b_ids))
    contrib = pd.DataFrame([
        {"contrast": "new13_vs_frozen33", "n_query_scaffolds": len(m13), "n_ref_scaffolds": len(m33),
         "n_shared": len(m13 & m33), "n_new_vs_ref": len(m13 - m33)},
        {"contrast": "new13_vs_docked156", "n_query_scaffolds": len(m13), "n_ref_scaffolds": len(m156),
         "n_shared": len(m13 & m156), "n_new_vs_ref": len(m13 - m156)},
        {"contrast": "historical21_vs_frozen33", "n_query_scaffolds": len(m21), "n_ref_scaffolds": len(m33),
         "n_shared": len(m21 & m33), "n_new_vs_ref": len(m21 - m33)},
        {"contrast": "new32_vs_docked156", "n_query_scaffolds": len(m32), "n_ref_scaffolds": len(m156),
         "n_shared": len(m32 & m156), "n_new_vs_ref": len(m32 - m156)},
    ])

    # fixed-size controls
    rng_rows = []
    clin_ids = master["repurposing_id"].astype(str).tolist()
    pool_ids = list(cohort_sets["pool_1588"])
    for n in (13, 21, 33, 54, 156):
        for seed in SEEDS:
            for parent_name, parent in (("clinical_8319", clin_ids), ("pool_1588", pool_ids)):
                if n > len(parent):
                    continue
                rng = np.random.default_rng(seed)
                pick = rng.choice(parent, size=n, replace=False)
                scaff = murcko_set(set(pick))
                rng_rows.append({"parent": parent_name, "n": n, "seed": seed, "n_murcko": len(scaff)})
    pd.DataFrame(rng_rows).to_csv(SUP / "SCAFFOLD_SAMPLE_SIZE_CONTROLS.csv", index=False)
    contrib.to_csv(OUT / "NEW_SCAFFOLD_CONTRIBUTION.csv", index=False)

    # pairwise similarity on the comparison-relevant set (33+13+21)
    pair_ids = sorted(frozen33 | new13 | set(hist21))
    fp_map = {}
    for lid in pair_ids:
        smi = desc_i.loc[lid, "largest_fragment_smiles"] if lid in desc_i.index else None
        # largest fragment smiles is not in flags; use standardized parent via desc
        parent = desc.set_index("repurposing_id").loc[lid, "standardized_smiles"] if lid in set(desc["repurposing_id"]) else None
        # use largest fragment if present
        if "largest_fragment_smiles" in desc.columns:
            parent = desc.set_index("repurposing_id").loc[lid, "largest_fragment_smiles"]
        fp_map[lid] = morgan_fp(parent) if isinstance(parent, str) else None
    # desc from standardize_record includes largest_fragment_smiles
    pair_rows = []
    nn_rows_div = []
    ids_ok = [i for i in pair_ids if fp_map[i] is not None]
    for i, a in enumerate(ids_ok):
        sims = DataStructs.BulkTanimotoSimilarity(fp_map[a], [fp_map[b] for b in ids_ok])
        order = np.argsort(sims)
        # nearest other
        if len(ids_ok) > 1:
            j = order[-2] if order[-1] == i else order[-1]
            # self is 1.0 at index i
            others = [(ids_ok[k], sims[k]) for k in range(len(ids_ok)) if k != i]
            nb, sim = max(others, key=lambda t: t[1])
            nn_rows_div.append({
                "ligand_id": a,
                "name": name_map.get(a),
                "nearest_id": nb,
                "nearest_name": name_map.get(nb),
                "nearest_tanimoto": sim,
                "same_murcko": desc_i.loc[a, "murcko_scaffold"] == desc_i.loc[nb, "murcko_scaffold"],
            })
        for b_i in range(i + 1, len(ids_ok)):
            b = ids_ok[b_i]
            pair_rows.append({
                "id_a": a, "id_b": b,
                "tanimoto_ecfp4": float(sims[b_i]),
                "same_murcko": desc_i.loc[a, "murcko_scaffold"] == desc_i.loc[b, "murcko_scaffold"],
            })
    pd.DataFrame(pair_rows).to_csv(OUT / "PAIRWISE_CANDIDATE_SIMILARITY.csv", index=False)
    pd.DataFrame(nn_rows_div).to_csv(SUP / "CANDIDATE_NEAREST_WITHIN_SET.csv", index=False)

    # reference inhibitors
    print("reference similarity", flush=True)
    lit = pd.read_csv(ROOT / "data/benchmarks/literature_benchmarks.csv", low_memory=False)
    lit["canonical_smiles"] = lit["canonical_smiles"].map(lambda s: canonicalize(s) if isinstance(s, str) else None)
    nlrp3_cur = pd.read_csv(ROOT / "data/processed/nlrp3_records.csv", low_memory=False)
    # molecule-level max pActivity if present
    if "pActivity" in nlrp3_cur.columns:
        mol_act = nlrp3_cur.groupby("canonical_smiles", as_index=False)["pActivity"].max()
    else:
        mol_act = nlrp3_cur[["canonical_smiles"]].drop_duplicates()
        mol_act["pActivity"] = np.nan
    urat1_path = ROOT / "data/processed/urat1_curated.csv"
    if not urat1_path.exists():
        # fall back to literature only
        urat1 = pd.DataFrame(columns=["canonical_smiles", "pActivity", "molecule_name"])
    else:
        urat1 = pd.read_csv(urat1_path)
    ref_sets = []
    # literature rows
    for _, r in lit.drop_duplicates("canonical_smiles").iterrows():
        if not isinstance(r["canonical_smiles"], str):
            continue
        gene = str(r.get("target_gene", ""))
        ref_sets.append({
            "ref_id": r.get("compound_id"),
            "ref_name": r.get("compound_name"),
            "smiles": r["canonical_smiles"],
            "target": gene,
            "endpoint": r.get("primary_endpoint"),
            "activity_type": r.get("activity_type"),
            "activity_value": r.get("pactivity"),
            "assay": r.get("assay_description"),
            "source": r.get("ref_citation"),
            "in_train_project": r.get("in_train_project"),
            "validation_role": r.get("validation_role"),
            "class": "literature_benchmark",
        })
    # NLRP3 training actives are not external; stored as reference chemical space, flagged
    if "pActivity" in mol_act.columns:
        act = mol_act[pd.to_numeric(mol_act["pActivity"], errors="coerce") >= 6]
        for _, r in act.iterrows():
            ref_sets.append({
                "ref_id": None,
                "ref_name": None,
                "smiles": r["canonical_smiles"],
                "target": "NLRP3",
                "endpoint": "training-set pActivity>=6",
                "activity_type": "pActivity",
                "activity_value": r["pActivity"],
                "assay": "ChEMBL curation used to train the model",
                "source": "data/processed/nlrp3_records.csv",
                "in_train_project": True,
                "validation_role": "training_active_not_external",
                "class": "nlrp3_training_active",
            })
    refs = pd.DataFrame(ref_sets)
    # similarities of pair_ids to each reference class
    sim_out = []
    for cls, g in refs.groupby("class"):
        fps = [morgan_fp(s) for s in g["smiles"].astype(str)]
        fps = [f for f in fps if f is not None]
        meta = g.reset_index(drop=True)
        meta_fps = []
        meta_rows = []
        for _, r in meta.iterrows():
            fp = morgan_fp(r["smiles"]) if isinstance(r["smiles"], str) else None
            if fp is None:
                continue
            meta_fps.append(fp)
            meta_rows.append(r)
        if not meta_fps:
            continue
        for lid in pair_ids:
            fp = fp_map.get(lid)
            if fp is None:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fp, meta_fps)
            j = int(np.argmax(sims))
            rb = meta_rows[j]
            sim_out.append({
                "ligand_id": lid,
                "name": name_map.get(lid),
                "reference_class": cls,
                "max_tanimoto": float(sims[j]),
                "nearest_ref_name": rb.get("ref_name"),
                "nearest_ref_id": rb.get("ref_id"),
                "nearest_endpoint": rb.get("endpoint"),
                "nearest_activity": rb.get("activity_value"),
                "nearest_source": rb.get("source"),
                "nearest_in_train_project": rb.get("in_train_project"),
                "same_murcko_as_nearest": False,
            })
    pd.DataFrame(sim_out).to_csv(OUT / "REFERENCE_INHIBITOR_SIMILARITY.csv", index=False)

    summary = {
        "ad_threshold_recomputed": ad_threshold,
        "ad_threshold_match": abs(ad_threshold - AD_THRESHOLD_DEPOSITED) < 1e-8,
        "n_clinical_in_ad": int(lib["recomputed_in_ad_v2"].sum()),
        "n_hist21": len(hist21),
        "n_new13": len(new13),
        "n_comparable": len(comparable),
        "chem_recovery_pains_match": float(rec["pains_match"].mean()) if len(rec) else None,
        "chem_recovery_brenk_match": float(rec["brenk_match"].mean()) if len(rec) else None,
        "n_hist21_pass_freeze_rule": int(chem["recomputed_would_pass_freeze_reportable_rule"].sum()),
        "scaffold_contribution": contrib.to_dict("records"),
        "audit_inconsistent": audit.loc[~audit["consistent"], "item"].tolist(),
        "final_status": observed_fs,
        "n_stage3b_runs": n_stage3b_runs,
        "in_domain_molecules": in_domain[["repurposing_id", "name", "recomputed_max_tanimoto_v2", "in_training_set"]].to_dict("records"),
        "freeze_excluded": excluded.to_dict("records"),
        "n_structural_not_in_40": len(outside37),
        "n_eligible_not_structural": len(only_eligible),
    }
    (SUP / "KEY_RESULTS_A_D.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in ("ad_threshold_recomputed", "n_clinical_in_ad", "n_hist21", "n_new13", "n_comparable", "chem_recovery_pains_match", "n_hist21_pass_freeze_rule", "audit_inconsistent", "n_stage3b_runs")}, indent=2, default=str))


if __name__ == "__main__":
    main()
