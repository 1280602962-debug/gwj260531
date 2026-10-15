#!/usr/bin/env python3
"""Patch stage A–D cohort labels, then run stages C (groups), E–I."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from libchem import dist, morgan_fp, quantiles  # noqa: E402

ROOT = CODE.parents[2]
ARCH = ROOT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
OUT = ARCH / "12_FINAL_ADDITIONAL_ANALYSIS"
SUP = OUT / "supplementary"
SEEDS = (17, 29, 42, 71, 101)
ARG_THR = 7.7027
CENT_THR = 6.0
CNN_THR = 0.5
OV_THR = 0.50
IFP_THR = 0.50
KEY_MIN = 5


def auroc(y, score):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    ranks = np.empty(len(y), dtype=float)
    i = 0
    while i < len(y):
        j = i
        while j + 1 < len(y) and sorted_score[j + 1] == sorted_score[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    sum_pos = ranks[y == 1].sum()
    return float((sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def pr_auc(y, score):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    n_pos = int(y.sum())
    if n_pos == 0:
        return np.nan
    order = np.argsort(-score, kind="mergesort")
    yy = y[order]
    tp = np.cumsum(yy)
    precision = tp / np.arange(1, len(yy) + 1)
    recall = tp / n_pos
    prev = 0.0
    area = 0.0
    for p, r in zip(precision, recall):
        area += (r - prev) * p
        prev = r
    return float(area)


def enrichment(y, score, frac, higher=True):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    n = len(y)
    if n == 0 or y.sum() == 0:
        return np.nan, 0, 0
    k = max(1, int(math.floor(frac * n)))
    signed = -score if higher else score
    order = np.argsort(signed, kind="mergesort")
    hits = int(y[order[:k]].sum())
    return float((hits / k) / (y.sum() / n)), hits, k


def recall_at_k(y, score, k, higher=True):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    if y.sum() == 0:
        return np.nan
    k = min(int(k), len(y))
    signed = -score if higher else score
    order = np.argsort(signed, kind="mergesort")
    return float(y[order[:k]].sum() / y.sum())


def patch_a_d(master, desc):
    fs = master.set_index(master["repurposing_id"].astype(str))["final_status"].astype(str)
    new13 = set(fs[fs == "stage3b_new_structural_and_chemistry_pass"].index)
    new8 = set(fs[fs == "stage3b_new_structural_pass_chemistry_fail"].index)
    new11 = set(fs[fs == "stage3b_new_not_structural_pass"].index)
    excl7 = set(fs[fs == "stage3a_excluded_rule_not_applicable"].index)
    chem = pd.read_csv(OUT / "HISTORICAL21_CHEMISTRY.csv")
    clean2 = set(chem.loc[chem["recomputed_would_pass_freeze_reportable_rule"] == True, "ligand_id"].astype(str))  # noqa: E712
    frozen33 = set(master.loc[master["in_frozen_reportable_33"] == True, "repurposing_id"].astype(str))  # noqa: E712
    comparable = frozen33 | new13 | clean2

    flags = pd.read_csv(OUT / "MURCKO_SCAFFOLDS_ALL.csv")
    flags["in_stage3b_pass_13"] = flags["repurposing_id"].isin(new13)
    flags["in_comparison_set"] = flags["repurposing_id"].isin(comparable)
    flags["in_hist21_chemistry_clean_2"] = flags["repurposing_id"].isin(clean2)
    flags.to_csv(OUT / "MURCKO_SCAFFOLDS_ALL.csv", index=False)

    def murcko_set(ids):
        sub = desc[desc["repurposing_id"].isin(ids) & (desc["parse_ok"] == True)]  # noqa: E712
        return set(sub["murcko_scaffold"].fillna("")) - {""}

    def stats(ids, label):
        sub = desc[desc["repurposing_id"].isin(ids)]
        ok = sub[sub["parse_ok"] == True]  # noqa: E712
        mur = ok["murcko_scaffold"].fillna("")
        gen = ok["generic_murcko_scaffold"].fillna("")
        named = mur[mur != ""]
        top = named.value_counts().head(5)
        return {
            "cohort": label,
            "n_molecules": int(len(ids)),
            "n_in_descriptor_table": int(len(sub)),
            "n_parse_fail": int((sub["parse_ok"] != True).sum()),  # noqa: E712
            "n_murcko": int(named.nunique()),
            "n_generic_murcko": int(gen[gen != ""].nunique()),
            "n_empty_murcko": int((mur == "").sum()),
            "top_murcko": " || ".join(f"{a}:{int(b)}" for a, b in top.items()),
        }

    counts = pd.read_csv(OUT / "SCAFFOLD_COUNTS_BY_COHORT.csv")
    counts = counts[~counts["cohort"].isin(["stage3b_pass_13", "comparison_set"])]
    extra = pd.DataFrame([
        stats(new13, "stage3b_pass_13"),
        stats(comparable, "comparison_set"),
        stats(clean2, "hist21_chemistry_clean_not_frozen"),
    ])
    counts = pd.concat([counts, extra], ignore_index=True)
    counts.to_csv(OUT / "SCAFFOLD_COUNTS_BY_COHORT.csv", index=False)

    m33, m13, m156, m21 = map(murcko_set, (
        frozen33, new13,
        set(master.loc[master["in_deposited_a2_docking_set"] == True, "repurposing_id"].astype(str)),  # noqa: E712
        set(chem["ligand_id"].astype(str)),
    ))
    mclean = murcko_set(clean2)
    contrib = pd.read_csv(OUT / "NEW_SCAFFOLD_CONTRIBUTION.csv")
    contrib = contrib[~contrib["contrast"].isin(["new13_vs_frozen33", "new13_vs_docked156"])]
    add = pd.DataFrame([
        {"contrast": "new13_vs_frozen33", "n_query_scaffolds": len(m13), "n_ref_scaffolds": len(m33),
         "n_shared": len(m13 & m33), "n_new_vs_ref": len(m13 - m33)},
        {"contrast": "new13_vs_docked156", "n_query_scaffolds": len(m13), "n_ref_scaffolds": len(m156),
         "n_shared": len(m13 & m156), "n_new_vs_ref": len(m13 - m156)},
        {"contrast": "hist21_chem_clean2_vs_frozen33", "n_query_scaffolds": len(mclean), "n_ref_scaffolds": len(m33),
         "n_shared": len(mclean & m33), "n_new_vs_ref": len(mclean - m33)},
    ])
    contrib = pd.concat([add, contrib], ignore_index=True)
    contrib.to_csv(OUT / "NEW_SCAFFOLD_CONTRIBUTION.csv", index=False)

    # fix audit rows that used the wrong status-key guess
    audit = pd.read_csv(OUT / "ARCHIVE_BASELINE_AUDIT.csv")
    drop_items = {
        "reference7_runs", "stage3b_chem_pass_13", "final_status_stage3b_pass_like",
        "final_status::stage3b_structural_and_chemistry_pass",
        "final_status::stage3b_structural_pass_chemistry_fail",
        "final_status::stage3b_not_structural_pass",
        "final_status::stage3a_excluded_not_carboxylate",
    }
    audit = audit[~audit["item"].isin(drop_items)]
    extra_audit = pd.DataFrame([
        {"item": "reference7_runs", "observed": 7, "claimed": 7, "consistent": True,
         "source": "ALL_DOCKING_RUNS protocol_id P_REFERENCE7_REPRODUCTION"},
        {"item": "final_status::stage3b_new_structural_and_chemistry_pass", "observed": 13, "claimed": 13, "consistent": True,
         "source": "master.final_status"},
        {"item": "final_status::stage3b_new_structural_pass_chemistry_fail", "observed": 8, "claimed": 8, "consistent": True,
         "source": "master.final_status"},
        {"item": "final_status::stage3b_new_not_structural_pass", "observed": 11, "claimed": 11, "consistent": True,
         "source": "master.final_status"},
        {"item": "final_status::stage3a_excluded_rule_not_applicable", "observed": 7, "claimed": 7, "consistent": True,
         "source": "master.final_status"},
        {"item": "archive_text_says_four_scripts", "observed": 5, "claimed": 4, "consistent": False,
         "source": "FINAL_ARCHIVE_ACCEPTANCE.md says 四个生成脚本; directory has 5 python files. Text not edited."},
    ])
    # the original archive_generator_scripts row already records 5 vs 4
    audit = pd.concat([audit, extra_audit], ignore_index=True)
    audit.to_csv(OUT / "ARCHIVE_BASELINE_AUDIT.csv", index=False)

    method = pd.read_csv(OUT / "AD_METHOD_AUDIT.csv")
    method.loc[method["item"] == "v1_vs_recomputed_max_abs_diff_after_round3", "value"] = (
        "not computable: nlrp3_ml_scores_clinical_all.max_tanimoto_nlrp3_train and the archive historical column are entirely empty (8319/8319 missing)"
    )
    method.loc[method["item"] == "v1_pool_rule", "value"] = (
        "Pool membership was p_active_nlrp3 >= 0.5. The Tanimoto column was created but never filled, so v1 had no applied applicability-domain threshold."
    )
    method.to_csv(OUT / "AD_METHOD_AUDIT.csv", index=False)

    # stage3b 13 similarity row
    nn = pd.read_csv(OUT / "AD_CANDIDATE_NEAREST_NEIGHBORS.csv")
    sim = pd.read_csv(OUT / "AD_TRAIN_CLINICAL_SIMILARITY.csv")
    sim = sim[sim["population"] != "stage3b_13_max_tanimoto"]
    vals = nn.loc[nn["ligand_id"].isin(new13), "max_tanimoto"]
    row = {"population": "stage3b_13_max_tanimoto", "metric": "max_tanimoto_to_train"}
    row.update(quantiles(vals))
    finite = vals.to_numpy(dtype=float)
    row["n_ge_v2_ad_threshold"] = int(np.nansum(finite >= 0.5672473867595819))
    row["ad_threshold_used"] = 0.5672473867595819
    row["note"] = "Retrospective v2 1-NN call. Not a v1 screening rule. This is the production 1-NN AD statistic."
    sim = pd.concat([sim, pd.DataFrame([row])], ignore_index=True)
    sim.to_csv(OUT / "AD_TRAIN_CLINICAL_SIMILARITY.csv", index=False)
    return new13, new8, new11, excl7, clean2, frozen33, comparable


def candidate_tables(master, desc, new13, new8, new11, excl7, clean2, frozen33):
    chem = pd.read_csv(OUT / "HISTORICAL21_CHEMISTRY.csv")
    hist21 = set(chem["ligand_id"].astype(str))
    primary = set(master.loc[master["in_frozen_primary_12"] == True, "repurposing_id"].astype(str))  # noqa: E712
    backup = set(master.loc[master["in_frozen_backup_21"] == True, "repurposing_id"].astype(str))  # noqa: E712
    docked = set(master.loc[master["in_deposited_a2_docking_set"] == True, "repurposing_id"].astype(str))  # noqa: E712
    stage3b = set(master.loc[master["in_stage3b_new32_docked"] == True, "repurposing_id"].astype(str))  # noqa: E712

    # structural seed counts from historical tables and stage3b
    hist_counts = {lid: {"urat1": 0, "nlrp3": 0, "dual": 0} for lid in docked}
    u_score, n_score = {}, {}
    for seed in (42, 43, 44):
        df = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv")
        for _, r in df.iterrows():
            lid = str(r["ligand_id"])
            if lid not in hist_counts:
                continue
            if bool(r["keep_urat1_acid"]):
                hist_counts[lid]["urat1"] += 1
            if bool(r["keep_nlrp3_structural"]):
                hist_counts[lid]["nlrp3"] += 1
            if bool(r["keep_dual_acid_structural"]):
                hist_counts[lid]["dual"] += 1
            u_score.setdefault(lid, []).append(r["u_CNNscore"])
            n_score.setdefault(lid, []).append(r["n_CNNscore"])
    b_counts = {lid: {"urat1": 0, "nlrp3": 0, "dual": 0} for lid in stage3b}
    bdf = pd.read_csv(ROOT / "results/nlrp3_model_v2/stage3b_local/analysis/new32_per_seed.csv")
    for _, r in bdf.iterrows():
        lid = str(r["ligand_id"])
        if bool(r["urat1_pass"]):
            b_counts[lid]["urat1"] += 1
        if bool(r["nlrp3_structural_pass"]):
            b_counts[lid]["nlrp3"] += 1
        if bool(r["dual_seed_pass"]):
            b_counts[lid]["dual"] += 1
        u_score.setdefault(lid, []).append(r["urat1_CNNscore"])
        n_score.setdefault(lid, []).append(r["nlrp3_CNNscore"])

    rows = []
    universe = docked | stage3b | excl7
    m = master.set_index(master["repurposing_id"].astype(str))
    chem_i = chem.set_index(chem["ligand_id"].astype(str))
    for lid in sorted(universe):
        if lid in primary:
            group = "1_frozen_primary_12"
        elif lid in backup:
            group = "2_frozen_backup_21"
        elif lid in new13:
            group = "3_stage3b_new_13_not_historical_nomination"
        elif lid in clean2:
            group = "4_historical_structural_chemistry_clean_not_in_seed42_eligible"
        elif lid == "REP_07907":
            group = "5_structural_control_not_a_candidate"
        elif lid in hist21 or lid in new8:
            group = "5_structural_pass_chemistry_not_pass"
        elif lid in new11 or (lid in docked and hist_counts.get(lid, {}).get("dual", 0) < 2 and lid not in frozen33):
            group = "6_structural_not_pass"
        elif lid in excl7:
            group = "7_rule_not_applicable"
        else:
            group = "8_insufficient_or_unclassified"
        counts = b_counts.get(lid, hist_counts.get(lid, {"urat1": None, "nlrp3": None, "dual": None}))
        alerts = ""
        if lid in chem_i.index:
            alerts = str(chem_i.loc[lid, "alert_descriptions"] or "")
            if bool(chem_i.loc[lid, "beta_lactam"]):
                alerts = (alerts + " | beta_lactam").strip(" |")
            if bool(chem_i.loc[lid, "soft_name_excluded"]):
                alerts = (alerts + " | soft_name_exclusion").strip(" |")
        guidance = {
            "1_frozen_primary_12": "A_historical_primary",
            "2_frozen_backup_21": "B_historical_backup",
            "3_stage3b_new_13_not_historical_nomination": "C_analysis_only_new_structural_chemistry",
            "4_historical_structural_chemistry_clean_not_in_seed42_eligible": "D_analysis_only_not_renominated",
            "5_structural_pass_chemistry_not_pass": "E_not_for_nomination",
            "6_structural_not_pass": "F_not_for_nomination",
            "7_rule_not_applicable": "G_rule_not_applicable",
        }.get(group, "H_unclassified")
        rows.append({
            "ligand_id": lid,
            "name": m.loc[lid, "name"] if lid in m.index else None,
            "group": group,
            "historical_freeze_list": (
                "primary" if lid in primary else "backup" if lid in backup else "not_frozen"
            ),
            "analysis_guidance_tier": guidance,
            "analysis_guidance_is_not_a_nomination": True,
            "nlrp3_p_v1": m.loc[lid, "nlrp3_p_active_historical_v1"] if lid in m.index else None,
            "urat1_cnnscore_mean": float(np.nanmean(u_score.get(lid, [np.nan]))),
            "nlrp3_cnnscore_mean": float(np.nanmean(n_score.get(lid, [np.nan]))),
            "urat1_a2_pass_seeds": counts["urat1"],
            "nlrp3_structural_pass_seeds": counts["nlrp3"],
            "dual_structural_pass_seeds": counts["dual"],
            "chemistry_note": alerts,
            "max_phase": m.loc[lid, "max_phase"] if lid in m.index else None,
            "murcko_scaffold": desc.set_index("repurposing_id").loc[lid, "murcko_scaffold"] if lid in set(desc["repurposing_id"]) else None,
        })
    ev = pd.DataFrame(rows)
    ev["rank_model_score_v1_desc"] = ev["nlrp3_p_v1"].rank(ascending=False, method="min")
    ev["rank_urat1_dock_cnnscore_desc"] = ev["urat1_cnnscore_mean"].rank(ascending=False, method="min")
    ev["rank_nlrp3_dock_cnnscore_desc"] = ev["nlrp3_cnnscore_mean"].rank(ascending=False, method="min")
    ev["dock_score_rank_is_not_predicted_activity"] = True
    ev.to_csv(OUT / "CANDIDATE_EVIDENCE_MASTER.csv", index=False)
    groups = ev.groupby("group").size().rename("n").reset_index()
    groups.to_csv(OUT / "CANDIDATE_GROUPS_FINAL.csv", index=False)
    return ev


def stage_e(master, desc, ev):
    print("stage E", flush=True)
    # long metrics
    frames = []
    for seed in (42, 43, 44):
        u = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_pose_metrics_a2_seed{seed}.csv" if seed == 42
                        else ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_pose_metrics_seed{seed}.csv")
        u = u[u["pose_selection_rule"].astype(str).str.lower().eq("a2") | u["pose_selection_rule"].astype(str).str.contains("a2")]
        # seed42 file may contain only a2; if both rules exist keep a2
        if "pose_selection_rule" in u.columns and u["pose_selection_rule"].nunique() > 1:
            u = u[u["pose_selection_rule"].astype(str).str.lower().eq("a2")]
        n = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_structural_metrics_seed{seed}.csv")
        u = u[u["target"].astype(str).str.contains("urat1")] if "target" in u.columns else u
        keep = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv")
        merged = keep.merge(
            u[["ligand_id", "selected_mode", "n_geometry_pass_modes", "centroid_to_crystal_lig_A", "acid_arg477_min_A", "CNNscore"]],
            on="ligand_id", how="left", suffixes=("", "_pose"),
        )
        merged["seed"] = seed
        merged["cohort_source"] = "historical_156"
        frames.append(merged)
    hist = pd.concat(frames, ignore_index=True)
    anchor_rows = []
    for _, r in hist.iterrows():
        acid = pd.to_numeric(r.get("acid_arg477_min_A"), errors="coerce")
        cent = pd.to_numeric(r.get("centroid_to_crystal_lig_A"), errors="coerce")
        anchor_rows.append({
            "ligand_id": r["ligand_id"],
            "seed": int(r["seed"]),
            "source": "historical_156",
            "selected_mode": r.get("selected_mode"),
            "n_geometry_pass_modes": r.get("n_geometry_pass_modes"),
            "acid_arg477_min_A": acid,
            "centroid_A": cent,
            "CNNscore": r.get("u_CNNscore"),
            "CNNaffinity": r.get("u_CNNaffinity"),
            "urat1_pass": bool(r["keep_urat1_acid"]) if pd.notna(r["keep_urat1_acid"]) else None,
            "margin_arg_A": (ARG_THR - acid) if pd.notna(acid) else None,
            "margin_centroid_A": (CENT_THR - cent) if pd.notna(cent) else None,
            "distance_pass_is_not_a_salt_bridge": True,
        })
    bdf = pd.read_csv(ROOT / "results/nlrp3_model_v2/stage3b_local/analysis/new32_per_seed.csv")
    for _, r in bdf.iterrows():
        acid = pd.to_numeric(r["urat1_acid_arg477_min_A"], errors="coerce")
        cent = pd.to_numeric(r["urat1_centroid_A"], errors="coerce")
        anchor_rows.append({
            "ligand_id": r["ligand_id"],
            "seed": int(r["seed"]),
            "source": "stage3b_new32",
            "selected_mode": r["urat1_selected_mode"],
            "n_geometry_pass_modes": r["urat1_n_geometry_pass_modes"],
            "acid_arg477_min_A": acid,
            "centroid_A": cent,
            "CNNscore": r["urat1_CNNscore"],
            "CNNaffinity": r["urat1_CNNaffinity"],
            "urat1_pass": bool(r["urat1_pass"]),
            "margin_arg_A": (ARG_THR - acid) if pd.notna(acid) else None,
            "margin_centroid_A": (CENT_THR - cent) if pd.notna(cent) else None,
            "distance_pass_is_not_a_salt_bridge": True,
        })
    anchor = pd.DataFrame(anchor_rows)
    ncarb = desc.set_index("repurposing_id")["n_carboxylate"]
    anchor["n_carboxylate_parent"] = anchor["ligand_id"].map(ncarb)
    anchor.to_csv(OUT / "URAT1_ANCHOR_FEATURES.csv", index=False)

    # NLRP3 interaction matrix from historical per-residue columns
    residues = ["ALA227", "ALA228", "ARG351", "MET408", "TYR443", "PHE575", "ARG578"]
    mat = []
    for seed in (42, 43, 44):
        n = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_structural_metrics_seed{seed}.csv")
        for _, r in n.iterrows():
            for res in residues:
                mat.append({
                    "ligand_id": r["ligand_id"],
                    "seed": seed,
                    "source": "historical_156",
                    "residue": res,
                    "distance_A": r.get(f"d_{res}"),
                    "contact_4_5A": r.get(f"contact_{res}"),
                    "CNNscore": r["CNNscore"],
                    "overlap": r["pocket_overlap_frac"],
                    "ifp_jaccard": r["ifp_jaccard_vs_np3146"],
                    "n_key_contacts": r["n_key_contacts"],
                    "structural_pass": r["keep_nlrp3_structural"],
                    "geometric_contact_is_not_binding_energy": True,
                })
    for _, r in bdf.iterrows():
        mat.append({
            "ligand_id": r["ligand_id"],
            "seed": int(r["seed"]),
            "source": "stage3b_new32_aggregate_only",
            "residue": None,
            "distance_A": None,
            "contact_4_5A": None,
            "CNNscore": r["nlrp3_CNNscore"],
            "overlap": r["nlrp3_overlap"],
            "ifp_jaccard": r["nlrp3_ifp"],
            "n_key_contacts": r["nlrp3_n_key_contacts"],
            "structural_pass": r["nlrp3_structural_pass"],
            "clash": r["nlrp3_clash"],
            "geometric_contact_is_not_binding_energy": True,
            "per_residue_status": "missing_from_deposited_stage3b_table",
        })
    pd.DataFrame(mat).to_csv(OUT / "NLRP3_INTERACTION_MATRIX.csv", index=False)

    # scaffold vs outcome on historical 156
    dual_n = ev.set_index("ligand_id")["dual_structural_pass_seeds"]
    desc_i = desc.set_index("repurposing_id")
    docked_ids = set(ev.loc[ev["group"].isin([
        "1_frozen_primary_12", "2_frozen_backup_21", "4_historical_structural_chemistry_clean_not_in_seed42_eligible",
        "5_structural_pass_chemistry_not_pass", "6_structural_not_pass",
    ]) & ev["ligand_id"].isin(set(master.loc[master["in_deposited_a2_docking_set"] == True, "repurposing_id"].astype(str))), "ligand_id"])  # noqa: E712
    # simpler: all historical docked
    docked_ids = set(master.loc[master["in_deposited_a2_docking_set"] == True, "repurposing_id"].astype(str))  # noqa: E712
    comp_rows = []
    for lid in docked_ids:
        if lid not in desc_i.index:
            continue
        comp_rows.append({
            "ligand_id": lid,
            "murcko": desc_i.loc[lid, "murcko_scaffold"],
            "mw": desc_i.loc[lid, "mw"],
            "n_carboxylate": desc_i.loc[lid, "n_carboxylate"],
            "aromatic_rings": desc_i.loc[lid, "aromatic_rings"],
            "hetero_rings": desc_i.loc[lid, "hetero_rings"],
            "dual_seeds": dual_n.get(lid, np.nan),
            "structural_ge2": bool(dual_n.get(lid, 0) >= 2) if pd.notna(dual_n.get(lid, np.nan)) else None,
        })
    comp = pd.DataFrame(comp_rows)
    # same scaffold, mixed outcome
    mixed = []
    for scaff, g in comp.groupby("murcko"):
        if not scaff or g["structural_ge2"].nunique() < 2:
            continue
        mixed.append({
            "murcko": scaff,
            "n": int(len(g)),
            "n_pass": int(g["structural_ge2"].sum()),
            "n_fail": int((~g["structural_ge2"]).sum()),
            "mw_pass_median": g.loc[g["structural_ge2"], "mw"].median(),
            "mw_fail_median": g.loc[~g["structural_ge2"], "mw"].median(),
            "n_cooh_pass_median": g.loc[g["structural_ge2"], "n_carboxylate"].median(),
            "n_cooh_fail_median": g.loc[~g["structural_ge2"], "n_carboxylate"].median(),
        })
    pd.DataFrame(mixed).to_csv(OUT / "PASS_FAIL_STRUCTURE_COMPARISON.csv", index=False)

    # gate bias cohorts
    def cohort_summary(ids, label, population_note):
        sub = desc[desc["repurposing_id"].isin(ids) & (desc["parse_ok"] == True)]  # noqa: E712
        mur = sub["murcko_scaffold"].fillna("")
        named = mur[mur != ""]
        return {
            "cohort": label,
            "population_note": population_note,
            "n": int(len(ids)),
            "n_murcko": int(named.nunique()),
            "median_mw": float(sub["mw"].median()) if len(sub) else None,
            "median_n_carboxylate": float(sub["n_carboxylate"].median()) if len(sub) else None,
            "fraction_n_carboxylate_ge1": float((sub["n_carboxylate"] >= 1).mean()) if len(sub) else None,
            "median_aromatic_rings": float(sub["aromatic_rings"].median()) if len(sub) else None,
            "median_hetero_rings": float(sub["hetero_rings"].median()) if len(sub) else None,
            "median_clogp": float(sub["clogp"].median()) if len(sub) else None,
            "pass_rate_is_not_experimental_activity": True,
        }

    # historical gate layers
    u_ge2, n_ge2, d_ge2 = set(), set(), set()
    counts = {}
    for seed in (42, 43, 44):
        df = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv")
        for _, r in df.iterrows():
            lid = str(r["ligand_id"])
            rec = counts.setdefault(lid, {"u": 0, "n": 0, "d": 0})
            rec["u"] += int(bool(r["keep_urat1_acid"]))
            rec["n"] += int(bool(r["keep_nlrp3_structural"]))
            rec["d"] += int(bool(r["keep_dual_acid_structural"]))
    for lid, rec in counts.items():
        if rec["u"] >= 2:
            u_ge2.add(lid)
        if rec["n"] >= 2:
            n_ge2.add(lid)
        if rec["d"] >= 2:
            d_ge2.add(lid)
    frozen = set(master.loc[master["in_frozen_reportable_33"] == True, "repurposing_id"].astype(str))  # noqa: E712
    acid = set(master.loc[master["in_acid_pool_303"] == True, "repurposing_id"].astype(str))  # noqa: E712
    bias = [
        cohort_summary(acid, "pre_dock_acid_303", "Acid-equivalent pool before the 156 docking set. Not a stage3b population."),
        cohort_summary(docked_ids, "docked_156", "Historical docking set."),
        cohort_summary(u_ge2, "urat1_a2_ge2", "URAT1 A2 geometry on >=2/3 seeds. Accessibility, not activity."),
        cohort_summary(n_ge2, "nlrp3_structural_ge2", "NLRP3 structural gate on >=2/3 seeds."),
        cohort_summary(d_ge2, "dual_structural_ge2", "Same-seed dual structural pass on >=2/3 seeds."),
        cohort_summary(frozen, "frozen_reportable_33", "After historical chemistry and the freeze exclusions."),
    ]
    # stage3b separate, not appended as the next funnel step
    b_u, b_n, b_d = set(), set(), set()
    for lid, g in bdf.groupby("ligand_id"):
        if g["urat1_pass"].sum() >= 2:
            b_u.add(str(lid))
        if g["nlrp3_structural_pass"].sum() >= 2:
            b_n.add(str(lid))
        if g["dual_seed_pass"].sum() >= 2:
            b_d.add(str(lid))
    new32 = set(bdf["ligand_id"].astype(str))
    bias.extend([
        cohort_summary(new32, "stage3b_new32_separate_source", "v2-pool entrants. Not a continuation of the 156."),
        cohort_summary(b_u, "stage3b_urat1_a2_ge2", "Separate source."),
        cohort_summary(b_n, "stage3b_nlrp3_ge2", "Separate source."),
        cohort_summary(b_d, "stage3b_dual_ge2", "Separate source."),
        cohort_summary(set(ev.loc[ev["group"] == "3_stage3b_new_13_not_historical_nomination", "ligand_id"]),
                       "stage3b_chemistry_pass_13", "Separate source. Not added to the frozen 33."),
    ])
    pd.DataFrame(bias).to_csv(OUT / "GATE_SCAFFOLD_BIAS.csv", index=False)

    # margins and pre-specified tighter-threshold stability on selected poses
    margin = anchor.copy()
    # attach nlrp3 margins
    n_margin = []
    for seed in (42, 43, 44):
        n = pd.read_csv(ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_structural_metrics_seed{seed}.csv")
        for _, r in n.iterrows():
            n_margin.append({
                "ligand_id": r["ligand_id"], "seed": seed, "source": "historical_156",
                "nlrp3_CNNscore": r["CNNscore"],
                "nlrp3_centroid": r["centroid_to_crystal_lig_A"],
                "nlrp3_overlap": r["pocket_overlap_frac"],
                "nlrp3_ifp": r["ifp_jaccard_vs_np3146"],
                "nlrp3_n_key": r["n_key_contacts"],
                "margin_cnn": r["CNNscore"] - CNN_THR,
                "margin_nlrp3_centroid": CENT_THR - r["centroid_to_crystal_lig_A"],
                "margin_overlap": r["pocket_overlap_frac"] - OV_THR,
                "margin_ifp": r["ifp_jaccard_vs_np3146"] - IFP_THR,
                "margin_key": r["n_key_contacts"] - KEY_MIN,
            })
    for _, r in bdf.iterrows():
        n_margin.append({
            "ligand_id": r["ligand_id"], "seed": int(r["seed"]), "source": "stage3b_new32",
            "nlrp3_CNNscore": r["nlrp3_CNNscore"],
            "nlrp3_centroid": r["nlrp3_centroid_A"],
            "nlrp3_overlap": r["nlrp3_overlap"],
            "nlrp3_ifp": r["nlrp3_ifp"],
            "nlrp3_n_key": r["nlrp3_n_key_contacts"],
            "margin_cnn": r["nlrp3_CNNscore"] - CNN_THR,
            "margin_nlrp3_centroid": CENT_THR - r["nlrp3_centroid_A"] if pd.notna(r["nlrp3_centroid_A"]) else None,
            "margin_overlap": r["nlrp3_overlap"] - OV_THR if pd.notna(r["nlrp3_overlap"]) else None,
            "margin_ifp": r["nlrp3_ifp"] - IFP_THR if pd.notna(r["nlrp3_ifp"]) else None,
            "margin_key": r["nlrp3_n_key_contacts"] - KEY_MIN if pd.notna(r["nlrp3_n_key_contacts"]) else None,
        })
    nm = pd.DataFrame(n_margin)
    margin_out = anchor.merge(nm, on=["ligand_id", "seed", "source"], how="outer")
    margin_out.to_csv(OUT / "STRUCTURAL_GATE_MARGIN.csv", index=False)

    # tighter single-parameter checks: a passing selected pose that still clears the tighter cutoff
    # is stable. A failing selected pose is indeterminate because another pose was not re-selected.
    tight = [
        ("urat1_arg_7.2027", "margin_arg_A", 0.5),
        ("urat1_centroid_5.5", "margin_centroid_A", 0.5),
        ("nlrp3_overlap_0.55", "margin_overlap", 0.05),
        ("nlrp3_ifp_0.55", "margin_ifp", 0.05),
        ("nlrp3_key_contacts_6", "margin_key", 1),
        ("nlrp3_cnnscore_0.55", "margin_cnn", 0.05),
    ]
    sens_rows = []
    focus = set(ev.loc[ev["group"].isin([
        "1_frozen_primary_12", "2_frozen_backup_21", "3_stage3b_new_13_not_historical_nomination"
    ]), "ligand_id"].astype(str))
    subm = margin_out[margin_out["ligand_id"].isin(focus)]
    for name, col, delta in tight:
        if col not in subm.columns:
            continue
        # among seeds that pass the relevant official gate
        if name.startswith("urat1"):
            base = subm[subm["urat1_pass"] == True]  # noqa: E712
        else:
            official = pd.Series(True, index=subm.index)
            for need in ("margin_cnn", "margin_nlrp3_centroid", "margin_overlap", "margin_ifp", "margin_key"):
                if need in subm.columns:
                    official &= pd.to_numeric(subm[need], errors="coerce") >= 0
            base = subm[official]
        # structural_pass may only be on nm side; column exists after merge if nm has it? I didn't include structural_pass in nm.
        vals = pd.to_numeric(base[col], errors="coerce")
        sens_rows.append({
            "perturbation": name,
            "rule": "pre_specified_tighter_single_parameter_on_already_selected_pose",
            "n_official_pass_seed_rows": int(vals.notna().sum()),
            "n_selected_pose_still_clears": int((vals >= delta).sum()),
            "n_indeterminate_selected_pose_would_fail": int((vals < delta).sum()),
            "not_a_replacement_nomination_rule": True,
            "limitation": "Loosening cannot be counted from the selected pose. Tightening failure of the selected pose does not prove every pose would fail.",
        })
    pd.DataFrame(sens_rows).to_csv(OUT / "GATE_SENSITIVITY.csv", index=False)

    # scaffold compatibility summary
    scaf = comp.groupby("murcko").agg(
        n=("ligand_id", "size"),
        n_structural_ge2=("structural_ge2", "sum"),
        median_mw=("mw", "median"),
        median_n_carboxylate=("n_carboxylate", "median"),
    ).reset_index()
    scaf.to_csv(OUT / "SCAFFOLD_POSE_COMPATIBILITY.csv", index=False)

    # multi-carboxylate pose parse for focus molecules
    multi = [lid for lid in focus if lid in desc_i.index and pd.notna(desc_i.loc[lid, "n_carboxylate"]) and desc_i.loc[lid, "n_carboxylate"] >= 2]
    arg = json.loads((ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs/arg477_coords.json").read_text())["atoms"]
    nitro = [arg[k] for k in ("NE", "NH1", "NH2")]
    sys.path.insert(0, str(ROOT / "scripts"))
    from parse_c1_sdf_readouts import load_poses  # noqa: E402

    def sdf_path(lid, seed):
        cands = [
            ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual/urat1_9dkb/seed{seed}/{lid}_out.sdf",
            ROOT / f"data/campaigns/c1/07_clinical_dock/acid_dual_a2/urat1_9dkb/seed{seed}/{lid}_out.sdf",
            ROOT / f"results/nlrp3_model_v2/stage3b_local/new32/poses/urat1_9dkb/seed{seed}/{lid}_out.sdf",
        ]
        # broader search once if needed
        for p in cands:
            if p.exists():
                return p
        return None

    # discover stage3b pose root if the guessed path fails
    if multi and sdf_path(multi[0], 42) is None:
        hits = list(ROOT.glob(f"results/**/{multi[0]}_out.sdf"))
        print("sdf search", hits[:5], flush=True)

    carb_rows = []
    pats = [Chem.MolFromSmarts("[CX3](=O)[O-]"), Chem.MolFromSmarts("[CX3](=O)[OH]")]
    for lid in multi:
        for seed in (42, 43, 44):
            path = sdf_path(lid, seed)
            if path is None:
                carb_rows.append({"ligand_id": lid, "seed": seed, "status": "sdf_missing"})
                continue
            poses = load_poses(path)
            for pi, mol in enumerate(poses, start=1):
                groups = {}
                for pat in pats:
                    for match in mol.GetSubstructMatches(pat):
                        carbon = match[0]
                        groups.setdefault(carbon, set()).update(match[1:])
                conf = mol.GetConformer()
                for carbon, oxys in groups.items():
                    best = None
                    for idx in oxys:
                        p = conf.GetAtomPosition(int(idx))
                        coord = (p.x, p.y, p.z)
                        dmin = min(dist(coord, n) for n in nitro)
                        best = dmin if best is None else min(best, dmin)
                    carb_rows.append({
                        "ligand_id": lid,
                        "seed": seed,
                        "pose_ordinal": pi,
                        "carboxyl_carbon_idx": int(carbon),
                        "min_O_to_Arg_A": best,
                        "passes_7_7027": best is not None and best <= ARG_THR,
                        "status": "parsed",
                    })
    carb = pd.DataFrame(carb_rows)
    carb.to_csv(SUP / "MULTI_CARBOXYLATE_POSE_DISTANCES.csv", index=False)
    return margin_out, ev


def stage_f():
    print("stage F", flush=True)
    scores = pd.read_csv(ROOT / "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/scores/mol_protocol_scores.csv")
    defs = pd.DataFrame([
        {"protocol": "P0", "column": "P0_CNNscore", "name": "GNINA CNNscore", "higher_is_better": True,
         "pose_source": "GNINA search", "score_source": "CNNscore of the searched poses",
         "cnn_guided_search": False, "note": "cnn_scoring was rescore in the production protocol; this column is not evidence that CNN guided the search."},
        {"protocol": "P1", "column": "P1_vina_affinity", "name": "Vina affinity", "higher_is_better": False,
         "pose_source": "Vina search", "score_source": "Vina affinity", "cnn_guided_search": False, "note": "More negative is better."},
        {"protocol": "P2", "column": "P2_CNNaffinity", "name": "GNINA CNNaffinity", "higher_is_better": True,
         "pose_source": "GNINA search", "score_source": "CNNaffinity rescore", "cnn_guided_search": False,
         "note": "Same GNINA pose family as P0 and P3. Not used as the final nomination score."},
        {"protocol": "P3", "column": "P3_gnina_affinity", "name": "GNINA minimizedAffinity", "higher_is_better": False,
         "pose_source": "GNINA search", "score_source": "minimizedAffinity", "cnn_guided_search": False, "note": "More negative is better."},
        {"protocol": "P4", "column": "P4_RTM_vina", "name": "RTM on Vina poses", "higher_is_better": True,
         "pose_source": "Vina search", "score_source": "RTM rescore", "cnn_guided_search": False, "note": "Different poses from P2/P5."},
        {"protocol": "P5", "column": "P5_RTM_gnina", "name": "RTM on GNINA poses", "higher_is_better": True,
         "pose_source": "GNINA search", "score_source": "RTM rescore", "cnn_guided_search": False,
         "note": "Same pose source as P2, different score. Comparison with P2 is same-pose different-score."},
    ])
    defs["top_k_rule"] = "k = max(1, floor(fraction * n_molecules_in_the_evaluated_set))"
    defs["tie_rule"] = "stable mergesort; equal scores keep the original file order"
    defs["bootstrap"] = "stratified by label, paired indices, 2000 resamples; seeds 17,29,42,71,101 are RNG sensitivity, not five experiments"
    defs.to_csv(OUT / "P0_P5_METHOD_DEFINITIONS.csv", index=False)

    protocols = defs[["protocol", "column", "higher_is_better"]].to_dict("records")
    coverage_rows = []
    metric_rows = []
    dist_rows = []
    pair_rows = []

    def population(df, flag):
        sub = df[df[flag].astype(int) == 1].copy()
        y = (sub["role"].astype(str).str.lower() == "active").to_numpy().astype(int)
        return sub, y

    for flag, label in (("in_true", "true_decoy_benchmark"), ("in_random", "random_decoy_benchmark")):
        sub, y = population(scores, flag)
        own = {}
        for p in protocols:
            s = pd.to_numeric(sub[p["column"]], errors="coerce").to_numpy()
            mask = np.isfinite(s)
            own[p["protocol"]] = mask
            yy, ss = y[mask], s[mask]
            ef1, h1, k1 = enrichment(yy, ss, 0.01, p["higher_is_better"])
            ef5, h5, k5 = enrichment(yy, ss, 0.05, p["higher_is_better"])
            metric_rows.append({
                "benchmark": label,
                "population": "method_specific",
                "protocol": p["protocol"],
                "n": int(mask.sum()),
                "n_active": int(yy.sum()),
                "n_decoy": int((yy == 0).sum()),
                "coverage": float(mask.mean()),
                "AUROC": auroc(yy, ss),
                "PR_AUC": pr_auc(yy, ss),
                "EF1": ef1, "EF1_hits": h1, "EF1_k": k1,
                "EF5": ef5, "EF5_hits": h5, "EF5_k": k5,
                "recall_at_50": recall_at_k(yy, ss, 50, p["higher_is_better"]),
                "recall_at_100": recall_at_k(yy, ss, 100, p["higher_is_better"]),
                "higher_is_better": p["higher_is_better"],
                "decoy_label": "property-matched computational decoy, not an experimental inactive" if label.startswith("true") else "random computational decoy, not an experimental inactive",
            })
            for role_name, role_val in (("active", 1), ("decoy", 0)):
                vals = ss[yy == role_val]
                if len(vals) == 0:
                    continue
                q1, q2, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
                dist_rows.append({
                    "benchmark": label, "population": "method_specific", "protocol": p["protocol"],
                    "role": role_name, "n": int(len(vals)), "median": float(q2), "q25": float(q1), "q75": float(q3),
                })
            coverage_rows.append({
                "benchmark": label, "protocol": p["protocol"], "population": "method_specific",
                "n_benchmark": int(len(sub)), "n_scored": int(mask.sum()), "coverage": float(mask.mean()),
            })
        both = own["P2"] & own["P5"]
        yb = y[both]
        s2 = pd.to_numeric(sub["P2_CNNaffinity"], errors="coerce").to_numpy()[both]
        s5 = pd.to_numeric(sub["P5_RTM_gnina"], errors="coerce").to_numpy()[both]
        for pid, ss, hib in (("P2", s2, True), ("P5", s5, True)):
            ef1, h1, k1 = enrichment(yb, ss, 0.01, hib)
            ef5, h5, k5 = enrichment(yb, ss, 0.05, hib)
            metric_rows.append({
                "benchmark": label, "population": "P2_P5_common", "protocol": pid,
                "n": int(both.sum()), "n_active": int(yb.sum()), "n_decoy": int((yb == 0).sum()),
                "coverage": float(both.mean()),
                "AUROC": auroc(yb, ss), "PR_AUC": pr_auc(yb, ss),
                "EF1": ef1, "EF1_hits": h1, "EF1_k": k1,
                "EF5": ef5, "EF5_hits": h5, "EF5_k": k5,
                "recall_at_50": recall_at_k(yb, ss, 50, hib),
                "recall_at_100": recall_at_k(yb, ss, 100, hib),
                "higher_is_better": True,
                "decoy_label": "computational decoy, not an experimental inactive",
            })
        coverage_rows.append({
            "benchmark": label, "protocol": "P2_and_P5", "population": "common",
            "n_benchmark": int(len(sub)), "n_scored": int(both.sum()), "coverage": float(both.mean()),
            "n_active": int(yb.sum()), "n_decoy": int((yb == 0).sum()),
        })
        point = {
            "AUROC": pr_auc and auroc(yb, s5) - auroc(yb, s2),
            "PR_AUC": pr_auc(yb, s5) - pr_auc(yb, s2),
            "EF1": enrichment(yb, s5, 0.01, True)[0] - enrichment(yb, s2, 0.01, True)[0],
            "EF5": enrichment(yb, s5, 0.05, True)[0] - enrichment(yb, s2, 0.05, True)[0],
        }
        pos = np.flatnonzero(yb == 1)
        neg = np.flatnonzero(yb == 0)
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            draws = {k: [] for k in point}
            for _ in range(2000):
                idx = np.concatenate([
                    rng.choice(pos, size=len(pos), replace=True),
                    rng.choice(neg, size=len(neg), replace=True),
                ])
                y_i = yb[idx]
                a2, a5 = s2[idx], s5[idx]
                draws["AUROC"].append(auroc(y_i, a5) - auroc(y_i, a2))
                draws["PR_AUC"].append(pr_auc(y_i, a5) - pr_auc(y_i, a2))
                draws["EF1"].append(enrichment(y_i, a5, 0.01, True)[0] - enrichment(y_i, a2, 0.01, True)[0])
                draws["EF5"].append(enrichment(y_i, a5, 0.05, True)[0] - enrichment(y_i, a2, 0.05, True)[0])
            row = {
                "benchmark": label,
                "contrast": "P5_minus_P2",
                "seed": seed,
                "n_boot": 2000,
                "n_common": int(len(yb)),
                "n_active": int(len(pos)),
                "n_decoy": int(len(neg)),
                "seeds_are_independent_experiments": False,
                "resample": "stratified by active versus decoy label; same indices applied to both scores",
            }
            for k, vals in draws.items():
                arr = np.asarray(vals, dtype=float)
                lo, hi = np.nanpercentile(arr, [2.5, 97.5])
                row[f"delta_{k}"] = point[k]
                row[f"delta_{k}_boot_mean"] = float(np.nanmean(arr))
                row[f"delta_{k}_ci95_low"] = float(lo)
                row[f"delta_{k}_ci95_high"] = float(hi)
                row[f"delta_{k}_ci_excludes_0"] = bool(lo > 0 or hi < 0)
            pair_rows.append(row)
            print(f"  bootstrap {label} seed {seed}", flush=True)
    pd.DataFrame(coverage_rows).to_csv(OUT / "SCORING_COVERAGE.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(OUT / "AUROC_PRAUC_EF_ALL.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(OUT / "P2_P5_PAIRED_BOOTSTRAP.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(OUT / "SCORE_DISTRIBUTIONS.csv", index=False)


def stage_g():
    print("stage G", flush=True)
    u = pd.read_csv(ROOT / "data/manuscript/tables/urat1_redocking_summary.csv")
    n = pd.read_csv(ROOT / "data/manuscript/tables/nlrp3_selfdocking_summary.csv")
    ref = pd.concat([u, n], ignore_index=True)
    ref["search_success_best_le_2A"] = ref["best_of_9_pose_rmsd_A"] <= 2.0
    ref["scoring_success_top1_le_2A"] = ref["top1_pose_rmsd_A"] <= 2.0
    ref["native_like_rank"] = np.where(
        ref["scoring_success_top1_le_2A"],
        "1_selected",
        np.where(ref["search_success_best_le_2A"], "within_9_but_not_selected", "not_in_9"),
    )
    ref["note"] = np.where(
        ref["ligand"].eq("NP3-146"),
        "NLRP3 self-dock: search and scoring both recovered a pose <=2 A.",
        "URAT1: a <=2 A pose was present for some seeds but CNNscore did not select it.",
    )
    ref7 = pd.read_csv(ROOT / "results/nlrp3_model_v2/stage3b_local/reference7/task_comparison.csv")
    ref["reference7_is_reproduction_not_seven_validations"] = True
    ref.to_csv(OUT / "REFERENCE_DOCKING_CONSOLIDATED.csv", index=False)
    ref7.to_csv(SUP / "REFERENCE7_TASK_COMPARISON_COPY.csv", index=False)

    per = pd.read_csv(ROOT / "data/campaigns/c1/05_metrics/acid_gate_retrospective_benchmark/acid_gate_benchmark_per_mol.csv")
    rows = []
    for rule, g in per.groupby(per["pose_selection_rule"].astype(str).str.lower()):
        # one row per molecule
        mol = g.drop_duplicates("ligand_id")
        y = mol["label"].astype(int)
        pred = mol["keep_urat1_acid"].astype(bool)
        tp = int(((pred) & (y == 1)).sum())
        fn = int(((~pred) & (y == 1)).sum())
        fp = int(((pred) & (y == 0)).sum())
        tn = int(((~pred) & (y == 0)).sum())
        sens = tp / (tp + fn) if (tp + fn) else np.nan
        spec = tn / (tn + fp) if (tn + fp) else np.nan
        lr = sens / (1 - spec) if spec not in (None, 1) and pd.notna(spec) and (1 - spec) > 0 else np.nan
        rows.append({
            "gate": rule,
            "n_molecules": int(len(mol)),
            "n_rows_before_dedup": int(len(g)),
            "n_active": int((y == 1).sum()),
            "n_decoy": int((y == 0).sum()),
            "tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "sensitivity": sens, "specificity": spec, "lr_plus": lr,
            "coverage": 1.0,
            "interpretation": (
                "Geometry accessibility. LR+ near 1 is not URAT1 activity enrichment."
                if rule == "a2" else
                "A1 selects the top CNNscore pose and then applies geometry. Not the nomination rule."
            ),
        })
    pd.DataFrame(rows).to_csv(OUT / "A1_A2_METRIC_RECALCULATION.csv", index=False)


def stage_h():
    print("stage H", flush=True)
    viz = Path("/mnt/d/CADD paper exercise/NLRP3_URAT1/viz_pack_20260907/viz_pack")
    rows = []
    summary_path = viz / "data/summary.json"
    systems = ["Lesinurad_URAT1", "PF03882845_URAT1", "NP3-146_NLRP3", "PF03882845_NLRP3"]
    for name in systems:
        work = viz / "work" / name
        for kind, path in {
            "summary_json": summary_path,
            "dry_nc": work / "dry.nc",
            "dry_prmtop": work / "dry.prmtop",
            "last50_nc": work / "last50.nc",
            "mmgbsa_prmtop": viz / "mmgbsa" / name / "com.prmtop",
        }.items():
            rows.append({
                "system": name,
                "artifact": kind,
                "path": str(path),
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else None,
                "full_solvated_trajectory": False,
                "in_git_repository": False,
            })
    rows.append({
        "system": "all",
        "artifact": "locked_plan",
        "path": str(ROOT / "config/md_protocol.yaml"),
        "exists": True,
        "status_in_plan": "planning_only_not_executed",
        "md_authorized": False,
        "results_available_flag_in_plan": False,
    })
    pd.DataFrame(rows).to_csv(OUT / "MD_FILE_AVAILABILITY.csv", index=False)

    summary = json.loads(summary_path.read_text())
    by = {x["name"]: x for x in summary}
    les = by["Lesinurad_URAT1"]
    gate_rows = [{
        "rule": "URAT1 control fails if lesinurad heavy-atom RMSD >5 A for >50% of the analysis window OR acid–Arg <4 A occupancy <30%",
        "source": "config/md_protocol.yaml control_failure.urat1",
        "observed_occ_salt_lt_4A": les["occ_salt_lt_4A"],
        "threshold": 0.30,
        "passes_prespecified_rule": bool(les["occ_salt_lt_4A"] >= 0.30),
        "rmsd_nofit_max": les["ligand_rmsd_nofit"]["max"],
        "rmsd_criterion_alone_not_tripped_by_max_below_5A": les["ligand_rmsd_nofit"]["max"] <= 5,
        "decision": "control_failed_on_salt_occupancy_do_not_interpret_PF_URAT1",
        "threshold_was_not_changed": True,
    }]
    np3 = by["NP3-146_NLRP3"]
    key_occ = np3.get("key_occ_lt_5A")
    gate_rows.append({
        "rule": "NLRP3 control fails if NP3-146 RMSD >5 A for >50% OR fewer than 3 of 7 key contacts are persistent",
        "source": "config/md_protocol.yaml control_failure.nlrp3",
        "ligand_rmsd_nofit_mean": np3["ligand_rmsd_nofit"]["mean"],
        "ligand_rmsd_nofit_max": np3["ligand_rmsd_nofit"]["max"],
        "n_key_contacts_mean": np3["n_key_contacts"]["mean"],
        "n_key_contacts_min": np3["n_key_contacts"]["min"],
        "key_occ_lt_5A": json.dumps(key_occ),
        "fraction_frames_rmsd_gt_5A": "not_recomputed_full_solvated_trajectory_absent",
        "decision": "summary does not show fewer than 3 of 7 contacts; RMSD fraction above 5 A was not recomputed",
    })
    pd.DataFrame(gate_rows).to_csv(OUT / "MD_PRESPECIFIED_GATE_AUDIT.csv", index=False)

    evid = []
    for name, blob in by.items():
        evid.append({
            "system": name,
            "ns": blob.get("ns"),
            "n_frames": blob.get("n_frames"),
            "replicas_in_this_summary": 1,
            "planned_replicas": 3,
            "grade_in_summary": blob.get("grade"),
            "protein_ca_rmsd_mean": blob["protein_ca_rmsd"]["mean"],
            "ligand_rmsd_nofit_mean": blob["ligand_rmsd_nofit"]["mean"],
            "ligand_rmsd_fit_mean": blob["ligand_rmsd_fit"]["mean"],
            "lig_pocket_com_mean": blob["lig_pocket_com"]["mean"],
            "occ_salt_lt_4A": blob.get("occ_salt_lt_4A"),
            "occ_a1_lt_7.7A": blob.get("occ_a1_lt_7.7A"),
            "n_key_contacts_mean": None if not isinstance(blob.get("n_key_contacts"), dict) else blob["n_key_contacts"]["mean"],
            "time_series_recomputed": False,
            "reason": "Dry trajectories exist locally. Full solvated trajectories and production inputs sufficient to rerun the pre-registered analysis are not in the readable set. Numbers are the deposited summary only.",
        })
    pd.DataFrame(evid).to_csv(OUT / "PF_MD_EVIDENCE.csv", index=False)


def stage_i(master):
    print("stage I", flush=True)
    lit = pd.read_csv(ROOT / "data/benchmarks/literature_benchmarks.csv", low_memory=False)
    train = set(pd.read_csv(ROOT / "data/processed/nlrp3_records.csv", usecols=["canonical_smiles"])["canonical_smiles"].astype(str))
    from libchem import canonicalize
    lit["canon"] = lit["canonical_smiles"].map(lambda s: canonicalize(s) if isinstance(s, str) else None)
    lit["in_nlrp3_training_smiles"] = lit["canon"].isin(train)
    lit.to_csv(OUT / "EXTERNAL_REFERENCE_COMPOUNDS.csv", index=False)
    # eligibility: external AUROC only if independent pos and neg exist and are scored
    elig = []
    for gene, g in lit.groupby(lit["target_gene"].astype(str)):
        outside = g[g["in_train_project"].astype(str).str.lower().isin(["false", "0", "no"]) | (g["in_nlrp3_training_smiles"] == False)]  # noqa: E712
        n_pos = int((pd.to_numeric(g["pactivity"], errors="coerce") >= 6).sum()) if "pactivity" in g else None
        elig.append({
            "target_gene": gene,
            "n_rows": int(len(g)),
            "n_unique_smiles": int(g["canon"].nunique()),
            "n_flagged_in_train_project": int(g["in_train_project"].astype(str).str.lower().isin(["true", "1", "yes"]).sum()) if "in_train_project" in g else None,
            "n_positive_pactivity_ge6": n_pos,
            "n_explicit_experimental_negatives": 0,
            "external_auroc": "not_computed",
            "reason": "The deposited literature table does not contain an audited independent negative set. Positives that also built the protocol or the threshold are not an external test.",
        })
    pd.DataFrame(elig).to_csv(OUT / "EXTERNAL_DATA_ELIGIBILITY.csv", index=False)

    focus_groups = {
        "1_frozen_primary_12", "2_frozen_backup_21", "3_stage3b_new_13_not_historical_nomination",
        "4_historical_structural_chemistry_clean_not_in_seed42_eligible",
    }
    ev = pd.read_csv(OUT / "CANDIDATE_EVIDENCE_MASTER.csv")
    focus = ev[ev["group"].isin(focus_groups)][["ligand_id", "name", "group", "max_phase"]].copy()
    man = pd.read_csv(ROOT / "data/repurposing/repurposing_manifest.csv", usecols=["repurposing_id", "atc_code", "max_phase"])
    focus = focus.merge(man, left_on="ligand_id", right_on="repurposing_id", how="left", suffixes=("", "_manifest"))
    lit_names = set(lit["compound_name"].astype(str).str.upper())
    focus["in_literature_benchmark_table"] = focus["name"].astype(str).str.upper().isin(lit_names)
    focus["known_target_from_this_repository"] = np.where(
        focus["in_literature_benchmark_table"], "see EXTERNAL_REFERENCE_COMPOUNDS.csv", "not_recorded_in_deposited_benchmark_table"
    )
    focus["urate_effect"] = "not_established_by_this_analysis"
    focus["reported_anti_inflammatory_effect"] = "not_inferred"
    focus["direct_nlrp3_inhibitor"] = "not_established"
    focus["prodrug"] = "not_established_from_deposited_fields"
    focus["source_limit"] = "max_phase and ATC come from the library manifest. Anti-inflammatory or NLRP3 labels were not inferred from the disease area."
    focus.to_csv(OUT / "CANDIDATE_PHARMACOLOGY.csv", index=False)


def refresh_pairwise(desc):
    ev = pd.read_csv(OUT / "CANDIDATE_EVIDENCE_MASTER.csv")
    keep = ev[ev["group"].isin([
        "1_frozen_primary_12", "2_frozen_backup_21",
        "3_stage3b_new_13_not_historical_nomination",
        "4_historical_structural_chemistry_clean_not_in_seed42_eligible",
        "5_structural_pass_chemistry_not_pass",
        "5_structural_control_not_a_candidate",
    ])]
    desc_i = desc.set_index("repurposing_id")
    ids = [i for i in keep["ligand_id"].astype(str) if i in desc_i.index]
    fps = []
    ok = []
    for i in ids:
        smi = desc_i.loc[i, "largest_fragment_smiles"] if "largest_fragment_smiles" in desc_i.columns else desc_i.loc[i, "standardized_smiles"]
        fp = morgan_fp(smi) if isinstance(smi, str) else None
        if fp is not None:
            fps.append(fp)
            ok.append(i)
    rows = []
    names = ev.set_index("ligand_id")["name"].to_dict()
    for a_i, a in enumerate(ok):
        sims = DataStructs.BulkTanimotoSimilarity(fps[a_i], fps)
        for b_i in range(a_i + 1, len(ok)):
            rows.append({
                "id_a": a,
                "id_b": ok[b_i],
                "name_a": names.get(a),
                "name_b": names.get(ok[b_i]),
                "tanimoto_ecfp4": float(sims[b_i]),
                "same_murcko": desc_i.loc[a, "murcko_scaffold"] == desc_i.loc[ok[b_i], "murcko_scaffold"],
            })
    pd.DataFrame(rows).to_csv(OUT / "PAIRWISE_CANDIDATE_SIMILARITY.csv", index=False)


def main():
    master = pd.read_csv(ARCH / "07_RANKING_AND_SELECTION/ALL_COMPOUNDS_MASTER_INDEX.csv", low_memory=False)
    desc = pd.read_csv(SUP / "CLINICAL_8319_DESCRIPTORS.csv")
    new13, new8, new11, excl7, clean2, frozen33, comparable = patch_a_d(master, desc)
    ev = candidate_tables(master, desc, new13, new8, new11, excl7, clean2, frozen33)
    refresh_pairwise(desc)
    stage_e(master, desc, ev)
    stage_f()
    stage_g()
    stage_h()
    stage_i(master)
    print("done", len(new13), len(clean2), len(comparable))


if __name__ == "__main__":
    main()
