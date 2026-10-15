#!/usr/bin/env python3
"""Stage 3A: pre-docking scientific audit. Does not launch gnina or MD.

1. Scaffold-CV performance of the uncalibrated structure-only model, binned by
   nearest-neighbour Tanimoto to the outer-training molecules.
2. One chemistry and geometry-scope annotation for every molecule that passes
   the production model (p >= 0.5) and the existing soft-chemistry gate.
3. A docking list for the previously undocked molecules, with a reason for
   every exclusion. Names, expected activity and the PF trajectory are not inputs.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score

from build_c1_acid_clinical_pool import acid_flags
from build_c5_tier_assignment import BETA_LACTAM_SMARTS
from nlrp3_model_v2_validate import classification_report, load_records, molecule_table
from nlrp3_v2_lib import StructureOnlyModel
from utils_ml import featurize_smiles, scaffold_cv_indices
from importlib.machinery import SourceFileLoader

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "nlrp3_model_v2" / "stage3a"
SCORES = PROJECT / "results" / "nlrp3_model_v2" / "stage2" / "clinical_library_scores.csv"
SEED = 42
THRESHOLD = 0.5
# Fixed from the similarity scale, not from which bin has the highest AUROC.
# 0.567 is the existing applicability threshold and is not a bin edge that we move.
BIN_EDGES = (0.0, 0.30, 0.40, 0.50, 0.567, 1.01)
LACTAM = Chem.MolFromSmarts(BETA_LACTAM_SMARTS)
filters = SourceFileLoader(
    "chem_filters", str(PROJECT / "scripts" / "09_cheminformatics_filters.py")
).load_module()


def _fps(smiles: list[str]):
    out = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        out.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None)
    return out


def _max_tanimoto(query_idx: np.ndarray, ref_idx: np.ndarray, fps) -> np.ndarray:
    ref = [fps[i] for i in ref_idx if fps[i] is not None]
    out = np.full(len(query_idx), np.nan)
    for k, i in enumerate(query_idx):
        q = fps[i]
        if q is None or not ref:
            continue
        out[k] = max(DataStructs.BulkTanimotoSimilarity(q, ref))
    return out


def _loo_nn(fps) -> np.ndarray:
    usable = [f for f in fps if f is not None]
    out = np.full(len(fps), np.nan)
    for i, q in enumerate(fps):
        if q is None:
            continue
        sims = sorted(DataStructs.BulkTanimotoSimilarity(q, usable), reverse=True)
        out[i] = sims[1] if len(sims) > 1 else sims[0]
    return out


def _bin_label(value: float) -> str:
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        if lo <= value < hi or (hi == BIN_EDGES[-1] and value <= hi):
            return f"[{lo:.3f}, {hi:.3f})"
    return "unbinned"


def similarity_cv(mol: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    smiles = mol["canonical_smiles"].tolist()
    y = mol["active"].to_numpy(dtype=int)
    folds = scaffold_cv_indices(smiles, n_splits=5)
    fps = _fps(smiles)
    loo = _loo_nn(fps)
    ad_threshold = float(np.nanpercentile(loo, 5))
    rows = []
    for fold_i, (tr, te) in enumerate(folds):
        fit = mol.iloc[tr][["canonical_smiles"]].copy()
        fit["active"] = y[tr]
        model = StructureOnlyModel().fit(fit, seed=SEED + fold_i)
        p = model.predict_molecules(mol.iloc[te]["canonical_smiles"].tolist())
        nn = _max_tanimoto(te, tr, fps)
        for j, idx in enumerate(te):
            rows.append(
                {
                    "fold": fold_i,
                    "canonical_smiles": smiles[idx],
                    "y_true": int(y[idx]),
                    "p_oof": float(p[j]),
                    "nn_tanimoto_to_outer_train": float(nn[j]),
                    "similarity_bin": _bin_label(float(nn[j])),
                    "above_existing_ad_threshold": bool(nn[j] >= ad_threshold),
                }
            )
    pred = pd.DataFrame(rows)
    bin_rows = []
    groups = [("all_outer_test", pred)] + [
        (b, pred[pred["similarity_bin"] == b]) for b in pred["similarity_bin"].unique()
    ]
    below = pred[pred["nn_tanimoto_to_outer_train"] < 0.30]
    groups.append(("[0.000, 0.300) clinical-like", below))
    for name, g in groups:
        if len(g) == 0:
            continue
        metrics = classification_report(g["y_true"].to_numpy(), g["p_oof"].to_numpy(), THRESHOLD)
        metrics["bin"] = name
        metrics["nn_min"] = float(g["nn_tanimoto_to_outer_train"].min())
        metrics["nn_median"] = float(g["nn_tanimoto_to_outer_train"].median())
        metrics["nn_max"] = float(g["nn_tanimoto_to_outer_train"].max())
        bin_rows.append(metrics)
    summary = {
        "protocol": "5-fold Murcko GroupKFold; each fold trains on the full outer-train set; no calibration",
        "similarity": "max Morgan radius-2 Tanimoto (2048 bits) to molecules in that fold's outer-training set",
        "existing_ad_rule": "5th percentile of leave-one-out nearest-neighbour Tanimoto inside the 513 training molecules",
        "ad_threshold_recomputed": ad_threshold,
        "ad_threshold_adjusted": False,
        "training_loo_nn": {
            "p05": ad_threshold,
            "p25": float(np.nanpercentile(loo, 25)),
            "median": float(np.nanmedian(loo)),
            "p75": float(np.nanpercentile(loo, 75)),
        },
        "oof_nn": {
            "min": float(pred["nn_tanimoto_to_outer_train"].min()),
            "p05": float(pred["nn_tanimoto_to_outer_train"].quantile(0.05)),
            "median": float(pred["nn_tanimoto_to_outer_train"].median()),
            "frac_below_0_30": float((pred["nn_tanimoto_to_outer_train"] < 0.30).mean()),
            "frac_below_ad_threshold": float((pred["nn_tanimoto_to_outer_train"] < ad_threshold).mean()),
            "n_below_0_30": int((pred["nn_tanimoto_to_outer_train"] < 0.30).sum()),
        },
        "why_threshold_excludes_clinical_candidates": (
            "The 0.567 line is the left tail of similarity among training-set analogues, "
            "not an estimate of where predictive performance collapses. A congeneric series "
            "still has a high 5th percentile. Clinical molecules sit near 0.2 because they "
            "are not members of that series. The threshold was not lowered."
        ),
    }
    return pred, pd.DataFrame(bin_rows), summary


def _alert_row(mol) -> dict:
    catalogs = _alert_row.catalogs
    pains_fam = []
    descriptions = []
    for label in ("PAINS_A", "PAINS_B", "PAINS_C"):
        match = catalogs[label].GetFirstMatch(mol)
        if match is not None:
            pains_fam.append(label)
            descriptions.append(f"{label}:{match.GetDescription()}")
    brenk = catalogs["BRENK"].GetFirstMatch(mol)
    nih = catalogs["NIH"].GetFirstMatch(mol)
    if brenk is not None:
        descriptions.append(f"BRENK:{brenk.GetDescription()}")
    if nih is not None:
        descriptions.append(f"NIH:{nih.GetDescription()}")
    return {
        "pains_any": len(pains_fam) > 0,
        "pains_families": ";".join(pains_fam),
        "brenk": brenk is not None,
        "nih": nih is not None,
        "alert_clean_pains_brenk": len(pains_fam) == 0 and brenk is None,
        "aggregation_risk_heuristic": bool(filters.aggregation_heuristic(mol)),
        "alert_descriptions": " | ".join(descriptions),
    }


_alert_row.catalogs = filters.build_catalogs()


def _acid_class(flags: dict) -> str:
    parts = []
    if flags.get("has_carboxylate"):
        parts.append("carboxylate")
    if flags.get("has_tetrazole"):
        parts.append("tetrazole")
    if flags.get("has_acyl_sulfonamide"):
        parts.append("acyl_sulfonamide")
    return "+".join(parts) if parts else "none"


def annotate_panel(scores: pd.DataFrame) -> pd.DataFrame:
    panel = scores[scores["in_v2_pool"] & scores["pass_chemistry_soft"]].copy()
    rows = []
    for _, r in panel.iterrows():
        mol = Chem.MolFromSmiles(str(r["canonical_smiles"]))
        flags = acid_flags(str(r["canonical_smiles"]))
        flags.pop("mol", None)
        beta = bool(mol and LACTAM and mol.HasSubstructMatch(LACTAM))
        alerts = _alert_row(mol) if mol is not None else {}
        carbox = bool(flags.get("has_carboxylate"))
        rows.append(
            {
                "repurposing_id": r["repurposing_id"],
                "name": r["name"],
                "chembl_id": r["chembl_id"],
                "max_phase": r["max_phase"],
                "canonical_smiles": r["canonical_smiles"],
                "p_n513": r["p_n513"],
                "max_tanimoto_train": r["max_tanimoto_train"],
                "in_training_set": bool(r["in_training_set"]),
                "in_historical_pool": bool(r["in_historical_pool"]),
                "in_deposited_a2_docking_set": bool(r["in_deposited_a2_docking_set"]),
                "cohort": "deposited_118" if r["in_deposited_a2_docking_set"] else "new_39",
                "acid_class": _acid_class(flags),
                "has_carboxylate": carbox,
                "has_tetrazole": bool(flags.get("has_tetrazole")),
                "has_acyl_sulfonamide": bool(flags.get("has_acyl_sulfonamide")),
                "beta_lactam": beta,
                "urat1_geometry_scope": (
                    "in_scope_existing_carboxylate_arg477_rule"
                    if carbox
                    else "outside_current_method_noncarboxylate_acid"
                ),
                **alerts,
            }
        )
    out = pd.DataFrame(rows)
    # Historical geometry outcome, so a non-carboxylate failure is not reread as inactivity.
    keep = pd.read_csv(
        PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/acid_pose_metrics_a2_seed42.csv"
    )
    u = keep[keep["target"] == "urat1_9dkb"][
        ["ligand_id", "pose_selection_status", "has_carboxylate_oxy_in_pose", "keep_urat1_acid"]
    ]
    out = out.merge(u, left_on="repurposing_id", right_on="ligand_id", how="left")
    out = out.drop(columns=["ligand_id"])
    return out


def assign_docking(panel: pd.DataFrame) -> pd.DataFrame:
    """Decide the new-molecule docking list without using names or activity."""
    rows = []
    new = panel[panel["cohort"] == "new_39"]
    for _, r in new.iterrows():
        if r["in_training_set"]:
            decision, reason = (
                "exclude",
                "canonical SMILES is already in the 513-molecule training set, so the classifier probability is a resubstitution",
            )
        elif not r["has_carboxylate"]:
            decision, reason = (
                "exclude",
                "no carboxylate; the locked URAT1 pose rule only measures carboxylate-oxygen to Arg477 distance, so a geometry miss would not be evidence of absent URAT1 activity. No separate anchor rule is defined in this phase",
            )
        else:
            decision, reason = (
                "dock",
                "carboxylate present, so the existing Arg477 geometry rule applies; beta-lactam status is recorded for the later nomination step and does not remove the molecule from docking",
            )
        rows.append(
            {
                "repurposing_id": r["repurposing_id"],
                "name": r["name"],
                "acid_class": r["acid_class"],
                "beta_lactam": bool(r["beta_lactam"]),
                "p_n513": r["p_n513"],
                "max_tanimoto_train": r["max_tanimoto_train"],
                "alert_clean_pains_brenk": bool(r["alert_clean_pains_brenk"]),
                "docking_decision": decision,
                "reason": reason,
            }
        )
    return pd.DataFrame(rows).sort_values(["docking_decision", "repurposing_id"])


def protocol_record() -> dict:
    import yaml

    eng = yaml.safe_load((PROJECT / "config/docking_c1_cpu.yaml").read_text())
    camp = yaml.safe_load((PROJECT / "config/campaign_c1.yaml").read_text())
    receptors = {}
    for key in ("urat1_9dkb", "nlrp3_7alv"):
        rel = eng["targets"][key]["prepared_receptor"]
        path = PROJECT / rel
        receptors[key] = {
            "prepared_receptor": rel,
            "exists": path.exists(),
            "bytes": path.stat().st_size if path.exists() else None,
            "center": eng["targets"][key]["center"],
            "size": eng["targets"][key]["size"],
            "qc_center": json.loads((PROJECT / f"data/structures/prepared/{'9DKB' if 'urat1' in key else '7ALV'}_receptor.qc.json").read_text())["center"],
            "qc_size": json.loads((PROJECT / f"data/structures/prepared/{'9DKB' if 'urat1' in key else '7ALV'}_receptor.qc.json").read_text())["size"],
        }
    sdf_counts = {
        "urat1_seed42": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual/urat1_9dkb/seed42").glob("*.sdf"))),
        "nlrp3_seed42": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual/nlrp3_7alv/seed42").glob("*.sdf"))),
        "urat1_seed43": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/urat1_9dkb/seed43").glob("*.sdf"))),
        "urat1_seed44": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/urat1_9dkb/seed44").glob("*.sdf"))),
        "nlrp3_seed43": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_7alv/seed43").glob("*.sdf"))),
        "nlrp3_seed44": len(list((PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_7alv/seed44").glob("*.sdf"))),
        "clinical_logs_deposited": 0,
        "gnina_binary_in_repo": (PROJECT / "tools/gnina").exists(),
    }
    return {
        "engine_config": "config/docking_c1_cpu.yaml",
        "gnina": eng["gnina"],
        "documented_selfdock_binary": "gnina v1.3.2 master:f23dd2b Built Jul 8 2025, from data/campaigns/c1/02_selfdock/urat1_9dkb/seed42/lesinurad_out_stdout.txt",
        "clinical_sdf_openbabel_headers": {
            "seed42": "OpenBabel08262613303D",
            "seed43_and_44": "OpenBabel08272609403D",
        },
        "sdf_counts": sdf_counts,
        "missing_nlrp3_sdf_seeds_43_44": ["REP_07837"],
        "ligand_prep": {
            "script": "scripts/prepare_ligands_c1.py",
            "protonation": "Dimorphite-DL pH 7.4",
            "embed_seed": camp["ligand_prep"]["embed_seed"] if "ligand_prep" in camp else "0xC0FFEE",
            "prefer_deprotonated_carboxylate": True,
            "writer": "Meeko PDBQTWriterLegacy",
            "existing_pdbqt_for_chemistry_pass": 156,
            "existing_pdbqt_for_the_39": 0,
        },
        "receptors": receptors,
        "locked_gates_not_to_be_retuned": {
            "urat1_carboxylate_arg477_A": camp["acid_track"]["urat1_acid_arg477_le_A"],
            "pocket_centroid_A": 6.0,
            "nlrp3_cnn_min": 0.5,
            "nlrp3_overlap_min": 0.50,
            "nlrp3_ifp_jaccard_min": 0.50,
            "nlrp3_key_recovery_min": "5/7",
            "seeds": camp["acid_track"]["clinical_seeds_a2"],
            "stability": camp["acid_track"]["clinical_stability_rule"],
            "urat1_pose_selection": camp["acid_track"]["urat1_pose_selection"],
        },
        "reference_redock_before_any_new_molecule": {
            "purpose": "Same seed should reproduce deposited CNNscore values if the binary, receptor and ligand PDBQT match. A mismatch means the new block cannot be merged with the historical scores.",
            "jobs": [
                "lesinurad at 9DKB, seed 42 (stdout already lists the nine CNN scores)",
                "NP3-146 at 7ALV, seed 42 (deposited self-dock stdout)",
                "lowest repurposing_id among deposited carboxylates with a2_geometry_then_cnnscore: REP_00015, both targets, seed 42",
                "lowest repurposing_id among deposited carboxylates with a2_no_geometry_pass: REP_00300, 9DKB, seed 42",
                "REP_00015 and REP_00300 at 9DKB, seed 43, because the seed 43/44 SDF headers are from the next calendar day",
            ],
            "accept_if": "all nine deposited CNNscore values match to 1e-4; otherwise keep the new block separate",
            "n_reference_jobs": 7,
        },
    }


def _figures(pred: pd.DataFrame, scores: pd.DataFrame, ad: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    bins = np.linspace(0, 1, 21)
    ax.hist(scores["max_tanimoto_train"], bins=bins, density=True, alpha=0.55, label="8319 clinical")
    ax.hist(pred["nn_tanimoto_to_outer_train"], bins=bins, density=True, alpha=0.55, label="513 CV test vs own fold")
    ax.axvline(ad, color="k", ls="--", lw=1, label=f"existing AD {ad:.3f}")
    ax.set_xlabel("max Morgan Tanimoto to the relevant training molecules")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "similarity_cv_vs_clinical.png", dpi=160)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_records()
    mol = molecule_table(df)
    print("similarity CV...")
    pred, bins, sim_summary = similarity_cv(mol)
    pred.to_csv(OUT / "oof_similarity.csv", index=False)
    bins.to_csv(OUT / "similarity_bin_metrics.csv", index=False)

    scores = pd.read_csv(SCORES)
    sim_summary["clinical_library_nn"] = {
        "n": int(len(scores)),
        "median": float(scores["max_tanimoto_train"].median()),
        "p95": float(scores["max_tanimoto_train"].quantile(0.95)),
        "frac_below_0_30": float((scores["max_tanimoto_train"] < 0.30).mean()),
        "frac_ge_ad": float((scores["max_tanimoto_train"] >= sim_summary["ad_threshold_recomputed"]).mean()),
        "n_ge_ad": int((scores["max_tanimoto_train"] >= sim_summary["ad_threshold_recomputed"]).sum()),
    }
    panel_scores = scores[scores["in_v2_pool"] & scores["pass_chemistry_soft"]]
    sim_summary["chemistry_pass_157_nn"] = {
        "n": int(len(panel_scores)),
        "median": float(panel_scores["max_tanimoto_train"].median()),
        "max": float(panel_scores["max_tanimoto_train"].max()),
        "n_ge_ad": int((panel_scores["max_tanimoto_train"] >= sim_summary["ad_threshold_recomputed"]).sum()),
        "frac_below_0_30": float((panel_scores["max_tanimoto_train"] < 0.30).mean()),
    }
    (OUT / "similarity_summary.json").write_text(json.dumps(sim_summary, indent=2))
    _figures(pred, scores, sim_summary["ad_threshold_recomputed"])

    print("annotating the 157...")
    panel = annotate_panel(scores)
    panel.to_csv(OUT / "chemistry_pass_157_annotation.csv", index=False)
    decisions = assign_docking(panel)
    decisions.to_csv(OUT / "new_molecule_docking_decisions.csv", index=False)

    frozen = pd.read_csv(PROJECT / "data/frozen/primary_candidates_frozen.csv")
    merged = panel.merge(
        frozen[["ligand_id", "pains_any", "brenk"]].rename(columns={"ligand_id": "repurposing_id", "pains_any": "frozen_pains", "brenk": "frozen_brenk"}),
        on="repurposing_id",
        how="inner",
    )
    alert_check = {
        "n_primary_in_panel": int(len(merged)),
        "pains_matches_frozen_table": bool((merged["pains_any"] == merged["frozen_pains"]).all()) if len(merged) else None,
        "brenk_matches_frozen_table": bool((merged["brenk"] == merged["frozen_brenk"]).all()) if len(merged) else None,
    }
    cohort = {}
    for name, g in panel.groupby("cohort"):
        cohort[name] = {
            "n": int(len(g)),
            "n_carboxylate": int(g["has_carboxylate"].sum()),
            "n_noncarboxylate": int((~g["has_carboxylate"]).sum()),
            "n_beta_lactam": int(g["beta_lactam"].sum()),
            "n_pains": int(g["pains_any"].sum()),
            "n_brenk": int(g["brenk"].sum()),
            "n_nih": int(g["nih"].sum()),
            "n_alert_clean_pains_brenk": int(g["alert_clean_pains_brenk"].sum()),
            "n_in_training_set": int(g["in_training_set"].sum()),
            "acid_class_counts": g["acid_class"].value_counts().to_dict(),
            "historical_urat1_status_among_noncarboxylate": (
                g.loc[~g["has_carboxylate"], "pose_selection_status"].fillna("not_docked").value_counts().to_dict()
            ),
        }
    chem_summary = {
        "n_panel": int(len(panel)),
        "cohorts": cohort,
        "alert_reproduction_on_frozen_primary": alert_check,
        "noncarboxylate_rule": (
            "Report as outside the current carboxylate-Arg477 method. Do not interpret a missing "
            "carboxylate match as absent URAT1 activity, and do not invent a tetrazole or acyl-sulfonamide "
            "distance threshold in this phase."
        ),
        "beta_lactam_rule": (
            "Flag with SMARTS [#6]1[#6][#7][#6]1=O. Carboxylate beta-lactams stay on the docking list "
            "because the historical 156 were docked before this flag removed them from the reportable reserve. "
            "The flag is applied at nomination, with the same threshold-free structural rule."
        ),
        "dock": {
            "n_dock": int((decisions["docking_decision"] == "dock").sum()),
            "n_exclude": int((decisions["docking_decision"] == "exclude").sum()),
            "n_dock_beta_lactam": int(((decisions["docking_decision"] == "dock") & decisions["beta_lactam"]).sum()),
            "gnina_jobs_if_three_seeds_two_targets": int((decisions["docking_decision"] == "dock").sum() * 6),
        },
    }
    (OUT / "chemistry_scope_summary.json").write_text(json.dumps(chem_summary, indent=2))
    protocol = protocol_record()
    protocol["software_present_in_this_environment"] = {
        "python": platform.python_version(),
        "rdkit_via_stage2": "used for this audit",
        "scikit_learn": sklearn.__version__,
        "xgboost": xgb.__version__,
        "dimorphite_dl": False,
        "meeko": False,
        "gnina_binary": False,
    }
    (OUT / "docking_protocol_audit.json").write_text(json.dumps(protocol, indent=2))
    print(json.dumps({"similarity_oof_nn": sim_summary["oof_nn"], "dock": chem_summary["dock"], "alert_check": alert_check}, indent=2))


if __name__ == "__main__":
    main()
