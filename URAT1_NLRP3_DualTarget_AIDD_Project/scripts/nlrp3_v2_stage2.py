#!/usr/bin/env python3
"""Stage 2 of the NLRP3 v2 rebuild.

1. Recount multi-record molecules and label conflicts from the training table.
2. Sensitivity of the historical max-activity aggregation, decided without any
   clinical candidate.
3. Fit the uncalibrated structure-only model on all 513 molecules with the
   validated hyperparameters, and keep the existing 392-molecule model.
4. Score the 8319-compound clinical library and compare the frozen pools.

Reads frozen files. Writes only under results/nlrp3_model_v2/.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from build_c1_acid_clinical_pool import acid_flags
from nlrp3_model_v2_validate import classification_report, load_records, molecule_table
from nlrp3_v2_lib import StructureOnlyModel, xgb_classifier
from utils_ml import canonicalize, featurize_smiles, scaffold_cv_indices

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "nlrp3_model_v2" / "stage2"
MODEL_DIR = PROJECT / "results" / "nlrp3_model_v2"
THRESHOLD = 0.5
ACTIVITY_CUTOFF = 6.0
SEED = 42
N_BOOT = 2000
# Fixed before the clinical library is opened. Median is the fallback because
# min lets one discordant assay veto a molecule and mean is outlier-sensitive.
AUROC_DROP_BAR = 0.02
RULES = ("max", "min", "median", "mean")


def recount(df: pd.DataFrame) -> dict:
    def _one(key: str) -> dict:
        sizes = df.groupby(key).size()
        flags = df.groupby(key)["active"].agg(["min", "max", "size"])
        discordant = (flags["size"] > 1) & (flags["min"] != flags["max"])
        return {
            "n_groups": int(len(sizes)),
            "n_multi_record": int((sizes > 1).sum()),
            "n_label_conflicts": int(discordant.sum()),
            "extra_records": int((sizes - 1).clip(lower=0).sum()),
        }

    by_smiles = df.groupby("canonical_smiles")
    sizes = by_smiles.size()
    flags = by_smiles["active"].agg(["min", "max"])
    discordant = sizes.index[(sizes > 1) & (flags["min"] != flags["max"])]
    n_assay = by_smiles["Assay ChEMBL ID"].nunique()
    within = df.groupby(["canonical_smiles", "Assay ChEMBL ID"])["active"].agg(["min", "max", "size"])
    within_conflict = (within["size"] > 1) & (within["min"] != within["max"])
    pact = by_smiles["pActivity"].agg(["min", "max"])
    gap = (pact.loc[discordant, "max"] - pact.loc[discordant, "min"]) if len(discordant) else pd.Series(dtype=float)

    # Definitions that were checked while looking for the prose figures 55 and 23.
    checked = {}
    multi = sizes[sizes > 1]
    ranges = (pact["max"] - pact["min"]).loc[multi.index]
    for thr in (0.0, 0.3, 0.5, 1.0, 1.5):
        checked[f"multi_and_pactivity_range_gt_{thr}"] = int((ranges > thr).sum())
        cross = pact.loc[multi.index]
        cross = cross[(cross["min"] < ACTIVITY_CUTOFF) & (cross["max"] >= ACTIVITY_CUTOFF) & ((cross["max"] - cross["min"]) > thr)]
        checked[f"conflict_and_range_gt_{thr}"] = int(len(cross))

    return {
        "source": "data/processed/nlrp3_records.csv",
        "n_records": int(len(df)),
        "n_molecules": int(df["canonical_smiles"].nunique()),
        "by_canonical_smiles": _one("canonical_smiles"),
        "by_molecule_chembl_id": _one("Molecule ChEMBL ID"),
        "by_raw_smiles": _one("Smiles"),
        "multi_record_size_distribution": {int(k): int(v) for k, v in sizes[sizes > 1].value_counts().sort_index().items()},
        "multi_record_molecules_spanning_more_than_one_assay": int(((sizes > 1) & (n_assay > 1)).sum()),
        "multi_record_molecules_with_a_single_assay": int(((sizes > 1) & (n_assay == 1)).sum()),
        "within_assay_label_conflicts": int(within_conflict.sum()),
        "conflict_pactivity_gap": {
            "n": int(len(gap)),
            "min": float(gap.min()) if len(gap) else None,
            "median": float(gap.median()) if len(gap) else None,
            "max": float(gap.max()) if len(gap) else None,
        },
        "definitions_checked_for_55_and_23": checked,
        "resolution": (
            "80 multi-record molecules and 32 cross-assay label conflicts, under canonical "
            "SMILES, Molecule ChEMBL ID and the raw SMILES column alike. No within-assay "
            "conflict. The figures 55 and 23 in the stage-1 narrative do not match any of "
            "these groupings or the pActivity-gap thresholds checked here; the machine-readable "
            "stage-1 audit already recorded 80 and 32."
        ),
    }


def aggregation_labels(df: pd.DataFrame, mol: pd.DataFrame) -> pd.DataFrame:
    agg = df.groupby("canonical_smiles")["pActivity"].agg(["min", "max", "median", "mean", "size"])
    agg = agg.reindex(mol["canonical_smiles"].tolist())
    out = pd.DataFrame({"canonical_smiles": mol["canonical_smiles"].values, "n_records": agg["size"].values})
    for rule in RULES:
        out[rule] = (agg[rule].values >= ACTIVITY_CUTOFF).astype(int)
    historical = mol["active"].astype(int).values
    if not np.array_equal(out["max"].values, historical):
        raise RuntimeError("max(pActivity)>=6 does not reproduce the stored molecule label")
    return out


def _cv_probabilities(mol: pd.DataFrame, y: np.ndarray, folds, seed: int) -> np.ndarray:
    p = np.full(len(mol), np.nan)
    for fold_i, (tr, te) in enumerate(folds):
        fit = mol.iloc[tr][["canonical_smiles"]].copy()
        fit["active"] = y[tr]
        model = StructureOnlyModel().fit(fit, seed=seed + fold_i)
        p[te] = model.predict_molecules(mol.iloc[te]["canonical_smiles"].tolist())
    if np.isnan(p).any():
        raise RuntimeError("a molecule received no outer-test prediction")
    return p


def _boot_delta(y_ref: np.ndarray, y_alt: np.ndarray, p: np.ndarray, n_boot: int, seed: int) -> dict:
    """AUROC under y_alt minus AUROC under y_ref, same scores, resampled molecules."""
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    y_ref = np.asarray(y_ref).astype(int)
    y_alt = np.asarray(y_alt).astype(int)
    observed = float(roc_auc_score(y_alt, p) - roc_auc_score(y_ref, p))
    deltas = []
    n = len(p)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_ref[idx])) < 2 or len(np.unique(y_alt[idx])) < 2:
            continue
        deltas.append(roc_auc_score(y_alt[idx], p[idx]) - roc_auc_score(y_ref[idx], p[idx]))
    deltas = np.asarray(deltas)
    return {
        "delta_auroc": observed,
        "ci95_low": float(np.percentile(deltas, 2.5)),
        "ci95_high": float(np.percentile(deltas, 97.5)),
        "n_bootstrap": int(len(deltas)),
    }


def sensitivity(df: pd.DataFrame, mol: pd.DataFrame) -> tuple[pd.DataFrame, dict, str]:
    labels = aggregation_labels(df, mol)
    folds = scaffold_cv_indices(mol["canonical_smiles"].tolist(), n_splits=5)
    probs = {}
    rows = []
    for rule in RULES:
        y = labels[rule].values.astype(int)
        p = _cv_probabilities(mol, y, folds, SEED)
        probs[rule] = p
        pooled = classification_report(y, p, THRESHOLD)
        rows.append({"evaluation": "trained_and_scored_under_this_rule", "rule": rule, "fold": "pooled", **pooled})
        for fold_i, (_, te) in enumerate(folds):
            m = classification_report(y[te], p[te], THRESHOLD)
            rows.append({"evaluation": "trained_and_scored_under_this_rule", "rule": rule, "fold": fold_i, **m})

    # Same max-trained scores, four label definitions. This is the comparison
    # that decides whether the historical label is load-bearing.
    p_max = probs["max"]
    cross = {}
    for rule in RULES:
        y = labels[rule].values.astype(int)
        pooled = classification_report(y, p_max, THRESHOLD)
        rows.append({"evaluation": "max_trained_scores_against_this_label", "rule": rule, "fold": "pooled", **pooled})
        if rule != "max":
            cross[rule] = _boot_delta(labels["max"].values, y, p_max, N_BOOT, SEED)

    flips = {}
    for rule in RULES:
        if rule == "max":
            continue
        differ = labels[rule].values != labels["max"].values
        flips[rule] = {
            "n_molecules_relabelled_vs_max": int(differ.sum()),
            "n_max_active_lost": int(((labels["max"] == 1) & (labels[rule] == 0)).sum()),
            "n_max_inactive_gained": int(((labels["max"] == 0) & (labels[rule] == 1)).sum()),
        }

    median_delta = cross["median"]
    # Switch only if the same scores discriminate the median label worse than
    # the max label by more than the pre-set bar, with the whole CI past it.
    reject_max = median_delta["ci95_high"] < -AUROC_DROP_BAR
    chosen = "median" if reject_max else "max"
    decision = {
        "criterion": (
            f"Keep max unless AUROC(max-trained scores, median labels) - "
            f"AUROC(max-trained scores, max labels) has its 95% CI entirely below "
            f"-{AUROC_DROP_BAR}. Median is the only fallback. Candidate identity is not an input."
        ),
        "auroc_drop_bar": AUROC_DROP_BAR,
        "cross_label_delta_vs_max": cross,
        "label_flips_vs_max": flips,
        "active_counts": {rule: int(labels[rule].sum()) for rule in RULES},
        "reject_max": bool(reject_max),
        "production_label_rule": chosen,
        "reason": (
            "The max-trained scores lose more than the pre-set AUROC bar on the median "
            "label, so the production model is fit on median labels."
            if reject_max
            else
            "The drop from the max label to the median label stays inside the pre-set "
            "bar, so the historical max rule remains the production label."
        ),
    }
    return pd.DataFrame(rows), decision, chosen, labels, probs


def fit_production(mol: pd.DataFrame, y: np.ndarray, rule: str) -> dict:
    fit = mol[["canonical_smiles"]].copy()
    fit["active"] = y
    model = StructureOnlyModel().fit(fit, seed=SEED)
    return {
        "model": model,
        "calibrators": {"none": None},
        "prediction_rule": (
            "single molecule-level prediction from Morgan2048+12 RDKit descriptors; "
            "no post-hoc calibration; threshold 0.5"
        ),
        "n_train_molecules": int(len(fit)),
        "label_rule": f"{rule}(pActivity) >= {ACTIVITY_CUTOFF}",
        "hyperparameters": xgb_classifier(SEED).get_params(),
        "seed": SEED,
        "threshold": THRESHOLD,
        "calibration": "none",
    }


def _fps(smiles: list[str]):
    out = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        out.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None)
    return out


def nearest_tanimoto(query, ref, exclude_self: bool = False) -> np.ndarray:
    ref = [f for f in ref if f is not None]
    out = np.full(len(query), np.nan)
    for i, q in enumerate(query):
        if q is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(q, ref)
        if not sims:
            continue
        if exclude_self:
            sims = sorted(sims, reverse=True)
            out[i] = sims[1] if len(sims) > 1 else sims[0]
        else:
            out[i] = max(sims)
    return out


def chemistry_table(smiles: list[str]) -> pd.DataFrame:
    rows = []
    for smi in smiles:
        flags = acid_flags(smi)
        flags.pop("mol", None)
        rows.append(
            {
                "acid_equivalent": bool(flags.get("ok")),
                "acid_reject_reason": None if flags.get("ok") else flags.get("reason"),
                "has_carboxylate": bool(flags.get("has_carboxylate", False)),
                "has_tetrazole": bool(flags.get("has_tetrazole", False)),
                "has_acyl_sulfonamide": bool(flags.get("has_acyl_sulfonamide", False)),
                "pass_veber": bool(flags.get("pass_veber", False)),
                "pass_ro5_hb": bool(flags.get("pass_ro5_hb", False)),
                "pass_mw_200_550": bool(flags.get("pass_mw_200_550", False)),
                "pass_chemistry_soft": bool(flags.get("pass_chemistry_soft", False)),
            }
        )
    return pd.DataFrame(rows)


def _a2_keep_tables() -> list[pd.DataFrame]:
    base = PROJECT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2"
    return [pd.read_csv(base / f"acid_dual_keep_a2_seed{seed}.csv") for seed in (42, 43, 44)]


def deposited_docking_ids() -> set[str]:
    """Every ligand the A2 campaign actually docked, on all three seeds.

    A null URAT1 score in the keep table means the geometry gate found no pose
    (``a2_no_geometry_pass``); the nine gnina poses are still deposited.
    """
    ids = None
    for d in _a2_keep_tables():
        present = set(d["ligand_id"].astype(str))
        ids = present if ids is None else ids & present
    return ids


def docking_complete_ids() -> set[str]:
    """Ligands with a non-null score for both targets on every A2 seed."""
    ids = None
    for d in _a2_keep_tables():
        ok = set(d.loc[d["u_CNNscore"].notna() & d["n_CNNscore"].notna(), "ligand_id"].astype(str))
        ids = ok if ids is None else ids & ok
    return ids


def score_library(production: dict, label_rule: str) -> None:
    manifest = pd.read_csv(PROJECT / "data/repurposing/repurposing_manifest.csv", low_memory=False)
    historical = pd.read_csv(PROJECT / "data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv", low_memory=False)
    old_pool = pd.read_csv(PROJECT / "data/repurposing/screening/docking_pool_p05.csv", low_memory=False)
    primary = pd.read_csv(PROJECT / "data/frozen/primary_candidates_frozen.csv")
    backup = pd.read_csv(PROJECT / "data/frozen/backup_candidates_frozen.csv")
    train = pd.read_csv(PROJECT / "data/processed/nlrp3_records.csv", usecols=["canonical_smiles"])
    train_smiles = set(train["canonical_smiles"].astype(str))

    lib = manifest.copy()
    lib["canonical_smiles_raw"] = lib["canonical_smiles"].astype(str)
    lib["canonical_smiles"] = lib["canonical_smiles_raw"].map(canonicalize)
    n_unparsed = int(lib["canonical_smiles"].isna().sum())
    lib = lib[lib["canonical_smiles"].notna()].drop_duplicates("canonical_smiles").reset_index(drop=True)

    smiles = lib["canonical_smiles"].tolist()
    x = featurize_smiles(smiles)
    p513 = production["model"].clf.predict_proba(production["model"].scaler.transform(x))[:, 1]

    held = joblib.load(MODEL_DIR / "nlrp3_v2_structure_only.joblib")
    p392 = held["model"].clf.predict_proba(held["model"].scaler.transform(x))[:, 1]

    train_fps = _fps(sorted(train_smiles))
    intra = nearest_tanimoto(train_fps, train_fps, exclude_self=True)
    ad_threshold = float(np.nanpercentile(intra, 5))
    query_fps = _fps(smiles)
    nn = nearest_tanimoto(query_fps, train_fps, exclude_self=False)

    chem = chemistry_table(smiles)
    docked_set = deposited_docking_ids()
    docked_complete = docking_complete_ids()
    hist = historical[["repurposing_id", "p_active_nlrp3", "nlrp3_pred_active"]].rename(
        columns={"p_active_nlrp3": "p_historical", "nlrp3_pred_active": "historical_pred_active"}
    )
    out = pd.concat([lib.reset_index(drop=True), chem], axis=1)
    out["p_historical"] = out["repurposing_id"].map(hist.set_index("repurposing_id")["p_historical"])
    out["p_n392"] = p392
    out["p_n513"] = p513
    out["in_historical_pool"] = out["repurposing_id"].isin(set(old_pool["repurposing_id"].astype(str)))
    out["in_v2_pool"] = out["p_n513"] >= THRESHOLD
    out["in_n392_pool"] = out["p_n392"] >= THRESHOLD
    out["max_tanimoto_train"] = nn
    out["ad_threshold"] = ad_threshold
    out["in_applicability_domain"] = nn >= ad_threshold
    out["in_training_set"] = out["canonical_smiles"].isin(train_smiles)
    out["in_deposited_a2_docking_set"] = out["repurposing_id"].astype(str).isin(docked_set)
    out["docking_complete_a2_three_seed"] = out["repurposing_id"].astype(str).isin(docked_complete)
    out["historical_name_exclusion_vecabrutinib"] = out["name"].astype(str).str.upper().str.contains("VECABRUTINIB")

    keep = [
        "repurposing_id", "chembl_id", "name", "canonical_smiles", "max_phase",
        "p_historical", "p_n392", "p_n513",
        "in_historical_pool", "in_v2_pool", "in_n392_pool",
        "max_tanimoto_train", "ad_threshold", "in_applicability_domain", "in_training_set",
        "acid_equivalent", "acid_reject_reason", "has_carboxylate", "has_tetrazole",
        "has_acyl_sulfonamide", "pass_veber", "pass_ro5_hb", "pass_mw_200_550",
        "pass_chemistry_soft", "in_deposited_a2_docking_set",
        "docking_complete_a2_three_seed",
        "historical_name_exclusion_vecabrutinib",
    ]
    scored = out[keep].sort_values("p_n513", ascending=False).reset_index(drop=True)
    scored.to_csv(OUT / "clinical_library_scores.csv", index=False)

    def _funnel(frame: pd.DataFrame) -> dict:
        return {
            "n": int(len(frame)),
            "n_acid_equivalent": int(frame["acid_equivalent"].sum()),
            "n_chemistry_soft_pass": int(frame["pass_chemistry_soft"].sum()),
            "n_in_applicability_domain": int(frame["in_applicability_domain"].sum()),
            "n_in_deposited_docking_set": int(frame["in_deposited_a2_docking_set"].sum()),
            "n_both_target_scores_all_three_seeds": int(frame["docking_complete_a2_three_seed"].sum()),
            "median_max_tanimoto": float(frame["max_tanimoto_train"].median()) if len(frame) else None,
        }

    v2 = scored[scored["in_v2_pool"]]
    old = scored[scored["in_historical_pool"]]
    both = scored[scored["in_v2_pool"] & scored["in_historical_pool"]]
    entered = scored[scored["in_v2_pool"] & ~scored["in_historical_pool"]]
    left = scored[scored["in_historical_pool"] & ~scored["in_v2_pool"]]
    qualified_new = entered[entered["pass_chemistry_soft"]]
    # New entrants are outside the frozen 1588, so none of them were in the
    # 156-ligand docking campaign. The gap is "never docked", not "geometry gate failed".
    gap = qualified_new[~qualified_new["in_deposited_a2_docking_set"]].copy()
    # Annotation only. Beta-lactam exclusion was applied at nomination, after
    # docking, so it does not remove anyone from this list.
    lactam = Chem.MolFromSmarts("[#6]1[#6][#7][#6]1=O")
    gap["beta_lactam_flag_not_used_as_filter"] = [
        bool(m is not None and m.HasSubstructMatch(lactam))
        for m in (Chem.MolFromSmiles(s) for s in gap["canonical_smiles"])
    ]

    summary = {
        "n_manifest_rows": int(len(manifest)),
        "n_unparsed_smiles": n_unparsed,
        "n_scored": int(len(scored)),
        "threshold": THRESHOLD,
        "formal_model": "nlrp3_v2_structure_only_n513_production.joblib",
        "label_rule": label_rule,
        "ad_rule": "in domain if max Morgan-r2 Tanimoto to the 513 training molecules >= the 5th percentile of leave-one-out nearest-neighbour Tanimoto inside that training set (the rule in scripts/12_ml_rigor_validation.py)",
        "ad_threshold": ad_threshold,
        "training_loo_nn_tanimoto": {
            "p05": ad_threshold,
            "median": float(np.nanmedian(intra)),
            "p25": float(np.nanpercentile(intra, 25)),
        },
        "v2_pool": _funnel(v2),
        "historical_pool_rescored": _funnel(old),
        "n_historical_pool_file": int(len(old_pool)),
        "overlap": {
            "n_both": int(len(both)),
            "n_entered": int(len(entered)),
            "n_left": int(len(left)),
            "jaccard": float(len(both) / len(set(v2["repurposing_id"]) | set(old["repurposing_id"]))),
            "historical_retained_fraction": float(len(both) / len(old)) if len(old) else None,
        },
        "n392_pool_n": int(scored["in_n392_pool"].sum()),
        "n392_vs_n513_jaccard": float(
            (scored["in_n392_pool"] & scored["in_v2_pool"]).sum()
            / (scored["in_n392_pool"] | scored["in_v2_pool"]).sum()
        ),
        "library_in_training_set": int(scored["in_training_set"].sum()),
        "v2_pool_tanimoto": {
            "median": float(v2["max_tanimoto_train"].median()),
            "p05": float(v2["max_tanimoto_train"].quantile(0.05)),
            "frac_below_0_4": float((v2["max_tanimoto_train"] < 0.4).mean()),
            "frac_in_domain": float(v2["in_applicability_domain"].mean()),
        },
        "new_entrants": _funnel(entered),
        "new_entrants_chemistry_pass_never_docked": int(len(gap)),
        "vecabrutinib_name_hits_in_v2_pool": int(v2["historical_name_exclusion_vecabrutinib"].sum()),
        "name_exclusion_applied": False,
        "supplementary_docking": {
            "needed": bool(len(gap) > 0),
            "n": int(len(gap)),
            "n_beta_lactam_annotated": int(gap["beta_lactam_flag_not_used_as_filter"].sum()),
            "reason": (
                "These molecules are outside the frozen 1588, pass the existing acid-equivalent "
                "and soft-chemistry gates, and are absent from the deposited 156-ligand A2 docking "
                "campaign. That campaign is what the historical funnel docked. No docking was run "
                "in this stage. A null score inside the deposited 156 is a geometry-gate miss "
                "(a2_no_geometry_pass), not a missing calculation, and those ligands are not in this list."
                if len(gap)
                else "Every newly entered chemistry-pass molecule is already in the deposited docking set."
            ),
        },
    }
    (OUT / "library_summary.json").write_text(json.dumps(summary, indent=2))
    entered.to_csv(OUT / "new_pool_entrants.csv", index=False)
    left.to_csv(OUT / "historical_pool_dropouts.csv", index=False)
    gap.sort_values("p_n513", ascending=False).to_csv(OUT / "supplementary_docking_candidates.csv", index=False)

    # Molecule-by-molecule comparison for the frozen sets. PF is reported both
    # inside the 12 and as its own row group; nothing is dropped by name.
    frames = []
    for group, table, id_col in (
        ("primary", primary, "ligand_id"),
        ("reserve", backup, "ligand_id"),
        ("historical_pool", old_pool, "repurposing_id"),
    ):
        ids = table[id_col].astype(str)
        part = scored[scored["repurposing_id"].astype(str).isin(set(ids))].copy()
        part.insert(0, "comparison_group", group)
        frames.append(part)
    comparison = pd.concat(frames, ignore_index=True)
    frozen_p = pd.concat(
        [
            primary[["ligand_id", "p_active_nlrp3"]].rename(columns={"ligand_id": "repurposing_id"}),
            backup[["ligand_id", "p_active_nlrp3"]].rename(columns={"ligand_id": "repurposing_id"}),
        ]
    ).drop_duplicates("repurposing_id")
    comparison["p_recorded_on_frozen_candidate_table"] = comparison["repurposing_id"].map(
        frozen_p.set_index("repurposing_id")["p_active_nlrp3"]
    )
    comparison["is_pf_03882845"] = comparison["name"].astype(str).eq("PF-03882845")
    comparison.to_csv(OUT / "frozen_sets_molecule_comparison.csv", index=False)

    candidate_groups = comparison[comparison["comparison_group"].isin(["primary", "reserve"])].copy()
    candidate_groups.to_csv(OUT / "candidates_12_21_pf_comparison.csv", index=False)

    _figures(scored, comparison)
    _candidate_summary(candidate_groups, summary)


def _candidate_summary(candidates: pd.DataFrame, library_summary: dict) -> None:
    rows = []
    for group, g in candidates.groupby("comparison_group"):
        rows.append(
            {
                "group": group,
                "n": int(len(g)),
                "n_retained_n513": int(g["in_v2_pool"].sum()),
                "n_retained_n392": int(g["in_n392_pool"].sum()),
                "n_in_domain": int(g["in_applicability_domain"].sum()),
                "n_in_training_set": int(g["in_training_set"].sum()),
                "min_p_n513": float(g["p_n513"].min()),
                "median_p_n513": float(g["p_n513"].median()),
            }
        )
    pf = candidates[candidates["is_pf_03882845"]]
    payload = {
        "groups": rows,
        "pf_03882845": pf.to_dict(orient="records"),
        "primary_not_retained": candidates[
            (candidates["comparison_group"] == "primary") & ~candidates["in_v2_pool"]
        ][["repurposing_id", "name", "p_historical", "p_n392", "p_n513"]].to_dict(orient="records"),
        "reserve_not_retained": candidates[
            (candidates["comparison_group"] == "reserve") & ~candidates["in_v2_pool"]
        ][["repurposing_id", "name", "p_historical", "p_n392", "p_n513"]].to_dict(orient="records"),
        "library": {
            "v2_pool_n": library_summary["v2_pool"]["n"],
            "jaccard_vs_historical": library_summary["overlap"]["jaccard"],
            "supplementary_docking_n": library_summary["supplementary_docking"]["n"],
        },
    }
    (OUT / "candidate_retention_summary.json").write_text(json.dumps(payload, indent=2))


def _figures(scored: pd.DataFrame, comparison: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(scored.loc[~scored["in_v2_pool"], "p_n513"], bins=40, range=(0, 1), alpha=0.7, label="below 0.5")
    ax.hist(scored.loc[scored["in_v2_pool"], "p_n513"], bins=40, range=(0, 1), alpha=0.7, label="v2 pool")
    ax.axvline(THRESHOLD, color="k", ls="--", lw=1)
    ax.set_xlabel("production-model probability")
    ax.set_ylabel("clinical molecules")
    ax.legend(fontsize=8)
    ax.set_title("8319 clinical compounds, uncalibrated structure-only model")
    fig.tight_layout()
    fig.savefig(OUT / "library_probability_distribution.png", dpi=160)
    plt.close(fig)

    old = comparison[comparison["comparison_group"] == "historical_pool"]
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    retained = old["in_v2_pool"].values
    ax.scatter(old.loc[~retained, "p_historical"], old.loc[~retained, "p_n513"], s=8, alpha=0.5, label="left the pool")
    ax.scatter(old.loc[retained, "p_historical"], old.loc[retained, "p_n513"], s=8, alpha=0.5, label="retained")
    ax.axhline(THRESHOLD, color="k", ls="--", lw=0.8)
    ax.axvline(THRESHOLD, color="k", ls="--", lw=0.8)
    ax.set_xlabel("historical calibrated probability")
    ax.set_ylabel("v2 production probability")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.set_title("Frozen 1588, rescored")
    fig.tight_layout()
    fig.savefig(OUT / "historical_pool_rescore.png", dpi=160)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_records()
    mol = molecule_table(df)

    audit = recount(df)
    (OUT / "label_recount.json").write_text(json.dumps(audit, indent=2))
    smiles_n = audit["by_canonical_smiles"]
    if smiles_n["n_multi_record"] != 80 or smiles_n["n_label_conflicts"] != 32:
        raise RuntimeError(f"recount is not 80/32: {smiles_n}")
    print(
        f"recount: {smiles_n['n_multi_record']} multi-record molecules, "
        f"{smiles_n['n_label_conflicts']} label conflicts"
    )

    metrics, decision, chosen, labels, _probs = sensitivity(df, mol)
    metrics.to_csv(OUT / "aggregation_sensitivity.csv", index=False)
    labels.to_csv(OUT / "molecule_labels_by_rule.csv", index=False)
    # Written before any clinical or candidate file is read.
    (OUT / "label_rule_decision.json").write_text(json.dumps(decision, indent=2))
    print("production label:", chosen, decision["reason"])

    y = labels[chosen].values.astype(int)
    bundle = fit_production(mol, y, chosen)
    prod_path = MODEL_DIR / "nlrp3_v2_structure_only_n513_production.joblib"
    joblib.dump(bundle, prod_path)

    index = {
        "formal_production_model": {
            "path": str(prod_path.relative_to(PROJECT)),
            "n_train_molecules": 513,
            "molecules": "all unique training molecules",
            "calibration": "none",
            "label_rule": bundle["label_rule"],
            "threshold": THRESHOLD,
            "hyperparameters_match_stage1": True,
        },
        "stage1_holdout_model": {
            "path": "results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib",
            "n_train_molecules": 392,
            "molecules": "scaffold holdout removed 121 molecules that stage 1 reserved for a calibration that was not used",
            "calibration": "none",
            "retained": "not overwritten",
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "xgboost": xgb.__version__,
        },
    }
    (OUT / "model_index.json").write_text(json.dumps(index, indent=2, default=str))

    print("scoring the clinical library...")
    score_library(bundle, bundle["label_rule"])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
