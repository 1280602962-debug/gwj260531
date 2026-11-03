#!/usr/bin/env python3
"""Analysis B: EGFR/HER2 S1 inversion diagnostics on the frozen D_vs_B set."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from posthoc_common import (
    DESCRIPTORS,
    INPUTS,
    PROJECT,
    auroc,
    formal_egfr_row,
    iqr,
    load_config,
    quartiles,
)

ACTIVITY_SOURCE = (
    "chembl37_dual_endpoint_crosswalk_805 historical_panel endpoint maxima; "
    "accepted only when dump_max, historical_panel and historical_aggregate_max "
    "agree exactly and class_matches_frozen is YES"
)


def _dist_row(cohort: str, endpoint: str, values: np.ndarray) -> dict:
    q1, q2, q3 = quartiles(values)
    return {
        "cohort": cohort,
        "endpoint": endpoint,
        "n": int(values.size),
        "median": q2,
        "q1": q1,
        "q3": q3,
        "iqr": iqr(values),
        "minimum": float(values.min()) if values.size else float("nan"),
        "maximum": float(values.max()) if values.size else float("nan"),
    }


def _activity_bin(value: float) -> str:
    if value < 5.0:
        return "<5.0"
    if value < 5.5:
        return "5.0_to_<5.5"
    if value < 6.0:
        return "5.5_to_<6.0"
    return ">=6.0"


def egfr_members() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    pop = pd.read_csv(PROJECT / INPUTS["population"])
    base = pop[(pop["pair"] == "EGFR/HER2") & (pop["activity_eligible"] == 1)].copy()
    primary = base[(base["class"].isin(["dual", "B_only"])) & (base["M0_A_valid"] == 1)].copy()
    neither = base[(base["class"] == "neither") & (base["M0_A_valid"] == 1)].copy()
    cfg = load_config()
    n_dual = int((primary["class"] == "dual").sum())
    n_b = int((primary["class"] == "B_only").sum())
    if n_dual != int(cfg["egfr_expected_n_dual"]) or n_b != int(cfg["egfr_expected_n_b_only"]):
        raise RuntimeError(f"EGFR D_vs_B membership changed: dual {n_dual}, B_only {n_b}")
    if primary.duplicated(["canonical_ligand_id"]).any():
        raise RuntimeError("duplicate EGFR ligand ids")
    return primary, neither, {"population_rows": int(len(pop)), "egfr_eligible": int(len(base))}


def build_analysis_b() -> dict[str, pd.DataFrame]:
    primary, neither, scope = egfr_members()
    formal = formal_egfr_row()
    dual_scores = primary.loc[primary["class"] == "dual", "M0_score_A"].to_numpy(dtype=float)
    b_scores = primary.loc[primary["class"] == "B_only", "M0_score_A"].to_numpy(dtype=float)
    neither_scores = neither["M0_score_A"].to_numpy(dtype=float)
    recomputed = auroc(dual_scores, b_scores)
    formal_auroc = float(formal["auroc"])
    # PRIMARY_DIRECTIONAL_METRICS stores this AUROC at 10 decimal places.
    # The exact Mann-Whitney value is retained for the diagnostic.
    if round(recomputed + 1e-15, 10) != round(formal_auroc + 1e-15, 10):
        raise RuntimeError(
            f"EGFR D_vs_B AUROC {recomputed} does not round to formal {formal_auroc}; diagnostics stopped"
        )
    diff = dual_scores[:, None] - b_scores[None, :]
    n_pairs = int(diff.size)
    n_forward = int((diff > 0).sum())
    n_reverse = int((diff < 0).sum())
    n_tie = int((diff == 0).sum())
    negative = pd.read_csv(PROJECT / INPUTS["negative_class"])
    stored = negative[(negative["pair"] == "EGFR/HER2") & (negative["direction"] == "D_vs_B")].iloc[0]
    neither_auroc = auroc(dual_scores, neither_scores)
    if abs(neither_auroc - float(stored["auc_neither"])) > 1e-12:
        raise RuntimeError("dual-neither AUROC does not match module_B")
    if abs(recomputed - float(stored["auc_single"])) > 1e-12:
        raise RuntimeError("dual-B-only AUROC does not match module_B auc_single")
    if int(stored["n_neither"]) != len(neither) or int(stored["n_dual"]) != len(dual_scores):
        raise RuntimeError("module_B EGFR membership does not match the frozen population")

    dist_rows = []
    for label, values in (("dual", dual_scores), ("B_only", b_scores), ("neither", neither_scores)):
        row = _dist_row(label, "M0_score_A_EGFR_S1", values)
        dist_rows.append(row)
    dist_rows.append({
        "cohort": "dual_vs_B_only_pairwise",
        "endpoint": "M0_score_A_EGFR_S1",
        "n": n_pairs,
        "median": n_forward / n_pairs,
        "q1": n_reverse / n_pairs,
        "q3": n_tie / n_pairs,
        "iqr": recomputed,
        "minimum": n_forward,
        "maximum": n_reverse,
    })
    # The pairwise row above overloads columns. Replace it with explicit fields in a second table
    # and keep distribution rows clean.
    distributions = pd.DataFrame([r for r in dist_rows if r["cohort"] in {"dual", "B_only", "neither"}])
    pairwise = pd.DataFrame([{
        "comparison": "dual_vs_B_only",
        "score": "M0_score_A",
        "publication_scheme": "S1",
        "n_dual": int(dual_scores.size),
        "n_B_only": int(b_scores.size),
        "n_pairwise": n_pairs,
        "n_dual_score_higher": n_forward,
        "n_B_only_score_higher": n_reverse,
        "n_tie": n_tie,
        "fraction_dual_higher": n_forward / n_pairs,
        "fraction_B_only_higher": n_reverse / n_pairs,
        "fraction_tie": n_tie / n_pairs,
        "auroc_ties_0_5": recomputed,
        "formal_auroc": formal_auroc,
        "auroc_minus_formal": recomputed - formal_auroc,
        "formal_ci_lo": float(formal["ci_lo"]),
        "formal_ci_hi": float(formal["ci_hi"]),
        "ci_recomputed_here": "NO",
        "stored_auc_neither": float(stored["auc_neither"]),
        "recomputed_auc_neither": neither_auroc,
        "stored_delta_negative": float(stored["delta_negative"]),
        "stored_delta_negative_ci_lo": float(stored["delta_negative_ci_lo"]),
        "stored_delta_negative_ci_hi": float(stored["delta_negative_ci_hi"]),
        "paired_delta_recomputed_here": "NO",
        "n_neither": int(len(neither)),
    }])

    d1 = pd.read_csv(PROJECT / INPUTS["composition"])
    d1_arm = d1[(d1["pair"] == "EGFR/HER2") & (d1["arm"] == "D_vs_B")].copy()
    joined = primary.merge(
        d1_arm,
        on=["pair", "canonical_ligand_id", "class"],
        how="left",
        indicator=True,
    )
    if joined.duplicated(["canonical_ligand_id"]).any():
        raise RuntimeError("descriptor join expanded EGFR rows")
    n_unmatched = int((joined["_merge"] != "both").sum())
    if n_unmatched:
        raise RuntimeError(f"{n_unmatched} EGFR M0 members have no D1 descriptor row")
    desc_rows = []
    assoc_rows = []
    for name in DESCRIPTORS:
        dual_v = joined.loc[joined["class"] == "dual", name].to_numpy(dtype=float)
        b_v = joined.loc[joined["class"] == "B_only", name].to_numpy(dtype=float)
        constant = int(pd.Series(np.concatenate([dual_v, b_v])).nunique(dropna=True) <= 1)
        q1d, med_d, q3d = quartiles(dual_v)
        q1b, med_b, q3b = quartiles(b_v)
        desc_rows.append({
            "descriptor": name,
            "n_dual": int(dual_v.size),
            "n_B_only": int(b_v.size),
            "dual_median": med_d,
            "dual_q1": q1d,
            "dual_q3": q3d,
            "dual_iqr": iqr(dual_v),
            "dual_min": float(dual_v.min()),
            "dual_max": float(dual_v.max()),
            "B_only_median": med_b,
            "B_only_q1": q1b,
            "B_only_q3": q3b,
            "B_only_iqr": iqr(b_v),
            "B_only_min": float(b_v.min()),
            "B_only_max": float(b_v.max()),
            "median_difference_dual_minus_B_only": med_d - med_b,
            "descriptor_auroc_dual_higher": auroc(dual_v, b_v),
            "status": "CONSTANT_VARIABLE" if constant else "OK",
        })
        for scope, frame in (
            ("pooled_dual_and_B_only", joined),
            ("dual", joined[joined["class"] == "dual"]),
            ("B_only", joined[joined["class"] == "B_only"]),
        ):
            x = frame[name].to_numpy(dtype=float)
            y = frame["M0_score_A"].to_numpy(dtype=float)
            if pd.Series(x).nunique(dropna=True) <= 1 or pd.Series(y).nunique(dropna=True) <= 1:
                rho, pvalue, status = float("nan"), float("nan"), "CONSTANT_VARIABLE"
            else:
                rho, pvalue = stats.spearmanr(x, y)
                rho, pvalue, status = float(rho), float(pvalue), "OK"
            assoc_rows.append({
                "descriptor": name,
                "scope": scope,
                "n": int(len(frame)),
                "spearman_rho": rho,
                "spearman_p_unadjusted": pvalue,
                "status": status,
                "interpretation": "exploratory_association_not_a_selection_threshold",
            })

    cross = pd.read_csv(PROJECT / INPUTS["activity_crosswalk"])
    act = primary.merge(cross, on=["pair", "panel_id"], how="left", suffixes=("", "_cross"), indicator=True)
    if int((act["_merge"] != "both").sum()) or act.duplicated(["canonical_ligand_id"]).any():
        raise RuntimeError("EGFR activity crosswalk did not match 1:1")
    class_ok_series = act["frozen_class"].astype(str) == act["class"].astype(str)
    value_ok = np.ones(len(act), dtype=bool)
    for left, right in (
        ("dump_max_A", "historical_panel_pA"),
        ("historical_panel_pA", "historical_aggregate_max_A"),
        ("dump_max_B", "historical_panel_pB"),
        ("historical_panel_pB", "historical_aggregate_max_B"),
    ):
        value_ok &= act[left].astype(float).to_numpy() == act[right].astype(float).to_numpy()
    act["activity_value_status"] = np.where(
        value_ok & class_ok_series.to_numpy() & act["class_matches_frozen"].eq("YES").to_numpy(),
        "VERIFIED_AGGREGATE",
        "NOT_VERIFIED",
    )
    if (act["activity_value_status"] != "VERIFIED_AGGREGATE").any():
        # Do not invent values. Continue with verified rows only and keep unverified rows blank.
        pass
    act["pchembl_EGFR"] = np.where(act["activity_value_status"] == "VERIFIED_AGGREGATE", act["historical_panel_pA"].astype(float), np.nan)
    act["pchembl_HER2"] = np.where(act["activity_value_status"] == "VERIFIED_AGGREGATE", act["historical_panel_pB"].astype(float), np.nan)
    activity_rows = []
    for cohort, frame in (
        ("B_only", act[act["class"] == "B_only"]),
        ("dual", act[act["class"] == "dual"]),
    ):
        for endpoint, column in (("EGFR", "pchembl_EGFR"), ("HER2", "pchembl_HER2")):
            values = frame[column].to_numpy(dtype=float)
            if not np.isfinite(values).all():
                activity_rows.append({
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "status": "NOT_COMPUTABLE",
                    "n": int(len(frame)),
                    "n_verified": int(np.isfinite(values).sum()),
                    "bin": "",
                    "n_in_bin": 0,
                    "median": float("nan"),
                    "q1": float("nan"),
                    "q3": float("nan"),
                    "minimum": float("nan"),
                    "maximum": float("nan"),
                    "value_source": ACTIVITY_SOURCE,
                })
                continue
            q1, med, q3 = quartiles(values)
            counts = pd.Series(values).map(_activity_bin).value_counts()
            for bin_name in ("<5.0", "5.0_to_<5.5", "5.5_to_<6.0", ">=6.0"):
                activity_rows.append({
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "status": "VERIFIED_AGGREGATE",
                    "n": int(values.size),
                    "n_verified": int(values.size),
                    "bin": bin_name,
                    "n_in_bin": int(counts.get(bin_name, 0)),
                    "median": med,
                    "q1": q1,
                    "q3": q3,
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                    "value_source": ACTIVITY_SOURCE,
                })
    # neither context, same verification rule, not part of the inversion numerator.
    neither_verified = False
    neither_act = neither.merge(cross, on=["pair", "panel_id"], how="left", suffixes=("", "_cross"), indicator=True)
    neither_ok = neither_act["_merge"].eq("both")
    if neither_ok.all() and len(neither_act) == len(neither) and not neither_act.duplicated(["panel_id"]).any():
        value_ok_n = np.ones(len(neither_act), dtype=bool)
        for left, right in (
            ("dump_max_A", "historical_panel_pA"),
            ("historical_panel_pA", "historical_aggregate_max_A"),
            ("dump_max_B", "historical_panel_pB"),
            ("historical_panel_pB", "historical_aggregate_max_B"),
        ):
            value_ok_n &= neither_act[left].astype(float).to_numpy() == neither_act[right].astype(float).to_numpy()
        class_ok_n = (
            neither_act["frozen_class"].astype(str).eq(neither_act["class"].astype(str))
            & neither_act["class_matches_frozen"].eq("YES")
        )
        if bool((value_ok_n & class_ok_n.to_numpy()).all()):
            neither_verified = True
            values = neither_act["historical_panel_pA"].astype(float).to_numpy()
            q1, med, q3 = quartiles(values)
            counts = pd.Series(values).map(_activity_bin).value_counts()
            for bin_name in ("<5.0", "5.0_to_<5.5", "5.5_to_<6.0", ">=6.0"):
                activity_rows.append({
                    "cohort": "neither_context",
                    "endpoint": "EGFR",
                    "status": "VERIFIED_AGGREGATE",
                    "n": int(values.size),
                    "n_verified": int(values.size),
                    "bin": bin_name,
                    "n_in_bin": int(counts.get(bin_name, 0)),
                    "median": med,
                    "q1": q1,
                    "q3": q3,
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                    "value_source": ACTIVITY_SOURCE,
                })

    provenance = pd.read_csv(PROJECT / INPUTS["provenance"])
    flag_cols = [
        "cross_source_or_endpoint_mix",
        "replicate_threshold_discordance",
        "value_or_label_mismatch",
        "document_identifier_missing",
        "historical_exception",
        "literature_complete",
    ]
    flagged = pd.concat([
        primary.assign(analysis_cohort="D_vs_B_primary"),
        neither.assign(analysis_cohort="neither_context"),
    ], ignore_index=True)
    prov = flagged.merge(
        provenance,
        on=["pair", "panel_id"],
        how="left",
        suffixes=("", "_prov"),
        indicator=True,
    )
    if prov.duplicated(["analysis_cohort", "canonical_ligand_id"]).any() or (prov["_merge"] != "both").any():
        raise RuntimeError("provenance join failed or expanded rows")
    prov_rows = []
    for cohort, frame in prov.groupby("analysis_cohort"):
        for cls, sub in frame.groupby(frame["class"]):
            for col in flag_cols:
                yes = int(sub[col].astype(str).str.upper().eq("YES").sum())
                prov_rows.append({
                    "analysis_cohort": cohort,
                    "class": cls,
                    "flag": col,
                    "n_class": int(len(sub)),
                    "n_yes": yes,
                    "fraction_yes": yes / len(sub) if len(sub) else float("nan"),
                })
    join_audit = pd.DataFrame([
        {
            "join": "egfr_primary_to_D1",
            "left_n": int(len(primary)),
            "right_n": int(len(d1_arm)),
            "matched_n": int((joined["_merge"] == "both").sum()),
            "unmatched_left_n": n_unmatched,
            "duplicate_key_n": 0,
            "one_to_many": "NO",
            "note": "pair + canonical_ligand_id + class; arm fixed to D_vs_B",
        },
        {
            "join": "egfr_primary_to_activity_crosswalk_805",
            "left_n": int(len(primary)),
            "right_n": int((cross["pair"] == "EGFR/HER2").sum()),
            "matched_n": int((act["_merge"] == "both").sum()),
            "unmatched_left_n": int((act["_merge"] != "both").sum()),
            "duplicate_key_n": int(act.duplicated(["panel_id"]).sum()),
            "one_to_many": "NO",
            "note": ACTIVITY_SOURCE,
        },
        {
            "join": "egfr_primary_and_neither_to_provenance_867",
            "left_n": int(len(flagged)),
            "right_n": int(len(provenance)),
            "matched_n": int((prov["_merge"] == "both").sum()),
            "unmatched_left_n": int((prov["_merge"] != "both").sum()),
            "duplicate_key_n": 0,
            "one_to_many": "NO",
            "note": "pair + panel_id; class columns were not used as the join key",
        },
    ])
    # Molecule-level activity retained for the figure and audit, without replacing labels.
    molecule_activity = act[[
        "pair", "panel_id", "canonical_ligand_id", "class", "M0_score_A",
        "activity_value_status", "pchembl_EGFR", "pchembl_HER2",
    ]].copy()
    if neither_verified:
        neither_rows = neither_act[["pair", "panel_id", "canonical_ligand_id", "class", "M0_score_A"]].copy()
        neither_rows["activity_value_status"] = "VERIFIED_AGGREGATE"
        neither_rows["pchembl_EGFR"] = neither_act["historical_panel_pA"].astype(float)
        neither_rows["pchembl_HER2"] = neither_act["historical_panel_pB"].astype(float)
        molecule_activity = pd.concat([molecule_activity, neither_rows], ignore_index=True)
    molecule_activity["egfr_bin"] = molecule_activity["pchembl_EGFR"].map(
        lambda v: _activity_bin(float(v)) if pd.notna(v) else ""
    )
    return {
        "B_EGFR_SCORE_DISTRIBUTIONS": distributions,
        "B_EGFR_PAIRWISE_RANKING": pairwise,
        "B_EGFR_DESCRIPTOR_COMPARISONS": pd.DataFrame(desc_rows),
        "B_EGFR_SCORE_DESCRIPTOR_ASSOCIATIONS": pd.DataFrame(assoc_rows),
        "B_EGFR_ACTIVITY_DISTRIBUTIONS": pd.DataFrame(activity_rows),
        "B_EGFR_ACTIVITY_MOLECULES": molecule_activity,
        "B_EGFR_PROVENANCE_FLAGS": pd.DataFrame(prov_rows),
        "B_JOIN_AUDIT": join_audit,
        "B_SCORE_POINTS": pd.concat([
            primary[["canonical_ligand_id", "class", "M0_score_A"]],
            neither[["canonical_ligand_id", "class", "M0_score_A"]],
        ], ignore_index=True),
        "B_DESCRIPTOR_POINTS": joined[[
            "canonical_ligand_id", "class", "M0_score_A", *DESCRIPTORS
        ]].copy(),
    }
