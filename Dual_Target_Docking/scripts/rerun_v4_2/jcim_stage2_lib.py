#!/usr/bin/env python3
"""Frozen Stage 2 helpers. Import must not train models or compute real AUROC.

Stage 2 reimplements Mann-Whitney AUROC / paired bootstrap / midrank here.
Do not import formal_metrics_lib locked functions.
Do not import analysis.bootstrap_metrics, analysis.uniform_protocol_config,
alt_score_master_and_stats, or alt_identity_test.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

BOOT_B = 10000
BOOT_SEED = 271828
CI_PERCENTILES = (2.5, 97.5)
CI_METHOD = "linear"
AUROC_TIE = 0.5
TOP_FRACTION = 0.10
LR_C = 1.0
LR_SOLVER = "liblinear"
LR_MAX_ITER = 1000
LR_RANDOM_STATE = 271828
ECFP_RADIUS = 2
ECFP_NBITS = 2048
SEEDS = (17, 29, 42, 71, 101)
METHODS = ("M0", "M1", "M1b", "M2", "M3")
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
ALT_PAIRS = ("PIK3CA/mTOR", "AChE/BChE", "F2/F10", "PPARA/PPARD")
ALT_PDBS = ("4L2Y", "4JT5", "4EY6", "1P0M", "3SHC", "2Y5F", "6KAX", "5U46")
BLOCKED_MODULES = ("G", "J", "K")
EXECUTABLE_MODULES = ("B", "C", "D", "E", "F", "H", "I")
PHYSCHEM_NAMES = (
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "FormalCharge",
)
FORBIDDEN_IMPORT_SUBSTRINGS = (
    "analysis.bootstrap_metrics",
    "analysis.uniform_protocol_config",
    "alt_score_master_and_stats",
    "alt_identity_test",
    "formal_metrics_lib",
)
NA = "NA"

RERUN = Path("reruns/UNIFORM_RERUN_V4_2_20260921")
PATHS = {
    "population": RERUN / "13_qa" / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv",
    "m0_master": RERUN / "13_qa" / "official_primary_seed42_score_master_8pair.csv",
    "mapping": RERUN / "00_protocol" / "pair_ligand_mapping.csv",
    "alias": RERUN / "00_protocol" / "parent_alias_registry.csv",
    "holdout": RERUN / "00_protocol" / "holdout_membership_freeze.csv",
    "folds": Path("results/canonical/model_fold_assignments.csv"),
    "fiveseed": RERUN / "13_qa" / "phase7_fiveseed_master.csv",
    "alt_independent": RERUN / "13_qa" / "alt_score_master_independent.csv",
    "primary_metrics": Path("results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv"),
    "primary_delta": Path("results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv"),
}

SCHEMA_B = (
    "module", "pair", "direction", "n_dual", "n_single", "n_neither",
    "auc_single", "auc_neither", "delta_negative",
    "auc_single_ci_lo", "auc_single_ci_hi",
    "auc_neither_ci_lo", "auc_neither_ci_hi",
    "delta_negative_ci_lo", "delta_negative_ci_hi",
    "reason", "bootstrap_started", "n_replicates",
)
SCHEMA_C = (
    "module", "pair", "arm", "n_dual", "n_negative_directional_m0",
    "n_negative_both_finite", "extra_drop_vs_directional_m0", "extra_drop_ids",
    "auc_correct", "auc_wrong", "delta_pocket",
    "auc_correct_ci_lo", "auc_correct_ci_hi",
    "auc_wrong_ci_lo", "auc_wrong_ci_hi",
    "delta_pocket_ci_lo", "delta_pocket_ci_hi",
    "reason", "bootstrap_started", "n_replicates",
)
SCHEMA_D_UNIT = (
    "module", "pair", "arm", "n_labeled_population", "n_m0_directional_population",
    "n_final_chemistry_oof", "fold_unassigned_in_labeled", "fold_unassigned_after_m0",
    "auroc_physchem", "auroc_ecfp4", "auroc_ecfp4_m0",
    "auroc_SUPPORTING_CHEMISTRY_MATCHED_M0", "delta_ecfp4_plus_m0",
    "test_single_class_folds", "blocking_reason", "evidence_status", "incomplete_reason",
)
SCHEMA_D_OOF = (
    "pair", "arm", "canonical_ligand_id", "class", "fold_id",
    "oof_physchem", "oof_ecfp4", "oof_ecfp4_m0", "relevant_m0",
)
SCHEMA_D1 = (
    "pair", "arm", "canonical_ligand_id", "class",
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "RingCount", "FormalCharge", "murcko",
    "max_ecfp4_tanimoto_other_class",
)
SCHEMA_E_SEED = (
    "module", "pair", "seed", "n_dual", "n_A_only", "n_B_only",
    "auc_B", "auc_A", "summary_min", "missing_jobs", "missing_ligands",
    "reason", "bootstrap_started",
)
SCHEMA_E_SUM = (
    "module", "pair", "metric", "seed42", "median", "min", "max", "iqr",
    "n_finite_seeds", "reason",
)
SCHEMA_F = (
    "module", "pair", "alt_id", "replaced_side",
    "n_dual_common", "n_A_only_common", "n_B_only_common",
    "common_dual_ids", "common_A_only_ids", "common_B_only_ids", "exclusion_reasons",
    "auc_B_primary", "auc_A_primary", "summary_min_primary",
    "auc_B_alt", "auc_A_alt", "summary_min_alt",
    "delta_AUC_affected", "delta_AUC_unreplaced", "delta_summary_min",
    "delta_AUC_affected_ci_lo", "delta_AUC_affected_ci_hi",
    "delta_AUC_unreplaced_ci_lo", "delta_AUC_unreplaced_ci_hi",
    "delta_summary_min_ci_lo", "delta_summary_min_ci_hi",
    "reason", "bootstrap_started", "n_replicates",
)
SCHEMA_H = (
    "module", "pair", "method", "class", "target",
    "expected", "valid", "missing", "missing_reason",
)
SCHEMA_I_DESC = (
    "module", "table", "pair", "method", "population_kind",
    "n_total", "n_dual", "n_A_only", "n_B_only", "k",
    "n_dual_top", "n_A_only_top", "n_B_only_top",
    "top10_dual_fraction", "top10_single_target_fraction", "dual_recall",
    "expected_dual_count", "expected_dual_fraction",
)
SCHEMA_I_LIG = (
    "module", "table", "pair", "method", "population_kind",
    "canonical_ligand_id", "global_ligand_entity_id", "class",
    "rank_A", "rank_B", "worst_rank", "joint_order", "in_topk",
)
SCHEMA_I_CMP = (
    "module", "table", "pair", "method",
    "n_rank_common", "n_directional_common_dual",
    "rank_vs_directional_member_diff",
    "delta_summary_min", "delta_top10_single_target_fraction",
    "concordance", "link_source", "link_reason",
    "supporting_only", "bootstrap_ci_added",
)


def resolve_project_root(explicit: str | Path | None = None) -> Path:
    if explicit:
        root = Path(explicit).resolve()
    else:
        root = Path(__file__).resolve().parents[2]
    if root.name != "Dual_Target_Docking":
        raise ValueError(f"PROJECT_ROOT must be Dual_Target_Docking, got {root}")
    return root


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def empty_record(fieldnames: tuple[str, ...], **values) -> dict:
    rec = {key: NA if key not in {"module", "table", "pair", "arm", "direction", "method", "seed"} else "" for key in fieldnames}
    rec.update(values)
    return rec


def finite(value) -> bool:
    try:
        return value not in (None, "", NA) and np.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + AUROC_TIE * (diff == 0).sum()) / (pos.size * neg.size))


def require_all_finite(values, label: str) -> list[float]:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        raise ValueError(f"{label}:empty_replicate_array")
    bad = int((~np.isfinite(arr)).sum())
    if bad:
        raise ValueError(f"{label}:nonfinite_replicates={bad}/{arr.size}")
    return [float(x) for x in arr]


def percentile_ci(values, n_boot: int = BOOT_B) -> tuple[float, float]:
    arr = np.asarray(require_all_finite(values, "percentile_ci"), dtype=float)
    if arr.size != n_boot:
        raise ValueError(f"percentile_ci:n_replicates={arr.size}!={n_boot}")
    lo, hi = np.percentile(arr, list(CI_PERCENTILES), method=CI_METHOD)
    if not (np.isfinite(lo) and np.isfinite(hi)):
        raise ValueError("percentile_ci:nonfinite_bounds")
    return float(lo), float(hi)


def five_seed_iqr(values) -> tuple[float, float, float, float, float]:
    arr = np.asarray(list(values), dtype=float)
    if arr.size != 5 or not np.isfinite(arr).all():
        raise ValueError("five_seed_iqr_requires_five_finite")
    p25, p75 = np.percentile(arr, [25, 75], method="linear")
    return float(np.median(arr)), float(np.min(arr)), float(np.max(arr)), float(p75 - p25)


def scores(rows: list[dict], key: str) -> list[float]:
    out = []
    for row in rows:
        if not finite(row.get(key)):
            raise ValueError(f"valid_row_nonfinite_score:{key}")
        out.append(float(row[key]))
    return out


def make_rng(seed: int = BOOT_SEED):
    return np.random.default_rng(seed)


def shared_dual_universe(rows: list[dict], method: str) -> dict:
    dual, a_only, b_only = [], [], []
    for row in rows:
        if row.get("activity_eligible") != "1":
            continue
        cls = row.get("class")
        va = row.get(f"{method}_A_valid") == "1"
        vb = row.get(f"{method}_B_valid") == "1"
        if cls == "dual" and va and vb:
            dual.append(row)
        elif cls == "A_only" and vb:
            a_only.append(row)
        elif cls == "B_only" and va:
            b_only.append(row)
    return {
        "method": method,
        "dual_shared": dual,
        "A_only_for_AUC_B": a_only,
        "B_only_for_AUC_A": b_only,
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
    }


def four_sided_common_complete(rows: list[dict], method_x: str, method_y: str) -> dict:
    dual, a_only, b_only = [], [], []
    for row in rows:
        if row.get("activity_eligible") != "1":
            continue
        xa = row.get(f"{method_x}_A_valid") == "1"
        xb = row.get(f"{method_x}_B_valid") == "1"
        ya = row.get(f"{method_y}_A_valid") == "1"
        yb = row.get(f"{method_y}_B_valid") == "1"
        cls = row.get("class")
        if cls == "dual" and xa and xb and ya and yb:
            dual.append(row)
        elif cls == "A_only" and xb and yb:
            a_only.append(row)
        elif cls == "B_only" and xa and ya:
            b_only.append(row)
    return {
        "method_x": method_x,
        "method_y": method_y,
        "dual_four_sided": dual,
        "A_only_both_score_B": a_only,
        "B_only_both_score_A": b_only,
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
    }


def member_keys(rows: list[dict]) -> set[tuple[str, str, str]]:
    return {(r.get("pair", ""), r["canonical_ligand_id"], r["class"]) for r in rows}


def four_sided_membership_sets(rows: list[dict], method_x: str, method_y: str) -> dict:
    pop = four_sided_common_complete(rows, method_x, method_y)
    return {
        "dual": member_keys(pop["dual_four_sided"]),
        "A_only": member_keys(pop["A_only_both_score_B"]),
        "B_only": member_keys(pop["B_only_both_score_A"]),
        "pop": pop,
    }


def four_sided_membership_identical(rows: list[dict], left: str, right: str, vs: str = "M0") -> bool:
    a = four_sided_membership_sets(rows, left, vs)
    b = four_sided_membership_sets(rows, right, vs)
    return a["dual"] == b["dual"] and a["A_only"] == b["A_only"] and a["B_only"] == b["B_only"]


def m1b_link_membership(rows: list[dict]) -> dict:
    """Compare M1b-M0, M1-M0, and M1b-M1 class sets separately.

    Adding PRIMARY(M1b-M1)+PRIMARY(M1-M0) is allowed only if all three
    four-sided populations have identical dual, A_only, and B_only members.
    """
    m1b_m0 = four_sided_membership_sets(rows, "M1b", "M0")
    m1_m0 = four_sided_membership_sets(rows, "M1", "M0")
    m1b_m1 = four_sided_membership_sets(rows, "M1b", "M1")
    diffs = []
    identical = True
    for grp in ("dual", "A_only", "B_only"):
        sets = (m1b_m0[grp], m1_m0[grp], m1b_m1[grp])
        if not (sets[0] == sets[1] == sets[2]):
            identical = False
            xor = (sets[0] ^ sets[1]) | (sets[0] ^ sets[2]) | (sets[1] ^ sets[2])
            diffs.append(f"{grp}:{sorted(xor)[:8]}")
    return {
        "identical": identical,
        "diffs": diffs,
        "M1b-M0": m1b_m0,
        "M1-M0": m1_m0,
        "M1b-M1": m1b_m1,
    }


def paired_negative_class_bootstrap(dual, single, neither, score_key: str,
                                    n_boot: int = BOOT_B, seed: int = BOOT_SEED,
                                    rng=None, record=None) -> dict:
    """delta = AUROC(dual vs neither) - AUROC(dual vs single). Shared dual index."""
    if not dual:
        return {"reason": "INSUFFICIENT_n_dual", "bootstrap_started": "NO"}
    if not single and not neither:
        return {"reason": "EMPTY_BOTH_NEGATIVE_CLASSES", "bootstrap_started": "NO"}
    auc_s = auroc(scores(dual, score_key), scores(single, score_key)) if single else float("nan")
    auc_n = auroc(scores(dual, score_key), scores(neither, score_key)) if neither else float("nan")
    out = {
        "auc_single": auc_s if single else NA,
        "auc_neither": auc_n if neither else NA,
        "delta_negative": (auc_n - auc_s) if (single and neither) else NA,
        "reason": "",
        "bootstrap_started": "NO",
        "n_replicates": 0,
    }
    if not single:
        out["reason"] = "EMPTY_SINGLE_TARGET_CLASS"
        return out
    if not neither:
        out["reason"] = "EMPTY_NEITHER_CLASS"
        return out
    rng = rng or make_rng(seed)
    s_reps, n_reps, d_reps = [], [], []
    for _ in range(n_boot):
        d_idx = rng.integers(0, len(dual), size=len(dual))
        s_idx = rng.integers(0, len(single), size=len(single))
        n_idx = rng.integers(0, len(neither), size=len(neither))
        if record is not None:
            record.append(("dual", tuple(int(i) for i in d_idx), "single", tuple(int(i) for i in s_idx), "neither", tuple(int(i) for i in n_idx)))
        d_s = [dual[i] for i in d_idx]
        a_s = auroc(scores(d_s, score_key), scores([single[i] for i in s_idx], score_key))
        a_n = auroc(scores(d_s, score_key), scores([neither[i] for i in n_idx], score_key))
        s_reps.append(a_s)
        n_reps.append(a_n)
        d_reps.append(a_n - a_s)
    slo, shi = percentile_ci(s_reps, n_boot=n_boot)
    nlo, nhi = percentile_ci(n_reps, n_boot=n_boot)
    dlo, dhi = percentile_ci(d_reps, n_boot=n_boot)
    out.update({
        "auc_single_ci_lo": slo, "auc_single_ci_hi": shi,
        "auc_neither_ci_lo": nlo, "auc_neither_ci_hi": nhi,
        "delta_negative_ci_lo": dlo, "delta_negative_ci_hi": dhi,
        "bootstrap_started": "YES",
        "n_replicates": n_boot,
        "delta_reps": d_reps,
    })
    return out


def paired_pocket_bootstrap(dual, negs, correct_key: str, wrong_key: str,
                            n_boot: int = BOOT_B, seed: int = BOOT_SEED,
                            rng=None, record=None) -> dict:
    if not dual or not negs:
        return {"reason": "INSUFFICIENT", "bootstrap_started": "NO"}
    auc_c = auroc(scores(dual, correct_key), scores(negs, correct_key))
    auc_w = auroc(scores(dual, wrong_key), scores(negs, wrong_key))
    rng = rng or make_rng(seed)
    c_reps, w_reps, d_reps = [], [], []
    for _ in range(n_boot):
        d_idx = rng.integers(0, len(dual), size=len(dual))
        n_idx = rng.integers(0, len(negs), size=len(negs))
        if record is not None:
            record.append((tuple(int(i) for i in d_idx), tuple(int(i) for i in n_idx)))
        d_s = [dual[i] for i in d_idx]
        n_s = [negs[i] for i in n_idx]
        a_c = auroc(scores(d_s, correct_key), scores(n_s, correct_key))
        a_w = auroc(scores(d_s, wrong_key), scores(n_s, wrong_key))
        c_reps.append(a_c)
        w_reps.append(a_w)
        d_reps.append(a_c - a_w)
    clo, chi = percentile_ci(c_reps, n_boot=n_boot)
    wlo, whi = percentile_ci(w_reps, n_boot=n_boot)
    dlo, dhi = percentile_ci(d_reps, n_boot=n_boot)
    return {
        "auc_correct": auc_c, "auc_wrong": auc_w, "delta_pocket": auc_c - auc_w,
        "auc_correct_ci_lo": clo, "auc_correct_ci_hi": chi,
        "auc_wrong_ci_lo": wlo, "auc_wrong_ci_hi": whi,
        "delta_pocket_ci_lo": dlo, "delta_pocket_ci_hi": dhi,
        "reason": "", "bootstrap_started": "YES", "n_replicates": n_boot,
        "delta_reps": d_reps,
    }


def paired_alt_bootstrap(pri_rows, alt_rows, replaced_side: str,
                         n_boot: int = BOOT_B, seed: int = BOOT_SEED,
                         rng=None) -> dict:
    pri = shared_dual_universe(pri_rows, "M0")
    alt = shared_dual_universe(alt_rows, "M0")
    if member_keys(pri["dual_shared"]) != member_keys(alt["dual_shared"]):
        raise ValueError("alt_primary_dual_membership_mismatch")
    if member_keys(pri["A_only_for_AUC_B"]) != member_keys(alt["A_only_for_AUC_B"]):
        raise ValueError("alt_primary_A_only_membership_mismatch")
    if member_keys(pri["B_only_for_AUC_A"]) != member_keys(alt["B_only_for_AUC_A"]):
        raise ValueError("alt_primary_B_only_membership_mismatch")
    dual, ao, bo = pri["dual_shared"], pri["A_only_for_AUC_B"], pri["B_only_for_AUC_A"]
    if not (dual and ao and bo):
        return {"reason": "INSUFFICIENT_COMMON_SIDES", "bootstrap_started": "NO"}

    def sm(rows_u):
        u = shared_dual_universe(rows_u, "M0")
        auc_b = auroc(scores(u["dual_shared"], "M0_score_B"), scores(u["A_only_for_AUC_B"], "M0_score_B"))
        auc_a = auroc(scores(u["dual_shared"], "M0_score_A"), scores(u["B_only_for_AUC_A"], "M0_score_A"))
        return auc_b, auc_a, min(auc_b, auc_a)

    pb, pa, psm = sm(pri_rows)
    ab, aa, asm = sm(alt_rows)
    d_aff_pt = (ab - pb) if replaced_side == "B" else (aa - pa)
    d_un_pt = (aa - pa) if replaced_side == "B" else (ab - pb)
    d_sm_pt = asm - psm
    rng = rng or make_rng(seed)
    aff_reps, un_reps, sm_reps = [], [], []
    alt_by_id = {r["canonical_ligand_id"]: r for r in alt_rows}

    def take(src, idx):
        return [src[i] for i in idx]

    def take_alt(src, idx):
        return [alt_by_id[src[i]["canonical_ligand_id"]] for i in idx]

    for _ in range(n_boot):
        d_idx = rng.integers(0, len(dual), size=len(dual))
        a_idx = rng.integers(0, len(ao), size=len(ao))
        b_idx = rng.integers(0, len(bo), size=len(bo))
        pd, pao, pbo = take(dual, d_idx), take(ao, a_idx), take(bo, b_idx)
        ad, aao, abo = take_alt(dual, d_idx), take_alt(ao, a_idx), take_alt(bo, b_idx)
        p_auc_b = auroc(scores(pd, "M0_score_B"), scores(pao, "M0_score_B"))
        p_auc_a = auroc(scores(pd, "M0_score_A"), scores(pbo, "M0_score_A"))
        a_auc_b = auroc(scores(ad, "M0_score_B"), scores(aao, "M0_score_B"))
        a_auc_a = auroc(scores(ad, "M0_score_A"), scores(abo, "M0_score_A"))
        d_b, d_a = a_auc_b - p_auc_b, a_auc_a - p_auc_a
        aff_reps.append(d_b if replaced_side == "B" else d_a)
        un_reps.append(d_a if replaced_side == "B" else d_b)
        sm_reps.append(min(a_auc_b, a_auc_a) - min(p_auc_b, p_auc_a))
    alo, ahi = percentile_ci(aff_reps, n_boot=n_boot)
    ulo, uhi = percentile_ci(un_reps, n_boot=n_boot)
    slo, shi = percentile_ci(sm_reps, n_boot=n_boot)
    return {
        "auc_B_primary": pb, "auc_A_primary": pa, "summary_min_primary": psm,
        "auc_B_alt": ab, "auc_A_alt": aa, "summary_min_alt": asm,
        "delta_AUC_affected": d_aff_pt, "delta_AUC_unreplaced": d_un_pt, "delta_summary_min": d_sm_pt,
        "delta_AUC_affected_ci_lo": alo, "delta_AUC_affected_ci_hi": ahi,
        "delta_AUC_unreplaced_ci_lo": ulo, "delta_AUC_unreplaced_ci_hi": uhi,
        "delta_summary_min_ci_lo": slo, "delta_summary_min_ci_hi": shi,
        "reason": "", "bootstrap_started": "YES", "n_replicates": n_boot,
    }


def midrank_desc(values) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(-x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    return ranks


def ranking_pool(rows: list[dict], method: str) -> list[dict]:
    return [
        r
        for r in rows
        if r.get("activity_eligible") == "1"
        and r.get("class") in {"dual", "A_only", "B_only"}
        and r.get(f"{method}_A_valid") == "1"
        and r.get(f"{method}_B_valid") == "1"
    ]


def pairwise_ranking_pool(rows: list[dict], method_x: str, method_y: str) -> list[dict]:
    return [
        r
        for r in rows
        if r.get("activity_eligible") == "1"
        and r.get("class") in {"dual", "A_only", "B_only"}
        and r.get(f"{method_x}_A_valid") == "1"
        and r.get(f"{method_x}_B_valid") == "1"
        and r.get(f"{method_y}_A_valid") == "1"
        and r.get(f"{method_y}_B_valid") == "1"
    ]


def joint_rank_rows(pool: list[dict], method: str) -> list[dict]:
    if not pool:
        return []
    rank_a = midrank_desc([float(r[f"{method}_score_A"]) for r in pool])
    rank_b = midrank_desc([float(r[f"{method}_score_B"]) for r in pool])
    enriched = []
    for row, ra, rb in zip(pool, rank_a, rank_b):
        item = dict(row)
        item["rank_A"] = float(ra)
        item["rank_B"] = float(rb)
        item["worst_rank"] = max(float(ra), float(rb))
        item["tie1"] = float(ra) + float(rb)
        item["tie2"] = row["global_ligand_entity_id"]
        enriched.append(item)
    enriched.sort(key=lambda r: (r["worst_rank"], r["tie1"], r["tie2"]))
    for i, row in enumerate(enriched, start=1):
        row["joint_order"] = i
    return enriched


def topk_composition(ranked: list[dict]) -> dict:
    n = len(ranked)
    k = int(math.ceil(TOP_FRACTION * n)) if n else 0
    top = ranked[:k]
    n_dual = sum(1 for r in top if r["class"] == "dual")
    n_a = sum(1 for r in top if r["class"] == "A_only")
    n_b = sum(1 for r in top if r["class"] == "B_only")
    n_dual_all = sum(1 for r in ranked if r["class"] == "dual")
    n_a_all = sum(1 for r in ranked if r["class"] == "A_only")
    n_b_all = sum(1 for r in ranked if r["class"] == "B_only")
    return {
        "n_total": n,
        "n_dual": n_dual_all,
        "n_A_only": n_a_all,
        "n_B_only": n_b_all,
        "k": k,
        "n_dual_top": n_dual,
        "n_A_only_top": n_a,
        "n_B_only_top": n_b,
        "top10_dual_fraction": n_dual / k if k else NA,
        "top10_single_target_fraction": (n_a + n_b) / k if k else NA,
        "dual_recall": n_dual / n_dual_all if n_dual_all else NA,
        "expected_dual_count": k * n_dual_all / n if n else NA,
        "expected_dual_fraction": n_dual_all / n if n else NA,
    }


def concordance_code(delta_summary_min: float, delta_top10_stf: float) -> str:
    if delta_summary_min == 0.0 or delta_top10_stf == 0.0:
        return "NO_CHANGE"
    if delta_summary_min > 0.0 and delta_top10_stf < 0.0:
        return "CONCORDANT_IMPROVE"
    if delta_summary_min < 0.0 and delta_top10_stf > 0.0:
        return "CONCORDANT_WORSEN"
    return "DISCORDANT"


def alias_ids(alias_rows: list[dict]) -> set[str]:
    return {r["alias_id"] for r in alias_rows}


def fold_keyset(fold_rows: list[dict]) -> set[tuple[str, str, str]]:
    return {(r["pair"], r["arm"], r["ligand_id"]) for r in fold_rows}


def directional_m0_ok(row: dict, arm: str) -> bool:
    if row["class"] == "dual":
        return row.get("M0_A_valid") == "1" and row.get("M0_B_valid") == "1"
    if arm == "D_vs_A":
        return row.get("M0_B_valid") == "1"
    return row.get("M0_A_valid") == "1"


def build_layer_members(pop_rows, alias_set, fold_keys, pair: str, arm: str) -> dict:
    classes = {"dual", "A_only"} if arm == "D_vs_A" else {"dual", "B_only"}
    labeled = [
        r
        for r in pop_rows
        if r["pair"] == pair
        and r.get("activity_eligible") == "1"
        and r.get("canonical_ligand_id") not in alias_set
        and r.get("class") in classes
    ]
    m0 = [r for r in labeled if directional_m0_ok(r, arm)]
    chem = [r for r in m0 if (pair, arm, r["canonical_ligand_id"]) in fold_keys]
    return {"labeled": labeled, "m0": m0, "chemistry_oof": chem}


def mapping_smiles(mapping_rows, alias_set) -> dict[tuple[str, str], str]:
    out = {}
    for row in mapping_rows:
        lig = row["canonical_ligand_id"]
        if lig in alias_set or str(row.get("parent_independent", "1")) in {"0", "false", "False"}:
            continue
        out[(row["pair"], lig)] = row.get("canonical_smiles") or ""
    return out


def load_stage2_tables(project_root: Path) -> dict:
    tables = {}
    for name, rel in PATHS.items():
        path = project_root / rel
        tables[name] = read_csv(path) if path.exists() else None
        tables[f"{name}_path"] = path
    return tables


def synthetic_self_check() -> list[str]:
    """Tiny synthetic checks for primitives only. Not module acceptance."""
    issues = []
    if abs(auroc([1, 2, 3], [0, 0, 0]) - 1.0) > 1e-12:
        issues.append("auroc_perfect")
    if abs(auroc([1, 1], [1, 1]) - 0.5) > 1e-12:
        issues.append("auroc_tie")
    ranks = midrank_desc([3.0, 3.0, 1.0])
    if list(ranks) != [1.5, 1.5, 3.0]:
        issues.append(f"midrank={list(ranks)}")
    if concordance_code(0.1, -0.2) != "CONCORDANT_IMPROVE":
        issues.append("improve")
    if concordance_code(-0.1, 0.2) != "CONCORDANT_WORSEN":
        issues.append("worsen")
    if concordance_code(0.1, 0.2) != "DISCORDANT":
        issues.append("discordant")
    if concordance_code(0.0, -0.2) != "NO_CHANGE":
        issues.append("zero")
    if int(math.ceil(TOP_FRACTION * 11)) != 2:
        issues.append("topk_ceil")
    med, mn, mx, iqr = five_seed_iqr([1.0, 2.0, 3.0, 4.0, 5.0])
    p25, p75 = np.percentile([1.0, 2.0, 3.0, 4.0, 5.0], [25, 75], method="linear")
    if abs(iqr - (p75 - p25)) > 1e-12 or med != 3.0 or mn != 1.0 or mx != 5.0:
        issues.append("iqr")
    return issues
