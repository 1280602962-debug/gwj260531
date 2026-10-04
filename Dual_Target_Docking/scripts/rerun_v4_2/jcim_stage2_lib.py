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


def finite(value) -> bool:
    try:
        return value not in (None, "") and np.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def auroc(pos, neg) -> float:
    """Mann-Whitney AUROC; tie = 0.5. For compute / synthetic tests only."""
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


def scores(rows: list[dict], key: str) -> list[float]:
    out = []
    for row in rows:
        if not finite(row.get(key)):
            raise ValueError(f"valid_row_nonfinite_score:{key}")
        out.append(float(row[key]))
    return out


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


def membership_ids(pop: dict, dual_key: str, a_key: str, b_key: str) -> tuple[set, set, set]:
    def keys(rows):
        return {(r["pair"], r["canonical_ligand_id"], r["class"]) for r in rows}

    return keys(pop[dual_key]), keys(pop[a_key]), keys(pop[b_key])


def four_sided_membership_identical(rows: list[dict], left: str, right: str, vs: str = "M0") -> bool:
    a = four_sided_common_complete(rows, left, vs)
    b = four_sided_common_complete(rows, right, vs)
    return membership_ids(a, "dual_four_sided", "A_only_both_score_B", "B_only_both_score_A") == membership_ids(
        b, "dual_four_sided", "A_only_both_score_B", "B_only_both_score_A"
    )


def single_method_point_and_bootstrap(rows: list[dict], method: str,
                                      n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    universe = shared_dual_universe(rows, method)
    dual = universe["dual_shared"]
    ao = universe["A_only_for_AUC_B"]
    bo = universe["B_only_for_AUC_A"]
    do_b = bool(dual and ao)
    do_a = bool(dual and bo)
    do_sm = do_b and do_a
    auc_b_pt = auroc(scores(dual, f"{method}_score_B"), scores(ao, f"{method}_score_B")) if do_b else None
    auc_a_pt = auroc(scores(dual, f"{method}_score_A"), scores(bo, f"{method}_score_A")) if do_a else None
    sm_pt = min(auc_b_pt, auc_a_pt) if do_sm else None
    auc_b_reps, auc_a_reps, sm_reps = [], [], []
    if do_b or do_a:
        rng = np.random.default_rng(seed)
        for _ in range(n_boot):
            d_s = [dual[i] for i in rng.integers(0, len(dual), size=len(dual))]
            auc_b = auc_a = None
            if do_b:
                a_s = [ao[i] for i in rng.integers(0, len(ao), size=len(ao))]
                auc_b = auroc(scores(d_s, f"{method}_score_B"), scores(a_s, f"{method}_score_B"))
                auc_b_reps.append(auc_b)
            if do_a:
                b_s = [bo[i] for i in rng.integers(0, len(bo), size=len(bo))]
                auc_a = auroc(scores(d_s, f"{method}_score_A"), scores(b_s, f"{method}_score_A"))
                auc_a_reps.append(auc_a)
            if do_sm:
                sm_reps.append(min(auc_b, auc_a))
    return {
        "method": method,
        "n_dual": len(dual),
        "n_A_only": len(ao),
        "n_B_only": len(bo),
        "dual_vs_A_only": (auc_b_pt, auc_b_reps) if do_b else None,
        "dual_vs_B_only": (auc_a_pt, auc_a_reps) if do_a else None,
        "summary_min": (sm_pt, sm_reps) if do_sm else None,
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


def joint_rank_rows(rows: list[dict], method: str) -> list[dict]:
    pool = [
        r
        for r in rows
        if r.get("activity_eligible") == "1"
        and r.get("class") in {"dual", "A_only", "B_only"}
        and r.get(f"{method}_A_valid") == "1"
        and r.get(f"{method}_B_valid") == "1"
    ]
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
    return enriched


def topk_composition(ranked: list[dict]) -> dict:
    n = len(ranked)
    k = int(math.ceil(TOP_FRACTION * n))
    top = ranked[:k]
    n_dual = sum(1 for r in top if r["class"] == "dual")
    n_single = sum(1 for r in top if r["class"] in {"A_only", "B_only"})
    n_dual_all = sum(1 for r in ranked if r["class"] == "dual")
    return {
        "n_total": n,
        "k": k,
        "n_dual_top": n_dual,
        "n_single_top": n_single,
        "top10_dual_fraction": n_dual / k if k else float("nan"),
        "top10_single_target_fraction": n_single / k if k else float("nan"),
        "dual_recall": n_dual / n_dual_all if n_dual_all else float("nan"),
        "expected_dual_count": k * n_dual_all / n if n else float("nan"),
        "expected_dual_fraction": n_dual_all / n if n else float("nan"),
    }


def concordance_code(delta_summary_min: float, delta_top10_stf: float) -> str:
    if delta_summary_min == 0.0 or delta_top10_stf == 0.0:
        return "NO_CHANGE"
    if delta_summary_min > 0.0 and delta_top10_stf < 0.0:
        return "CONCORDANT"
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


def required_m0_key(row: dict, arm: str) -> str:
    if arm == "D_vs_A":
        return "M0_score_B"
    return "M0_score_A"


def build_layer_members(pop_rows, alias_set, fold_keys, pair: str, arm: str) -> dict:
    """Compute-time layer builder. Audit must not treat this as independent proof."""
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


def load_stage2_tables(project_root: Path) -> dict:
    tables = {}
    for name, rel in PATHS.items():
        path = project_root / rel
        tables[name] = read_csv(path) if path.exists() else None
        tables[f"{name}_path"] = path
    return tables


def synthetic_self_check() -> list[str]:
    """Tiny synthetic checks only. No project ligands."""
    issues = []
    if abs(auroc([1, 2, 3], [0, 0, 0]) - 1.0) > 1e-12:
        issues.append("auroc_perfect")
    if abs(auroc([1, 1], [1, 1]) - 0.5) > 1e-12:
        issues.append("auroc_tie")
    if abs(auroc([0, 1], [0, 1]) - 0.5) > 1e-12:
        issues.append("auroc_identical")
    ranks = midrank_desc([3.0, 3.0, 1.0])
    if list(ranks) != [1.5, 1.5, 3.0]:
        issues.append(f"midrank={list(ranks)}")
    if concordance_code(0.1, -0.2) != "CONCORDANT":
        issues.append("concordant")
    if concordance_code(0.0, -0.2) != "NO_CHANGE":
        issues.append("zero_not_no_change")
    if concordance_code(-0.1, -0.2) != "DISCORDANT":
        issues.append("discordant")
    if int(math.ceil(TOP_FRACTION * 11)) != 2:
        issues.append("topk_ceil")
    rng = np.random.default_rng(BOOT_SEED)
    reps = [float(rng.random()) for _ in range(BOOT_B)]
    lo, hi = percentile_ci(reps)
    if not (0.0 <= lo <= hi <= 1.0):
        issues.append("percentile_ci")
    return issues
