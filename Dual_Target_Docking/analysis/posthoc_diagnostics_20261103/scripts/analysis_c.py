#!/usr/bin/env python3
"""Analysis C: Top-5/10/20% joint-list sensitivity after reproducing Top-10%."""
from __future__ import annotations

import math

import pandas as pd

from posthoc_common import (
    COMPARATORS,
    INPUTS,
    METHOD_TO_SCHEME,
    PROJECT,
    as_number,
    composition,
    concordance_code,
    joint_rank_rows,
    load_config,
    load_population_strings,
    near,
    pairwise_ranking_pool,
    ranking_pool,
    read_csv,
)
from jcim_stage2_compute import module_i
from jcim_stage2_lib import METHODS


def _pool_map(pop: list[dict]) -> dict[tuple[str, str, str], list[dict]]:
    """Return ranked lists keyed by pair, method, population_kind."""
    out = {}
    pairs = sorted({r["pair"] for r in pop})
    for pair in pairs:
        rows = [r for r in pop if r["pair"] == pair]
        for method in METHODS:
            out[(pair, method, "method_specific")] = joint_rank_rows(ranking_pool(rows, method), method)
        for method in COMPARATORS:
            common = pairwise_ranking_pool(rows, "M0", method)
            out[(pair, "M0", f"pairwise_common_vs_{method}")] = joint_rank_rows(common, "M0")
            out[(pair, method, "pairwise_common_vs_M0")] = joint_rank_rows(common, method)
    return out


def reproduce_topk() -> tuple[pd.DataFrame, dict]:
    pop = load_population_strings()
    primary = read_csv(PROJECT / INPUTS["method_delta"])
    desc, ligs, cmps = module_i(pop, primary)
    stored_desc = read_csv(PROJECT / INPUTS["descriptive"])
    stored_cmp = read_csv(PROJECT / INPUTS["method_compare"])
    stored_lig = read_csv(PROJECT / INPUTS["ligand_ranks"])
    pools = _pool_map(pop)
    checks = []

    def add_check(name: str, n_compared: int, n_mismatch: int, max_abs: float):
        checks.append({
            "check": name,
            "n_compared": n_compared,
            "n_mismatch": n_mismatch,
            "max_abs_numeric_diff": max_abs,
            "status": "PASS" if n_mismatch == 0 else "FAIL",
        })

    # Descriptive counts and fractions, including method-specific and pairwise pools.
    desc_key = lambda r: (r["pair"], r["method"], r["population_kind"])
    stored_desc_map = {desc_key(r): r for r in stored_desc}
    computed_desc_map = {desc_key(r): r for r in desc}
    if set(stored_desc_map) != set(computed_desc_map):
        add_check("descriptive_keys", len(stored_desc_map), 1, float("nan"))
    else:
        mismatches = 0
        max_abs = 0.0
        fields = (
            "n_total", "n_dual", "n_A_only", "n_B_only", "k",
            "n_dual_top", "n_A_only_top", "n_B_only_top",
            "top10_dual_fraction", "top10_single_target_fraction",
            "dual_recall", "expected_dual_fraction",
        )
        for key, stored in stored_desc_map.items():
            got = computed_desc_map[key]
            for field in fields:
                if not near(stored[field], got[field], tol=1e-12):
                    mismatches += 1
                else:
                    a, b = as_number(stored[field]), as_number(got[field])
                    if a is not None and b is not None:
                        max_abs = max(max_abs, abs(a - b))
        add_check("module_I_descriptive_vs_module_i", len(stored_desc_map) * len(fields), mismatches, max_abs)

    # Independent composition at 10% must match module_i / stored k and class counts.
    mismatches = 0
    compared = 0
    max_abs = 0.0
    for key, ranked in pools.items():
        stats = composition(ranked, 0.10)
        stored = stored_desc_map[key]
        compared += 1
        for field, value in (
            ("n_total", stats["n_total"]),
            ("k", stats["k"]),
            ("n_dual_top", stats["n_dual_top"]),
            ("n_A_only_top", stats["n_A_only_top"]),
            ("n_B_only_top", stats["n_B_only_top"]),
            ("top10_single_target_fraction", stats["single_target_fraction"]),
        ):
            if not near(stored[field], value, tol=1e-12):
                mismatches += 1
            else:
                max_abs = max(max_abs, abs(float(stored[field]) - float(value)))
    add_check("independent_ceil_10pct_vs_stored_descriptive", compared, mismatches, max_abs)

    stored_cmp_map = {(r["pair"], r["method"]): r for r in stored_cmp}
    computed_cmp_map = {(r["pair"], r["method"]): r for r in cmps}
    mismatches = 0
    max_abs = 0.0
    if set(stored_cmp_map) != set(computed_cmp_map) or len(stored_cmp_map) != 32:
        add_check("method_compare_keys", len(stored_cmp_map), 1, float("nan"))
    else:
        for key, stored in stored_cmp_map.items():
            got = computed_cmp_map[key]
            if stored["concordance"] != got["concordance"]:
                mismatches += 1
            for field in ("delta_summary_min", "delta_top10_single_target_fraction", "n_rank_common"):
                if not near(stored[field], got[field], tol=1e-12):
                    mismatches += 1
                else:
                    max_abs = max(max_abs, abs(float(stored[field]) - float(got[field])))
        add_check("module_I_method_compare_32", 32, mismatches, max_abs)

    # Four-state totals by internal method.
    order = ["CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"]
    expected = {
        "M1": (2, 2, 3, 1),
        "M1b": (3, 3, 1, 1),
        "M2": (5, 3, 0, 0),
        "M3": (4, 2, 2, 0),
    }
    mismatches = 0
    for method, counts in expected.items():
        got = tuple(sum(1 for r in stored_cmp if r["method"] == method and r["concordance"] == label) for label in order)
        if got != counts:
            mismatches += 1
    add_check("formal_four_state_counts", 4, mismatches, 0.0)

    # Member identity of stored ligand ranks versus recomputed pools.
    mismatches = 0
    compared = 0
    for key, ranked in pools.items():
        pair, method, kind = key
        stored_ids = [
            r["global_ligand_entity_id"]
            for r in stored_lig
            if r["pair"] == pair and r["method"] == method and r["population_kind"] == kind
        ]
        got_ids = [r["global_ligand_entity_id"] for r in ranked]
        compared += 1
        if stored_ids != got_ids:
            mismatches += 1
    add_check("ligand_rank_order_and_membership", compared, mismatches, 0.0)

    # Within each pairwise comparison, both schemes contain the same molecules.
    mismatches = 0
    compared = 0
    for pair in {r["pair"] for r in pop}:
        for method in COMPARATORS:
            left = {r["global_ligand_entity_id"] for r in pools[(pair, "M0", f"pairwise_common_vs_{method}")]}
            right = {r["global_ligand_entity_id"] for r in pools[(pair, method, "pairwise_common_vs_M0")]}
            compared += 1
            if left != right:
                mismatches += 1
    add_check("pairwise_common_member_identity", compared, mismatches, 0.0)
    report = pd.DataFrame(checks)
    if (report["status"] != "PASS").any():
        return report, {}
    return report, {"pop": pop, "pools": pools, "stored_cmp": stored_cmp_map, "desc": desc, "cmps": cmps}


def build_analysis_c() -> dict[str, pd.DataFrame]:
    reproduction, payload = reproduce_topk()
    if payload == {}:
        return {"C_TOPK_BASELINE_REPRODUCTION": reproduction}
    cfg = load_config()
    pools = payload["pools"]
    stored_cmp = payload["stored_cmp"]
    all_rows = []
    cmp_rows = []
    for fraction in cfg["topk_fractions"]:
        label = f"top_{int(round(fraction * 100))}pct"
        for pair in cfg["pair_order"]:
            for method in COMPARATORS:
                left = pools[(pair, "M0", f"pairwise_common_vs_{method}")]
                right = pools[(pair, method, "pairwise_common_vs_M0")]
                ids_left = [r["global_ligand_entity_id"] for r in left]
                ids_right = [r["global_ligand_entity_id"] for r in right]
                identical = ids_left != [] and set(ids_left) == set(ids_right) and len(ids_left) == len(ids_right)
                if not identical:
                    raise RuntimeError(f"common-pool mismatch {pair} {method}")
                stats_s1 = composition(left, fraction)
                stats_new = composition(right, fraction)
                if stats_s1["n_total"] != stats_new["n_total"] or stats_s1["k"] != stats_new["k"]:
                    raise RuntimeError("N or k differs within a comparison")
                if stats_s1["k"] != math.ceil(fraction * stats_s1["n_total"]):
                    raise RuntimeError("k rule violated")
                official = stored_cmp[(pair, method)]
                delta_summary = float(official["delta_summary_min"])
                delta_single = stats_new["single_target_fraction"] - stats_s1["single_target_fraction"]
                state = concordance_code(delta_summary, delta_single)
                if abs(fraction - 0.10) < 1e-15:
                    if not near(delta_single, official["delta_top10_single_target_fraction"], tol=1e-12):
                        raise RuntimeError(f"Top-10% single-target delta drifted for {pair} {method}")
                    if state != official["concordance"]:
                        raise RuntimeError(f"Top-10% concordance drifted for {pair} {method}")
                for role, stats in (("S1", stats_s1), (METHOD_TO_SCHEME[method], stats_new)):
                    all_rows.append({
                        "pair": pair,
                        "comparison": f"{METHOD_TO_SCHEME[method]}_vs_S1",
                        "internal_method": method if role != "S1" else "M0",
                        "list_scheme": role,
                        "cutoff": label,
                        "fraction": fraction,
                        "supporting_only": official["supporting_only"],
                        **stats,
                    })
                cmp_rows.append({
                    "pair": pair,
                    "comparison": f"{METHOD_TO_SCHEME[method]}_vs_S1",
                    "internal_method": method,
                    "cutoff": label,
                    "fraction": fraction,
                    "N": stats_s1["n_total"],
                    "k": stats_s1["k"],
                    "member_ids_identical": "YES",
                    "full_n_dual": stats_s1["n_dual"],
                    "full_n_A_only": stats_s1["n_A_only"],
                    "full_n_B_only": stats_s1["n_B_only"],
                    "baseline_dual_fraction": stats_s1["baseline_dual_fraction"],
                    "S1_n_dual_top": stats_s1["n_dual_top"],
                    "S1_n_A_only_top": stats_s1["n_A_only_top"],
                    "S1_n_B_only_top": stats_s1["n_B_only_top"],
                    "S1_dual_fraction": stats_s1["dual_fraction"],
                    "S1_single_target_fraction": stats_s1["single_target_fraction"],
                    "S1_dual_recall": stats_s1["dual_recall"],
                    "S1_enrichment_factor": stats_s1["enrichment_factor"],
                    "new_n_dual_top": stats_new["n_dual_top"],
                    "new_n_A_only_top": stats_new["n_A_only_top"],
                    "new_n_B_only_top": stats_new["n_B_only_top"],
                    "new_dual_fraction": stats_new["dual_fraction"],
                    "new_single_target_fraction": stats_new["single_target_fraction"],
                    "new_dual_recall": stats_new["dual_recall"],
                    "new_enrichment_factor": stats_new["enrichment_factor"],
                    "delta_dual_fraction": stats_new["dual_fraction"] - stats_s1["dual_fraction"],
                    "delta_single_target_fraction": delta_single,
                    "delta_dual_count": stats_new["n_dual_top"] - stats_s1["n_dual_top"],
                    "delta_summary_min": delta_summary,
                    "concordance": state,
                    "supporting_only": official["supporting_only"],
                    "summary_min_source": official["link_source"],
                    "new_confidence_interval": "NOT_COMPUTED",
                })
    comparisons = pd.DataFrame(cmp_rows)
    count_rows = []
    for cutoff, sub in comparisons.groupby("cutoff", sort=False):
        for scope, frame in (
            ("all_32", sub),
            ("exclude_S4", sub[sub["internal_method"] != "M1b"]),
        ):
            for state in ("CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"):
                count_rows.append({
                    "cutoff": cutoff,
                    "scope": scope,
                    "concordance": state,
                    "n": int((frame["concordance"] == state).sum()),
                    "n_comparisons": int(len(frame)),
                })
    reproduction = pd.concat([
        reproduction,
        pd.DataFrame([{
            "check": "sensitivity_accepted",
            "n_compared": int(len(comparisons)),
            "n_mismatch": 0,
            "max_abs_numeric_diff": 0.0,
            "status": "PASS",
        }]),
    ], ignore_index=True)
    return {
        "C_TOPK_SENSITIVITY_ALL": pd.DataFrame(all_rows),
        "C_TOPK_METHOD_COMPARISONS": comparisons,
        "C_TOPK_CONCORDANCE_COUNTS": pd.DataFrame(count_rows),
        "C_TOPK_BASELINE_REPRODUCTION": reproduction,
    }
