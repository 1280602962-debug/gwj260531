#!/usr/bin/env python3
"""Formal metrics implementation for V4.2.

Source METRICS_EXECUTION_UNLOCKED must remain False.
Only the PRIMARY runner / POST QA may set the module flag True in-process.

This module must not import scripts/analysis/*.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np

METRICS_EXECUTION_UNLOCKED = False
BOOT_B = 10000
BOOT_SEED = 271828
THETA_PRIMARY = 6.0
SCORE_DIRECTION = "higher_better"
D_VS_A_USES = "score_B"  # dual vs A_only
D_VS_B_USES = "score_A"  # dual vs B_only
PRIMARY_SEED = 42
FORMAL_METHODS = ("M0", "M1", "M1b", "M2", "M3")
FORBIDDEN_METHODS = ("OLD_M1", "OLD_M3", "LEGACY")
POINT_ESTIMATE = "original_sample_not_bootstrap_mean"
CI = "percentile_2.5_97.5"
CI_PERCENTILES = (2.5, 97.5)
CI_METHOD = "linear"
AUROC_TIE = 0.5
DELTA_SIGN = "summary_min(A)-summary_min(B)"
RNG_POLICY = "reset_default_rng_per_analysis_unit"
PRIMARY_DELTAS = (("M1", "M0"), ("M1b", "M1"), ("M2", "M0"), ("M3", "M0"))
SECONDARY_DELTAS = (("M3", "M1b"), ("M3", "M1"))  # not computed in PRIMARY
NAN_POLICY = "no_silent_drop"

# Merge keys (pair observation level)
PAIR_OBS_KEY = ("pair", "global_ligand_entity_id")
PHYSICAL_KEY = ("pdb_id", "global_ligand_entity_id")
POSE_KEY = ("pdb_id", "global_ligand_entity_id", "vina_mode")

RECEPTORS = {
    "EGFR/HER2": ("3POZ", "3RCD"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "PIK3CA/mTOR": ("4L23", "4JT6"),
    "AChE/BChE": ("4EY7", "4BDS"),
    "F2/F10": ("4UDW", "2JKH"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}


def _locked(name: str) -> None:
    if not METRICS_EXECUTION_UNLOCKED:
        raise RuntimeError(
            f"{name} is locked until METRICS_IMPLEMENTATION_GATE=PASS "
            "and a human unlock of METRICS_EXECUTION_UNLOCKED"
        )


def auroc(pos, neg) -> float:
    """Mann-Whitney AUROC. LOCKED — do not call from readiness audit."""
    _locked("auroc")
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + AUROC_TIE * (diff == 0).sum()) / (pos.size * neg.size))


def finite(x) -> bool:
    try:
        return x not in (None, "") and np.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def method_valid(method: str, rec: dict) -> bool:
    if method == "M0":
        return rec.get("vina_status") == "SUCCESS" and finite(rec.get("vina_score"))
    if method == "M1":
        return rec.get("job_status") == "SUCCESS" and finite(rec.get("M1_CNNscore"))
    if method == "M1b":
        return rec.get("job_status") == "SUCCESS" and finite(rec.get("M1b_CNNscore"))
    if method == "M2":
        return rec.get("job_status") == "SUCCESS" and finite(rec.get("M2_RTMScore"))
    if method == "M3":
        return (
            rec.get("status") == "SUCCESS"
            and rec.get("official_m3", "YES") != "NO"
            and finite(rec.get("CNNscore"))
        )
    return False


def method_score(method: str, rec: dict) -> str:
    if method == "M0":
        return rec.get("vina_score", "")
    if method == "M1":
        return rec.get("M1_CNNscore", "")
    if method == "M1b":
        return rec.get("M1b_CNNscore", "")
    if method == "M2":
        return rec.get("M2_RTMScore", "")
    if method == "M3":
        return rec.get("CNNscore", "")
    return ""


def method_missing_reason(method: str, rec: dict | None) -> str:
    if rec is None:
        return "MISSING_NO_ROW"
    if method == "M0":
        if rec.get("vina_status") == "TIMEOUT":
            return "TIMEOUT"
        return rec.get("vina_status") or "MISSING"
    if method in {"M1", "M1b"}:
        return rec.get("job_status") or rec.get("reason") or "MISSING"
    if method == "M2":
        return rec.get("job_status") or "MISSING"
    if method == "M3":
        if rec.get("status") == "TIMEOUT":
            return "TIMEOUT"
        if rec.get("official_m3") == "NO":
            return rec.get("score_source") or "NOT_OFFICIAL"
        return rec.get("status") or "MISSING"
    return "MISSING"


def assert_one_to_one(left_n: int, right_n: int, merged_n: int, how: str, keys: tuple[str, ...]) -> None:
    """Refuse unexplained observation inflation."""
    if how == "inner":
        if merged_n > min(left_n, right_n):
            raise ValueError(f"many-to-many inflation on {keys}: {left_n}x{right_n}->{merged_n}")
    elif how == "left":
        if merged_n > left_n:
            raise ValueError(f"left merge inflated {left_n}->{merged_n} on {keys}")


def shared_dual_universe(rows: list[dict], method: str) -> dict:
    """One-method summary_min universe.

    Dual ligands must have valid score_A AND score_B (both pockets).
    A_only ligands need valid score_B (the B-pocket score used for dual vs A_only).
    B_only ligands need valid score_A.
    Dual resample is drawn once from dual_both and applied to both directions.
    """
    dual, a_only, b_only = [], [], []
    for r in rows:
        if r.get("activity_eligible") != "1" or r.get("class") not in {"dual", "A_only", "B_only", "neither"}:
            continue
        va = r.get(f"{method}_A_valid") == "1"
        vb = r.get(f"{method}_B_valid") == "1"
        cls = r["class"]
        if cls == "dual" and va and vb:
            dual.append(r)
        elif cls == "A_only" and vb:
            a_only.append(r)
        elif cls == "B_only" and va:
            b_only.append(r)
    return {
        "method": method,
        "level": "pair_observation",
        "dual_shared": dual,
        "A_only_for_AUC_B": a_only,
        "B_only_for_AUC_A": b_only,
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "AUC_B_definition": "dual vs A_only using score_B",
        "AUC_A_definition": "dual vs B_only using score_A",
        "summary_min_definition": "min(AUC_B, AUC_A) inside each bootstrap replicate",
        "forbidden": "AUROC of min(score_A, score_B)",
    }


def four_sided_common_complete(rows: list[dict], method_x: str, method_y: str) -> dict:
    """Inter-method Δsummary_min population.

    Four sides for each dual ligand: score_A_X, score_B_X, score_A_Y, score_B_Y.
    A_only: score_B on both methods. B_only: score_A on both methods.
    """
    dual, a_only, b_only = [], [], []
    for r in rows:
        if r.get("activity_eligible") != "1":
            continue
        xa, xb = r.get(f"{method_x}_A_valid") == "1", r.get(f"{method_x}_B_valid") == "1"
        ya, yb = r.get(f"{method_y}_A_valid") == "1", r.get(f"{method_y}_B_valid") == "1"
        cls = r.get("class")
        if cls == "dual" and xa and xb and ya and yb:
            dual.append(r)
        elif cls == "A_only" and xb and yb:
            a_only.append(r)
        elif cls == "B_only" and xa and ya:
            b_only.append(r)
    return {
        "method_x": method_x,
        "method_y": method_y,
        "delta": f"{method_x}-{method_y}",
        "level": "pair_observation",
        "dual_four_sided": dual,
        "A_only_both_score_B": a_only,
        "B_only_both_score_A": b_only,
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "bootstrap": (
            "One dual resample shared by both methods and both directions; "
            "one A_only resample; one B_only resample; "
            "inside replicate sm_x=min(AUC_B_x,AUC_A_x), sm_y=min(AUC_B_y,AUC_A_y), "
            "delta=sm_x-sm_y. Forbidden: independent method bootstraps then subtract CIs."
        ),
    }


def paired_delta_summary_min_bootstrap(pop: dict, n_boot: int = BOOT_B, seed: int = BOOT_SEED):
    """Paired bootstrap of Δsummary_min. LOCKED. Caller must pre-check sufficiency."""
    _locked("paired_delta_summary_min_bootstrap")
    rng = np.random.default_rng(seed)
    dual = pop["dual_four_sided"]
    ao = pop["A_only_both_score_B"]
    bo = pop["B_only_both_score_A"]
    mx, my = pop["method_x"], pop["method_y"]
    deltas = []
    for _ in range(n_boot):
        d_idx = rng.integers(0, len(dual), size=len(dual))
        a_idx = rng.integers(0, len(ao), size=len(ao))
        b_idx = rng.integers(0, len(bo), size=len(bo))
        d_s = [dual[i] for i in d_idx]
        a_s = [ao[i] for i in a_idx]
        b_s = [bo[i] for i in b_idx]
        auc_b_x = auroc(_finite_scores(d_s, f"{mx}_score_B"), _finite_scores(a_s, f"{mx}_score_B"))
        auc_a_x = auroc(_finite_scores(d_s, f"{mx}_score_A"), _finite_scores(b_s, f"{mx}_score_A"))
        auc_b_y = auroc(_finite_scores(d_s, f"{my}_score_B"), _finite_scores(a_s, f"{my}_score_B"))
        auc_a_y = auroc(_finite_scores(d_s, f"{my}_score_A"), _finite_scores(b_s, f"{my}_score_A"))
        deltas.append(min(auc_b_x, auc_a_x) - min(auc_b_y, auc_a_y))
    return deltas


def _finite_scores(rows: list[dict], key: str) -> list[float]:
    out = []
    for r in rows:
        if not finite(r.get(key)):
            raise ValueError(f"valid_row_nonfinite_score:{key}")
        out.append(float(r[key]))
    return out


def require_all_finite(values, label: str) -> list[float]:
    """Refuse silent NaN/Inf drop. Legal bootstrap must be all finite."""
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


def _na_metric(reason: str) -> dict:
    return {
        "value": "",
        "ci_lo": "",
        "ci_hi": "",
        "reason": reason,
        "bootstrap_started": "NO",
        "n_replicates": 0,
        "n_finite": 0,
    }


def _ok_metric(value: float, reps: list[float], n_boot: int = BOOT_B) -> dict:
    finite_reps = require_all_finite(reps, "bootstrap")
    lo, hi = percentile_ci(finite_reps, n_boot=n_boot)
    return {
        "value": f"{value:.10f}",
        "ci_lo": f"{lo:.10f}",
        "ci_hi": f"{hi:.10f}",
        "reason": "",
        "bootstrap_started": "YES",
        "n_replicates": len(finite_reps),
        "n_finite": len(finite_reps),
    }


def single_method_point_and_bootstrap(rows: list[dict], method: str,
                                      n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    """One pair×method analysis unit. Resets default_rng(seed). LOCKED via auroc."""
    _locked("single_method_point_and_bootstrap")
    u = shared_dual_universe(rows, method)
    dual, ao, bo = u["dual_shared"], u["A_only_for_AUC_B"], u["B_only_for_AUC_A"]
    do_b = len(dual) > 0 and len(ao) > 0
    do_a = len(dual) > 0 and len(bo) > 0
    do_sm = do_b and do_a

    auc_b_pt = auc_a_pt = sm_pt = None
    if do_b:
        auc_b_pt = auroc(_finite_scores(dual, f"{method}_score_B"),
                         _finite_scores(ao, f"{method}_score_B"))
    if do_a:
        auc_a_pt = auroc(_finite_scores(dual, f"{method}_score_A"),
                         _finite_scores(bo, f"{method}_score_A"))
    if do_sm:
        if not (np.isfinite(auc_b_pt) and np.isfinite(auc_a_pt)):
            raise ValueError(f"{method}:nonfinite_point_estimate")
        sm_pt = min(auc_b_pt, auc_a_pt)

    auc_b_reps, auc_a_reps, sm_reps = [], [], []
    if do_b or do_a:
        rng = np.random.default_rng(seed)
        for _ in range(n_boot):
            d_idx = rng.integers(0, len(dual), size=len(dual))
            d_s = [dual[i] for i in d_idx]
            auc_b = auc_a = None
            if do_b:
                a_idx = rng.integers(0, len(ao), size=len(ao))
                a_s = [ao[i] for i in a_idx]
                auc_b = auroc(_finite_scores(d_s, f"{method}_score_B"),
                              _finite_scores(a_s, f"{method}_score_B"))
                auc_b_reps.append(auc_b)
            if do_a:
                b_idx = rng.integers(0, len(bo), size=len(bo))
                b_s = [bo[i] for i in b_idx]
                auc_a = auroc(_finite_scores(d_s, f"{method}_score_A"),
                              _finite_scores(b_s, f"{method}_score_A"))
                auc_a_reps.append(auc_a)
            if do_sm:
                sm_reps.append(min(auc_b, auc_a))

    def pack(ok, pt, reps, insuff):
        if not ok:
            return _na_metric(insuff)
        if not np.isfinite(pt):
            raise ValueError(f"{method}:nonfinite_point_estimate")
        return _ok_metric(pt, reps, n_boot=n_boot)

    return {
        "method": method,
        "n_dual": len(dual),
        "n_A_only": len(ao),
        "n_B_only": len(bo),
        "universe": u,
        "dual_vs_A_only": pack(
            do_b, auc_b_pt, auc_b_reps,
            "INSUFFICIENT_n_dual" if len(dual) == 0 else "INSUFFICIENT_n_A_only",
        ),
        "dual_vs_B_only": pack(
            do_a, auc_a_pt, auc_a_reps,
            "INSUFFICIENT_n_dual" if len(dual) == 0 else "INSUFFICIENT_n_B_only",
        ),
        "summary_min": pack(
            do_sm, sm_pt, sm_reps,
            "INSUFFICIENT_for_summary_min",
        ),
    }


def paired_delta_point_and_bootstrap(rows: list[dict], method_x: str, method_y: str,
                                     n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    """One pair×primary-delta analysis unit. Δ(X−Y)=sm(X)−sm(Y). LOCKED via auroc."""
    _locked("paired_delta_point_and_bootstrap")
    pop = four_sided_common_complete(rows, method_x, method_y)
    dual, ao, bo = pop["dual_four_sided"], pop["A_only_both_score_B"], pop["B_only_both_score_A"]
    reasons = []
    if len(dual) == 0:
        reasons.append("INSUFFICIENT_n_dual")
    if len(ao) == 0:
        reasons.append("INSUFFICIENT_n_A_only")
    if len(bo) == 0:
        reasons.append("INSUFFICIENT_n_B_only")
    if reasons:
        out = _na_metric("+".join(reasons))
        out.update({
            "method_x": method_x,
            "method_y": method_y,
            "delta": f"{method_x}-{method_y}",
            "n_dual": len(dual),
            "n_A_only": len(ao),
            "n_B_only": len(bo),
            "universe": pop,
            "summary_min_x": "",
            "summary_min_y": "",
        })
        return out

    auc_b_x = auroc(_finite_scores(dual, f"{method_x}_score_B"),
                    _finite_scores(ao, f"{method_x}_score_B"))
    auc_a_x = auroc(_finite_scores(dual, f"{method_x}_score_A"),
                    _finite_scores(bo, f"{method_x}_score_A"))
    auc_b_y = auroc(_finite_scores(dual, f"{method_y}_score_B"),
                    _finite_scores(ao, f"{method_y}_score_B"))
    auc_a_y = auroc(_finite_scores(dual, f"{method_y}_score_A"),
                    _finite_scores(bo, f"{method_y}_score_A"))
    if not all(np.isfinite(v) for v in (auc_b_x, auc_a_x, auc_b_y, auc_a_y)):
        raise ValueError(f"{method_x}-{method_y}:nonfinite_point_estimate")
    sm_x = min(auc_b_x, auc_a_x)
    sm_y = min(auc_b_y, auc_a_y)
    sm_x_s = f"{sm_x:.10f}"
    sm_y_s = f"{sm_y:.10f}"
    point = float(sm_x_s) - float(sm_y_s)
    reps = paired_delta_summary_min_bootstrap(pop, n_boot=n_boot, seed=seed)
    packed = _ok_metric(point, reps, n_boot=n_boot)
    packed.update({
        "method_x": method_x,
        "method_y": method_y,
        "delta": f"{method_x}-{method_y}",
        "n_dual": len(dual),
        "n_A_only": len(ao),
        "n_B_only": len(bo),
        "universe": pop,
        "summary_min_x": sm_x_s,
        "summary_min_y": sm_y_s,
    })
    return packed


def static_self_audit() -> dict:
    """AST checks only. Does not compute AUROC. Ignores docstring mentions of forbidden patterns."""
    src_path = Path(__file__).resolve()
    src = src_path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    issues = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "analysis" or alias.name.startswith("analysis."):
                    issues.append("imports_old_analysis_tree")
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "analysis" or node.module.startswith("analysis."):
                issues.append("imports_old_analysis_tree")
    unlocked = None
    boot_b = boot_seed = None
    auroc_tie = ci_method = None
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if "METRICS_EXECUTION_UNLOCKED" in names:
                unlocked = ast.literal_eval(node.value)
            if "BOOT_B" in names:
                boot_b = ast.literal_eval(node.value)
            if "BOOT_SEED" in names:
                boot_seed = ast.literal_eval(node.value)
            if "AUROC_TIE" in names:
                auroc_tie = ast.literal_eval(node.value)
            if "CI_METHOD" in names:
                ci_method = ast.literal_eval(node.value)
    if unlocked is not False:
        issues.append("METRICS_EXECUTION_UNLOCKED_not_false")
    if boot_b != 10000:
        issues.append("BOOT_B_ne_10000")
    if boot_seed != 271828:
        issues.append("BOOT_SEED_ne_271828")
    if auroc_tie != 0.5:
        issues.append("AUROC_TIE_ne_0.5")
    if ci_method != "linear":
        issues.append("CI_METHOD_ne_linear")
    required = {
        "auroc", "shared_dual_universe", "four_sided_common_complete",
        "paired_delta_summary_min_bootstrap", "assert_one_to_one",
        "percentile_ci", "require_all_finite",
        "single_method_point_and_bootstrap", "paired_delta_point_and_bootstrap",
    }
    funs = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
    missing = sorted(required - funs)
    if missing:
        issues.append("missing_functions:" + ",".join(missing))
    uses_default_rng = any(
        isinstance(n, ast.Attribute) and n.attr == "default_rng"
        for n in ast.walk(tree)
    )
    if not uses_default_rng:
        issues.append("default_rng_missing")
    return {
        "implementation_file": str(src_path),
        "METRICS_EXECUTION_UNLOCKED": METRICS_EXECUTION_UNLOCKED,
        "BOOT_B": BOOT_B,
        "BOOT_SEED": BOOT_SEED,
        "issues": issues,
        "ok": not issues,
    }


def main() -> int:
    audit = static_self_audit()
    print(audit)
    print("AUROC_NOT_COMPUTED")
    return 0 if audit["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
