#!/usr/bin/env python3
"""Q110 statistics runner.

Uses jcim_stage2_lib and the accepted Stage 2 module functions.
Does not import or unlock formal_metrics_lib.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "rerun_v4_2"
sys.path.insert(0, str(SCRIPTS))

from jcim_stage2_lib import (  # noqa: E402
    BOOT_B,
    BOOT_SEED,
    SEEDS,
    auroc,
    five_seed_iqr,
    four_sided_common_complete,
    load_stage2_tables,
    percentile_ci,
    scores,
    shared_dual_universe,
)
from jcim_stage2_compute import (  # noqa: E402
    module_b,
    module_c,
    module_d,
    module_d1,
    module_e,
    module_h,
    module_i,
)

METHODS = ("M0", "M1", "M1b", "M2", "M3")
PRIMARY_DELTAS = (("M1", "M0"), ("M1b", "M1"), ("M2", "M0"), ("M3", "M0"))
PAIR = "PIK3CA/mTOR"


def fmt10(value: float) -> str:
    return f"{float(value):.10f}"


def single_method(rows: list[dict], method: str, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    """Same analysis unit as the locked primary metric, using the Stage 2 AUROC."""
    universe = shared_dual_universe(rows, method)
    dual = universe["dual_shared"]
    a_only = universe["A_only_for_AUC_B"]
    b_only = universe["B_only_for_AUC_A"]
    do_b = len(dual) > 0 and len(a_only) > 0
    do_a = len(dual) > 0 and len(b_only) > 0
    do_sm = do_b and do_a
    auc_b = auroc(scores(dual, f"{method}_score_B"), scores(a_only, f"{method}_score_B")) if do_b else None
    auc_a = auroc(scores(dual, f"{method}_score_A"), scores(b_only, f"{method}_score_A")) if do_a else None
    summary = min(auc_b, auc_a) if do_sm else None
    reps_b, reps_a, reps_s = [], [], []
    if do_b or do_a:
        rng = np.random.default_rng(seed)
        for _ in range(n_boot):
            drawn_dual = [dual[i] for i in rng.integers(0, len(dual), size=len(dual))]
            one_b = one_a = None
            if do_b:
                drawn_a = [a_only[i] for i in rng.integers(0, len(a_only), size=len(a_only))]
                one_b = auroc(scores(drawn_dual, f"{method}_score_B"), scores(drawn_a, f"{method}_score_B"))
                reps_b.append(one_b)
            if do_a:
                drawn_b = [b_only[i] for i in rng.integers(0, len(b_only), size=len(b_only))]
                one_a = auroc(scores(drawn_dual, f"{method}_score_A"), scores(drawn_b, f"{method}_score_A"))
                reps_a.append(one_a)
            if do_sm:
                reps_s.append(min(one_b, one_a))

    def pack(ok, point, reps):
        if not ok:
            return {"value": "", "ci_lo": "", "ci_hi": "", "n_replicates": 0}
        lo, hi = percentile_ci(reps, n_boot=n_boot)
        return {"value": fmt10(point), "ci_lo": fmt10(lo), "ci_hi": fmt10(hi), "n_replicates": len(reps)}

    return {
        "method": method,
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "dual_vs_A_only": pack(do_b, auc_b, reps_b),
        "dual_vs_B_only": pack(do_a, auc_a, reps_a),
        "summary_min": pack(do_sm, summary, reps_s),
    }


def paired_delta(rows: list[dict], method_x: str, method_y: str, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    pop = four_sided_common_complete(rows, method_x, method_y)
    dual = pop["dual_four_sided"]
    a_only = pop["A_only_both_score_B"]
    b_only = pop["B_only_both_score_A"]
    out = {
        "delta": f"{method_x}-{method_y}",
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "n_dual_four_sided": len(dual),
        "n_A_only_both_score_B": len(a_only),
        "n_B_only_both_score_A": len(b_only),
    }
    if not dual or not a_only or not b_only:
        out.update(value="", ci_lo="", ci_hi="")
        return out
    def sm(sample_d, sample_a, sample_b, method):
        auc_b = auroc(scores(sample_d, f"{method}_score_B"), scores(sample_a, f"{method}_score_B"))
        auc_a = auroc(scores(sample_d, f"{method}_score_A"), scores(sample_b, f"{method}_score_A"))
        return min(auc_b, auc_a)
    sm_x = float(fmt10(sm(dual, a_only, b_only, method_x)))
    sm_y = float(fmt10(sm(dual, a_only, b_only, method_y)))
    rng = np.random.default_rng(seed)
    reps = []
    for _ in range(n_boot):
        drawn_d = [dual[i] for i in rng.integers(0, len(dual), size=len(dual))]
        drawn_a = [a_only[i] for i in rng.integers(0, len(a_only), size=len(a_only))]
        drawn_b = [b_only[i] for i in rng.integers(0, len(b_only), size=len(b_only))]
        reps.append(sm(drawn_d, drawn_a, drawn_b, method_x) - sm(drawn_d, drawn_a, drawn_b, method_y))
    lo, hi = percentile_ci(reps, n_boot=n_boot)
    out.update(delta_summary_min=fmt10(sm_x - sm_y), ci_lo=fmt10(lo), ci_hi=fmt10(hi), summary_min_x=fmt10(sm_x), summary_min_y=fmt10(sm_y))
    return out


def primary_rows(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    metrics, deltas = [], []
    for method in METHODS:
        result = single_method(rows, method)
        for contrast, key, n_neg in (
            ("dual_vs_A_only", "dual_vs_A_only", result["n_A_only"]),
            ("dual_vs_B_only", "dual_vs_B_only", result["n_B_only"]),
            ("summary_min", "summary_min", ""),
        ):
            packed = result[key]
            metrics.append({
                "pair": rows[0]["pair"],
                "method": method,
                "contrast": contrast,
                "n_dual": result["n_dual"],
                "n_negative": n_neg,
                "auroc": packed["value"],
                "ci_lo": packed["ci_lo"],
                "ci_hi": packed["ci_hi"],
                "n_replicates": packed["n_replicates"],
            })
    for left, right in PRIMARY_DELTAS:
        result = paired_delta(rows, left, right)
        deltas.append({"pair": rows[0]["pair"], **result})
    return metrics, deltas


def five_seed_fixed_membership(seed_rows: list[dict]) -> list[dict]:
    """Recompute is supplied by module E's per-seed rows only when n matches.

    The production fixed-membership table is built by q110_five_seed.py from
    raw scores. This helper checks the published full-sample case.
    """
    by_seed = {int(r["seed"]): r for r in seed_rows}
    counts = {(r["n_dual"], r["n_A_only"], r["n_B_only"]) for r in seed_rows}
    same = len(counts) == 1
    out = []
    for seed in SEEDS:
        row = by_seed[seed]
        out.append({
            "population": "fixed_membership" if same else "not_identical_to_per_seed",
            "seed": seed,
            "n_dual": row["n_dual"] if same else "",
            "n_A_only": row["n_A_only"] if same else "",
            "n_B_only": row["n_B_only"] if same else "",
            "auc_B": row["auc_B"] if same else "",
            "auc_A": row["auc_A"] if same else "",
            "summary_min": row["summary_min"] if same else "",
        })
    return out


def _norm(value: str) -> str:
    text = "" if value is None else str(value)
    if text in {"", "NA", "None"}:
        return ""
    try:
        number = float(text)
    except ValueError:
        return text
    if not np.isfinite(number):
        return text
    return fmt10(number)


def _key_part(value) -> str:
    """Join keys across int and str forms. Does not rescale measured values."""
    text = "" if value is None else str(value).strip()
    if text in {"", "NA", "None", "nan"}:
        return ""
    try:
        number = float(text)
    except ValueError:
        return text
    if np.isfinite(number) and number.is_integer():
        return str(int(number))
    return text


def _row_key(row: dict, keys: list[str]) -> tuple[str, ...]:
    return tuple(_key_part(row[k]) for k in keys)


def compare_rows(name: str, got: list[dict], ref: list[dict], keys: list[str], fields: list[str]) -> list[dict]:
    ref_map: dict[tuple[str, ...], dict] = {}
    checks = []
    for old in ref:
        key = _row_key(old, keys)
        if key in ref_map:
            checks.append(_check(name, "|".join(key), "", "FAIL", "", "", "duplicate reference key"))
        ref_map[key] = old
    seen = set()
    for row in got:
        key = _row_key(row, keys)
        label = "|".join(key)
        if key in seen:
            checks.append(_check(name, label, "", "FAIL", "", "", "duplicate computed key"))
        seen.add(key)
        old = ref_map.get(key)
        if old is None:
            checks.append(_check(name, label, "", "FAIL", "", "", "missing reference"))
            continue
        for field in fields:
            if field not in row or field not in old:
                checks.append(_check(name, label, field, "FAIL", row.get(field, ""), old.get(field, ""), "field absent"))
                continue
            status = "PASS" if _norm(row[field]) == _norm(old[field]) else "FAIL"
            checks.append(_check(name, label, field, status, row[field], old[field], ""))
    for key in ref_map:
        if key not in seen:
            checks.append(_check(name, "|".join(key), "", "FAIL", "", "", "missing computed row"))
    return checks


def _check(name: str, key: str, field: str, status: str, got, ref, detail: str) -> dict:
    return {
        "module": name,
        "key": key,
        "field": field,
        "status": status,
        "got": got,
        "ref": ref,
        "detail": detail,
    }


def load_pair(project: Path, pair: str = PAIR) -> list[dict]:
    tables = load_stage2_tables(project)
    rows = [r for r in tables["population"] if r["pair"] == pair]
    if len(rows) != 48:
        raise SystemExit(f"expected 48 PM48 rows, found {len(rows)}")
    return rows, tables


def run_pm48_regression(project: Path, outdir: Path) -> dict:
    rows, tables = load_pair(project)
    pik_tables = dict(tables)
    pik_tables["population"] = rows
    outdir.mkdir(parents=True, exist_ok=True)
    metrics, deltas = primary_rows(rows)
    write(outdir / "pm48_primary_metrics.csv", metrics)
    write(outdir / "pm48_primary_deltas.csv", deltas)
    published = list(csv.DictReader((project / "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv").open()))
    published_d = list(csv.DictReader((project / "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv").open()))
    checks = compare_rows(
        "primary",
        metrics,
        [r for r in published if r["pair"] == PAIR],
        ["method", "contrast"],
        ["auroc", "ci_lo", "ci_hi", "n_dual", "n_negative"],
    )
    checks += compare_rows(
        "delta",
        deltas,
        [r for r in published_d if r["pair"] == PAIR],
        ["delta"],
        ["delta_summary_min", "ci_lo", "ci_hi", "n_dual_four_sided", "n_A_only_both_score_B", "n_B_only_both_score_A", "summary_min_x", "summary_min_y"],
    )
    # Stage 2 accepted modules on the PM48 population.
    from jcim_stage2_lib import (
        SCHEMA_B, SCHEMA_C, SCHEMA_D1, SCHEMA_D_OOF, SCHEMA_D_UNIT,
        SCHEMA_E_SEED, SCHEMA_E_SUM, SCHEMA_H, SCHEMA_I_CMP, SCHEMA_I_DESC, SCHEMA_I_LIG,
        write_csv,
    )
    b_rows = module_b(rows)
    c_rows = module_c(rows)
    d_units, d_oof = module_d(pik_tables)
    d1 = module_d1(pik_tables)
    e_seed, e_sum = module_e(pik_tables)
    h_rows = module_h(rows)
    i_desc, i_lig, i_cmp = module_i(rows, [r for r in published_d if r["pair"] == PAIR])
    write_csv(outdir / "module_B_negative_class.csv", b_rows, SCHEMA_B)
    write_csv(outdir / "module_C_wrong_pocket.csv", c_rows, SCHEMA_C)
    write_csv(outdir / "module_D_chemistry_units.csv", d_units, SCHEMA_D_UNIT)
    write_csv(outdir / "module_D_chemistry_oof_ligands.csv", d_oof, SCHEMA_D_OOF)
    write_csv(outdir / "module_D1_composition.csv", d1, SCHEMA_D1)
    write_csv(outdir / "module_E_five_seed.csv", e_seed, SCHEMA_E_SEED)
    write_csv(outdir / "module_E_five_seed_summary.csv", e_sum, SCHEMA_E_SUM)
    write_csv(outdir / "module_H_missingness.csv", h_rows, SCHEMA_H)
    write_csv(outdir / "module_I_descriptive.csv", i_desc, SCHEMA_I_DESC)
    write_csv(outdir / "module_I_ligand_ranks.csv", i_lig, SCHEMA_I_LIG)
    write_csv(outdir / "module_I_method_compare.csv", i_cmp, SCHEMA_I_CMP)
    stage = project / "results/jcim_stage2"
    checks += compare_rows("B", b_rows, [r for r in csv.DictReader((stage / "module_B_negative_class.csv").open()) if r["pair"] == PAIR], ["direction"], ["auc_single", "auc_neither", "delta_negative", "auc_single_ci_lo", "auc_single_ci_hi", "auc_neither_ci_lo", "auc_neither_ci_hi", "delta_negative_ci_lo", "delta_negative_ci_hi", "n_dual", "n_single", "n_neither"])
    checks += compare_rows("C", c_rows, [r for r in csv.DictReader((stage / "module_C_wrong_pocket.csv").open()) if r["pair"] == PAIR], ["arm"], ["auc_correct", "auc_wrong", "delta_pocket", "auc_correct_ci_lo", "auc_correct_ci_hi", "auc_wrong_ci_lo", "auc_wrong_ci_hi", "delta_pocket_ci_lo", "delta_pocket_ci_hi", "n_dual", "n_negative_directional_m0"])
    checks += compare_rows("D", d_units, [r for r in csv.DictReader((stage / "module_D_chemistry_units.csv").open()) if r["pair"] == PAIR], ["arm"], ["n_labeled_population", "n_m0_directional_population", "n_final_chemistry_oof", "auroc_physchem", "auroc_ecfp4", "auroc_ecfp4_m0", "delta_ecfp4_plus_m0", "evidence_status"])
    checks += compare_rows("Doof", d_oof, [r for r in csv.DictReader((stage / "module_D_chemistry_oof_ligands.csv").open()) if r["pair"] == PAIR], ["arm", "canonical_ligand_id"], ["class", "fold_id", "oof_physchem", "oof_ecfp4", "oof_ecfp4_m0", "relevant_m0"])
    checks += compare_rows("D1", d1, [r for r in csv.DictReader((stage / "module_D1_composition.csv").open()) if r["pair"] == PAIR], ["arm", "canonical_ligand_id"], ["class", "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "RingCount", "FormalCharge", "murcko", "max_ecfp4_tanimoto_other_class"])
    checks += compare_rows("E", e_seed, [r for r in csv.DictReader((stage / "module_E_five_seed.csv").open()) if r["pair"] == PAIR], ["seed"], ["n_dual", "n_A_only", "n_B_only", "auc_B", "auc_A", "summary_min"])
    checks += compare_rows("Esum", e_sum, [r for r in csv.DictReader((stage / "module_E_five_seed_summary.csv").open()) if r["pair"] == PAIR], ["metric"], ["seed42", "median", "min", "max", "iqr", "n_finite_seeds"])
    checks += compare_rows("H", h_rows, [r for r in csv.DictReader((stage / "module_H_missingness.csv").open()) if r["pair"] == PAIR], ["method", "class", "target"], ["expected", "valid", "missing", "missing_reason"])
    checks += compare_rows("I", i_desc, [r for r in csv.DictReader((stage / "module_I_descriptive.csv").open()) if r["pair"] == PAIR], ["method", "population_kind"], ["k", "n_total", "n_dual", "n_A_only", "n_B_only", "n_dual_top", "n_A_only_top", "n_B_only_top", "top10_dual_fraction", "top10_single_target_fraction", "dual_recall"])
    checks += compare_rows("Ilig", i_lig, [r for r in csv.DictReader((stage / "module_I_ligand_ranks.csv").open()) if r["pair"] == PAIR], ["method", "population_kind", "canonical_ligand_id"], ["global_ligand_entity_id", "class", "rank_A", "rank_B", "worst_rank", "joint_order", "in_topk"])
    checks += compare_rows("Icmp", i_cmp, [r for r in csv.DictReader((stage / "module_I_method_compare.csv").open()) if r["pair"] == PAIR], ["method"], ["n_rank_common", "n_directional_common_dual", "rank_vs_directional_member_diff", "delta_summary_min", "delta_top10_single_target_fraction", "concordance", "link_source", "link_reason", "supporting_only", "bootstrap_ci_added"])
    fixed = five_seed_fixed_membership(e_seed)
    write(outdir / "pm48_five_seed_fixed_membership.csv", fixed)
    same_population = all(r["population"] == "fixed_membership" for r in fixed)
    checks.append(_check(
        "Efixed_population",
        PAIR,
        "population",
        "PASS" if same_population else "FAIL",
        "fixed_membership" if same_population else "not_identical_to_per_seed",
        "fixed_membership",
        "" if same_population else "five-seed populations differ; fixed-membership check needs raw rescoring",
    ))
    if same_population:
        checks += compare_rows("Efixed", fixed, e_seed, ["seed"], ["n_dual", "n_A_only", "n_B_only", "auc_B", "auc_A", "summary_min"])
    failures = [row for row in checks if row["status"] != "PASS"]
    write(outdir / "PM48_REGRESSION_PASSFAIL.csv", checks)
    iqr_check = five_seed_iqr([float(r["summary_min"]) for r in e_seed])
    by_module: dict[str, dict[str, int]] = {}
    for row in checks:
        slot = by_module.setdefault(row["module"], {"PASS": 0, "FAIL": 0})
        slot[row["status"]] = slot.get(row["status"], 0) + 1
    report = {
        "verdict": "PASS" if not failures else "FAIL",
        "n_checks": len(checks),
        "n_pass": len(checks) - len(failures),
        "n_fail": len(failures),
        "by_module": by_module,
        "failures": failures,
        "pm48_summary_min_iqr": iqr_check,
        "d1_rows": len(d1),
        "oof_rows": len(d_oof),
        "i_ligands": len(i_lig),
        "i_compare": len(i_cmp),
    }
    (outdir / "PM48_REGRESSION_REPORT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("verdict", "n_checks", "n_pass", "n_fail", "by_module")}, indent=2))
    for item in failures:
        print(f"FAIL {item['module']} {item['key']} {item['field']}: {item['detail']} got={item['got']} ref={item['ref']}")
    return report


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    project = ROOT
    outdir = Path(__file__).resolve().parent / "regression_pm48"
    report = run_pm48_regression(project, outdir)
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
