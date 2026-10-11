#!/usr/bin/env python3
"""Activity-label sensitivity for the fixed 805-observation set.

Reads restored ChEMBL 37 rows and the formal M0 population.
Does not import scripts/analysis or the locked formal_metrics_lib.auroc.
Does not write mapping, population, PRIMARY, YAML, or score masters.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
OUT = Path(__file__).resolve().parent
SRC = ROOT / "activity_label_source_audit_20261008"
SEM = ROOT / "activity_semantics_audit_20261009"

POP_PATH = RUN / "13_qa" / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
MAP_PATH = RUN / "00_protocol" / "pair_ligand_mapping.csv"
ALIAS_PATH = RUN / "00_protocol" / "parent_alias_registry.csv"
FOLD_PATH = ROOT / "results" / "canonical" / "model_fold_assignments.csv"
PRIMARY_PATH = ROOT / "results" / "formal_metrics" / "PRIMARY_DIRECTIONAL_METRICS.csv"
LONG_PATH = SRC / "chembl37_activity_rows_long.csv"
CROSS_PATH = SRC / "chembl37_dual_endpoint_crosswalk_805.csv"
DOC_PATH = SRC / "document_identifier_coverage_805.csv"
HIST_PATH = SRC / "historical_adjudication_differences.csv"
PPAR_PATH = SEM / "PPAR_activity_semantics.csv"
ADME_PATH = SEM / "ADME_assay_review.csv"
ENRICHED = Path("/tmp/semantics_work/enriched_3555.csv")

BOOT_B = 10000
BOOT_SEED = 271828
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
FOUR = ("dual", "A_only", "B_only", "neither")
SCENARIOS = ("S0", "S1", "S2", "S3", "S4")
ROLE = {
    "S0": "formal_reference_max_theta_6.0",
    "S1": "prespecified_threshold_sensitivity",
    "S2": "prespecified_threshold_sensitivity",
    "S3": "prespecified_threshold_sensitivity",
    "S4": "additional_median_sensitivity_not_prespecified",
}
DEFINITION = {
    "S0": "最大 pChEMBL，θ=6.0。正式参照，不是新阈值。",
    "S1": "预先规定的阈值敏感性：最大 pChEMBL，θ=5.5。有活性为 ≥5.5。",
    "S2": "预先规定的阈值敏感性：最大 pChEMBL，θ=6.5。有活性为 ≥6.5。",
    "S3": "预先规定的严格分离：最大 pChEMBL。dual 两端≥6.5；A-only 为 A≥6.5 且 B≤5.5；B-only 为 A≤5.5 且 B≥6.5；neither 两端≤5.5。其余为 strict_gray_excluded，不是 neither。",
    "S4": "附加汇总敏感性，不是早期预先规定的实验：同一批原始记录的中位 pChEMBL，θ=6.0。中位数不比最大值更正确。",
}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def fnum(x: float) -> str:
    return f"{x:.10f}"


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def scores(rows: list[dict], key: str) -> list[float]:
    out = []
    for row in rows:
        val = float(row[key])
        if not np.isfinite(val):
            raise ValueError(f"nonfinite {key} {row.get('panel_id')}")
        out.append(val)
    return out


def four_class(pa: float, pb: float, theta: float) -> str:
    a_on = pa >= theta
    b_on = pb >= theta
    if a_on and b_on:
        return "dual"
    if a_on:
        return "A_only"
    if b_on:
        return "B_only"
    return "neither"


def strict_class(pa: float, pb: float) -> str:
    def band(p: float) -> str:
        if p >= 6.5:
            return "H"
        if p <= 5.5:
            return "L"
        return "G"

    sa, sb = band(pa), band(pb)
    if sa == "G" or sb == "G":
        return "strict_gray_excluded"
    if sa == "H" and sb == "H":
        return "dual"
    if sa == "H" and sb == "L":
        return "A_only"
    if sa == "L" and sb == "H":
        return "B_only"
    return "neither"


def universe(rows: list[dict], class_key: str):
    dual, a_only, b_only = [], [], []
    for row in rows:
        if row.get("activity_eligible") != "1":
            continue
        cls = row.get(class_key, "")
        if cls not in FOUR:
            continue
        va = row.get("M0_A_valid") == "1"
        vb = row.get("M0_B_valid") == "1"
        if cls == "dual" and va and vb:
            dual.append(row)
        elif cls == "A_only" and vb:
            a_only.append(row)
        elif cls == "B_only" and va:
            b_only.append(row)
    return dual, a_only, b_only


def analyze(rows: list[dict], class_key: str) -> dict:
    """Match formal_metrics_lib.single_method_point_and_bootstrap for M0.

    Population order is the caller order. Each replicate draws dual, then
    A-only, then B-only. An empty required class is not drawn, matching the
    formal branch that skips that integers() call.
    """
    dual, a_only, b_only = universe(rows, class_key)
    do_b = len(dual) > 0 and len(a_only) > 0
    do_a = len(dual) > 0 and len(b_only) > 0
    do_sm = do_b and do_a
    auc_b = auroc(scores(dual, "M0_score_B"), scores(a_only, "M0_score_B")) if do_b else None
    auc_a = auroc(scores(dual, "M0_score_A"), scores(b_only, "M0_score_A")) if do_a else None
    sm = min(auc_b, auc_a) if do_sm else None
    b_reps, a_reps, sm_reps = [], [], []
    if do_b or do_a:
        rng = np.random.default_rng(BOOT_SEED)
        for _ in range(BOOT_B):
            d_idx = rng.integers(0, len(dual), size=len(dual))
            d_s = [dual[i] for i in d_idx]
            auc_b_i = auc_a_i = None
            if do_b:
                a_idx = rng.integers(0, len(a_only), size=len(a_only))
                a_s = [a_only[i] for i in a_idx]
                auc_b_i = auroc(scores(d_s, "M0_score_B"), scores(a_s, "M0_score_B"))
                b_reps.append(auc_b_i)
            if do_a:
                b_idx = rng.integers(0, len(b_only), size=len(b_only))
                b_s = [b_only[i] for i in b_idx]
                auc_a_i = auroc(scores(d_s, "M0_score_A"), scores(b_s, "M0_score_A"))
                a_reps.append(auc_a_i)
            if do_sm:
                sm_reps.append(min(auc_b_i, auc_a_i))

    def pack(ok, point, reps, reason):
        if not ok:
            return {
                "auroc": "",
                "ci_lo": "",
                "ci_hi": "",
                "reason": reason,
                "bootstrap_started": "NO",
                "n_replicates": "0",
            }
        arr = np.asarray(reps, dtype=float)
        if arr.size != BOOT_B or not np.all(np.isfinite(arr)) or not np.isfinite(point):
            raise ValueError("bootstrap replicates were not all finite")
        lo, hi = np.percentile(arr, [2.5, 97.5], method="linear")
        return {
            "auroc": fnum(point),
            "ci_lo": fnum(float(lo)),
            "ci_hi": fnum(float(hi)),
            "reason": "",
            "bootstrap_started": "YES",
            "n_replicates": str(BOOT_B),
        }

    return {
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "dual_ids": [r["canonical_ligand_id"] for r in dual],
        "a_ids": [r["canonical_ligand_id"] for r in a_only],
        "b_ids": [r["canonical_ligand_id"] for r in b_only],
        "dual_vs_A_only": pack(
            do_b, auc_b, b_reps,
            "INSUFFICIENT_n_dual" if len(dual) == 0 else "INSUFFICIENT_n_A_only",
        ),
        "dual_vs_B_only": pack(
            do_a, auc_a, a_reps,
            "INSUFFICIENT_n_dual" if len(dual) == 0 else "INSUFFICIENT_n_B_only",
        ),
        "summary_min": pack(do_sm, sm, sm_reps, "INSUFFICIENT_for_summary_min"),
    }


def layer_counts(pop, alias_ids, fold_keys):
    arms = (("D_vs_A", {"dual", "A_only"}), ("D_vs_B", {"dual", "B_only"}))
    n1 = n2 = n3 = 0
    for pair in PAIRS:
        for arm, classes in arms:
            labeled = [
                r for r in pop
                if r["pair"] == pair
                and r.get("activity_eligible") == "1"
                and r.get("canonical_ligand_id") not in alias_ids
                and r.get("class") in classes
            ]
            def m0_ok(row):
                if row["class"] == "dual":
                    return row.get("M0_A_valid") == "1" and row.get("M0_B_valid") == "1"
                if arm == "D_vs_A":
                    return row.get("M0_B_valid") == "1"
                return row.get("M0_A_valid") == "1"
            m0 = [r for r in labeled if m0_ok(r)]
            chem = [r for r in m0 if (pair, arm, r["canonical_ligand_id"]) in fold_keys]
            n1 += len(labeled)
            n2 += len(m0)
            n3 += len(chem)
    return n1, n2, n3


def change_reason(scenario, old, new, pa, pb, med_a, med_b) -> str:
    if scenario == "S0":
        return "最大 pChEMBL 与 θ=6.0，类别与正式标签相同" if old == new else "最大 pChEMBL 与正式标签不一致"
    if scenario == "S3":
        if new == "strict_gray_excluded":
            bits = []
            if 5.5 < pa < 6.5:
                bits.append(f"A={pa} 处于 (5.5, 6.5)")
            if 5.5 < pb < 6.5:
                bits.append(f"B={pb} 处于 (5.5, 6.5)")
            if not bits:
                bits.append("两端组合不属于四条严格定义")
            return "严格灰区，未分类：" + "；".join(bits)
        if old == new:
            return "严格规则下仍为原类别"
        return f"严格 6.5/5.5 将类别从 {old} 改为 {new}"
    if scenario in {"S1", "S2"}:
        theta = 5.5 if scenario == "S1" else 6.5
        if old == new:
            return f"最大 pChEMBL 在 θ={theta} 下仍为 {old}"
        ends = []
        if (pa >= 6.0) != (pa >= theta):
            ends.append(f"A={pa} 相对 {theta} 与相对 6.0 的活性侧不同")
        if (pb >= 6.0) != (pb >= theta):
            ends.append(f"B={pb} 相对 {theta} 与相对 6.0 的活性侧不同")
        detail = "；".join(ends) if ends else "两端组合改变了四类"
        return f"阈值改为 {theta}：" + detail
    if old == new:
        return "中位 pChEMBL 在 θ=6.0 下仍为原类别；这不表示中位数更正确"
    ends = []
    if (pa >= 6.0) != (med_a >= 6.0):
        ends.append(f"A 最大 {pa}、中位 {med_a}")
    if (pb >= 6.0) != (med_b >= 6.0):
        ends.append(f"B 最大 {pb}、中位 {med_b}")
    detail = "；".join(ends) if ends else "两端中位数的组合改变了四类"
    return "附加中位数汇总：" + detail


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    pop = read_csv(POP_PATH)
    mapping = read_csv(MAP_PATH)
    alias = read_csv(ALIAS_PATH)
    folds = read_csv(FOLD_PATH)
    primary = read_csv(PRIMARY_PATH)
    long_rows = read_csv(LONG_PATH)
    cross = read_csv(CROSS_PATH)
    docs = read_csv(DOC_PATH)
    hist = read_csv(HIST_PATH)
    ppar = read_csv(PPAR_PATH) if PPAR_PATH.exists() else []
    adme = read_csv(ADME_PATH) if ADME_PATH.exists() else []
    enriched = read_csv(ENRICHED) if ENRICHED.exists() else []

    alias_ids = {r["alias_id"] for r in alias}
    fold_keys = {(r["pair"], r["arm"], r["ligand_id"]) for r in folds}
    mol_by = {(r["pair"], r["panel_id"]): r["molecule_chembl_id"] for r in cross}
    for row in pop:
        row["molecule_chembl_id"] = mol_by.get((row["pair"], row["panel_id"]), "")
    by_pair = defaultdict(list)
    for row in pop:
        by_pair[row["pair"]].append(row)

    eligible = [r for r in pop if r.get("activity_eligible") == "1" and r.get("class") in FOUR]
    ineligible = [r for r in pop if r.get("activity_eligible") != "1" or r.get("class") not in FOUR]
    molecules = {r["molecule_chembl_id"] for r in cross}
    n1, n2, n3 = layer_counts(pop, alias_ids, fold_keys)

    values = defaultdict(list)
    for row in long_rows:
        values[(row["pair"], row["panel_id"], row["endpoint"])].append(float(row["pchembl_value"]))
    agg = {}
    for key, vals in values.items():
        arr = np.asarray(vals, dtype=float)
        agg[key] = (float(arr.max()), float(np.median(arr)), len(vals))

    recon = []

    def add(check, ok, detail):
        recon.append({"check": check, "result": "PASS" if ok else "FAIL", "detail": detail})

    add("population_rows", len(pop) == 807, str(len(pop)))
    add("mapping_rows", len(mapping) == 808, str(len(mapping)))
    add("alias_only_AB_040", alias_ids == {"AB_040"}, ",".join(sorted(alias_ids)))
    add("AB_040_absent_from_population", all(r["canonical_ligand_id"] != "AB_040" for r in pop), "")
    add("eligible_fourclass_805", len(eligible) == 805, str(len(eligible)))
    add("ineligible_not_called_neither", len(ineligible) == 2 and all(r["class"] == "" for r in ineligible),
        ",".join(r["canonical_ligand_id"] + ":" + r["class"] for r in ineligible))
    add("unique_molecules_783", len(molecules) == 783, str(len(molecules)))
    add("activity_rows_3555", len(long_rows) == 3555, str(len(long_rows)))
    add("direction_layers_938_934_928", (n1, n2, n3) == (938, 934, 928), f"{n1}/{n2}/{n3}")
    add("crosswalk_805", len(cross) == 805, str(len(cross)))

    missing_end = []
    class_mismatch = []
    panel_mismatch = []
    for row in eligible:
        key = (row["pair"], row["panel_id"])
        if (key + ("A",)) not in agg or (key + ("B",)) not in agg:
            missing_end.append(row["panel_id"])
            continue
        pa, med_a, _ = agg[key + ("A",)]
        pb, med_b, _ = agg[key + ("B",)]
        got = four_class(pa, pb, 6.0)
        if got != row["class"]:
            class_mismatch.append(f"{row['panel_id']}:{row['class']}->{got}")
        row["pA_max"] = pa
        row["pB_max"] = pb
        row["pA_median"] = med_a
        row["pB_median"] = med_b
        row["S0"] = got
        row["S1"] = four_class(pa, pb, 5.5)
        row["S2"] = four_class(pa, pb, 6.5)
        row["S3"] = strict_class(pa, pb)
        row["S4"] = four_class(med_a, med_b, 6.0)
    cross_by = {(r["pair"], r["panel_id"]): r for r in cross}
    for row in eligible:
        src = cross_by[(row["pair"], row["panel_id"])]
        if src["class_matches_frozen"] != "YES" or src["dump_theta6_class"] != row["class"]:
            panel_mismatch.append(row["panel_id"])
        if abs(float(src["dump_max_A"]) - row["pA_max"]) > 0 or abs(float(src["dump_max_B"]) - row["pB_max"]) > 0:
            panel_mismatch.append(row["panel_id"] + ":max")
    add("both_endpoints_present", not missing_end, ",".join(missing_end[:10]))
    add("max_theta6_matches_formal_class_805", not class_mismatch, ",".join(class_mismatch[:10]))
    add("crosswalk_max_and_class", not panel_mismatch, ",".join(panel_mismatch[:10]))
    add("fixed_set_no_resampling", True, "S1–S4 只对这 805 条重新分类，不从 ChEMBL 重新抽样")

    write_csv(OUT / "01_input_reconciliation.csv", recon, ["check", "result", "detail"])
    if any(r["result"] == "FAIL" for r in recon):
        (OUT / "PROBLEM_REPORT.md").write_text(
            "输入核对未通过，未计算 S1–S4，未合并。\n\n"
            + "\n".join(f"- {r['check']}: {r['result']} {r['detail']}" for r in recon)
            + "\n",
            encoding="utf-8",
        )
        print("INPUT_RECONCILIATION_FAIL")
        return 2

    primary_m0 = [r for r in primary if r["method"] == "M0"]
    acceptance = []
    s0_blocks = {}
    point_fail = ci_fail = False
    for pair in PAIRS:
        block = analyze(by_pair[pair], "class")
        s0_blocks[pair] = block
        for contrast, n_neg_key in (
            ("dual_vs_A_only", "n_A_only"),
            ("dual_vs_B_only", "n_B_only"),
            ("summary_min", None),
        ):
            got = block[contrast]
            exp = next(r for r in primary_m0 if r["pair"] == pair and r["contrast"] == contrast)
            n_neg = "" if n_neg_key is None else str(block[n_neg_key])
            checks = {
                "n_dual": (str(block["n_dual"]), exp["n_dual"]),
                "n_negative": (n_neg, exp["n_negative"]),
                "n_A_only": (str(block["n_A_only"]), exp["n_A_only"]),
                "n_B_only": (str(block["n_B_only"]), exp["n_B_only"]),
                "auroc": (got["auroc"], exp["auroc"]),
                "ci_lo": (got["ci_lo"], exp["ci_lo"]),
                "ci_hi": (got["ci_hi"], exp["ci_hi"]),
            }
            for field, (left, right) in checks.items():
                ok = left == right
                if not ok and field == "auroc":
                    point_fail = True
                if not ok and field in {"ci_lo", "ci_hi"}:
                    ci_fail = True
                acceptance.append({
                    "pair": pair,
                    "contrast": contrast,
                    "field": field,
                    "computed": left,
                    "primary_m0": right,
                    "match": "YES" if ok else "NO",
                })
    write_csv(
        OUT / "00_s0_acceptance.csv",
        acceptance,
        ["pair", "contrast", "field", "computed", "primary_m0", "match"],
    )
    if point_fail or ci_fail:
        kind = "点估计" if point_fail else "置信区间"
        bad = [r for r in acceptance if r["match"] == "NO"]
        lines = [
            f"S0 验收失败：{kind}与 PRIMARY 的 M0 不一致。未计算 S1–S4，未合并。",
            "",
            "抽样实现：人群表原顺序；每次重复依次 default_rng(271828).integers 抽取 dual、A-only、B-only；dual 两个方向共用；np.percentile(..., method='linear')。",
            "",
        ]
        for row in bad[:30]:
            lines.append(f"- {row['pair']} {row['contrast']} {row['field']}: {row['computed']} != {row['primary_m0']}")
        (OUT / "PROBLEM_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("S0_ACCEPTANCE_FAIL", "point" if point_fail else "ci", len(bad))
        return 2

    blocks = {pair: {"S0": s0_blocks[pair]} for pair in PAIRS}
    for pair in PAIRS:
        rows = by_pair[pair]
        for scenario in ("S1", "S2", "S3", "S4"):
            blocks[pair][scenario] = analyze(rows, scenario)

    transitions = []
    for pair in PAIRS:
        for row in by_pair[pair]:
            if row.get("activity_eligible") != "1" or row.get("class") not in FOUR:
                continue
            for scenario in SCENARIOS:
                new = row[scenario]
                transitions.append({
                    "pair": pair,
                    "panel_id": row["panel_id"],
                    "molecule_chembl_id": row.get("molecule_chembl_id", ""),
                    "global_ligand_entity_id": row["global_ligand_entity_id"],
                    "scenario": scenario,
                    "scenario_role": ROLE[scenario],
                    "original_class": row["class"],
                    "new_class": new,
                    "class_changed": "YES" if new != row["class"] else "NO",
                    "change_reason": change_reason(
                        scenario, row["class"], new,
                        row["pA_max"], row["pB_max"], row["pA_median"], row["pB_median"],
                    ),
                    "strict_gray_excluded": "YES" if new == "strict_gray_excluded" else "NO",
                    "pA_max": row["pA_max"],
                    "pB_max": row["pB_max"],
                    "pA_median": row["pA_median"],
                    "pB_median": row["pB_median"],
                    "M0_A_valid": row["M0_A_valid"],
                    "M0_B_valid": row["M0_B_valid"],
                    "M0_score_A": row["M0_score_A"],
                    "M0_score_B": row["M0_score_B"],
                    "score_missing_not_imputed": "YES",
                    "evaluation_set": "fixed_805_not_an_independent_validation",
                })
    write_csv(OUT / "03_label_transitions.csv", transitions, list(transitions[0].keys()))

    metric_rows = []
    for pair in PAIRS:
        obs = [r for r in by_pair[pair] if r.get("activity_eligible") == "1" and r.get("class") in FOUR]
        for scenario in SCENARIOS:
            block = blocks[pair][scenario]
            counts = {k: sum(r[scenario] == k for r in obs) for k in (*FOUR, "strict_gray_excluded")}
            changed = sum(r[scenario] != r["class"] for r in obs)
            for contrast, n_neg in (
                ("dual_vs_A_only", block["n_A_only"]),
                ("dual_vs_B_only", block["n_B_only"]),
                ("summary_min", ""),
            ):
                packed = block[contrast]
                if contrast == "summary_min":
                    unstable = (
                        block["n_dual"] < 10 or block["n_A_only"] < 10 or block["n_B_only"] < 10
                    )
                    n_negative = ""
                else:
                    unstable = block["n_dual"] < 10 or n_neg < 10
                    n_negative = str(n_neg)
                metric_rows.append({
                    "pair": pair,
                    "scenario": scenario,
                    "scenario_role": ROLE[scenario],
                    "definition": DEFINITION[scenario],
                    "contrast": contrast,
                    "score_used": {"dual_vs_A_only": "score_B", "dual_vs_B_only": "score_A", "summary_min": "min(AUC_B,AUC_A)"}[contrast],
                    "n_observations": len(obs),
                    "n_dual_class": counts["dual"],
                    "n_A_only_class": counts["A_only"],
                    "n_B_only_class": counts["B_only"],
                    "n_neither_class": counts["neither"],
                    "n_strict_gray": counts["strict_gray_excluded"],
                    "n_class_changed": changed,
                    "n_dual_scored": block["n_dual"],
                    "n_negative": n_negative,
                    "n_A_only_scored": block["n_A_only"],
                    "n_B_only_scored": block["n_B_only"],
                    "auroc": packed["auroc"],
                    "ci_lo": packed["ci_lo"],
                    "ci_hi": packed["ci_hi"],
                    "unstable_n_below_10": "YES" if unstable and packed["bootstrap_started"] == "YES" else "NO",
                    "bootstrap_started": packed["bootstrap_started"],
                    "reason": packed["reason"],
                    "n_replicates": packed["n_replicates"],
                    "boot_B": BOOT_B,
                    "boot_seed": BOOT_SEED,
                    "rng_policy": "reset_default_rng_per_pair_scenario; draw dual then A-only then B-only",
                    "ci_method": "numpy.percentile linear 2.5 97.5",
                    "member_order": "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv appearance order",
                    "scores": "existing M0 seed42 mode1 higher-better; no redocking; no imputation",
                    "scope": "fixed 805 observations; not an independent validation; not a rebuilt panel",
                })
    write_csv(OUT / "04_sensitivity_metrics.csv", metric_rows, list(metric_rows[0].keys()))

    by_metric = {(r["pair"], r["scenario"], r["contrast"]): r for r in metric_rows}
    comparisons = []
    for pair in PAIRS:
        s0_a = by_metric[(pair, "S0", "dual_vs_A_only")]
        s0_b = by_metric[(pair, "S0", "dual_vs_B_only")]
        weak_s0 = "dual_vs_A_only" if float(s0_a["auroc"]) <= float(s0_b["auroc"]) else "dual_vs_B_only"
        if float(s0_a["auroc"]) == float(s0_b["auroc"]):
            weak_s0 = "tie"
        for scenario in ("S1", "S2", "S3", "S4"):
            a = by_metric[(pair, scenario, "dual_vs_A_only")]
            b = by_metric[(pair, scenario, "dual_vs_B_only")]
            if a["auroc"] and b["auroc"]:
                weak = "dual_vs_A_only" if float(a["auroc"]) <= float(b["auroc"]) else "dual_vs_B_only"
                if float(a["auroc"]) == float(b["auroc"]):
                    weak = "tie"
            else:
                weak = "not_computable"
            for contrast in ("dual_vs_A_only", "dual_vs_B_only", "summary_min"):
                base = by_metric[(pair, "S0", contrast)]
                alt = by_metric[(pair, scenario, contrast)]
                if base["auroc"] and alt["auroc"]:
                    delta = fnum(float(alt["auroc"]) - float(base["auroc"]))
                    below0 = "YES" if float(base["auroc"]) < 0.5 else "NO"
                    below1 = "YES" if float(alt["auroc"]) < 0.5 else "NO"
                    cross05 = "YES" if below0 != below1 else "NO"
                else:
                    delta, below0, below1, cross05 = "", "", "", ""
                comparisons.append({
                    "pair": pair,
                    "scenario": scenario,
                    "scenario_role": ROLE[scenario],
                    "contrast": contrast,
                    "auroc_s0": base["auroc"],
                    "auroc_scenario": alt["auroc"],
                    "delta_auroc_point_only": delta,
                    "difference_has_confidence_interval": "NO",
                    "ci_lo_s0": base["ci_lo"],
                    "ci_hi_s0": base["ci_hi"],
                    "ci_lo_scenario": alt["ci_lo"],
                    "ci_hi_scenario": alt["ci_hi"],
                    "n_dual_scored_s0": base["n_dual_scored"],
                    "n_negative_s0": base["n_negative"],
                    "n_dual_scored_scenario": alt["n_dual_scored"],
                    "n_negative_scenario": alt["n_negative"],
                    "n_class_changed": alt["n_class_changed"],
                    "n_strict_gray": alt["n_strict_gray"],
                    "weaker_direction_s0": weak_s0,
                    "weaker_direction_scenario": weak,
                    "weaker_direction_changed": "YES" if weak != weak_s0 else "NO",
                    "s0_below_0.5": below0,
                    "scenario_below_0.5": below1,
                    "crossed_0.5_descriptive_only": cross05,
                    "unstable_n_below_10": alt["unstable_n_below_10"],
                    "label_change_is_not_a_binding_mechanism_change": "YES",
                })
    write_csv(OUT / "05_sensitivity_comparison.csv", comparisons, list(comparisons[0].keys()))

    enr = {r["activity_id"]: r for r in enriched} if enriched and "activity_id" in enriched[0] else {}
    ppar_by = {(r["pair"], r["panel_id"], r["activity_id"]): r for r in ppar}
    adme_by = {(r["pair"], r["panel_id"], r["activity_id"]): r for r in adme}
    named = {
        "EH40_27", "PAPD_015", "PAPD_082", "PGPA_067", "PGPA_068", "PGPA_069",
        "PM48_04", "PM48_05", "PM48_22",
    }
    doc_gap = {(r["pair"], r["panel_id"]) for r in docs if r.get("all_tied_maxima_both_ends") == "NO"}
    hist_ids = {r["activity_id"] for r in hist}
    review = []
    for row in long_rows:
        keep = (
            row["panel_id"] in named
            or (row["pair"], row["panel_id"], row["activity_id"]) in adme_by
            or ((row["pair"], row["panel_id"]) in doc_gap and row["is_tied_maximum"] == "YES")
            or row["activity_id"] in hist_ids
        )
        if not keep or row["is_tied_maximum"] != "YES" and row["activity_id"] not in hist_ids and (row["pair"], row["panel_id"], row["activity_id"]) not in adme_by:
            if not (row["panel_id"] in named and row["is_tied_maximum"] == "YES"):
                continue
        meta = enr.get(row["activity_id"], {})
        ph = ppar_by.get((row["pair"], row["panel_id"], row["activity_id"]), {})
        am = adme_by.get((row["pair"], row["panel_id"], row["activity_id"]), {})
        action = ph.get("pharmacology_class") or am.get("semantic_judgment") or meta.get("action_type") or ""
        variant = meta.get("mutation") or ""
        species = row.get("assay_organism") or meta.get("assay_organism") or ph.get("assay_organism") or am.get("assay_organism") or ""
        status = "已核实，分类依据明确"
        effect = "不改变 θ=6.0 正式类别"
        evidence = meta.get("assay_description") or ph.get("assay_description") or am.get("assay_description") or ""
        if row["panel_id"] == "EH40_27" and row["endpoint"] == "A":
            status = "已核实数值，但测定变体与对接用野生型结构不是同一对象"
            effect = "若只保留未标明突变的记录，A 端最高约 5.43，假设类别可由 A_only 变为 neither；未标明突变不是已核实野生型"
            if not variant:
                variant = "决定不低于 6 的记录为 G719C 与 L858R/T790M"
        elif row["panel_id"] == "PAPD_015" and row["endpoint"] == "A":
            status = "已核实数值，但药理学作用方式不同"
            effect = "A 端是 PPAR-alpha 拮抗，不把该观察写成激动剂；θ=6.0 类别仍为 A_only"
            action = action or "PPAR拮抗剂相关功能实验"
        elif row["panel_id"] == "PAPD_082" and row["endpoint"] == "A":
            status = "已核实数值，但最大值来自仓鼠测定"
            effect = "仓鼠记录决定最大值；同端人源记录仍不低于 6，dual 不变"
        elif row["panel_id"] in {"PGPA_067", "PGPA_068", "PGPA_069"}:
            status = "已核实来源，但作用方式仍不明确"
            effect = "FRET 作用方式未确认。没有其他活性行。不把记录判为无效，也不自行归为激动或拮抗"
            action = action or "作用方式未确认"
        elif row["activity_id"] in hist_ids:
            status = "已核实数值；历史排除的实质理由未证实"
            effect = "原库最大值与后续聚合不同，但冻结类别不变。不把原始记录写成错误"
            evidence = (evidence + " " + next(h["substantive_reason"] for h in hist if h["activity_id"] == row["activity_id"])).strip()
        elif (row["pair"], row["panel_id"], row["activity_id"]) in adme_by:
            status = "已核实，分类依据明确"
            effect = "描述是靶点酶活或激酶活。assay_type=A 与描述不一致，不据此删除"
            action = am.get("semantic_judgment", action)
            evidence = am.get("assay_description", evidence)
        elif (row["pair"], row["panel_id"]) in doc_gap:
            status = "已核实来源，但指定文献标识缺失"
            effect = "原始活性存在。缺少 DOI、PubMed 或专利号不等于没有实验，也不改变类别"
            srcdoc = next(d for d in docs if d["pair"] == row["pair"] and d["panel_id"] == row["panel_id"])
            evidence = srcdoc.get("dataset_titles", "") + " | " + srcdoc.get("dataset_sources", "")
        if ph.get("pharmacology_class") == "不能明确确定药理作用的实验" and row["panel_id"] not in {"PGPA_067", "PGPA_068", "PGPA_069"}:
            if status == "已核实，分类依据明确":
                status = "已核实来源，但作用方式仍不明确"
        obs = next(r for r in eligible if r["pair"] == row["pair"] and r["panel_id"] == row["panel_id"])
        review.append({
            "pair": row["pair"],
            "panel_id": row["panel_id"],
            "endpoint": row["endpoint"],
            "activity_id": row["activity_id"],
            "assay_id": row["assay_chembl_id"],
            "document_id": row["document_chembl_id"],
            "standard_type": row["standard_type"],
            "pChEMBL": row["pchembl_value"],
            "species": species,
            "variant": variant,
            "pharmacological_action": action,
            "original_class": obs["class"],
            "evidence": evidence[:500],
            "verification_status": status,
            "possible_label_effect": effect,
            "is_tied_maximum": row["is_tied_maximum"],
        })
    write_csv(OUT / "02_activity_record_review.csv", review, list(review[0].keys()))

    named_rows = []
    median_changed = sorted({
        (r["pair"], r["panel_id"], r["original_class"], r["new_class"])
        for r in transitions
        if r["scenario"] == "S4" and r["class_changed"] == "YES"
    })
    watch = named | {pid for _, pid, _, _ in median_changed}
    for pair in PAIRS:
        for row in by_pair[pair]:
            if row["panel_id"] not in watch:
                continue
            if row.get("activity_eligible") != "1":
                continue
            for scenario in SCENARIOS:
                new = row[scenario]
                in_b = new == "dual" and row["M0_A_valid"] == "1" and row["M0_B_valid"] == "1" or new == "A_only" and row["M0_B_valid"] == "1"
                in_a = new == "dual" and row["M0_A_valid"] == "1" and row["M0_B_valid"] == "1" or new == "B_only" and row["M0_A_valid"] == "1"
                named_rows.append({
                    "pair": pair,
                    "panel_id": row["panel_id"],
                    "scenario": scenario,
                    "original_class": row["class"],
                    "new_class": new,
                    "class_changed": "YES" if new != row["class"] else "NO",
                    "in_dual_vs_A_only": "YES" if in_b and new in {"dual", "A_only"} else "NO",
                    "in_dual_vs_B_only": "YES" if in_a and new in {"dual", "B_only"} else "NO",
                    "exited_because_gray": "YES" if new == "strict_gray_excluded" else "NO",
                    "score_missing": "YES" if row["M0_A_valid"] != "1" or row["M0_B_valid"] != "1" else "NO",
                    "median_changed_class": "YES" if row["S4"] != row["class"] else "NO",
                    "named_special": "YES" if row["panel_id"] in named else "NO",
                })
    write_csv(OUT / "10_special_observation_scenarios.csv", named_rows, list(named_rows[0].keys()))

    summary = {
        "s0_acceptance": "PASS",
        "n_acceptance_fields": len(acceptance),
        "n_review_rows": len(review),
        "median_class_changes": [
            {"pair": a, "panel_id": b, "original": c, "median_class": d}
            for a, b, c, d in median_changed
        ],
        "pairs": {},
    }
    for pair in PAIRS:
        summary["pairs"][pair] = {}
        for scenario in SCENARIOS:
            summary["pairs"][pair][scenario] = {
                contrast: {
                    k: by_metric[(pair, scenario, contrast)][k]
                    for k in (
                        "auroc", "ci_lo", "ci_hi", "n_dual_scored", "n_negative",
                        "n_dual_class", "n_A_only_class", "n_B_only_class", "n_neither_class",
                        "n_strict_gray", "n_class_changed", "unstable_n_below_10", "bootstrap_started",
                    )
                }
                for contrast in ("dual_vs_A_only", "dual_vs_B_only", "summary_min")
            }
    (OUT / "summary_for_report.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print("S0_ACCEPTANCE_PASS")
    print("median_class_changes", len(median_changed))
    print("review_rows", len(review))
    return 0


if __name__ == "__main__":
    sys.exit(main())
