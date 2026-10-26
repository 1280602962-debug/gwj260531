#!/usr/bin/env python3
"""CAP5 sensitivity statistics. Reads the frozen membership and does not rewrite it."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

ROOT = Path("/home/gwj/q110_v42")
DTD = ROOT / "Dual_Target_Docking"
OUT = DTD / "results/current/chemical_diversity_cap5_20261023"
sys.path.insert(0, str(DTD / "scripts/rerun_v4_2"))
sys.path.insert(0, str(DTD / "extension_jcim_validation/pik3ca_mtor_q110_v42"))

from jcim_stage2_compute import module_b, module_c, module_d, module_i  # noqa: E402
from jcim_stage2_lib import alias_ids, build_layer_members, fold_keyset  # noqa: E402
from q110_stats import paired_delta, single_method  # noqa: E402

METHODS = ("M0", "M1", "M1b", "M2", "M3")
DELTAS = (("M1", "M0"), ("M1b", "M1"), ("M2", "M0"), ("M3", "M0"))
PAIRS_ORDER = None


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        path.write_text("status,reason\nBLOCKED,no rows\n", encoding="utf-8")
        return
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def blank(value) -> bool:
    return value is None or str(value).strip() in {"", "NA", "None", "nan"}


def as_float(value):
    if blank(value):
        return None
    return float(value)


def same(a, b) -> bool:
    if blank(a) and blank(b):
        return True
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)


def point_class(original, cap, null: float, band: float) -> str:
    if original is None or cap is None:
        return "NOT_ESTIMABLE"
    effect_o = original - null
    effect_c = cap - null
    if abs(effect_o) < band and abs(effect_c) < band:
        return "NEAR_NULL"
    if effect_o * effect_c < 0 and abs(effect_o) >= band and abs(effect_c) >= band:
        return "REVERSED"
    if abs(cap - original) < band:
        return "STABLE"
    if abs(effect_c) < abs(effect_o):
        return "ATTENUATED"
    if abs(effect_c) > abs(effect_o):
        return "AMPLIFIED"
    return "STABLE"


def evidence(olo, ohi, clo, chi, null: float) -> str:
    if None in (olo, ohi, clo, chi):
        return "NOT_ESTIMABLE"
    original = olo <= null <= ohi
    current = clo <= null <= chi
    if original and current:
        return "BOTH_INTERVALS_INCLUDE_NULL"
    if original:
        return "ORIGINAL_INTERVAL_INCLUDES_NULL"
    if current:
        return "CAP5_INTERVAL_INCLUDES_NULL"
    return "NEITHER_INTERVAL_INCLUDES_NULL"


def compare_row(layer, key, fields, left, right) -> dict:
    mismatches = [name for name in fields if not same(left.get(name, ""), right.get(name, ""))]
    return {
        "layer": layer,
        "key": key,
        "result": "PASS" if not mismatches else "FAIL",
        "fields": ",".join(mismatches),
        "detail": "" if not mismatches else "; ".join(
            f"{name}:{left.get(name, '')}|{right.get(name, '')}" for name in mismatches[:6]
        ),
    }


def rows_for(pop, pair):
    return [row for row in pop if row["pair"] == pair]


def primary_bundle(pop) -> tuple[list[dict], list[dict]]:
    metrics, deltas = [], []
    for pair in sorted({row["pair"] for row in pop}):
        subset = rows_for(pop, pair)
        for method in METHODS:
            result = single_method(subset, method)
            for contrast, key, n_neg in (
                ("dual_vs_A_only", "dual_vs_A_only", result["n_A_only"]),
                ("dual_vs_B_only", "dual_vs_B_only", result["n_B_only"]),
                ("summary_min", "summary_min", ""),
            ):
                packed = result[key]
                metrics.append({
                    "pair": pair,
                    "method": method,
                    "contrast": contrast,
                    "n_dual": result["n_dual"],
                    "n_negative": n_neg,
                    "n_A_only": result["n_A_only"],
                    "n_B_only": result["n_B_only"],
                    "auroc": packed["value"],
                    "ci_lo": packed["ci_lo"],
                    "ci_hi": packed["ci_hi"],
                    "n_replicates": packed["n_replicates"],
                    "status": "ESTIMABLE" if packed["value"] != "" else "NOT_ESTIMABLE",
                })
        for left, right in DELTAS:
            result = paired_delta(subset, left, right)
            result["pair"] = pair
            result["status"] = "ESTIMABLE" if result.get("delta_summary_min", "") != "" else "NOT_ESTIMABLE"
            deltas.append(result)
    return metrics, deltas


def index(rows, *keys):
    return {"|".join(str(row[key]) for key in keys): row for row in rows}


def top_ids(ligs) -> dict[str, list[str]]:
    grouped = defaultdict(list)
    for row in ligs:
        if row.get("in_topk") == "YES":
            grouped[f"{row['pair']}|{row['method']}|{row['population_kind']}"].append(row["canonical_ligand_id"])
    return dict(grouped)


def fold_audit(pop, folds, alias_rows, reference_unassigned) -> list[dict]:
    keys = fold_keyset(folds)
    aliases = alias_ids(alias_rows)
    out = []
    for pair in sorted({row["pair"] for row in pop}):
        for arm in ("D_vs_A", "D_vs_B"):
            layers = build_layer_members(pop, aliases, keys, pair, arm)
            missing = [row["canonical_ligand_id"] for row in layers["labeled"] if (pair, arm, row["canonical_ligand_id"]) not in keys]
            original = set(reference_unassigned.get((pair, arm), []))
            new_fail = [item for item in missing if item not in original]
            out.append({
                "pair": pair,
                "arm": arm,
                "n_labeled": len(layers["labeled"]),
                "n_chemistry_oof": len(layers["chemistry_oof"]),
                "originally_missing_retained": ";".join(item for item in missing if item in original),
                "new_fold_join_failures": ";".join(new_fail),
                "qa_status": "QA_FAIL" if new_fail else "PASS",
            })
    return out


def reference_missing(pop, folds, alias_rows) -> dict:
    keys = fold_keyset(folds)
    aliases = alias_ids(alias_rows)
    out = {}
    for pair in sorted({row["pair"] for row in pop}):
        for arm in ("D_vs_A", "D_vs_B"):
            layers = build_layer_members(pop, aliases, keys, pair, arm)
            out[(pair, arm)] = [row["canonical_ligand_id"] for row in layers["labeled"] if (pair, arm, row["canonical_ligand_id"]) not in keys]
    return out


def chemistry_status(unit, qa) -> str:
    if qa == "QA_FAIL":
        return "QA_FAIL"
    reason = unit.get("blocking_reason") or ""
    if reason in {"TRAIN_SINGLE_CLASS", "OOF_SINGLE_CLASS", "EMPTY_CHEMISTRY_OOF", "CHEMISTRY_SMILES_PARSE_FAILURE"} or str(unit.get("incomplete_reason") or "").startswith("INCOMPLETE"):
        return "NOT_ESTIMABLE"
    if unit.get("test_single_class_folds"):
        return "TEST_SINGLE_CLASS_WARNING"
    return "ESTIMABLE"


def comparison_record(layer, pair, method, contrast, original, cap, null, band, n_orig, n_cap) -> dict:
    return {
        "layer": layer,
        "pair": pair,
        "method": method,
        "contrast": contrast,
        "original_value": "" if original[0] is None else f"{original[0]:.10f}",
        "original_ci_lo": "" if original[1] is None else f"{original[1]:.10f}",
        "original_ci_hi": "" if original[2] is None else f"{original[2]:.10f}",
        "original_n": n_orig,
        "cap5_value": "" if cap[0] is None else f"{cap[0]:.10f}",
        "cap5_ci_lo": "" if cap[1] is None else f"{cap[1]:.10f}",
        "cap5_ci_hi": "" if cap[2] is None else f"{cap[2]:.10f}",
        "cap5_n": n_cap,
        "point_change_class": point_class(original[0], cap[0], null, band),
        "evidence_status": evidence(original[1], original[2], cap[1], cap[2], null),
    }


def fmt_triplet(row, value_key, lo_key, hi_key):
    return (as_float(row.get(value_key)), as_float(row.get(lo_key)), as_float(row.get(hi_key)))


def write_reports(summary: dict) -> None:
    excluded = summary["excluded"]
    lines = []
    lines.append("# CAP5 方法与选择\n")
    lines.append("这是正式八对评价集冻结之后的事后化学多样性敏感性分析。它不是预先登记的抽样规则，也不是八对分子当初共同使用的选样算法。\n")
    lines.append("骨架调用现行 `jcim_stage2_compute._murcko`：`GetScaffoldForMol` 之后 `MolToSmiles`。空骨架且环数为 0 记为无环，不参与非空骨架的 5 个名额。上限只在同一 pair、同一实验类别、同一非空骨架内计算。超额组内按 UTF-8、`|` 分隔的 `BM5_SENSITIVITY_20261023_V1|pair|class|global_ligand_entity_id|panel_id` SHA-256 升序保留前 5 个，哈希相同再按 `panel_id`。哈希不改变统计输入的原行顺序。\n")
    lines.append(f"父提交 `{summary['head']}`。合格观察 867，保留 {summary['retained']}，排除 {summary['excluded_n']}，无环 {summary['acyclic']}。排除名单：{excluded}。neither 的最大同骨架重复均不超过 5，因此 neither 没有排除。\n")
    lines.append("描述性阈值在敏感性统计之前写入 `CAP5_SELECTION_FREEZE.json`：AUROC 的零效应是 0.5，差值的零效应是 0，点变化带宽是 0.02。`point_change_class` 与 `evidence_status` 分开。这些标签只描述点估计移动和区间是否包含零效应，不是显著性或等效性结论。\n")
    lines.append("PRIMARY 调用 `q110_stats.single_method` 和 `paired_delta`。Stage 2 只调用 `module_b`、`module_c`、`module_d`、`module_i`。CAP5 的 Module I 使用 CAP5 八对重新计算的 `paired_delta`。没有调用 `execute()`、`run_stage2_modules()` 或 `run_primary_metrics_v42.py`。\n")
    (OUT / "CAP5_METHOD_AND_SELECTION_zh.md").write_text("\n".join(lines), encoding="utf-8")

    comp = ["# CAP5 与正式结果的对照\n", "完整点估计、区间和样本量见 `CAP5_PRIMARY_VS_ORIGINAL.csv` 及各模块表。下表只摘录 M0 `summary_min`。\n", "| pair | 正式 | CAP5 | 类别 | 区间状态 |", "|---|---:|---:|---|---|"]
    for row in summary["m0_rows"]:
        comp.append(f"| {row['pair']} | {row['original_value']} | {row['cap5_value']} | {row['point_change_class']} | {row['evidence_status']} |")
    comp.append("\n变化类别计数：" + json.dumps(summary["class_counts"], ensure_ascii=False))
    comp.append("\n反向变化：" + (summary["reversed"] or "无"))
    (OUT / "CAP5_RESULTS_COMPARISON_zh.md").write_text("\n".join(comp) + "\n", encoding="utf-8")

    assess = ["# CAP5 科学评估\n"]
    assess.append(f"无过滤回归：{summary['regression_text']}。\n")
    assess.append("CAP5 只回答一件事：在已经冻结的八对评价集里，把同一 pair、同一类别、同一非空 Bemis–Murcko 骨架限制到 5 个分子后，方向性对接评价、化学基线和双靶排序是否改变。它不消除活性阈值、命名偏好、文献来源或跨 pair 化学重叠带来的其他偏差。\n")
    assess.append(summary["assessment"] + "\n")
    assess.append("五随机种子、独立 GNINA、RTMScore 重算、共晶再对接、替代受体、外部验证、活性重新聚合以及其他 Stage 2 模块本轮未重评，不能把旧数字当成 CAP5 结果。\n")
    (OUT / "CAP5_SCIENTIFIC_ASSESSMENT_zh.md").write_text("\n".join(assess), encoding="utf-8")

    guide = ["# 稿件修改提示\n"]
    guide.append("本轮没有修改稿件或补充材料。若正文提到化学多样性，应把它写成冻结评价集上的事后敏感性分析，而不是八对分子原有的统一抽样规则。EGFR/HER2 历史构建和 AChE/BChE 的弱替代键都不能改写成当前 CAP5 规则。\n")
    guide.append("正文可以保留正式人群的主结果。CAP5 宜放入补充材料，并列 `point_change_class` 和 `evidence_status`，并保留点估计、区间和样本量。0.02 带宽不能写成显著性或等效性界值。\n")
    guide.append("需要说明：排除由骨架重复和冻结哈希决定，不按对接分数保留；无环分子不占非空骨架名额；Module D 在 CAP5 上重拟合，排除成员既不训练也不测试；原本缺 fold 的成员沿用原规则退出化学 OOF。\n")
    guide.append("本轮未重评的模块应写“未在本轮重评”，不要补写 CAP5 数字。\n")
    (OUT / "CAP5_MANUSCRIPT_REVISION_GUIDE_zh.md").write_text("\n".join(guide), encoding="utf-8")


def main() -> int:
    freeze = json.loads((OUT / "CAP5_SELECTION_FREEZE.json").read_text(encoding="utf-8"))
    membership_path = OUT / "CAP5_MEMBERSHIP_867.csv"
    if sha256_file(membership_path) != freeze["membership_sha256"]:
        raise SystemExit("MEMBERSHIP_HASH_MISMATCH")
    band = float(freeze["descriptive_thresholds"]["point_band"])
    null_auc = float(freeze["descriptive_thresholds"]["auroc_null"])
    null_delta = float(freeze["descriptive_thresholds"]["delta_null"])
    membership = read_csv(membership_path)
    retained = {(row["pair"], row["panel_id"]) for row in membership if row["retained_cap5"] == "YES"}
    population = read_csv(DTD / "results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")
    mapping = read_csv(DTD / "results/current/pair_ligand_mapping.csv")
    folds = read_csv(DTD / "results/canonical/model_fold_assignments.csv")
    alias_rows = read_csv(DTD / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/parent_alias_registry.csv")
    formal_metrics = read_csv(DTD / "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv")
    formal_deltas = read_csv(DTD / "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv")
    formal_b = read_csv(DTD / "results/jcim_stage2/module_B_negative_class.csv")
    formal_c = read_csv(DTD / "results/jcim_stage2/module_C_wrong_pocket.csv")
    formal_d = read_csv(DTD / "results/jcim_stage2/module_D_chemistry_units.csv")
    formal_oof = read_csv(DTD / "results/jcim_stage2/module_D_chemistry_oof_ligands.csv")
    formal_i = read_csv(DTD / "results/jcim_stage2/module_I_ligand_ranks.csv")
    formal_cmp = read_csv(DTD / "results/jcim_stage2/module_I_method_compare.csv")
    cap_pop = [row for row in population if (row["pair"], row["panel_id"]) in retained]
    if [row["panel_id"] for row in cap_pop] != [row["panel_id"] for row in membership if row["retained_cap5"] == "YES"]:
        # Pair-aware order check follows; panel IDs can repeat across pairs.
        pass
    cap_keys = [(row["pair"], row["panel_id"]) for row in cap_pop]
    mem_keys = [(row["pair"], row["panel_id"]) for row in membership if row["retained_cap5"] == "YES"]
    if cap_keys != mem_keys:
        raise SystemExit("CAP5_ORDER_MISMATCH")
    checks = []
    print("primary regression", flush=True)
    base_metrics, base_deltas = primary_bundle(population)
    metric_fields = ["n_dual", "n_negative", "n_A_only", "n_B_only", "auroc", "ci_lo", "ci_hi", "n_replicates"]
    base_m = index(base_metrics, "pair", "method", "contrast")
    for row in formal_metrics:
        checks.append(compare_row("PRIMARY", f"{row['pair']}|{row['method']}|{row['contrast']}", metric_fields, row, base_m[f"{row['pair']}|{row['method']}|{row['contrast']}"]))
    delta_fields = ["n_dual_four_sided", "n_A_only_both_score_B", "n_B_only_both_score_A", "delta_summary_min", "summary_min_x", "summary_min_y", "ci_lo", "ci_hi"]
    base_d = index(base_deltas, "pair", "delta")
    for row in formal_deltas:
        checks.append(compare_row("PRIMARY_DELTA", f"{row['pair']}|{row['delta']}", delta_fields, row, base_d[f"{row['pair']}|{row['delta']}"]))
    print("module B/C regression", flush=True)
    base_b = module_b(population)
    base_c = module_c(population)
    b_fields = ["n_dual", "n_single", "n_neither", "auc_single", "auc_neither", "delta_negative", "auc_single_ci_lo", "auc_single_ci_hi", "auc_neither_ci_lo", "auc_neither_ci_hi", "delta_negative_ci_lo", "delta_negative_ci_hi", "n_replicates"]
    got_b = index(base_b, "pair", "direction")
    for row in formal_b:
        checks.append(compare_row("MODULE_B", f"{row['pair']}|{row['direction']}", b_fields, row, got_b[f"{row['pair']}|{row['direction']}"]))
    c_fields = ["n_dual", "n_negative_directional_m0", "n_negative_both_finite", "auc_correct", "auc_wrong", "delta_pocket", "auc_correct_ci_lo", "auc_correct_ci_hi", "auc_wrong_ci_lo", "auc_wrong_ci_hi", "delta_pocket_ci_lo", "delta_pocket_ci_hi", "n_replicates"]
    got_c = index(base_c, "pair", "arm")
    for row in formal_c:
        checks.append(compare_row("MODULE_C", f"{row['pair']}|{row['arm']}", c_fields, row, got_c[f"{row['pair']}|{row['arm']}"]))
    print("module D/I regression", flush=True)
    tables = {"population": population, "alias": alias_rows, "folds": folds, "mapping": mapping}
    base_units, base_oof = module_d(tables)
    d_fields = ["n_labeled_population", "n_m0_directional_population", "n_final_chemistry_oof", "fold_unassigned_in_labeled", "fold_unassigned_after_m0", "auroc_physchem", "auroc_ecfp4", "auroc_ecfp4_m0", "auroc_SUPPORTING_CHEMISTRY_MATCHED_M0", "delta_ecfp4_plus_m0", "test_single_class_folds", "blocking_reason"]
    got_d = index(base_units, "pair", "arm")
    for row in formal_d:
        checks.append(compare_row("MODULE_D", f"{row['pair']}|{row['arm']}", d_fields, row, got_d[f"{row['pair']}|{row['arm']}"]))
    oof_fields = ["class", "fold_id", "oof_physchem", "oof_ecfp4", "oof_ecfp4_m0", "relevant_m0"]
    got_oof = index(base_oof, "pair", "arm", "canonical_ligand_id")
    formal_oof_index = index(formal_oof, "pair", "arm", "canonical_ligand_id")
    checks.append({
        "layer": "MODULE_D_OOF_IDS",
        "key": "all",
        "result": "PASS" if set(got_oof) == set(formal_oof_index) else "FAIL",
        "fields": "canonical_ligand_id",
        "detail": f"computed={len(got_oof)} formal={len(formal_oof_index)}",
    })
    for key, row in formal_oof_index.items():
        if key not in got_oof:
            continue
        item = compare_row("MODULE_D_OOF", key, oof_fields, row, got_oof[key])
        if item["result"] != "PASS":
            checks.append(item)
    if not any(row["layer"] == "MODULE_D_OOF" and row["result"] == "FAIL" for row in checks):
        checks.append({"layer": "MODULE_D_OOF", "key": "all_values", "result": "PASS", "fields": "", "detail": str(len(formal_oof_index))})
    base_desc, base_ligs, base_cmps = module_i(population, formal_deltas)
    got_top = top_ids(base_ligs)
    formal_top = top_ids(formal_i)
    top_keys = set(got_top) | set(formal_top)
    top_fail = [key for key in sorted(top_keys) if got_top.get(key) != formal_top.get(key)]
    checks.append({
        "layer": "MODULE_I_TOPK",
        "key": "all",
        "result": "PASS" if not top_fail else "FAIL",
        "fields": "canonical_ligand_id",
        "detail": "" if not top_fail else ";".join(top_fail[:8]),
    })
    cmp_fields = ["n_rank_common", "n_directional_common_dual", "delta_summary_min", "delta_top10_single_target_fraction", "concordance", "link_source"]
    got_cmp = index(base_cmps, "pair", "method")
    for row in formal_cmp:
        checks.append(compare_row("MODULE_I_COMPARE", f"{row['pair']}|{row['method']}", cmp_fields, row, got_cmp[f"{row['pair']}|{row['method']}"]))
    write_csv(OUT / "UNFILTERED_REGRESSION.csv", checks, ["layer", "key", "result", "fields", "detail"])
    layers = sorted({row["layer"] for row in checks})
    layer_ok = {layer: all(row["result"] == "PASS" for row in checks if row["layer"] == layer) for layer in layers}
    print(json.dumps(layer_ok), flush=True)

    def blocked(name, reason):
        return [{"status": "BLOCKED", "reason": reason, "layer": name}]

    if layer_ok.get("PRIMARY") and layer_ok.get("PRIMARY_DELTA"):
        print("CAP5 primary", flush=True)
        cap_metrics, cap_deltas = primary_bundle(cap_pop)
    else:
        cap_metrics, cap_deltas = blocked("PRIMARY", "unfiltered primary regression failed"), []
    if layer_ok.get("MODULE_B"):
        print("CAP5 B", flush=True)
        cap_b = module_b(cap_pop)
    else:
        cap_b = blocked("MODULE_B", "unfiltered module B regression failed")
    if layer_ok.get("MODULE_C"):
        print("CAP5 C", flush=True)
        cap_c = module_c(cap_pop)
    else:
        cap_c = blocked("MODULE_C", "unfiltered module C regression failed")
    original_missing = reference_missing(population, folds, alias_rows)
    qa_rows = fold_audit(cap_pop, folds, alias_rows, original_missing)
    qa_by = {(row["pair"], row["arm"]): row["qa_status"] for row in qa_rows}
    if layer_ok.get("MODULE_D") and layer_ok.get("MODULE_D_OOF") and layer_ok.get("MODULE_D_OOF_IDS"):
        print("CAP5 D", flush=True)
        cap_tables = {"population": cap_pop, "alias": alias_rows, "folds": folds, "mapping": mapping}
        cap_units, cap_oof = module_d(cap_tables)
    else:
        cap_units, cap_oof = blocked("MODULE_D", "unfiltered module D regression failed"), []
    if layer_ok.get("MODULE_I_TOPK") and layer_ok.get("MODULE_I_COMPARE") and layer_ok.get("PRIMARY_DELTA"):
        print("CAP5 I", flush=True)
        cap_desc, cap_ligs, cap_cmps = module_i(cap_pop, cap_deltas)
    else:
        cap_desc, cap_ligs, cap_cmps = blocked("MODULE_I", "unfiltered module I or delta regression failed"), [], []

    write_csv(OUT / "CAP5_PRIMARY_DIRECTIONAL_METRICS.csv", cap_metrics)
    write_csv(OUT / "CAP5_PRIMARY_METHOD_DELTAS.csv", cap_deltas)
    write_csv(OUT / "CAP5_MODULE_B_NEGATIVE_CLASS.csv", cap_b)
    write_csv(OUT / "CAP5_MODULE_C_POCKET_CORRESPONDENCE.csv", cap_c)
    for row in cap_units:
        if "pair" in row and "arm" in row:
            row["cap5_status"] = chemistry_status(row, qa_by.get((row["pair"], row["arm"]), "QA_FAIL"))
            row["new_fold_join_failures"] = next(item["new_fold_join_failures"] for item in qa_rows if item["pair"] == row["pair"] and item["arm"] == row["arm"])
            row["originally_missing_retained"] = next(item["originally_missing_retained"] for item in qa_rows if item["pair"] == row["pair"] and item["arm"] == row["arm"])
    write_csv(OUT / "CAP5_MODULE_D_CHEMISTRY_UNITS.csv", cap_units)
    write_csv(OUT / "CAP5_MODULE_D_CHEMISTRY_OOF.csv", cap_oof)
    for row in qa_rows:
        unit = next((item for item in cap_units if item.get("pair") == row["pair"] and item.get("arm") == row["arm"]), {})
        row["test_single_class_folds"] = unit.get("test_single_class_folds", "")
        row["blocking_reason"] = unit.get("blocking_reason", "")
        row["cap5_status"] = unit.get("cap5_status", "BLOCKED")
    write_csv(OUT / "CAP5_MODULE_D_FOLD_AUDIT.csv", qa_rows)
    write_csv(OUT / "CAP5_MODULE_I_TOPK.csv", cap_desc)
    write_csv(OUT / "CAP5_MODULE_I_METHOD_COMPARISON.csv", cap_cmps)
    write_csv(OUT / "CAP5_MODULE_I_RANKS.csv", cap_ligs)

    versus = []
    if cap_metrics and "auroc" in cap_metrics[0]:
        cap_m = index(cap_metrics, "pair", "method", "contrast")
        for row in formal_metrics:
            key = f"{row['pair']}|{row['method']}|{row['contrast']}"
            got = cap_m[key]
            versus.append(comparison_record(
                "PRIMARY", row["pair"], row["method"], row["contrast"],
                fmt_triplet(row, "auroc", "ci_lo", "ci_hi"),
                fmt_triplet(got, "auroc", "ci_lo", "ci_hi"),
                null_auc, band,
                f"{row['n_dual']}/{row['n_negative']}",
                f"{got['n_dual']}/{got['n_negative']}",
            ))
        cap_delta_index = index(cap_deltas, "pair", "delta")
        for row in formal_deltas:
            got = cap_delta_index[f"{row['pair']}|{row['delta']}"]
            versus.append(comparison_record(
                "PRIMARY_DELTA", row["pair"], row["delta"], "summary_min",
                fmt_triplet(row, "delta_summary_min", "ci_lo", "ci_hi"),
                fmt_triplet(got, "delta_summary_min", "ci_lo", "ci_hi"),
                null_delta, band,
                row["n_dual_four_sided"], got.get("n_dual_four_sided", ""),
            ))
    if cap_b and "auc_single" in cap_b[0]:
        cap_b_index = index(cap_b, "pair", "direction")
        for row in formal_b:
            got = cap_b_index[f"{row['pair']}|{row['direction']}"]
            for contrast, value_key, lo_key, hi_key, null in (
                ("auc_single", "auc_single", "auc_single_ci_lo", "auc_single_ci_hi", null_auc),
                ("auc_neither", "auc_neither", "auc_neither_ci_lo", "auc_neither_ci_hi", null_auc),
                ("delta_negative", "delta_negative", "delta_negative_ci_lo", "delta_negative_ci_hi", null_delta),
            ):
                versus.append(comparison_record("MODULE_B", row["pair"], "M0", f"{row['direction']}:{contrast}", fmt_triplet(row, value_key, lo_key, hi_key), fmt_triplet(got, value_key, lo_key, hi_key), null, band, f"{row['n_dual']}/{row['n_single']}/{row['n_neither']}", f"{got['n_dual']}/{got['n_single']}/{got['n_neither']}"))
    if cap_c and "auc_correct" in cap_c[0]:
        cap_c_index = index(cap_c, "pair", "arm")
        for row in formal_c:
            got = cap_c_index[f"{row['pair']}|{row['arm']}"]
            for contrast, value_key, lo_key, hi_key, null in (
                ("auc_correct", "auc_correct", "auc_correct_ci_lo", "auc_correct_ci_hi", null_auc),
                ("auc_wrong", "auc_wrong", "auc_wrong_ci_lo", "auc_wrong_ci_hi", null_auc),
                ("delta_pocket", "delta_pocket", "delta_pocket_ci_lo", "delta_pocket_ci_hi", null_delta),
            ):
                versus.append(comparison_record("MODULE_C", row["pair"], "M0", f"{row['arm']}:{contrast}", fmt_triplet(row, value_key, lo_key, hi_key), fmt_triplet(got, value_key, lo_key, hi_key), null, band, row["n_dual"], got["n_dual"]))
    if cap_units and "auroc_ecfp4" in cap_units[0]:
        cap_d_index = index(cap_units, "pair", "arm")
        for row in formal_d:
            got = cap_d_index[f"{row['pair']}|{row['arm']}"]
            if got.get("cap5_status") in {"QA_FAIL", "NOT_ESTIMABLE", "BLOCKED"}:
                versus.append({
                    "layer": "MODULE_D", "pair": row["pair"], "method": got.get("cap5_status"), "contrast": row["arm"],
                    "original_value": row.get("auroc_ecfp4", ""), "original_ci_lo": "", "original_ci_hi": "", "original_n": row.get("n_final_chemistry_oof", ""),
                    "cap5_value": got.get("auroc_ecfp4", ""), "cap5_ci_lo": "", "cap5_ci_hi": "", "cap5_n": got.get("n_final_chemistry_oof", ""),
                    "point_change_class": "NOT_ESTIMABLE" if got.get("cap5_status") == "NOT_ESTIMABLE" else got.get("cap5_status"),
                    "evidence_status": "NOT_ESTIMABLE" if got.get("cap5_status") != "QA_FAIL" else "QA_FAIL",
                })
                continue
            for contrast, key, null in (("auroc_ecfp4", "auroc_ecfp4", null_auc), ("auroc_ecfp4_m0", "auroc_ecfp4_m0", null_auc), ("delta_ecfp4_plus_m0", "delta_ecfp4_plus_m0", null_delta)):
                record = comparison_record("MODULE_D", row["pair"], "chemistry", f"{row['arm']}:{contrast}", (as_float(row.get(key)), None, None), (as_float(got.get(key)), None, None), null, band, row.get("n_final_chemistry_oof", ""), got.get("n_final_chemistry_oof", ""))
                record["evidence_status"] = "NOT_ESTIMABLE"
                versus.append(record)
    if cap_desc and "top10_dual_fraction" in cap_desc[0]:
        formal_desc = read_csv(DTD / "results/jcim_stage2/module_I_descriptive.csv")
        cap_desc_index = index(cap_desc, "pair", "method", "population_kind")
        for row in formal_desc:
            if row["population_kind"] != "method_specific":
                continue
            got = cap_desc_index.get(f"{row['pair']}|{row['method']}|method_specific")
            if not got:
                continue
            original_effect = as_float(row["top10_dual_fraction"]) - as_float(row["expected_dual_fraction"])
            cap_effect = as_float(got["top10_dual_fraction"]) - as_float(got["expected_dual_fraction"])
            versus.append(comparison_record(
                "MODULE_I", row["pair"], row["method"], "method_specific_enrichment",
                (original_effect, None, None), (cap_effect, None, None), null_delta, band,
                f"{row['k']}/{row['n_total']};top={row['n_dual_top']}/{row['n_A_only_top']}/{row['n_B_only_top']}",
                f"{got['k']}/{got['n_total']};top={got['n_dual_top']}/{got['n_A_only_top']}/{got['n_B_only_top']}",
            ))
    write_csv(OUT / "CAP5_PRIMARY_VS_ORIGINAL.csv", versus)

    counts = defaultdict(int)
    reversed_rows = []
    for row in versus:
        counts[row["point_change_class"]] += 1
        if row["point_change_class"] == "REVERSED":
            reversed_rows.append(f"{row['layer']}:{row['pair']}:{row['method']}:{row['contrast']}")
    m0_rows = [row for row in versus if row["layer"] == "PRIMARY" and row["method"] == "M0" and row["contrast"] == "summary_min"]
    excluded_txt = ", ".join(f"{row['pair']} {row['panel_id']} {row['class']}" for row in membership if row["retained_cap5"] == "NO")
    fail_layers = [layer for layer, ok in layer_ok.items() if not ok]
    if not fail_layers and all(row["qa_status"] == "PASS" for row in qa_rows):
        regression_text = "PRIMARY、方法差值、Module B、Module C、Module D OOF 和 Module I Top-k 成员均与正式表一致"
        audit_result = "PASS"
    elif fail_layers:
        regression_text = "未复现的层：" + ", ".join(fail_layers)
        audit_result = "FAIL"
    else:
        regression_text = "统计层复现，但 Module D 出现新的 fold 连接失败"
        audit_result = "FAIL"
    notable = [row for row in m0_rows if row["point_change_class"] in {"REVERSED", "AMPLIFIED", "ATTENUATED"}]
    assessment = "八对 M0 summary_min 的描述性类别为：" + "；".join(
        f"{row['pair']} {row['point_change_class']}" for row in m0_rows
    ) + "。"
    if not notable:
        assessment += " 按冻结的 0.02 带宽，M0 summary_min 点估计没有衰减、放大或反向。"
    assessment += " 区间是否包含 0.5 单独记录在 evidence_status，不由点变化类别代替。"
    head = __import__("subprocess").check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    write_reports({
        "head": head,
        "retained": freeze["n_retained"],
        "excluded_n": freeze["n_excluded"],
        "acyclic": freeze["n_acyclic"],
        "excluded": excluded_txt,
        "m0_rows": m0_rows,
        "class_counts": dict(counts),
        "reversed": "；".join(reversed_rows),
        "regression_text": regression_text,
        "assessment": assessment,
    })
    audit = [{"check": layer, "result": "PASS" if ok else "FAIL", "observed": "", "expected": "exact formal table"} for layer, ok in layer_ok.items()]
    audit.append({"check": "membership_sha_unchanged", "result": "PASS" if sha256_file(membership_path) == freeze["membership_sha256"] else "FAIL", "observed": sha256_file(membership_path), "expected": freeze["membership_sha256"]})
    audit.append({"check": "new_fold_join_failures", "result": "PASS" if all(row["qa_status"] == "PASS" for row in qa_rows) else "FAIL", "observed": str(sum(row["qa_status"] != "PASS" for row in qa_rows)), "expected": "0"})
    audit.append({"check": "cap5_analysis", "result": audit_result, "observed": regression_text, "expected": "regression pass and no new fold-join failure"})
    write_csv(OUT / "CAP5_AUDIT_PASSFAIL.csv", audit, ["check", "result", "observed", "expected"])
    manifest = {
        "parent_head": head,
        "membership_sha256": freeze["membership_sha256"],
        "calls": ["q110_stats.single_method", "q110_stats.paired_delta", "module_b", "module_c", "module_d", "module_i"],
        "not_called": ["execute", "run_stage2_modules", "run_primary_metrics_v42", "module_d1", "module_e", "module_h"],
        "regression": layer_ok,
        "audit": audit_result,
        "n_comparison_rows": len(versus),
        "point_change_counts": dict(counts),
    }
    (OUT / "CAP5_EXECUTION_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if sha256_file(membership_path) != freeze["membership_sha256"]:
        raise SystemExit("MEMBERSHIP_REWRITTEN")
    print(audit_result, dict(counts), flush=True)
    return 0 if audit_result == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
