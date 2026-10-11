# 本次授权与数据来源

日期：2026-10-11。

工作分支：`cursor/activity-label-sensitivity-20261011`。起点是正式分支 `cursor/v4-2-final-audited-results-20260930` 的 `4ebef6d8a765ba8b5872c699ed4af68dc88ceda4`。

本次任务指令授权只读使用已经恢复的 ChEMBL 37 活性记录，完成 S0–S4。授权只覆盖这一次敏感性分析。没有改写 `analysis_plan_jcim/THRESHOLD_PROVENANCE_GATE.md`。该历史文件里的 `GATE_ACTIVITY_SOURCE_READ=NOT_GRANTED` 仍然是当时的记录。

实际读入：

- `activity_label_source_audit_20261008/chembl37_activity_rows_long.csv`：3555 行，STANDARD_OK，未取整 pChEMBL。最大值和中位数都来自这批行。
- `chembl37_dual_endpoint_crosswalk_805.csv`：对照正式类别和两端最大值。
- `document_identifier_coverage_805.csv`、`historical_adjudication_differences.csv`。
- `activity_semantics_audit_20261009/` 中已完成的药理核查。
- `PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv`：配体身份、正式类别和已有 M0 分数。M0 是 seed 42、mode 1，较高分更有利。没有重新对接，没有填补缺失分数。
- `PRIMARY_DIRECTIONAL_METRICS.csv`：只用于 S0 对照，不覆盖。
- `pair_ligand_mapping.csv`、`parent_alias_registry.csv`、`model_fold_assignments.csv`：只读核对 808/807/805 和 938/934/928。

没有读取历史 `ligand_activity_aggregate_v1.csv`，没有调用 `scripts/analysis/compute_label_sensitivity`。Mann–Whitney 和 bootstrap 按正式定义在本目录脚本中重写，配体顺序使用人群表出现顺序。

S1、S2、S3 是冻结 YAML 里预先写出的 `sensitivity_only: [5.5, 6.5, strict_6.5_5.5]`。S4 中位数是本次附加汇总，不是那份预先规定清单里的实验。中位数不比最大值更正确。

全部情景只对固定的 805 条合格观察重新分类。这不是独立验证，也不能推断重新构建评价集后会得到相同结果。
