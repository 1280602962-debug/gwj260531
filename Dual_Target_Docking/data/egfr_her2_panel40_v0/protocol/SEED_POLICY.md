# SEED POLICY（egfr_her2_panel40_v0）

## 结论
- `v0_pilot = as_run_per_job_random`
- `v0_1_forward = fixed_global`
- 正向默认固定 seed：`20260727`
- 历史结果仅做归档，不以新 seed 重解释，不覆盖原始 tables / poses。

## as-run 归档说明
本次对 `logs/vina_seeds.json` 与原始 Vina 运行日志做了交叉核对。

- 期望作业数：80（`3POZ/3RCD × EH40_01..40`）
- `vina_seeds.json` 直接记录：41 条
- 从原始运行日志恢复到的 seed：80/80
- `seeds_as_run.csv` 最终覆盖率：80/80
- 缺失作业：无

说明：
- `vina_seeds.json` 只是一份**部分捕获**，不能单独作为完整 as-run seed 档案。
- 完整归档以 `protocol/seeds_as_run.csv` 为准；其 seed 值来自 Vina 日志中的 `Performing docking (random seed: ...)`。
- Vina 历史 seed 允许为负整数；这是运行时输出格式的一部分，不做人为改写。

## 前向协议（v0.1）
后续 exhaustiveness 敏感性与任何正向补充实验，统一采用：

- seed policy: `fixed_global`
- `seed_fixed_global = 20260727`
- `n_modes = 9`
- `energy_range = 3`
- `exhaustiveness` 先不预设，待 `analysis/exhaustiveness_sensitivity_v1/SENSITIVITY_VERDICT.md` 给出建议后回写。

## 科学边界
1. as-run（多 seed）只回答“旧结果能否大致复现”，**不能**用于选择 exhaustiveness。
2. exhaustiveness 对照必须在**同一配体输入 + 同一盒子 + 同一主 seed**下进行。
3. 历史 v0 pilot 的原始排名、分数、poses 不因固定 seed 政策而被追认或重算。
