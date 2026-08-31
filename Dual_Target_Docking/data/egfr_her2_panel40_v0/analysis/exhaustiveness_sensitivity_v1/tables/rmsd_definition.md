# RMSD 定义（exhaustiveness sensitivity v1）

## 目的
本任务的 RMSD 只用于回答“提高 exhaustiveness 是否更容易找到近晶构象”。

## 对照对象
- EGFR / 3POZ 共晶参考：`analysis/exhaustiveness_sensitivity_v1/tables/3POZ_cocrystal_03P.pdb`
- HER2 / 3RCD 共晶参考：`analysis/exhaustiveness_sensitivity_v1/tables/3RCD_cocrystal_03P.pdb`
- 共晶配体残基名：`03P`
- 面板对应分子：`EH40_01`（TAK-285）

## 计算定义
- 原子集合：**重原子 only**（排除所有氢；EH40_01 / 03P = 38 重原子）
- 化学键级：两端均用 `ligands_sdf/EH40_01.sdf` 经 RDKit `AssignBondOrdersFromTemplate` 统一
- 对称处理：在模板子结构匹配集合上取 **min CalcRMS**（模板约束的图自同构/对称校正）
  - 该定义与 panel40 as-run 历史数一致（3POZ as-run mode1 ≈ 9.51 Å，mode2 ≈ 1.02 Å）
  - 注意：裸 `GetBestRMS` 在本化学型上可能因键级歧义而系统性偏低，本任务不采用
- 不做蛋白叠合微调；直接使用 prepared 结构坐标系下的共晶坐标作为参考

## 必报指标
对每个 target、每个 exhaustiveness（固定主 seed = `20260727`）分别报告：
1. **Vina mode_01 RMSD**
2. **9 个 mode 中的 min RMSD（best_of_9）**
3. **RTM-best mode RMSD**（仅当 RTMScore 可复现时报告）

## 概念澄清
- `exhaustiveness` = 搜索强度，影响“找不找得到更接近共晶的构象”
- `n_modes` = 输出条数，本任务固定为 9，**不是比较因子**
