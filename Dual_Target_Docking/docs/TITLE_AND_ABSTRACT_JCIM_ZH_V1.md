# 题名与摘要（JCIM Articles 中文工作稿）

## 题名

**基于对接的双靶识别：四靶对的评价设定审计**

## 摘要

两端有利的对接分数能否作为双靶活性的证据，取决于评价时采用的负类。我们用 ChEMBL 衍生的 dual、A-only、B-only 与 neither 操作性状态，对四个靶对进行评价设定审计。两条口袋匹配方向 AUROC 分别将 dual 与对应单靶选择性配体比较；二者最小值 `summary_min` 仅作保守的描述性摘要。在 EGFR/HER2 上，Vina 的 Dual versus neither AUROC 为 0.756，而方向性 `summary_min` 为 0.430；独立 GNINA 姿态生成得到的对应数值为 0.783 与 0.220。其他靶对没有出现同样的设定差距，且 PIK3CA/mTOR 的 neither 对照只有 4 个分子。在支架分组模型中，把 docking 加入 ECFP4 后交叉验证 AUROC 的最大绝对变化为 0.020。替代受体使 PIK3CA/mTOR 的 `summary_min` 从 0.692 降至 0.486/0.505，但使 PIK3CA/PIK3CB 升高，说明结果依赖受体实现而不是具有结构稳健性。四个主 `summary_min` 的 95% CI 均包含 0.5。同一套冻结 EGFR/HER2 分数上，以 Dual 中位 `vina_worst` 做 AND 式双口袋过滤时，通过者多数仍是实验选择性配体（precision 0.298；硬负比例 0.702）。在四对完整 ChEMBL 图上，仅用配体 ECFP4 仍比 Dual versus 选择性更容易分开 Dual versus neither，说明设定问题不是 n ≈ 28 对接面板的抽样伪影。按文献 `document_id` 阻断后，EGFR/HER2 的弱方向臂仍为 0.430；预先冻结的 2018 文献年份分割没有两个可评估靶对，因此不作为外部验证。因此，这一受数据供给约束的案例面板支持把选择性硬负与混淆感知对照作为评价要求，但不能建立靶标通用的对接性能或生物学识别结论。

## 关键词

双靶对接；基准设定；选择性硬负样本；化学混淆；受体实现；虚拟筛选
