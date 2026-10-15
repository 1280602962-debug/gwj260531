# NLRP3 数据审计

| 项目 | 已核实内容 | 文件 |
| --- | --- | --- |
| 来源 | ChEMBL 导出，版本号无法恢复 | `data/raw/NLRP3_CHEMBL_4807.csv`，1334 行 |
| 靶点 | CHEMBL1741208 | 同上 |
| 活性类型 | 存档只有 IC50 | 原始表 |
| 清洗 | 测定描述含 IL-1、测定类型 B、等号关系、pChEMBL 在 4 到 10、测定至少 5 个化合物 | `curate_nlrp3_records` |
| 记录 | 609 条，513 个分子，39 个测定 | `data/processed/nlrp3_records.csv` |
| 标签 | 记录 pActivity ≥ 6.0；分子取最大 | 同上 |
| 分子活性比 | 310 活性 / 203 非活性，约 0.604 | `data/processed/data_summary.json` |
| one-hot | 训练用频次最高的 25 个测定；部署时对频次最高的 5 个测定上下文取最大概率 | v1 训练脚本与筛选摘要 |

多测定冲突有两个不同数字：`data_summary.json` 的极差大于 1 个对数单位是 37 个分子；v2 审计的标签不一致是 32 个。它们不是同一个定义，不要合并。

0.8934 对应 v1 的五折平均值，记录数 609、分子数 513、使用测定 25。见 `03_MACHINE_LEARNING/MODEL_METHODS_zh.md`。
