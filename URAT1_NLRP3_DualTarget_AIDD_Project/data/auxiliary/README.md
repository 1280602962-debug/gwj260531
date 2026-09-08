# SLC22 辅助库（OAT / OCT）

数据保留供对照，**当前投稿漏斗不使用 OAT 迁移作为主结果**。配置：`config/targets.yaml` → `auxiliary_targets`。

URAT1（SLC22A12）是 OAT 亚家族尿酸交换体，不是 OCT。

| 库 | 基因 | 清洗后 SMILES | 角色 |
|----|------|---------------|------|
| OAT1 | SLC22A6 | 63 | 阴离子邻域对照（非主贡献） |
| OAT3 | SLC22A8 | 41 | 同上 |
| OAT 合并 | — | 73 | `oat_combined_transfer.csv` |
| OCT1 | SLC22A1 | 108 | 阳离子脱靶讨论 |
| OCT2 | SLC22A2 | 105 | 同上 |

清洗规则与 URAT1 主任务相同（见 `docs/DATA_FACT_CHECK.md`）。错误 ChEMBL ID 黑名单亦见该文档。
