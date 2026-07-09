# NLRP3 共价抑制剂主表

本目录数据文件：`NLRP3_covalent_inhibitors_master_table.csv`

## 纳入标准

实验明确共价结合位点，至少满足以下两项：

1. 位点突变挽救（Cys/Lys → Ala 后结合或抑制消失）
2. LC-MS/MS 肽段直接鉴定修饰位点
3. 不可逆性对照（洗脱实验或还原弹头失活）

## 字段说明

| 字段 | 说明 |
|------|------|
| `compound_type` | 天然产物 / 天然产物类似物 / 合成药物 / 合成类似物 / 内源代谢物 |
| `mechanism_route` | `A_NEK7_block`（NEK7阻断）/ `B_ATPase_oligomerization`（ATP酶活/寡聚）/ `C_crosslinking`（多Cys交联） |
| `af3_calibration_priority` | AF3 共价回顾校准推荐优先级 |

## 统计（2026-07）

- 总条目：10
- 天然产物：4；天然产物类似物：2；合成：3；内源代谢物：1
- 明确单 Cys 位点：5（279/280/409/548/598）

## 关联文档

- `NLRP3_共价抑制剂双路线筛选技术报告.md`
