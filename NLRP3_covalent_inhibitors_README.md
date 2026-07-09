# NLRP3 共价抑制剂主表

本目录数据文件：`NLRP3_covalent_inhibitors_master_table.csv`

## 纳入标准

实验明确共价结合位点，至少满足以下两项：

1. 位点突变挽救（Cys/Lys → Ala 后结合或抑制消失）
2. LC-MS/MS 肽段直接鉴定修饰位点
3. 不可逆性对照（洗脱实验或还原弹头失活）

## 物种编号说明（重要）

人源 NLRP3（UniProt Q96P20，1036 aa）与鼠源 NLRP3（1033 aa）在 NACHT 域存在 **1 个氨基酸偏移**，导致部分 Cys 编号不同。**这不是录入错误，而是原文使用的物种编号体系不同。**

| 人源残基 | 鼠源残基 | 代表分子 | 说明 |
|----------|----------|----------|------|
| **Cys279** | **Cys280** | Oridonin / 化合物49（人源） vs DCL / E1（鼠源） | 同一 NEK7 界面同源半胱氨酸；机制均为 Michael 加成后阻断 NEK7–NLRP3 |
| Cys409 | Cys409 | RRx-001 / 149-01 | 人源实验体系；鼠源同编号（此位点两侧保守） |
| Cys598 | Cys598 | Costunolide | 人源 rhNLRP3-NACHT 实验；ATP 区邻域，与 279/280/409 位置不同 |
| Cys548 | Cys548 | Itaconate / 4-OI | 原文基于 **murine** NLRP3 HD2；人源引用需核对 |
| Lys377 | Lys377 | Manoalide | 人源 HD1；**非 Cys**，共价机制差异化 |

**统一写法建议**：表格中保留原文 `covalent_site` + `site_species`；跨物种比较时使用 `human_ortholog_site` 列。

> 英文备注模板：`Mouse Cys280 corresponds to human Cys279 (species numbering difference; same orthologous NEK7-interface cysteine).`

## 字段说明

| 字段 | 说明 |
|------|------|
| `compound_type` | 天然产物 / 天然产物类似物 / 合成药物 / 合成类似物 / 内源代谢物 |
| `site_species` | 原文位点编号所依据的物种（Human / Mouse / Murine / Human/mouse） |
| `human_ortholog_site` | 人源同源残基编号（便于跨文献比较） |
| `mechanism_route` | `A_NEK7_block` / `B_ATPase_oligomerization` / `C_crosslinking` |
| `pubmed_status` | PubMed 收录状态（引用前建议复核） |
| `af3_calibration_priority` | AF3 共价回顾校准推荐优先级 |

## E1 文献状态（2026-07 复核）

| 项目 | 状态 |
|------|------|
| DOI | `10.1021/acs.jmedchem.5c01663` — **可正常解析** |
| 期刊 | *J. Med. Chem.* **2025**, 68(20), 21534–21559 |
| PubMed | **已收录**（检索 DOI 或标题可定位） |
| 位点 | 原文 **鼠源 Cys280** = 人源 **Cys279** |

## 机制路线与验证读数

| 路线 | 代表分子 | 必要验证实验 |
|------|----------|-------------|
| NEK7 阻断（A） | Oridonin, DCL, RRx-001, Manoalide, 4-OI | NEK7 Co-IP / Co-IP 挽救 |
| ATPase/寡聚（B） | Costunolide | **NLRP3 ATPase assay** + 寡聚（SDD-AGE） |
| 蛋白交联（C） | VLX1570 | HMW 交联；阴性对照 |

## 统计（2026-07）

- 总条目：10
- 天然产物：4；天然产物类似物：2；合成：3；内源代谢物：1
- 人源明确单 Cys 位点：C279, C409, C598
- 鼠源报道、人源同源 Cys279：DCL, E1（鼠 Cys280）

## 关联文档

- `NLRP3_共价抑制剂双路线筛选技术报告.md`
