# NLRP3 抑制剂化合物信息表

经文献核对整理的 NLRP3 直接抑制剂数据集，并扩展 **BAL Glu-switch 位点抑制剂发现课题** 相关数据与方案。

## 文件说明

| 文件 | 行数/规模 | 用途 |
|------|-----------|------|
| [presentation/BAL_PROJECT_PRESENTATION_GUIDE.md](./presentation/BAL_PROJECT_PRESENTATION_GUIDE.md) | — | **课题汇报叙事指南（调研→决策→方案）** |
| [presentation/SLIDES_OUTLINE.md](./presentation/SLIDES_OUTLINE.md) | 20 页 | **逐页幻灯片文字稿** |
| [presentation/IMAGE_SOURCES.md](./presentation/IMAGE_SOURCES.md) | — | **PPT 图片清单与引用规范** |
| [presentation/images/](./presentation/images/) | 11 张 PNG | **可直接插入 PPT 的本地图片** |
| [BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md](./BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md) | — | **课题可行性分析、四阶段实施方案、资源清单** |
| [BAL_PUBLICATION_TIMELINE.md](./BAL_PUBLICATION_TIMELINE.md) | — | **BAL 新位点分子发表时间线** |
| [bal_publication_timeline.csv](./bal_publication_timeline.csv) | 19 事件 | 时间线机器可读 CSV |
| [REFERENCES.md](./REFERENCES.md) | — | **完整参考文献链接（按主题分类）** |
| [patent_bal_compounds_merged.csv](./patent_bal_compounds_merged.csv) | 1087 行 / 1039 独特 SMILES | 五篇 BAL 专利合并清洗数据 |
| [patent_data_summary.json](./patent_data_summary.json) | — | 专利数据统计摘要 |
| [non_mcc950_site_inhibitors.csv](./non_mcc950_site_inhibitors.csv) | 18 | 非 MCC950 口袋抑制剂（含 BAL 系列） |
| [mcc950_pocket_inhibitors_reference.csv](./mcc950_pocket_inhibitors_reference.csv) | 14 | MCC950 口袋及临床对照 |
| [nlrp3_inhibitor_structures_timeline.csv](./nlrp3_inhibitor_structures_timeline.csv) | 19 | PDB 结构时间线 |
| [CORRECTIONS.md](./CORRECTIONS.md) | — | 勘误与数据修正记录 |

## 课题概要（BAL 位点新抑制剂发现）

- **目标**：在原 indazole 骨架上发现活性类似物 + 拓展多骨架候选
- **可行性**：★★★★☆（配体驱动可立即启动；结构驱动需约束验证）
- **核心数据**：1039 独特分子，893 有活性标签，166 Murcko 骨架
- **关键约束**：BAL 位点为变构沟槽，AI 共折叠需指定 Y258/H260；9IHN 仍 HPUB

详见 [BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md](./BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md)。

## 字段定义（patent_bal_compounds_merged.csv）

| 字段 | 说明 |
|------|------|
| `patent_id` | 专利号 |
| `compound_number` | 化合物编号 |
| `smiles` | 规范化 SMILES |
| `activity_raw` | 原始活性标签 |
| `activity_score` | 统一评分 0–3 |
| `murcko_scaffold` | Murcko 骨架 |
| `mw`, `logp`, `tpsa`, `hbd`, `hba` | 理化描述符 |

## 课题汇报

1. **汇报准备**：先读 `presentation/BAL_PROJECT_PRESENTATION_GUIDE.md`（叙事逻辑），再按 `SLIDES_OUTLINE.md` 制作 PPT。
2. **图片资源**：`presentation/images/` 含结构图、流程图、数据图；论文原图索引见 `IMAGE_SOURCES.md`。

## 使用建议

1. **课题实施**：先读 `BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md`，按 Phase 1→4 推进。
2. **活性建模**：按 `patent_id` 或 `murcko_scaffold` 分层，不宜五篇简单合并回归。
3. **结构建模**：用 7PZC 单体 + AF3/Boltz（约束 Y258/H260），不用 MCC950 口袋对接 BAL 化合物。
4. **验证**：专利 +++ 分子重对接 + THP-1 实验闭环。

## 关键参考文献

完整列表见 [REFERENCES.md](./REFERENCES.md)。核心链接：

| 主题 | 链接 |
|------|------|
| BAL-0028 发现 | https://doi.org/10.1016/j.bmcl.2024.129675 |
| BAL-0028 机制 | https://doi.org/10.1084/jem.20242403 |
| BAL-1516 结构 | https://doi.org/10.1101/2025.07.01.662566 |
| 7PZC 结构 | https://doi.org/10.1038/s41586-022-04467-w |
| FoldBench 基准 | https://doi.org/10.1038/s41467-025-67127-3 |

## 许可与引用

数据整理自公开专利与文献，仅供研究参考。引用时请核对原始来源。
