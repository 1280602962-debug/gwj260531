# NLRP3 抑制剂化合物信息表

经文献核对整理的 NLRP3 直接抑制剂数据集，用于区分 **MCC950/CRID3 口袋抑制剂** 与 **其他结合位点/机制** 的化合物。

## 文件说明

| 文件 | 行数 | 用途 |
|------|------|------|
| [non_mcc950_site_inhibitors.csv](./non_mcc950_site_inhibitors.csv) | 18 | 非 MCC950 口袋或不同机制的 tool/clinical 化合物 |
| [mcc950_pocket_inhibitors_reference.csv](./mcc950_pocket_inhibitors_reference.csv) | 14 | MCC950 口袋结合剂及临床管线对照 |
| [nlrp3_inhibitor_structures_timeline.csv](./nlrp3_inhibitor_structures_timeline.csv) | 19 | PDB 结构条目与尚无结构的化合物标注 |
| [CORRECTIONS.md](./CORRECTIONS.md) | — | 相对先前版本的勘误说明 |

## 字段定义

| 字段 | 说明 |
|------|------|
| `compound_id` | 主要化合物名称 |
| `aliases` | 别名、代号 |
| `scaffold` | 化学骨架类型 |
| `binding_site_category` | 结合位点分类（见 CORRECTIONS.md 图示） |
| `binding_site_detail` | 残基/结构域细节 |
| `mechanism` | 作用机制简述 |
| `covalent_reversible` | 共价/可逆性 |
| `assay_cells` | 主要细胞体系 |
| `stimuli` | 激活刺激 |
| `readout` | 读出指标 |
| `IC50_KD` | 活性数据（注明体系） |
| `ATPase_inhibition` | 是否抑制 NLRP3 ATP 酶 |
| `PDB` | 共晶 PDB ID（如有） |
| `structure_year` | 结构发布年 |
| `key_residues` | 关键氨基酸 |
| `species_selectivity` | 种属选择性 |
| `clinical_status` | 临床阶段 |
| `DOI` | 主要参考文献 DOI |
| `notes` | 勘误、争议、注意事项 |

## 使用建议

1. **虚拟筛选 / 对接**：按 `binding_site_category` 分口袋建模，不宜将五类位点混合训练单一 QSAR 模型。
2. **活性比较**：`IC50_KD` 列中不同细胞系、刺激条件、预孵育时间的结果**不可横向比较**。
3. **结构药物设计**：优先参考 `nlrp3_inhibitor_structures_timeline.csv` 中有 PDB 的 MCC950 口袋条目（7PZC, 7ALV, 8RI2 等）。
4. **双口袋策略**：CRID3 口袋（MCC950 类）+ BAL-0028 变构位点可作为互补筛选假设（见 Wilhelmsen 2025 JEM）。

## 关键参考文献速查

| 主题 | DOI |
|------|-----|
| MCC950 结构 | [10.1038/s41586-022-04467-w](https://doi.org/10.1038/s41586-022-04467-w) |
| NP3-146 结构 | [10.1016/j.jmb.2021.167309](https://doi.org/10.1016/j.jmb.2021.167309) |
| CY-09 | [10.1084/jem.20171419](https://doi.org/10.1084/jem.20171419) |
| Oridonin | [10.1038/s41467-018-04947-6](https://doi.org/10.1038/s41467-018-04947-6) |
| Tranilast | [10.15252/emmm.201708689](https://doi.org/10.15252/emmm.201708689) |
| BAL-0028 机制 | [10.1084/jem.20242403](https://doi.org/10.1084/jem.20242403) |
| BAL-1516 结构 | [10.1101/2025.07.01.662566](https://doi.org/10.1101/2025.07.01.662566) (PDB 9IHN/9Q8V) |
| 4-OI / Itaconate | [10.1016/j.cmet.2020.07.016](https://doi.org/10.1016/j.cmet.2020.07.016) |
| RRx-001 | [10.1038/s41423-021-00683-y](https://doi.org/10.1038/s41423-021-00683-y) |
| VLX1570 | [10.1021/acschembio.3c00330](https://doi.org/10.1021/acschembio.3c00330) |

## 许可与引用

数据整理自公开文献，仅供研究参考。引用时请核对原始论文中的活性数据与实验条件。
