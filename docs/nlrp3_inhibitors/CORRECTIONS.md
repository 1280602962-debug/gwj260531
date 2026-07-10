# NLRP3 抑制剂信息勘误说明

本文档列出相对先前对话整理版本的**具体勘误**，并说明数据来源。

## 文献归属错误

| 化合物 | 先前错误 | 更正 |
|--------|----------|------|
| **CY-09** | 误写为 *Nature Medicine* | 正确期刊：*Journal of Experimental Medicine*，2017；DOI: [10.1084/jem.20171419](https://doi.org/10.1084/jem.20171419) |
| **OLT1177 / Dapansutrile** | 误写为 *Nature Medicine* (Marchetti 2018) | 正确期刊：*PNAS*，2018；DOI: [10.1073/pnas.1716095115](https://doi.org/10.1073/pnas.1716095115) |
| **BAL-0028** | 机制文献未区分 | 发现：Hartman et al. *Bioorg Med Chem Lett* 2024 ([10.1016/j.bmcl.2024.129675](https://doi.org/10.1016/j.bmcl.2024.129675))；机制：Wilhelmsen et al. *JEM* 2025 ([10.1084/jem.20242403](https://doi.org/10.1084/jem.20242403)) |

## 结合位点分类错误

| 化合物 | 先前错误 | 更正 |
|--------|----------|------|
| **NIC-12 / NP3-562 / compound 32 (8WSM)** | 列为「非 MCC950 位点」 | **与 MCC950 同一口袋**（CRID3/sulfonylurea 口袋），仅化学骨架与 linker 锚定不同；见 EMBO Mol Med 2024 与 PDB 8RI2/8WSM |
| **Oridonin** | 部分综述写 Cys273 | 人源 NLRP3 共价位点为 **Cys279**；鼠源功能验证使用 **C275A** 突变（He et al. 2018） |
| **4-OI** | 未说明种属差异 | 主要生化证据为**鼠源 Cys548**；人源是否有直接同源修饰位点仍有争议 |

## 遗漏的重要化合物（已补入表格）

| 化合物 | 机制要点 | 关键文献 |
|--------|----------|----------|
| **4-octyl itaconate (4-OI)** | 代谢物 itaconate 衍生物；C548 修饰阻断 NLRP3-NEK7 | Hooftman et al. *Cell Metab* 2020 |
| **INF58** | 丙烯酰胺共价抑制剂；ATP 酶 IC50 74 µM；推测 Cys419 | Cocco et al. *J Med Chem* 2016 |
| **VLX1570** | PYD 域 Cys108/Cys130 交联；全新机制 | Stanton et al. *ACS Chem Biol* 2024 |
| **MNS** | NACHT+LRR 共价修饰；ATP 酶抑制 | He et al. *J Biol Chem* 2013 |
| **BOT-4-one** | NLRP3 烷基化 + 泛素化增强 | *Cell Death & Disease* 2017 |

## 活性数据注意事项

| 化合物 | 注意 |
|--------|------|
| **BAL-0028** | DEL 筛选 IC50=25 nM 与 THP-1+nigericin IC50=57.5 nM 来自**不同实验体系**，不可直接比较 |
| **CY-09** | Kd=0.5 µM (MST) 与细胞 IL-1β IC50~6 µM 并存，反映结合 vs 细胞功能差异 |
| **OLT1177** | Teske et al. 2024 在多种细胞实验中质疑其直接抑制 NLRP3 的效力；临床安全性数据仍来自 Marchetti 2018 Phase 1 |
| **INF58 Cys419** | 来自同源建模假设，**尚未有共晶或突变验证实验** |

## PDB 结构时间线补充

先前版本遗漏或需更新的条目：

- **8ETR** (2022): GDC-2394 + 人源 NACHT
- **8WSM** (2023): compound 32 / NIC-12 系列
- **8RI2** (2024): NP3-562 三环骨架
- **8ZEM** (2024): SN3-1 / NP3-1
- **9GU4, 9SFG** (2025): NP3-253, NP3-742 脑渗透系列

**尚无共晶结构**：Oridonin, Tranilast, CY-09, RRx-001, 4-OI, VLX1570

**2025-07 补充（此前遗漏）**：
- **BAL-1516**：BAL-0028 的 CNS 优化衍生物；首个解析结合位点的 BAL 系列共晶
- **PDB 9IHN**（单体 C1, 3.06 Å）和 **9Q8V**（十聚体 D5, 3.06 Å）
- 来源：Torp et al. bioRxiv 2025 ([10.1101/2025.07.01.662566](https://doi.org/10.1101/2025.07.01.662566))
- 结合位点：NBD 表面沟槽（FISNA+WHD+β-sheet Glu-switch），与 MCC950 **可加合**（非重叠）

## 结合位点分类（修正版）

```
┌─────────────────────────────────────────────────────────────────┐
│ CRID3/MCC950 口袋（Walker B 区域，稳定 ADP 结合封闭构象）        │
│   MCC950, NP3-146, GDC-2394, SN3-1, NP3-562, NIC-12, Inzomelid  │
├─────────────────────────────────────────────────────────────────┤
│ Walker A / ATP 袋（阻断 ATP 结合或 ATP 酶）                      │
│   CY-09, MNS, INF58, INF39, OLT1177, Bay11-7082, Parthenolide   │
├─────────────────────────────────────────────────────────────────┤
│ NEK7 界面（共价修饰阻断 NLRP3-NEK7）                             │
│   Oridonin (C279), RRx-001 (C409), 4-OI (mouse C548)            │
├─────────────────────────────────────────────────────────────────┤
│ 寡聚化界面（ATP 酶非依赖）                                       │
│   Tranilast                                                     │
├─────────────────────────────────────────────────────────────────┤
│ 新变构位点（不抑制 ATP 酶，与 MCC950 不重叠）                     │
│   BAL-0028 / BAL-1516 / BAL-0598（FISNA Glu-switch 沟槽；9IHN/9Q8V）│
├─────────────────────────────────────────────────────────────────┤
│ 共价交联（PYD 域分子间交联）                                     │
│   VLX1570                                                       │
└─────────────────────────────────────────────────────────────────┘
```

## 数据文件

| 文件 | 内容 |
|------|------|
| `non_mcc950_site_inhibitors.csv` | 非 MCC950 口袋 / 不同机制抑制剂（主攻表） |
| `mcc950_pocket_inhibitors_reference.csv` | MCC950 口袋及临床对照化合物 |
| `nlrp3_inhibitor_structures_timeline.csv` | PDB 结构时间线 |

*最后更新：2026-07-10*
