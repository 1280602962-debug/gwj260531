# 课题汇报幻灯片大纲（逐页文字稿）

> 共 20 页主内容 + 可选附录。每页含：标题、要点、讲稿提示、推荐图片。  
> 叙事逻辑详见 [BAL_PROJECT_PRESENTATION_GUIDE.md](./BAL_PROJECT_PRESENTATION_GUIDE.md)。

---

## Slide 1 | 封面

**标题**：基于 NLRP3 BAL Glu-switch 变构位点的抑制剂发现

**副标题**：原骨架活性类似物 + 多骨架候选拓展

**要点**
- 汇报人 / 单位 / 日期
- 数据来源：5 篇 WO 专利 + 文献结构生物学

**讲稿提示**：一句话概括——「利用专利大数据和结构生物学，在 BAL 新位点上发现下一代 NLRP3 抑制剂。」

---

## Slide 2 | 研究背景：NLRP3 为什么是重要靶点？

**要点**
- NLRP3 炎性小体：感知危险信号 → caspase-1 → IL-1β/IL-18 + 焦亡
- 疾病：痛风、CAPS、阿尔茨海默、动脉粥样硬化等
- 现有疗法局限：抗体只阻断下游细胞因子
- 小分子直接抑制 NLRP3 可阻断 IL-18 释放和焦亡，且可能穿透 BBB

**推荐图**：自制流程图或 NLRP3 激活示意图（可参考 [Broz & Dixit 2016](https://doi.org/10.1038/nri.2016.58) Figure 1）

**讲稿提示**：「炎症是多种老年疾病的共同通路，NLRP3 是这条通路上最可成药的关键节点之一。」

---

## Slide 3 | NLRP3 抑制剂有多种结合位点

**要点**
- 并非所有「NLRP3 抑制剂」结合同一口袋
- MCC950 口袋 vs Walker A vs NEK7 界面 vs **BAL Glu-switch**
- 本课题聚焦：**BAL 变构位点**（Y258/H260 β2 链）

**推荐图**：`images/binding_sites_schematic.png`

**讲稿提示**：「如果不对位点做区分，计算和实验都会走错方向——这是我们调研的第一个重要结论。」

---

## Slide 4 | 调研方法论

**要点**
1. 文献检索 → 抑制剂分类表（18 个非 MCC950 位点化合物）
2. PDB 结构时间线 → 19 个结构条目
3. 5 篇 WO 专利数据挖掘 → 1087 条化合物
4. AI 方法基准（FoldBench）→ 共折叠策略选择
5. 勘误与数据清洗 → 可复现数据集

**推荐图**：调研流程 mermaid 图（见 GUIDE 文档）

**讲稿提示**：「课题不是从零假设开始，而是建立了可核对、可复现的数据基础。」

---

## Slide 5 | 调研中的关键纠偏

**要点**
| 纠偏项 | 意义 |
|--------|------|
| 五篇专利 = 同一 BAL indazole 系列 | 数据可合并用于 ML |
| BAL 不对接 MCC950 口袋 | 结构建模方向正确 |
| WO1468 SMILES 修复 | 回收 75 条高活性分子 |
| DEL IC50 ≠ 细胞 IC50 | 避免跨体系比较 |

**推荐图**：无（表格为主）

**讲稿提示**：「主动展示纠错过程，说明调研是迭代深化的，不是一次性检索。」

---

## Slide 6 | BAL 系列：从 DEL 到 cryo-EM 的证据链

**要点**
- Hartman 2024：DEL 筛选 → BAL-0028（IC50=25 nM, KD=104–123 nM）
- Wilhelmsen 2025：THP-1 IC50=57.5 nM；不抑制 ATP 酶；与 MCC950 不同位点
- Torp 2025：BAL-1516 cryo-EM → PDB 9IHN/9Q8V（HPUB）
- 机制：3 个 H 键到 β2 链，稳定 inactive decamer

**推荐图**：
- 本地：`images/7PZC_assembly.png`（NLRP3 背景）
- **论文原图**（需从 bioRxiv 截取）：Torp 2025 **Fig. 2**（BAL-1516 结合模式）

**讲稿提示**：「BAL 位点有完整的发现-验证-结构解析链条，不是计算假设出来的。」

---

## Slide 7 | 结构生物学资源与限制

**要点**
- 7PZC（MCC950 decamer）：公开，BAL-1516 精修起点
- 9IHN/9Q8V（BAL-1516）：HPUB，暂不可下载
- **决策**：Phase 1 不依赖共晶；Phase 2 用 7PZC + 约束建模

**推荐图**：`images/7PZC_assembly.png`、`images/7ALV_assembly.png`

**讲稿提示**：「结构未公开不是课题 blocker——配体驱动路线可立即启动。」

---

## Slide 8 | 专利数据资产

**要点**
- 5 篇 WO 专利，1087 行 → 1039 独特 SMILES
- 893 个有活性标签（86%）
- 166 个 Murcko 骨架
- 活性分布：高 479 / 中 177 / 低 194

**推荐图**：`images/patent_data_summary_chart.png`

**讲稿提示**：「这是本课题相对一般虚拟筛选项目的核心优势——有近 900 个标注分子。」

---

## Slide 9 | 化学空间与药效团

**要点**
- 主骨架：indazole-酰胺-芳基（248+170 条）
- 共同药效团：indazole HBD/HBA + 酰胺 linker + 疏水芳环
- 双轨基础：骨架内（Top 3 Murcko）+ 骨架跃迁（166 骨架聚类）

**推荐图**：`images/murcko_top1.png`、`images/murcko_top2.png`

**讲稿提示**：「166 个骨架看起来多，但 89% 共享 indazole 核心——是同一药效团的取代基变体。」

---

## Slide 10 | 可行性评估

**要点**
| 维度 | ★ | 说明 |
|------|---|------|
| 科学 | ★★★★☆ | 证据链完整 |
| 数据 | ★★★★☆ | 893 标注分子 |
| 计算 | ★★★★☆ | Phase 1 立即可做 |
| 实验 | ★★★☆☆ | 人源细胞必须 |
| IP | ★★★☆☆ | 专利密集 |

**结论**：可行，三轨并行

**讲稿提示**：「整体四星，实验和 IP 是主要风险点，但有明确对策。」

---

## Slide 11 | 决策 1：为何选 BAL 位点？

**要点**
- MCC950 口袋：肝毒性失败、专利拥挤
- BAL Glu-switch：新 modality、与 MCC950 可加合、CNS 潜力
- 专利数据丰富，可立即 ML

**推荐图**：Wilhelmsen 2025 **Fig. 4D**（nanoDSF 竞争实验，证明不同位点）

**讲稿提示**：「选 BAL 不是因为它『更新』，而是因为有数据、有结构、有差异化机制。」

---

## Slide 12 | 决策 2：配体驱动为主

**要点**
- 不选纯对接：9IHN 未公开 + 变构位点假阳性高
- 不选纯 ML：缺结构过滤
- **选配体 ML + 约束结构**：Phase 1 立即启动，Phase 2 结构验证

**推荐图**：`images/four_phase_workflow.png`

**讲稿提示**：「这是风险最低的推进顺序——先用数据产生假设，再用结构验证。」

---

## Slide 13 | 决策 3：AI 共折叠策略

**要点**
- FoldBench：AF3 64.9% > Boltz 55% > Chai 51%
- 变构位点是短板（Allosteric Paradox 2026）
- **必须加 Y258/H260 pocket 约束**；SiteAF3 有约束时 ~72%
- 多方法交叉验证，不一致的 pose 不采纳

**推荐图**：FoldBench 2025 **Figure 2**（方法比较柱状图）

**讲稿提示**：「AI 共折叠不是不能用，但不能盲用——约束是必要条件。」

---

## Slide 14 | 决策 4：双轨候选生成

**要点**
| 轨道 | 方法 | 规模 | 产出 |
|------|------|------|------|
| A 骨架内 | R-group 枚举 + ML | 500–2000/系列 | 30–50 个 |
| B 骨架跃迁 | 药效团 + REAL 库 | 3–5 系列 | 20–40 个 |

**推荐图**：BAL 药效团示意图（Torp 2025 Fig. 2c 氢键模式）

**讲稿提示**：「A 轨保成功率，B 轨保创新性——并行推进。」

---

## Slide 15 | 决策 5：实验验证体系

**要点**
- 金标准：THP-1 LPS(100 ng/mL) → 化合物 30 min → nigericin 10 µM 1 h → IL-1β ELISA
- BAL 不抑制 ATP 酶（机制分型必做）
- 种属：BAL-0028 在鼠细胞 IC50 >6 µM，**必须人源**
- 可选：ASC speck 流式、SPR、nanoDSF

**推荐图**：`images/assay_pyramid.png`

**讲稿提示**：「实验设计直接沿用 Wilhelmsen 和 Torp 的验证体系，确保可比性。」

---

## Slide 16 | 四阶段实施方案

**要点**
- Phase 1：RDKit 清洗 → Murcko 聚类 → XGBoost → Top 50
- Phase 2：7PZC 单体 → AF3 约束 → GNINA 重对接验证
- Phase 3：双轨枚举 → ML+对接双阈值 → ADMET 过滤
- Phase 4：THP-1 测活 → 迭代

**推荐图**：`images/four_phase_workflow.png`

---

## Slide 17 | 里程碑

**要点**
| 里程碑 | 标准 | 时间 |
|--------|------|------|
| M1 | 数据+SAR 报告 | Phase 1 末 |
| M2 | ML AUC>0.75 | Phase 1 末 |
| M3 | 结构模型验证 | Phase 2 末 |
| M4 | 多骨架库 20–40 | Phase 3 末 |
| M5 | ≥3 hit IC50<100 nM | Phase 4 末 |

---

## Slide 18 | 风险与对策

**要点**（表格，见 GUIDE Slide 18）

**讲稿提示**：「风险都已识别并有对策，不是未知未知。」

---

## Slide 19 | 创新点总结

**要点**
1. 新变构位点（BAL Glu-switch）
2. 专利大数据驱动（1039 分子）
3. 配体+结构双轨计算
4. 骨架内+骨架跃迁双目标
5. 实验闭环（THP-1 金标准）

---

## Slide 20 | 下一步 & Q&A

**要点**
- 本周：启动 Phase 1 ML pipeline
- 并行：BAL-1516 SMILES + AF3 输入准备
- 采购：THP-1 细胞、LPS、nigericin、IL-1β ELISA kit

**推荐图**：GitHub 仓库 QR 码或链接

---

## 附录 A（可选）| 抑制剂完整分类表

- 引用 `non_mcc950_site_inhibitors.csv` 前 10 行

## 附录 B（可选）| PDB 结构时间线

- 引用 `nlrp3_inhibitor_structures_timeline.csv`

## 附录 C（可选）| 参考文献速查

- 见 [REFERENCES.md](../REFERENCES.md)

---

*配合 [IMAGE_SOURCES.md](./IMAGE_SOURCES.md) 使用，确保图片引用规范。*
