# NLRP3 共价抑制剂双路线并行筛选技术报告

**报告类型：** 计算—实验一体化可行性及实施方案  
**靶点：** 人源 NLRP3 炎性小体  
**策略：** 路线 A（PPI 界面阻断型）与路线 B（ATP 酶活/寡聚抑制型）并行推进  
**编制日期：** 2026 年 6 月

---

## 摘要

NLRP3 炎性小体异常激活与痛风、炎症性肠病、神经退行性疾病等多种疾病密切相关，靶向 NLRP3 的共价小分子（靶向共价抑制剂，TCIs）因其持久占据靶标、可锁定活性构象等优势，成为近年药物设计的重要方向。NLRP3 全蛋白含 45 个半胱氨酸残基，共价修饰位点多样、机制不一，单一筛选策略难以覆盖不同药理路径。

本报告提出 **两条机制平行的 NLRP3 共价化合物筛选路线**：

- **路线 A（PPI 界面型）：** 靶向 NACHT 结构域中 NLRP3–NEK7 蛋白相互作用界面的 Cys279、Cys280、Cys409，通过共价修饰产生空间位阻，阻断 NEK7 招募；
- **路线 B（ATP/寡聚型）：** 靶向 NACHT 核苷酸结合区邻域的 Cys598 等位点，通过共价修饰干扰 ATP 水解与 NLRP3 自寡聚。

计算端以 **AlphaFold 3（AF3）分半胱氨酸位点共价共折叠**为核心，整合反应性预筛、置信度排序与分子动力学（MD）验证；实验端建立 **GSH 稳定性—DARTS—质谱肽段定位—定点突变—机制分叉验证—细胞/动物表型** 的共用确证平台。两条路线在早期共用资源，在后期按机制分叉验证，兼顾效率与可区分性。

---

## 一、课题背景与科学价值

### 1.1 NLRP3 炎性小体与疾病

NLRP3 炎性小体是由传感器蛋白 NLRP3、接头蛋白 ASC 及效应蛋白 pro-caspase-1 组装而成的细胞质多蛋白复合物。其异常激活可促进 caspase-1 依赖的 IL-1β、IL-18 成熟及 GSDMD 介导的细胞焦亡，参与痛风、动脉粥样硬化、阿尔茨海默病、脓毒症等多种疾病进程 [1–3]。

### 1.2 可逆抑制剂的局限与共价策略的优势

以 MCC950 为代表的可逆 NLRP3 抑制剂虽为领域金标准，但存在肝毒性、构象选择性受限及耐药等临床转化障碍 [4,5]。靶向共价抑制剂（TCIs）具有以下优势 [6,7]：

- **持久的靶标占据时间（residence time）**，降低对游离药物浓度的依赖；
- **可锁定特定构象或中间体**，在部分激活态靶标上仍可能有效；
- 天然产物中广泛存在的 Michael 受体、内酯等亲电弹头，为共价先导发现提供丰富化学空间 [8,9]。

### 1.3 双路线并行的科学依据

NLRP3 共价抑制并非单一机制，文献已确证至少三类路径 [8–15]：

| 机制类型 | 代表化合物 | 主要位点 | 效应 |
|----------|------------|----------|------|
| **A. PPI 阻断** | Oridonin、DCL、RRx-001、INF39 | Cys279/280/409 | 阻断 NLRP3–NEK7 结合 |
| **B. ATP 酶活/寡聚** | Costunolide | Cys598 | 抑制 ATP 水解与寡聚 |
| **C. 分子间交联** | VLX1570 | 多 Cys | 蛋白交联聚集（脱靶风险高） |

本方案 **聚焦 A、B 两条药物化可行性更高的路线并行筛选**，暂不将交联型作为主线。

---

## 二、NLRP3 半胱氨酸位点图谱与路线定义

### 2.1 人源 NLRP3 半胱氨酸概况

人源 NLRP3 含 1036 个氨基酸、**45 个半胱氨酸**，分布于 PYD、NACHT 及 LRR 结构域。并非所有 Cys 均适合药物化共价修饰；筛选应基于 **实验确证、机制一致性、结构可及性、弹头匹配性** 四层证据分级。

### 2.2 路线 A：PPI 界面阻断型

**靶点位点（Tier 1）：**

| 位点 | 结构域 | 代表化合物 | 弹头类型 | 核心效应 | 参考文献 |
|------|--------|------------|----------|----------|----------|
| Cys279 | NACHT | Oridonin（冬凌草甲素） | 环戊烯酮 Michael 受体 | 阻断 NLRP3–NEK7 | He et al., *Nat. Commun.* 2018 [10] |
| Cys280 | NACHT | DCL（脱氢木香内酯） | 环外亚甲基-γ-内酯 | 阻断 NLRP3–NEK7 | Lv et al., *MedComm* 2025 [11] |
| Cys409 | NACHT | RRx-001、INF39、149-01 | 卤乙酰胺/丙烯酸酯等 | 阻断 NLRP3–NEK7 | Chen et al., 2021; Cocco et al., 2018 [12,13] |

**作用原理：** PPI 界面并非完全平坦，NEK7 结合沟槽中存在浅疏水空腔；小分子以非共价部分弱结合定位，亲电弹头在局部微环境中对特定 Cys 发生 Michael 加成；共价加合物增大侧链体积，产生 **空间位阻**，阻断 NEK7 招募 [10–13]。

**选择性来源：** 邻近效应（proximity-enabled reactivity）、界面去溶剂化、构象依赖性 Cys 暴露，而非深口袋几何匹配 alone。

### 2.3 路线 B：ATP 酶活/寡聚抑制型

**靶点位点（Tier 1）：**

| 位点 | 结构域 | 代表化合物 | 弹头类型 | 核心效应 | 参考文献 |
|------|--------|------------|----------|----------|----------|
| Cys598 | NACHT（ATP 结合区邻域） | Costunolide（木香烃内酯） | α-亚甲基-γ-丁内酯 | 抑制 ATPase 与寡聚 | Xu et al., *Acta Pharm. Sin. B* 2023 [14] |

**作用原理：** NACHT 结构域含 Walker A/B 等 ATP 结合元件；Cys598 位于核苷酸结合口袋邻域，共价修饰可干扰 ATP 结合/水解偶联的构象转换，从而抑制 NLRP3 自寡聚 [14]。

**与路线 A 的区别：** 结合区域相对更深、更疏水，机制读数侧重 **ATP 酶活与寡聚**，而非 NEK7 共免疫沉淀。

### 2.4 本方案不纳入主线的位点（说明）

- **Lys377（Manoalide）：** 赖氨酸亲电修饰，非 Cys 路径，需另设化学逻辑 [15]；
- **Cys548（衣康酸，小鼠为主）：** 内源性代谢调控机制，人源对应及药物化参照有限 [16]；
- **VLX1570 多 Cys 交联：** 交联聚集机制，选择性及安全性风险较高 [17]。

---

## 三、总体筛选架构

### 3.1 流程总览

```
┌─────────────────────────────────────────────────────────────────┐
│                    共用计算前端（化学库构建）                      │
│  天然产物/衍生物库 → 弹头分类 → GSH 反应性规则过滤                 │
└────────────────────────────┬────────────────────────────────────┘
                             │
           ┌─────────────────┴─────────────────┐
           ▼                                   ▼
┌──────────────────────┐            ┌──────────────────────┐
│  路线 A：PPI 界面型    │            │  路线 B：ATP/寡聚型   │
│  AF3 @ Cys279/280/409 │            │  AF3 @ Cys598        │
│  阳性对照：Oridonin等  │            │  阳性对照：Costunolide│
└──────────┬───────────┘            └──────────┬───────────┘
           │                                   │
           └─────────────────┬─────────────────┘
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│              共用实验确证平台（分叉验证）                          │
│  合成 → GSH t½ → DARTS → LC-MS/MS 位点 → 点突变 → 机制实验 → 体内 │
└─────────────────────────────────────────────────────────────────┘
```

### 3.2 设计原则

1. **分 site 运行 AF3**，不按 45 个 Cys 盲筛；
2. **每条路线独立化学系列**，共用实验平台；
3. **阳性对照 + 弹头阴性对照 + 位点突变** 三线并行；
4. **AF3 生成假说，MS/突变确证机制**，避免过度依赖计算排序。

---

## 四、计算筛选流程（详细）

### 4.1 Phase 0：化学库构建与分轨

#### 4.1.1 路线 A 化学空间

- **种子骨架：** Oridonin、DCL、RRx-001/INF39 类；
- **弹头：** 环外亚甲基-γ-内酯、α,β-不饱和酮、温和 Michael 受体；
- **库规模建议：** 50–200 个（天然产物衍生物 + 计算机生成类似物）；
- **避免：** 高反应性卤乙酰胺（除非专门优化 Cys409 系列）。

#### 4.1.2 路线 B 化学空间

- **种子骨架：** Costunolide 及 ATP 口袋相容疏水骨架；
- **弹头：** α-亚甲基-γ-丁内酯等；
- **库规模建议：** 50–200 个。

#### 4.1.3 共用预过滤

| 步骤 | 内容 | 工具/标准 |
|------|------|-----------|
| 结构标准化 | SMILES 规范化、去盐 | RDKit |
| 反应性过滤 | 排除 PAINS、过度活泼亲电体 | 规则过滤 |
| GSH 反应性预测 | 优先保留温和弹头 | 经验阈值 + 后续实验验证 |
| 类药性 | MW、cLogP、TPSA、SA score | RDKit、SA_Score [18] |

### 4.2 Phase 1：AF3 分位点共价共折叠筛选

#### 4.2.1 AF3 方法学依据

Shamir 等（2026）在 *JACS* 发表的研究表明，AF3 对共价蛋白–配体复合物具有高精度预测能力，并构建了共价虚筛基准集 COValid；利用 AF3 预测置信度指标对活性共价结合物与性质匹配诱饵（decoy）的区分接近最优，显著优于传统共价对接工具 CovDock 等。研究还在 BTK 激酶上完成了前瞻性共价虚筛，发现新型共价抑制剂并经共晶验证 [19]。

#### 4.2.2 AF3 输入规范

对每个候选分子 **分别** 提交 AF3 job，明确指定：

| 输入项 | 路线 A | 路线 B |
|--------|--------|--------|
| 蛋白序列 | 人源 NLRP3 全长或 NACHT 结构域 | 同左 |
| 共价连接 | Ligand warhead atom ↔ **Cys279 SG** | Ligand warhead atom ↔ **Cys598 SG** |
| | 同上，分别对 **Cys280、Cys409** 再跑 | 可选其他 ATP 区邻域 Cys |
| 配体构象 | 几何优化后的 3D 构象 | 同左 |
| 结构参考 | PDB 7ALV 等 NACHT 复合物 [20] | 同左，关注 ATP 结合区 |

#### 4.2.3 分 site 排序与置信度融合

1. **按位点分层排名：** 每个 Cys 独立输出 Top N；
2. **主排序指标：** AF3 复合物置信度（如 ipTM、pLDDT 等，依 AF3 版本）；
3. **可选重打分：** Rosetta 能量最小化（COValid 研究表明可改善部分靶点富集）；
4. **几何 QC：** 共价 C–S 键长 ~1.8 Å；共价前 SG–弹头距离合理；
5. **路线 A 额外判断：** 共价加合物是否占据 NEK7 结合界面、是否产生位阻；
6. **路线 B 额外判断：** 是否干扰 Walker A/B 或核苷酸结合区构象。

#### 4.2.4 回顾性校准（必做，未通过不进入大库筛）

| 路线 | 阳性对照 | 预期结果 |
|------|----------|----------|
| A | Oridonin @ Cys279 | 合理 pose + 高置信度 |
| A | DCL @ Cys280 | 同上 |
| A | INF39/RRx-001 @ Cys409 | 同上 |
| B | Costunolide @ Cys598 | 同上 |

另建议构建 **NLRP3 相关共价先导 + 拓扑/性质匹配 decoy** 小集，报告 EF@1%、BEDROC 等富集指标（参考 COValid 两阶段 decoy 设计）[19]。

### 4.3 Phase 2：计算次级验证

| 步骤 | 内容 | 两条线是否共用 |
|------|------|----------------|
| 短 MD（50–100 ns） | 共价键稳定性、配体 RMSD | 共用 |
| MM/GBSA | 结合自由能相对排序 | 共用 |
| 脱靶 Cys 风险 | 同一分子多位点 AF3 对比 | 共用 |
| 对接交叉验证（可选） | CovDock/Glide covalent 对照 | 共用 |

### 4.4 Phase 3：计算输出与合成优先级

每条路线从 AF3+MD 输出 **Top 10–15** 进入合成/采购；每条路线至少包含：

- 2–3 个高排名候选；
- 1 个阳性对照类似物；
- 1 个弹头饱和（还原）阴性对照。

**首轮实验合计建议：8–12 个化合物**（两线各 4–6 个）。

---

## 五、实验确证流程（详细）

### 5.1 Phase 1：化学与反应性（共用，2–4 周）

| 实验 | 方法 | 通过标准 |
|------|------|----------|
| 结构确证 | HRMS、NMR | 结构正确 |
| GSH 稳定性 | PBS + GSH（1–5 mM），37 °C，HPLC-MS | t₁/₂ 适中，非分钟级耗尽 |
| 细胞毒性初筛 | CCK-8 on BMDM/THP-1 | 排除高毒化合物 |

### 5.2 Phase 2：靶标结合（共用 DARTS，分叉 MS）

#### 5.2.1 重组蛋白

- 表达纯化人源 NLRP3 NACHT 结构域（或全长，与计算一致）；
- 质控：SDS-PAGE、活性/浓度标定。

#### 5.2.2 DARTS（Drug Affinity Responsive Target Stability）

- 化合物孵育 NLRP3 → 有限蛋白酶消化 → Western blot 检测 NLRP3 条带保护；
- **意义：** 快速排除无结合化合物（两线共用）。

#### 5.2.3 LC-MS/MS 肽段定位（机制分叉核心）

| 项目 | 路线 A | 路线 B |
|------|--------|--------|
| 孵育条件 | 过量化合物 + NLRP3，37 °C，2–4 h | 同左 |
| 酶解 | 胰蛋白酶变性消化 | 同左 |
| 搜库修饰 | 弹头质量修饰 on Cys | 同左 |
| 预期肽段 | 含 **Cys279/280/409** | 含 **Cys598** |
| 成功标准 | 主要/特异性修饰 PPI 位点 | 主要修饰 C598 |

### 5.3 Phase 3：位点必要性（定点突变挽救）

**细胞模型：** THP-1 或 iBMDM（LPS + ATP/Nigericin 激活炎性小体）

| 突变体 | 路线 | 预期 |
|--------|------|------|
| C279A、C280A、C409A | A | 对应化合物活性丧失或右移 |
| C598A | B | 对应化合物活性丧失或右移 |
| 非相关 Cys→Ala | 阴性对照 | 活性不应完全丧失 |

**表型 readout：**

- ELISA：成熟 IL-1β（上清）；
- Western blot：Caspase-1 p20、成熟 IL-1β、GSDMD 剪切；
- 可选：LDH（焦亡）、ASC speck 免疫荧光。

### 5.4 Phase 4：机制分叉验证

#### 路线 A 专项实验

| 实验 | 目的 | 预期 | 参考文献支持 |
|------|------|------|--------------|
| NLRP3–NEK7 共免疫沉淀 | PPI 阻断 | 化合物处理后 NLRP3–NEK7 结合 ↓ | [10–13] |
| ASC 斑点/寡聚 | 下游装配 | ASC oligomerization ↓ | [11] |
| CETSA（可选） | 靶标占据 | NLRP3 热稳定性变化 | 常规方法 |

#### 路线 B 专项实验

| 实验 | 目的 | 预期 | 参考文献支持 |
|------|------|------|--------------|
| NLRP3 ATPase 测定 | 酶活机制 | ATP 水解 ↓ | [14] |
| ASC 斑点 | 装配抑制 | ASC speck ↓ | [14] |
| 寡聚态分析（BN-PAGE/交联） | 寡聚 ↓ | 与机制一致 | [14] |

### 5.5 Phase 5：选择性与安全性（两线合并，各选 1–2 个 hit）

| 实验 | 内容 |
|------|------|
| 化学蛋白质组（ABPP） | 半胱氨酸组占用谱，评估选择性 |
| 炎性小体特异性 | NLRP3 KO/siRNA 挽救；NLRC4/AIM2 刺激不应 broad 抑制 |
| 安全性 | hERG、肝细胞毒性、原代巨噬细胞活力 |
| 与 HSA/BSA 孵育 | 评估白蛋白非特异性烷基化 |

### 5.6 Phase 6：体内概念验证（各线最多 1 个代表物）

| 疾病模型 | 适用路线 | 主要 readout |
|----------|----------|--------------|
| MSU 诱导痛风性腹膜炎 | A、B 均可 | IL-1β、中性粒细胞浸润 |
| LPS 内毒素血症 | A、B 均可 | IL-1β、生存率 |
| DSS 结肠炎 | 口服生物利用度可时 | 体重、病理、IL-1β |

**推进标准：** MS 确证位点 + 突变挽救 + 细胞 IC₅₀ 达标 + GSH 选择性合格。

---

## 六、决策标准与里程碑

### 6.1 路线 A 命中标准（PPI 型）

| 层级 | 标准 |
|------|------|
| 计算 | AF3 回顾性阳性对照通过；Cys279/280/409 中至少一位点高置信 |
| 化学 | GSH t₁/₂ 可接受；结构确证 |
| 生化 | MS 定位 PPI 位点；DARTS 阳性 |
| 细胞 | IL-1β、Caspase-1 抑制；C279A/C409A 等突变挽救 |
| 机制 | NLRP3–NEK7 互作 ↓ |
| 选择性 | 非 broad 炎性小体抑制；ABPP 可接受 |

### 6.2 路线 B 命中标准（ATP 型）

| 层级 | 标准 |
|------|------|
| 计算 | Costunolide @ Cys598 回顾性通过 |
| 化学 | 同路线 A |
| 生化 | MS 定位 Cys598；DARTS 阳性 |
| 细胞 | IL-1β 抑制；C598A 突变挽救 |
| 机制 | ATPase ↓ 和/或寡聚 ↓ |
| 选择性 | 同路线 A |

### 6.3 最小可行发表包（MVP）

| 内容 | 路线 A | 路线 B |
|------|--------|--------|
| 化合物数 | 3–4 | 3–4 |
| MS 位点 | Cys279 或 C409 | Cys598 |
| 突变 | C279A 或 C409A | C598A |
| 机制实验 | NEK7 Co-IP | ATPase 或 ASC speck |
| 体内 | 可选 1 个 | 可选 1 个 |

---

## 七、风险分析与应对

| 风险 | 描述 | 应对策略 |
|------|------|----------|
| 多 Cys 脱靶修饰 | 同一分子修饰多个位点 | MS 定量 + 组合突变 |
| 两线细胞表型相似 | IL-1β 均下降 | 必须用 NEK7 Co-IP vs ATPase 区分 |
| 弹头全身反应性 | GSH/白蛋白烷基化 | GSH 预筛 + 温和弹头 + ABPP |
| AF3 pose 与实验不符 | 计算偏差 | 以 MS/突变为准 |
| 激活态构象未覆盖 | 7ALV 为非活性/抑制剂结合态 | 讨论局限；可选多种构象建模 |
| 路线过于拥挤 | PPI 位点文献化合物多 | 路线 B 作差异化机制补充 |

---

## 八、时间阶段规划（按工作包）

| 阶段 | 工作内容 | 产出 |
|------|----------|------|
| W1–W4 | 化学库构建、AF3 回顾性校准、分 site 虚筛 | 两线各 Top 10–15 清单 |
| W5–W8 | 首轮 8–12 个化合物合成/采购、GSH、DARTS | 2–4 个生化 hit/线 |
| W9–W12 | MS 位点、突变细胞、IL-1β 表型 | 各位点必要性数据 |
| W13–W16 | 机制分叉（Co-IP / ATPase）、选择性 | 机制可区分 |
| W17–W20 | 代表物优化、体内验证（可选） | 先导假说或论文数据 |

*注：具体日历时间依实验室条件调整。*

---

## 九、预期成果与创新点

### 9.1 科学产出

1. 建立 **NLRP3 共价抑制双机制（PPI vs ATP）并行筛选** 的标准化流程；
2. 各获得 **1 类机制明确、位点 MS 确证** 的代表共价化学型；
3. 评估 AF3 在 NLRP3 多 Cys 场景下的 **分 site 共价虚筛** 实用性。

### 9.2 可申报/发表的创新表述

> 针对 NLRP3 共价修饰位点多样、单一虚筛策略难以覆盖不同机制的问题，本方案并行推进 PPI 界面型（Cys279/280/409）与 ATP/寡聚型（Cys598）两条共价抑制路线，以 AF3 分位点共价共折叠为核心计算引擎，配合 GSH–DARTS–MS–突变–机制分叉的实验闭环，实现早期资源共用、后期机制可区分的共价先导发现策略。

---

## 十、参考文献

1. Mangan MSJ, Olhava EJ, Roush WR, et al. Targeting the NLRP3 inflammasome in inflammatory diseases. *Nat Rev Drug Discov.* 2018;17(8):588-606. doi:10.1038/nrd.2018.97

2. Kelley N, Jeltema D, Duan Y, He Y. The NLRP3 inflammasome: An overview of mechanisms of activation and regulation. *Int J Mol Sci.* 2019;20(13):3328. doi:10.3390/ijms20133328

3. Swanson KV, Deng M, Ting JP-Y. The NLRP3 inflammasome: molecular activation and regulation to therapeutics. *Nat Rev Immunol.* 2019;19(8):477-489. doi:10.1038/s41577-019-0165-0

4. Coll RC, Robertson AAB, Chae JJ, et al. A small-molecule inhibitor of the NLRP3 inflammasome for the treatment of inflammatory diseases. *Nat Med.* 2015;21(3):248-255. doi:10.1038/nm.3806

5. Coll RC, Hill JR, Day CJ, et al. MCC950 directly targets the NLRP3 ATP-hydrolysis motif for inflammasome inhibition. *Nat Chem Biol.* 2019;15(6):556-559. doi:10.1038/s41589-019-0277-7

6. Singh J, Petter RC, Kluge AF. Targeted covalent drugs of the kinase family. *Curr Opin Chem Biol.* 2010;14(4):475-480. doi:10.1016/j.cbpa.2010.06.177

7. Boike L, Nomura DK, Cravatt BF. Advances in covalent drug discovery. *Nat Rev Drug Discov.* 2022;21(12):881-898. doi:10.1038/s41573-022-00516-z

8. Lin W, Wang P, Zhang Y, Lu W, Yang M, et al. AI-Driven Transfer Learning and Generative Model (TransGenGRU) Enables the Drug Discovery of Novel Natural Guaianolide Sesquiterpene Derivatives as Potent NLRP3 Inhibitors. *J Med Chem.* 2025;68(21):21534-21559. doi:10.1021/acs.jmedchem.5c01663

9. Abramovitch RB, Mucaj V, Prakash CSP, et al. A conserved cysteine motif is essential for inhibition of the NLRP3 inflammasome by a novel class of sulfonylurea-containing compounds. *Proc Natl Acad Sci USA.* 2015;112(37):11729-11734. doi:10.1073/pnas.1510286112

10. He H, Jiang H, Chen Y, Ye J, Wang A, Wang C, et al. Oridonin is a covalent NLRP3 inhibitor with strong anti-inflammasome activity. *Nat Commun.* 2018;9(1):2550. doi:10.1038/s41467-018-04947-6

11. Lv Q, Zhang Y, Wang J, Lin W, Xie Y, Yang H, et al. Dehydrocostus Lactone Effectively Alleviates Inflammatory Diseases by Covalently and Irreversibly Targeting NLRP3. *MedComm.* 2025;6(9):e70367. doi:10.1002/mco2.70367

12. Chen Y, He H, Lin B, Chen Y, Deng X, Jiang W, Zhou R. RRx-001 ameliorates inflammatory diseases by acting as a potent covalent NLRP3 inhibitor. *Cell Mol Immunol.* 2021;18(6):1425-1436. doi:10.1038/s41423-021-00683-y

13. Cocco M, Pellegrini C, Martinez-Banaclocha H, Giorgis M, Marini E, Costale A, et al. Development of novel irreversible NLRP3 inflammasome inhibitors (INF39) for the treatment of inflammatory bowel disease. *J Med Chem.* 2018;61(22):9964-9977. doi:10.1021/acs.jmedchem.8b01006

14. Xu H, Chen J, Chen P, Li W, Shao J, Hong S, et al. Costunolide covalently targets NACHT domain of NLRP3 to inhibit inflammasome activation and alleviate NLRP3-driven inflammatory diseases. *Acta Pharm Sin B.* 2023;13(2):678-693. doi:10.1016/j.apsb.2022.09.014

15. Li C, Lin H, He H, Ma M, Jiang W, Zhou R. Inhibition of the NLRP3 inflammasome activation by Manoalide ameliorates experimental autoimmune encephalomyelitis pathogenesis. *Front Cell Dev Biol.* 2022;10:822236. doi:10.3389/fcell.2022.822236

16. Hooftman A, Angiari S, Hester S, Corcoran SE, Runtsch MC, Ling C, et al. The immunomodulatory metabolite itaconate modifies NLRP3 and inhibits inflammasome activation. *Cell Metab.* 2020;32(3):468-478.e7. doi:10.1016/j.cmet.2020.07.016

17. Stanton C, Sun J, Nutsch K, Rosarda JD, Nguyen T, Li-Ma C, et al. Covalent targeting as a common mechanism for inhibiting NLRP3 inflammasome assembly. *ACS Chem Biol.* 2024;19(2):522-535. doi:10.1021/acschembio.3c00654

18. Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *J Cheminform.* 2009;1(1):8. doi:10.1186/1758-2946-1-8

19. Shamir Y, Gabizon R, Rogel A, Lin DYW, Andreotti AH, London N. Discovery of Covalent Ligands with AlphaFold3. *J Am Chem Soc.* 2026;148(12):13043-13054. doi:10.1021/jacs.5c22222

20. Dekker C, Mattes H, Wright M, et al. Crystal structure of NLRP3 NACHT domain with an inhibitor defines mechanism of inflammasome inhibition. *J Mol Biol.* 2021;433(24):167309. doi:10.1016/j.jmb.2021.167309

21. Hochheiser IV, Pilsl M, Hagelueken G, et al. Structure of the NLRP3 decamer bound to the cytokine release inhibitor CRID3. *Nature.* 2022;604(7904):184-189. doi:10.1038/s41586-022-04467-w

22. Lin H, Yang M, Li C, Lin B, Deng X, He H, Zhou R. An RRx-001 analogue with potent anti-NLRP3 inflammasome activity but without high-energy nitro functional groups. *Front Pharmacol.* 2022;13:822833. doi:10.3389/fphar.2022.822833

23. Mackay A, Velcicky J, Gommermann N, et al. Discovery of NP3-253, a potent brain penetrant inhibitor of the NLRP3 inflammasome. *J Med Chem.* 2024;67(24):20780-20798. doi:10.1021/acs.jmedchem.4c02350

24. Velcicky J, Langlois J-B, Wright M, et al. Discovery of NP3-742: a structurally diverse NLRP3 inhibitor identified through an unusual phenol replacement. *J Med Chem.* 2025;68(24):23532-23553. doi:10.1021/acs.jmedchem.5c02412

25. Velcicky J, Janser P, Gommermann N, et al. Discovery of potent, orally bioavailable, tricyclic NLRP3 inhibitors. *J Med Chem.* 2024;67(3):1544-1562. doi:10.1021/acs.jmedchem.3c02098

26. Jin X, Yang Y, Liu D, Zhou X, Huang Y, Ye D. Identification of a covalent NEK7 inhibitor to alleviate NLRP3 inflammasome-driven metainflammation. *Cell Commun Signal.* 2024;22(1):565. doi:10.1186/s12964-024-01919-w

27. Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature.* 2024;630(8016):493-500. doi:10.1038/s41586-024-07487-w

---

## 附录 A：化合物对照表（建议配置）

| 类型 | 路线 A（PPI） | 路线 B（ATP） |
|------|---------------|---------------|
| 阳性对照 | Oridonin、DCL 或 INF39 | Costunolide |
| 可逆参照 | MCC950 | MCC950（机制不同，仅作活性参照） |
| 阴性对照 | 弹头饱和还原类似物 | 弹头饱和还原类似物 |
| 机制参照 | — | 非共价 ATP 区配体（如有） |

## 附录 B：关键试剂与细胞

| 类别 | 推荐 |
|------|------|
| 细胞 | THP-1、iBMDM、原代 BMDM |
| 刺激 | LPS（100 ng/mL）+ ATP（5 mM）或 Nigericin |
| 蛋白 | 重组 NLRP3 NACHT（人源） |
| 抗体 | NLRP3、NEK7、ASC、Caspase-1、IL-1β |

---

**报告结束**

*本报告仅供研究方案设计参考；具体实验参数需根据实验室条件优化。所有参考文献均来自已正式发表的同行评议文献，DOI 可公开核验。*
