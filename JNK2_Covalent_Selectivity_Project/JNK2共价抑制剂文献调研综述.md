# JNK2 抑制剂文献调研综述（共价 + 可逆）

**编制日期：** 2026 年 7 月  
**关联项目：** `JNK2_Covalent_Selectivity_Project/`  
**版本：** v2.0

---

## 执行摘要

本综述覆盖 **JNK2 相关小分子抑制剂** 两条主线：**共价（Cys116）** 与 **可逆 ATP 竞争**，**不含 JNK3 共价（Cys154）** 文献。截至 2026 年 7 月：

| 主线 | 分类 | 代表化合物 | 亚型选择性 | 文献 |
|------|------|-----------|-----------|------|
| **共价** | JNK2 > JNK1 | YL5084 / YL2056 | kinact/KI ~21× | Lu 2023 [R1] |
| **共价** | JNK2/3 > JNK1 | 56d | kinact/KI(JNK2)=38,200 M⁻¹s⁻¹ | Wydra 2025 [R20] |
| **共价** | Pan-JNK 起点 | JNK-IN-8 | 三亚型近等效 | Zhang 2012 [R2] |
| **共价** | 可逆共价 pan-JNK | 1aR-IN-8 | 环己烯酮；GSH 稳定 | Tóth 2024 [R21] |
| **可逆** | JNK2/3 > JNK1 | 21b, 51d, 45a | **148–340×** vs JNK1 | Wydra/Park/Zheng [R10,R20,R29] |
| **可逆** | Pan-JNK 临床/工具 | CC-930, SP600125 | 无 JNK1/2 选择性 | [R11,R12] |
| **可逆** | JNK1 > JNK2（反向） | CC-90001 | 细胞 **12.9×** JNK1 偏好 | Nagy 2021 [R36] |
| **可逆** | JNK2 激酶 + PPI | **6l** | 双靶 MKK7–JNK2 + JNK2 激酶 | Chen 2026 [R37] |
| **工具** | 化学遗传 / 探针 | JNK2C116S, Qian 探针 | 功能验证 | Du [R25]; Qian [R26] |

**关键结论：**

1. **共价 JNK2 > JNK1** 仅 **YL5084** 充分表征；**Wydra 56d** 为 JNK2/3 共价第二独立验证（kinome clean）。
2. **可逆 JNK2/3 > JNK1** 以氨基吡唑 **21b/51d** 为最强表征（Wydra 2025 可逆系列）；**尚无** 可逆 **严格 JNK2> JNK1** 临床级 probe（CC-90001 反向偏 JNK1）。
3. **Ligand-first 逻辑贯穿共价与可逆**：先优化可逆 JNK2/3 fit（Leu106/Ile106）→ 再共价化（56d）或作可逆对照（21b）。
4. 共价 readout：**kinact/KI + 预孵育 IC50 + C116S/washout**；可逆 readout：**三亚型 IC50 + 细胞 p-c-Jun + JNK1/JNK2 KO 纤维母细胞**（CC-90001 范式）。
5. **JNK-IN-8 / Soleimani TNBC** 等应用文献存在 **JNK 非依赖 off-target**；表型须用 C116S 或 isoform-selective 化合物拆分。
6. **2026 最新可逆 JNK2 先导 6l**（Chen）：JNK2 激酶 + MKK7–JNK2 PPI 双机制，代表 upstream+PPI 新策略，**非共价**。

---

## 一、文献全景与 JNK2 相关性分级

### 1.1 为何 JNK2 亚型选择性极难

JNK1 与 JNK2 ATP 口袋 **98% 序列同一**，仅 **Met77/Leu77** 与 **Ile106/Leu106** 两处差异；JNK2 与 JNK3 亦仅 **Leu77/Met115** 一处差异 [R1, R8]。传统可逆抑制剂（如 CC-930、SP600125）几乎无法实现 JNK1 vs JNK2 选择性 [R9, R11]。

共价策略利用 **JNK 家族独有的 Cys116（JNK1/2）/ Cys154（JNK3）**（MAPK 激酶组中罕见），通过 **非共价预定位（KI）× 成键速率（kinact）** 两个维度实现亚型差异 [R1, R2, R13]。

### 1.2 本综述纳入标准

| 纳入 | 排除 |
|------|------|
| **共价：** Cys116 共价机制 + JNK2 活性/选择性数据 | **JNK3 共价**（Cys154；Muth、JC16I、Reynders 等） |
| **可逆：** 报告 JNK2 生化/细胞活性或 JNK2 功能依赖 | 纯 JNK3 选择性可逆（indazole 25c 等） |
| JNK1/2/3 至少两亚型对比，或 JNK2 特异性功能验证 | 仅 JNK 通路抑制剂、未测 JNK2 |
| 工具/化学遗传/计算（服务 JNK2 项目） | 仅专利摘要、无活性数据 |

### 1.3 完整文献清单（按 JNK2 相关性分级）

| Tier | 类型 | 文献 | JNK2 关联 |
|------|------|------|-----------|
| **T0 共价** | JNK2> JNK1 | Lu 2023 [R1] | YL5084；唯一严格 JNK2 共价先导 |
| **T0 共价** | JNK2/3> JNK1 | Wydra 2025 [R20] | 56d；kinome clean |
| **T1 共价** | Pan-JNK | Zhang 2012 [R2]; Tóth 2024 [R21] | JNK-IN-8 起点；可逆共价 |
| **T0 可逆** | JNK2/3> JNK1 | Wydra 2025 可逆系列 [R20]; Park [R10]; Zheng [R29] | 21b 148×；51d >340× |
| **T1 可逆** | Pan-JNK 临床/工具 | CC-930 [R12]; SP600125 [R11] | 可逆对照；无亚型选择性 |
| **T1 可逆** | JNK1> JNK2（反向） | CC-90001 [R36] | 纤维化 JNK1 偏好参照 |
| **T2 可逆** | JNK2 + PPI | Chen 2026 [R37] | 6l；MKK7–JNK2 双靶；ALI |
| **T3 工具** | 化学遗传/探针/计算 | Du [R25]; Qian [R26]; Liu [R27,R28] | C116S；live-cell；CpHMD |
| **T3 应用** | JNK-IN-8 下游 | Soleimani 2022 [R33] | off-target 警示 |
| **T4 综述** | 领域背景 | Wang [R9]; Koch [R34]; Feng 2024 [R38] | JNK 抑制剂综述/专利 |

---

## 二、共价 JNK2 抑制剂 — 分篇深度复盘

---

### [P1] Zhang et al. 2012 — JNK-IN-8：首个 irreversible pan-JNK 共价抑制剂

**文献：** Discovery of Potent and Selective Covalent Inhibitors of JNK. *Chemistry & Biology* 19(1):140–154.  
**DOI:** [10.1016/j.chembiol.2011.11.010](https://doi.org/10.1016/j.chembiol.2011.11.010) | **PMC:** [PMC3270411](https://pmc.ncbi.nlm.nih.gov/articles/PMC3270411/) | **PDB:** 3V6R (JNK3–JNK-IN-7)

#### 2.1.1 发现路径（结构导向，imatinib → acrylamide）

1. **起点：** imatinib 苯胺氨基嘧啶 **type-2** 骨架；在 c-Kit/PDGFR 晶体中发现 DFG 前 **保守 Cys** 可被共价弹头利用。
2. **JNK 序列比对：** JNK1/2 **Cys116**、JNK3 **Cys154** 为 MAPK 中独特位点 → 设计 **JNK-IN-1**（imatinib + 丙烯酰胺）。
3. **MedChem 迭代（Table 1）：**
   - 去 flag methyl → **JNK-IN-2**（IC50 改善 5–10×）
   - 优化 linker（1,4-苯酰胺 + 1,3-苯胺）→ **JNK-IN-5/7**（~500×  potency 提升）
   - 非共价对照 **JNK-IN-6**（丙酰胺替换丙烯酰胺）→ ~100× 失活 → 确证共价贡献
   - 引入 flag methyl → **JNK-IN-8**；替换吡啶为 2-苯基吡唑并[1,5-a]吡啶 → **JNK-IN-11**

#### 2.1.2 优化策略

- **Linker 几何：** 使 acrylamide 与 Cys154（JNK3 编号）形成更优 Michael 加成角度（JNK-IN-7 共晶验证）。
- **Kinome 选择性：** flag methyl 消除 IRAK1、PIK3C3、PIP4K2C、PIP5K3 等 off-target（与 imatinib 经验一致）。
- **末端替换：** 大体积芳环（JNK-IN-11）提高细胞 potency 但拓宽 kinome；苯并噻唑乙腈（JNK-IN-12）提高特异性。

#### 2.1.3 共价机制确证

| 实验 | 结果 |
|------|------|
| **ESI-MS 全蛋白** | JNK-IN-2 + JNK1 → +493 Da（单分子加成） |
| **LC-MS/MS 肽段** | 肽段 LMDANLC*QVIQME（JNK1 110–122）**Cys116** 独家修饰 |
| **JNK3 共晶** | JNK-IN-2 (2.60 Å)、JNK-IN-7 (2.97 Å)；Cys154 连续电子密度 |
| **JNK-IN-6 对照** | 无 acrylamide → ~100× 生化 IC50 右移 |
| **JNK2 C116S 突变** | JNK-IN-7/8 IC50 **≥100×** 右移；JNK-IN-11 仅 ~10× |
| **细胞 electrophoretic mobility** | JNK-IN-8 处理后 JNK 条带迁移率改变（共价加合物 marker） |

#### 2.1.4 亚型选择性（JNK1 vs JNK2 vs JNK3）

**JNK-IN-8 对三亚型近等效 — 非 JNK2 选择性化合物：**

| 指标 | JNK1 | JNK2 | JNK3 |
|------|------|------|------|
| 生化 IC50 (nM) | 4.67 | 18.7 | 0.98 |
| p-c-Jun EC50 Hela (nM) | 486 | — | — |
| p-c-Jun EC50 A375 (nM) | 338 | — | — |

- **JNK-IN-11** 对三亚型更均衡（IC50 0.50 nM 级），细胞 EC50 更低（A375 8.6 nM）。
- 选择性体现在 **JNK 家族 vs 其他激酶**（KINOMEscan S(10)=0.031 @ 1 μM），**非** JNK2> JNK1。

#### 2.1.5 测活指标与湿实验清单

**生化：** 放射性激酶 assay（JNK1/2/3 IC50）；Km/Vmax（C116S 突变体与 WT 相当）。

**细胞：**
- TR-FRET p-c-Jun（GFP-c-Jun + Tb-anti-pSer73）
- Western blot p-c-Jun（anisomycin 刺激 A375）
- ATP-biotin pulse-chase（细胞内共价标记动力学，~3 h 饱和）
- KiNativ（A375，1 μM，~200 激酶）
- KINOMEscan（442 激酶 @ 1 μM）+ Dundee 121 激酶酶学 panel
- 通路 sentinel 显微成像（Erk/p38/Akt/Stat/NF-κB/Rsk）

**选择性 gate：** C116S ≥100×；S(10) < 0.05；仅 on-pathway（c-Jun）抑制。

---

### [P2] Lu et al. 2023 — YL5084：首个 JNK2 > JNK1 共价抑制剂 ★

**文献：** Development of a Covalent Inhibitor of c-Jun N-Terminal Protein Kinase (JNK) 2/3 with Selectivity over JNK1. *J. Med. Chem.* 66(5):3356–3371.  
**DOI:** [10.1021/acs.jmedchem.2c01834](https://doi.org/10.1021/acs.jmedchem.2c01834) | **PMC:** [PMC11190964](https://pmc.ncbi.nlm.nih.gov/articles/PMC11190964/) | **PDB:** 8ELC, 7N8T

#### 2.2.1 发现路径（JNK-IN-8 骨架 MedChem + 意外亚型偏好）

1. **起点：** JNK-IN-8 / JNK-IN-11 的 2-苯基吡唑并[1,5-a]吡啶 scaffold。
2. **关键线索：** JNK-IN-11 区域异构体 **THZ-3-60-1** 在固定时间生化 assay 中显示 **JNK2 > JNK1** 趋势，但 kinome 差（EGFR、CK1、DDR1、CDK7、Aurora A/B）[R1]。
3. **THZ531 类比：** CDK7/12/13 共价抑制剂中，用 **(R)-3-氨基哌啶** 替换苯环可改善 kinome → 迁移至 JNK scaffold。
4. **饱和环扫描：** 哌啶、**吡咯烷**、氮杂环庚烷、环己烷 → **YL2056**（(R)-3-氨基吡咯烷）JNK2/JNK1 选择性最佳；对映体 YL2012 保留 JNK1 活性、选择性消失。

#### 2.2.2 优化策略

| 阶段 | 修饰 | 目的 | 代表数据 |
|------|------|------|---------|
| Scaffold | (R)-3-aminopyrrolidine 替换苯环 | JNK2 选择性 + kinome | YL2056: IC50 JNK1/JNK2 = 166/5 nM (33×) |
| 取代基扫描 | C6-Cl、C4-Me、C6-Me 等 | SAR 边界 | YL5189 C6-Me: 18/1 nM (18×) |
| Off-target 剔除 | 吡咯烷 **flag methyl**（JNK-IN-8 经验） | 消除 PIKFYVE 共价 off-target | YL5084: PIKFYVE Kd 5000 vs 270 nM (YL2056) |
| 动力学优化 | 保留 Cys116 共价，微调非共价 fit | kinact/KI 最大化 | 见 §2.2.4 |

#### 2.2.3 共价机制确证

| 实验 | 结果 |
|------|------|
| **LC-MS 全蛋白** | YL2056 + JNK2 (60 min, RT) → 单分子质量位移 |
| **CE-MS/MS 肽段** | **Cys116 独家修饰**（无其他 Cys 位点） |
| **8ELC 共晶 (2.0 Å)** | Cys116 SG 连续密度至 acrylamide β-碳；DFG-in |
| **YL2056R 对照** | 还原 acrylamide → 共价 off-target 表型消失（PIKFYVE 空泡） |

#### 2.2.4 活性数据（核心表）

**固定时间 IC50（1 h 孵育，生化）：**

| 化合物 | IC50 JNK1 (nM) | IC50 JNK2 (nM) | 选择性 (JNK1/JNK2) |
|--------|---------------|---------------|-------------------|
| YL2056 | 166 ± 35 | 5 ± 1 | **33×** |
| YL5084 | 2173 ± 90 | 70 ± 1 | **31×** |
| JNK-IN-8 | ~等效 | ~等效 | ~1× |

**kinact/KI（Sox 荧光肽底物 chemosensor）[R1]：**

| 化合物 | JNK1 (M⁻¹s⁻¹) | JNK2 (M⁻¹s⁻¹) | 倍数 |
|--------|--------------|--------------|------|
| YL5084 | 335 | **7166** | **~21×** |
| JNK-IN-8 | 近等效 | 近等效 | ~1× |

**JNK3：** Z'-LYTE IC50 = 84 ± 10 nM，但 **最大抑制百分比低于 JNK2**（偏好 JNK2）。

**细胞靶向：**

| 实验 | 结果 |
|------|------|
| biotin-JNK-IN-7 竞争 pull-down (MM.1S, 6 h) | YL5084 **~500 nM** 占据 JNK2；**2 μM** 仍难占据 JNK1 |
| NanoBRET (HEK293T) | YL5084 JNK2 活性 ≈ JNK-IN-8；JNK1 显著弱 |
| KiNativ (MM.1S, 1 μM) | JNK2 为主要保护靶点 |

**选择性结构机制 [R1]：**
- **Leu106 (JNK2/3) vs Ile106 (JNK1)：** 后口袋芳环占据；JNK1 back pocket clash
- **Val54 vs Ile54：** 主链柔性；JNK2 可诱导 fit ~0.6 Å
- **Arg50–Glu109 盐桥 (JNK1)：** P-loop 稳定，不利于结合；JNK2 Ile50 无此桥
- **500 ns MD + Glide 同源建模：** JNK2/3 ΔG_bind ≈ −8.6 kcal/mol vs JNK1 −6.1 kcal/mol

#### 2.2.5 局限性（原文承认）

- MM.1S 抗增殖 **不完全依赖 JNK2**（可能存在 p38α 等非共价 off-target；YL5084 vs YL5084R 对 p38α IC50 ~15 nM）
- 尚未完成全细胞共价靶标组学（chemoproteomics follow-up ongoing）

#### 2.2.6 湿实验完整清单

结构：X-ray (7N8T, 8ELC)；同源建模 (Prime)；分子对接 (Glide)；MD (500 ns)。

生化：固定时间 IC50；**kinact/KI (Sox sensor)**；Z'-LYTE JNK3；SelectScreen p38α；KdELECT PIKFYVE。

组学：KINOMEscan (>400)；KiNativ (237 激酶)。

细胞：biotin 探针竞争 WB；NanoBRET JNK1/2；PIKFYVE 空泡表型；MM.1S 增殖。

ADME：人/小鼠肝微粒体 t½ (16/11 min)；GSH t½ (46 min)。

---

### [P3] Wydra et al. 2025 — 56d：ligand-first 氨基吡唑 JNK2/3 共价抑制剂 ★

**文献：** A "Ligand First" Approach toward Selective, Covalent JNK2/3 Inhibitors. *J. Med. Chem.* 2025.  
**DOI:** [10.1021/acs.jmedchem.5c00884](https://doi.org/10.1021/acs.jmedchem.5c00884)

#### 2.3.1 发现路径（可逆 JNK2/3 先导 → 共价化）

1. **起点：** 已报道 **氨基吡唑 reversible JNK2/3 选择性** scaffold [R10, Park 2015]；Tübingen 组 pyridinylimidazole / aminopyrazole 长期积累。
2. **Ligand-first：** 在 **所有三亚型** 上并行测活，结构导向迭代 → 可逆先导 **>100× vs JNK1**。
3. **共价化：** 可逆先导 + **丙烯酰胺 warhead** → 靶向 JNK 保守 Cys（JNK1/2 Cys116，JNK3 Cys154）。
4. **Lead：** **56d** — JNK2/3 选择性共价先导，kinome clean。

#### 2.3.2 优化策略

- 利用 **Leu144/Ile106** 等 JNK2/3 vs JNK1 差异（与 Park 2015 氨基吡唑机制一致）[R10]
- 可逆先导充分优化 **代谢稳定性、细胞活性、PK** 后再共价化（降低 warhead 引入的 off-target 风险）
- 博士论文补充 [Wydra 2024]：肝再生 JNK2 探针需求 → JNK2 选择性 + JNK3 共抑制可接受

#### 2.3.3 可逆先导关键数据（共价化前）

| 化合物 | JNK1/JNK2 选择性 | 备注 |
|--------|-----------------|------|
| **16a** | **114×** (vs JNK2) | 去溶剂区后的 naphthyl 系列 |
| **16b** | 系列最高选择性 | 2-pyridyl 核心 |
| **21b** | **148×** | 苯甲酰胺 HR-II 优化；细胞 p-c-Jun ~1.11 μM |
| **21h** | 83× | JNK2 IC50=73 nM；JNK3=45 nM |
| **45a** | 12× | 最强可逆双 JNK2/3：JNK2 **4 nM**, JNK3 **6 nM** |
| **51d** | **>340×** | 倒置酰胺；完全丧失 JNK1 亲和力 |

**21b 细胞验证：** 小鼠 HCC 细胞（Nras G12V; Cdkn2a ARF−/−），sorbitol 渗透压刺激 + 2 h 预孵育；**21b @ 1.11 μM** 抑制 p-c-Jun 与 SP600125 相当。

**ADME/PK（可逆系列）：** 21b/21h/51d 小鼠肝微粒体 2 h **>90%** 剩余；21h i.v. t½=**3.33 h**；21b 细胞渗透优于带正电荷 piperidine 的 A-1。

#### 2.3.4 共价化设计

- **叠加策略：** 16a 倒置酰胺 scaffold + JNK-IN-7 共晶（PDB 3V6S）与 A-2/4WHZ 叠合 → **meta-aminobenzamide + acrylamide** linker 可达 Cys116。
- **56 系列：** acrylamide（56b, **56d**）vs 丙酰胺对照（56a, 56c）；**56d（meta 取代）** 为 lead。

#### 2.3.5 共价与选择性确证（全文级）

**33PanQuinase IC50（60 min, 无 ATP 预孵育）：** 56d 对 JNK2/3 有 modest 活性，**JNK1 无活性**；但此格式 **不足以** 评价共价效率。

**PhosphoSens TDI（AssayQuant, 60 min 预孵育, 无 ATP 竞争）：**

| 化合物 | JNK1 (0/60 min) | JNK2 (0/60 min) | JNK3 (0/60 min) |
|--------|----------------|----------------|----------------|
| 56b | >10,000 / >10,000 | >10,000 / >10,000 | >10,000 / >10,000 |
| **56d** | >10,000 / >10,000 | >10,000 / **25 nM** | >10,000 / **40 nM** |

**kinact/KI（PhosphoSens 全局拟合）：**

| 亚型 | kinact/KI (M⁻¹s⁻¹) |
|------|-------------------|
| JNK2 | **38,200** |
| JNK3 | **70,100** |

**IPMS：** 重组 JNK2 + 56d（5× 过量, 4.25 h, 20 °C）→ **单分子共价标记**（Mr=516 位移）。

**NanoBRET washout（JNK2, 2 h 预孵育）：** 56d **不可 washout**（vs 可逆参照 CTX-0294885）；56b 保留弱于 56d。

**GSH 稳定性：** 56d  extrapolated t½ = **88.9 h**（afatinib 8.64 h 为参照）；提示 warhead 选择性良好。

**NanoBRET TE（细胞 EC50）：**

| 化合物 | EC50 JNK2 | vs JNK1 |
|--------|----------|---------|
| **56d** | **883 nM** | **>11.3×** |
| 45c（可逆） | 1555 nM | >6.4× |
| 56b | 有活性 | 无亚型区分 |

**KINOMEscan（97 激酶 @ 500 nM）：** JNK2 POC=**1.4%**, JNK3 POC=**0.3%**；其余激酶 POC **>35%** → **clean profile**。Kd 剂量–反应 orthogonally 确认。

> **与 YL5084 对比：** 56d 为 **JNK2/3 双选择性**（非严格 JNK2> JNK1 单轴），但 kinact/KI(JNK2) 数值更高；独立 scaffold 验证 ligand-first 策略。全文/SI：[PMC12169684](https://pmc.ncbi.nlm.nih.gov/articles/PMC12169684/) | [SI](https://doi.org/10.1021/acs.jmedchem.5c00884.s001)。**尚无 56d 共晶 PDB**（叠合设计基于 4WHZ + 3V6S）。

#### 2.3.6 湿实验完整清单

33PanQuinase IC50；PhosphoSens TDI + kinact/KI；IPMS；NanoBRET TE + washout；GSH t½；KINOMEscan scanEDGE (97) + Kd；Western blot p-c-Jun（可逆系列）；MLM；小鼠 i.v. PK cassette；化学合成（SI Schemes S1–S6）。

---

### [P4] Tóth et al. 2024 — 环己烯酮可逆共价 JNK 抑制剂（JNK-IN-8 衍生）

**文献：** Reversible covalent c-Jun N-terminal kinase inhibitors targeting a specific cysteine by precision-guided Michael-acceptor warheads. *Nature Communications* 15:8269.  
**DOI:** [10.1038/s41467-024-52573-2](https://doi.org/10.1038/s41467-024-52573-2) | **PDB:** 8PTA, 8PT9, 8PT8 (JNK1–环己烯酮复合物)

#### 2.4.1 发现路径（warhead-first → JNK-IN-8 骨架嫁接）

1. **弹头发现：** electrophile-first 筛选 → **环己烯酮/戊烯酮** 可逆 Michael 受体（类萜立体 frustration 设计）。
2. **嫁接：** JNK-IN-8 的 IN-8 ATP 结合 moiety + 环己烯酮 → **1aR-IN-8**（主要对比物）。
3. **目标：** 解决 acrylamide 不可逆 + GSH 消耗 / 脱靶担忧。

#### 2.4.2 共价确证

| 实验 | 结果 |
|------|------|
| **JNK1 共晶** | 1aR/S-IN-8, 1a'R-IN-8；Cys116 共价加合物 Fo−Fc omit map |
| **SPR** | WT vs **C116S**：C116S 使 k_off 回升至非共价水平 |
| **2-step kinetic model** | k3/k4 可拟合；JNK-IN-8 k4=0（不可逆） |
| **Dialysis (5 d)** | 1aR-IN-8 可恢复 ~50% JNK1 活性；JNK-IN-8 不可 |
| **GSH 挑战 (10 mM, 18 h)** | JNK-IN-8 PhALC IC50 暴跌；1aR-IN-8 不变 |
| **SDS-PAGE mobility** | JNK-IN-8 改变 JNK 迁移率；1aR-IN-8 不改变（可逆） |
| **Photocage** | 笼锁吡啶 N → 蓝光释放 → 恢复 c-Jun 抑制（证明 ATP 口袋结合必要性） |

#### 2.4.3 活性数据

| 化合物 | PhALC IC50 (JNK1, +10 mM GSH) | NanoBRET TE (HEK293T) | p-c-Jun EC50 (SH-SY5Y MKK7 ACT) |
|--------|-------------------------------|----------------------|--------------------------------|
| JNK-IN-8 | 强（GSH 预孵育后失效） | ~10 nM | ~0.7 μM (WB) / ~100 nM (HTRF) |
| 1aR-IN-8 | 相当 | ~10 nM | 相当 |
| IN-8 (无 warhead) | 弱 | ~500 nM | — |
| CA-IN-8 (cyanoacrylamide) | 较弱 | ~200 nM | — |

**Kinome（1 μM，340 激酶，Reaction Biology）：** 1aR-IN-8 强抑制 JNK1/2/3；仅 LIMK1 (38% 剩余)、TNK1 (44%) 为 off-target；**20 个 JNK-IN-8 off-target 激酶均未被 1aR-IN-8 抑制**。

#### 2.4.4 亚型选择性（可调，非默认 JNK2）

**环己烯酮 C4 延伸可「编程」亚型偏好 [R21]：**

| 化合物 | 选择性趋势 |
|--------|-----------|
| 1aR-IN-8 | Pan-JNK（三亚型近等） |
| **1bR-IN-8** | **JNK1 > JNK2 ~10×** |
| **1cR-IN-8** | Pan-JNK 或 JNK3 偏好（Tóth C4 延伸；**非 JNK2 主线**） |

机制：C4 取代基投射至 **exon 6 编码的底物结合沟**（JNK1/2/3 序列差异区），环己烯酮 3D 形状感知 subtle 差异。

**对本项目的意义：** 证明 **共价 + 弹头侧链** 可实现亚型 tuning；若目标 JNK2> JNK1，需避免 1bR 型 C4 延伸方向，可参考 1aR/1cR SAR。

#### 2.4.5 湿实验清单

PhALC（自研 MAPK 活性 luciferase 互补）；NanoBRET TE + washout；SPR（WT/C116S）；X-ray；kinome scan ×2 平台交叉验证；AP-1 reporter；MKK7 ACT 工程细胞系；PROTAC (PRT_1) 降解 JNK1；大鼠 hepatocyte / plasma PK。

---

## 三、可逆 JNK2 抑制剂 — 分篇深度复盘

> 本章聚焦 **ATP 竞争性可逆** 小分子；JNK3 选择性可逆（indazole 25c 等）及 JNK3 共价文献 **不纳入**（见 §1.2）。

---

### [B1] Park 2015 + Zheng 2014 — 氨基吡唑可逆 JNK2/3 选择性（Wydra 前体）★

**Park 2015 [R10]：** Structural Basis and Biological Consequences for JNK2/3 Isoform Selective Aminopyrazoles. *Sci. Rep.* 5:8047.  
**DOI:** [10.1038/srep08047](https://doi.org/10.1038/srep08047) | **PDB:** 4WHZ (JNK3–A-2 衍生物)

**Zheng 2014 [R29]：** Design and Synthesis of Highly Potent and Isoform Selective JNK3 Inhibitors: SAR Studies on Aminopyrazole Derivatives. *J. Med. Chem.* 57(23):10013–10030.  
**DOI:** [10.1021/jm501256y](https://doi.org/10.1021/jm501256y) | **PDB:** 4WHZ (26k, 1.8 Å)

#### 3.1.1 发现与选择性机制

- **Scaffold：** 中央氨基吡唑 + HR-I 芳环 + 溶剂暴露 pyrazole/piperidine（A-1/A-2 系列）。
- **结构基础：** JNK3 **Leu144** / JNK2 **Leu106** vs JNK1 **Ile106** — HR-I 小甲基/亮氨酸差异允许 JNK2/3 后口袋更大芳环占据 [R10]。
- **Zheng SAR：** 26n @ 464 kinome 仅 7 个激酶 >80% @ 10 μM；26k 共晶确认 hinge H-bond（pyrazole N2 + amide NH → Met149）。

#### 3.1.2 关键活性

| 化合物 | 选择性声明 | 备注 |
|--------|-----------|------|
| A-1 (inhibitor) | ~30× vs JNK1 | Wydra 系列起点 |
| 26n | JNK3 高选择性 | 464 kinome clean |
| 16a/21b (Wydra) | 114–148× vs JNK1 | 见 [P3] / [B2] |

**对本项目：** Wydra **ligand-first** 的直接化学前体；Leu106/Ile106 机制与 YL5084 [R1] **一致但 scaffold 正交** — 可作为非 JNK-IN-8 骨架的 MedChem 参照。

---

### [B2] Wydra 2025 — 可逆氨基吡唑系列（21b/51d/45a；共价化前体）

**交叉引用：** 共价 lead **56d** 见 **[P3]**；本节汇总 **可逆** 数据（全文 [R20]）。

#### 3.2.1 策略

- **Ligand-first：** 在 JNK1/2/3 并行测活 → 可逆先导 **>100× vs JNK1** → 再引入 acrylamide（56d）。
- **Scaffold 来源：** Park/Zheng 氨基吡唑 [R10, R29]；倒置酰胺（51d）可 **完全消除 JNK1 亲和力**。

#### 3.2.2 关键可逆活性

| 化合物 | JNK1/JNK2 选择性 | JNK2 IC50 量级 | 细胞 / ADME |
|--------|-----------------|---------------|------------|
| **16a** | 114× | — | naphthyl 系列 |
| **21b** | **148×** | 双位数 nM | p-c-Jun ~**1.11 μM**（HCC, sorbitol） |
| **21h** | 83× | JNK2 **73 nM** | i.v. t½=**3.33 h** |
| **45a** | 12× (JNK2/3) | JNK2 **4 nM**, JNK3 **6 nM** | 最强可逆双 JNK2/3 |
| **51d** | **>340×** | 活性 | 倒置酰胺；JNK1 无结合 |

**湿实验：** 33PanQuinase 三亚型 IC50；Western blot p-c-Jun；小鼠肝微粒体稳定性；i.v. PK cassette。

**对本项目：** 可逆 JNK2/3> JNK1 的 **benchmark 系列**；共价筛选中 **21b/51d** 应作为可逆对照 Tier-1（与 CC-930 互补）。

---

### [B3] Bennett et al. 2001 — SP600125（pan-JNK 历史工具化合物）

**文献：** SP600125, an anthrapyrazolone inhibitor of Jun N-terminal kinase. *PNAS* 98(24):13681–13686.  
**DOI:** [10.1073/pnas.251194298](https://doi.org/10.1073/pnas.251194298) | **[R11]**

#### 3.3.1 活性与局限

| 亚型 | IC50 |
|------|------|
| JNK1 | **40 nM** |
| JNK2 | **40 nM** |
| JNK3 | **90 nM** |

- **机制：** ATP 竞争性、**可逆** anthrapyrazolone；广泛用于 p-c-Jun / AP-1 readout。
- **选择性 caveat：** 对 JNK 以外激酶 off-target 显著；**不可单独** 归因表型至 JNK2 [R9]。
- **在本综述中的角色：** Wydra 21b 细胞对照 [B2]；共价/可逆实验的 **pan-JNK 历史参照**。

---

### [B4] Plantevin Krenitsky et al. 2012 — CC-930（Tanzisertib，可逆 pan-JNK 临床化合物）

**文献：** Discovery of CC-930, an orally active anti-fibrotic JNK inhibitor. *Bioorg. Med. Chem. Lett.* 22(3):1433–1438.  
**DOI:** [10.1016/j.bmcl.2011.12.010](https://doi.org/10.1016/j.bmcl.2011.12.010) | **[R12]**

#### 3.4.1 活性

| 指标 | 数据 |
|------|------|
| Ki JNK1/2/3 | **5–61 nM**（三亚型近等） |
| 细胞 | c-Jun 磷酸化抑制 |
| 临床 | Phase II IPF（后终止；获益/风险 profile 不足） |

- **Scaffold：** 2,4-dialkylamino-pyrimidine-5-carboxamide 系列；**无 JNK1 vs JNK2 亚型选择性**。
- **对本项目：** Phase 0 **可逆 pan-JNK 阴性对照**；与 YL5084 共价选择性形成 readout 对比（§7.2）。

---

### [B5] Nagy et al. 2021 — CC-90001（JNK1 > JNK2 可逆；纤维化反向参照）★

**文献：** Discovery of the c-Jun N-Terminal Kinase Inhibitor CC-90001. *J. Med. Chem.* 64(24):18193–18208.  
**DOI:** [10.1021/acs.jmedchem.1c01716](https://doi.org/10.1021/acs.jmedchem.1c01716) | **PMID:** [34894681](https://pubmed.ncbi.nlm.nih.gov/34894681/) | **[R36]**

#### 3.5.1 发现路径

- **生物学动机：** 纤维化病理中 **JNK1 而非 JNK2/JNK3** 可能为主要驱动 → 在 CC-930 基础上优化 **JNK1 偏好**。
- **MedChem：** 2,4-dialkylamino-pyrimidine-5-carboxamide SAR；优化 PK 后选定 **CC-90001** 进入临床（Phase II IPF, NCT03142191）。

#### 3.5.2 活性数据

| 指标 | 数据 |
|------|------|
| Ki JNK1 / JNK2 / JNK3 | **0.3 / 0.4 / 0.6 nM**（生化近等） |
| 细胞亚型选择性 | **JNK1 > JNK2 12.9×**（JNK1/JNK2 KO 纤维母细胞） |
| LPS c-Jun pEC50 | ~480 nM |
| Kinome | >1000× vs p38/ERK/Akt 等 |
| 体内 | 口服暴露优；IPF 纤维母细胞 TGF-β 胶原抑制 |

**对本项目：** **反向亚型选择性** 参照 — 证明可逆系列可实现 **JNK1> JNK2**（与 YL5084 共价 JNK2> JNK1 对偶）；细胞 readout 须用 **isoform KO** 而非仅生化 Ki。

---

### [B6] Chen et al. 2026 — 6l：JNK2 激酶 + MKK7–JNK2 PPI 双靶可逆先导 ★

**文献：** Discovery of 6l with a New Skeleton of 6-Oxo-*N*-(4-(Quinolin-4-yloxy)phenyl)-1,6-dihydropyridine-3-carboxamide as a Dual-Action Inhibitor against JNK2 and the MKK7–JNK2 Protein–Protein Interaction for Acute Lung Injury. *J. Med. Chem.* 69(5):5648–5676.  
**DOI:** [10.1021/acs.jmedchem.5c02860](https://doi.org/10.1021/acs.jmedchem.5c02860) | **PMID:** [41729167](https://pubmed.ncbi.nlm.nih.gov/41729167/) | **[R37]**

#### 3.6.1 发现路径

1. **靶点：** ALI 中 **JNK2** 驱动 MAPK 炎症；除激酶域外，**MKK7–JNK2 PPI** 为 upstream 调控节点。
2. **方法：** 结构导向虚拟筛选 + SAR → 新骨架 **6-oxo-1,6-dihydropyridine-3-carboxamide**。
3. **Lead 6l：** **非共价**；同时抑制 JNK2 激酶活性 **与** MKK7–JNK2 结合。

#### 3.6.2 活性与机制

| 指标 | 数据 |
|------|------|
| JNK2 激酶 IC50 | **0.99 μM** |
| MKK7–JNK2 PPI Kd | **81.6 μM**（双靶） |
| IL-6（THP-1, LPS） | IC50 = **0.14 μM** |
| TNF-α（J774A, LPS） | IC50 = **0.55 μM** |
| 机制 | ↓ c-Jun 磷酸化；↓ 炎症因子 |
| 体内 | LPS / CLP 诱导 ALI 小鼠有效；PK / 安全性良好 |

**对本项目：** 2026 最新 **可逆 JNK2 导向** 先导；代表 **激酶 + PPI 双机制** 新策略（与 Cys116 共价路线正交）。**非共价**，不应与 YL5084/56d 共价 gate 混用。

---

## 四、工具、化学遗传与计算辅助

---

### [P8] Du et al. 2019 — JNK2C116S 化学遗传 + JNK-IN-8 功能验证

**文献：** JNK2 Is Required for the Tumorigenic Properties of Melanoma Cells. *ACS Chem. Biol.* 14(5):994–1003.  
**DOI:** [10.1021/acschembio.9b00083](https://doi.org/10.1021/acschembio.9b00083)

#### 2.8.1 方法

- 在黑色素瘤细胞中表达 **WT 或 C116S** JNK1/JNK2；JNK-IN-8 仅抑制 WT，C116S **≥100× 耐药**。
- **化学遗传 readout：** 在共表达背景下，JNK-IN-8 处理可 **分离** JNK1 vs JNK2 依赖表型。

#### 2.8.2 主要结论

| 表型 | JNK1 | JNK2 |
|------|------|------|
| 黑色素瘤增殖/侵袭 | 次要 | **必需** |
| BRAFi 耐药 | — | JNK2 依赖 |

**对本项目：** 确立 **C116S 化学遗传** 作为 JNK2 功能验证 gold standard；与 YL5084 biotin 竞争 [R1] 互补。注意：Du 2019 使用 **pan-JNK JNK-IN-8**，非 JNK2 选择性化合物。

---

### [P9] Qian et al. 2019 — JNK-IN-8 衍生 live-cell 探针

**文献：** Live-cell imaging and profiling of c-Jun N-terminal kinases using covalent inhibitor-derived probes. *Chem. Commun.* 55:1092–1095.  
**DOI:** [10.1039/C8CC09558B](https://doi.org/10.1039/C8CC09558B)

#### 2.9.1 内容

- 基于 JNK-IN-8 scaffold 设计 **共价探针**（Michael 受体 + 双光子 Turn-ON 荧光 / chemoproteomic handle）。
- 联合 **TBET 荧光 turn-on** 实现免洗 live-cell JNK 成像；chemoproteomics 验证靶标 engagement。

#### 2.9.2 对本项目

- 提供 **JNK-IN-8 探针化** 先例（biotin-JNK-IN-7 竞争实验 [R1] 的同源思路）。
- 探针 **不区分 JNK1/2/3** — 需配合 isoform-specific 抑制剂（YL5084/56d）做竞争实验方可推断 JNK2 占据。

---

### [P10] Liu et al. 2021/2022 — CpHMD 预测 JNK Cys116/154 反应性

**Liu 2021 [R27]：** Profiling MAP kinase cysteines for targeted covalent inhibitor design. *RSC Med. Chem.* 13:773–781.  
**DOI:** [10.1039/D1MD00277E](https://doi.org/10.1039/D1MD00277E)

**Liu 2022 [R28]：** Reactivities of the Front Pocket N-Terminal Cap Cysteines in Human Kinases. *J. Med. Chem.* 65(2):1525–1535.  
**DOI:** [10.1021/acs.jmedchem.1c01186](https://doi.org/10.1021/acs.jmedchem.1c01186) | **PMC:** [PMC8812259](https://pmc.ncbi.nlm.nih.gov/articles/PMC8812259/)

#### 2.10.1 核心预测（JNK 相关）

| 位点 | 亚型 | CpHMD pKa | 反应性 |
|------|------|-----------|--------|
| Cys116 | JNK1 | 7.5 | reactive |
| Cys116 | JNK2 | 8.0 | reactive |
| Cys154 | JNK3 | **6.3** | **hyper-reactive** |

- 机制：JNK3 Cys154 额外受 **Ser72** H-bond 稳定 → pKa 更低；解释 JNK 家族内 **Cys154 反应性高于 Cys116** [R28]。
- 14 个 MAPK 系统扫描：FP Ncap (Cys116)、EFP、DFG−1 三处为可药化共价位点。

#### 4.3.2 对本项目

- 支持 Cys116 acrylamide 策略的 **理论合理性**；但 **不能替代** Phase 0 中 YL5084 C116S gate [R19]。
- JNK3 Cys154 更高反应性 → 共价 MedChem 中若需 **严格 JNK2> JNK3**，须依赖非共价 fit（Leu106 等）而非 Cys 位点本身 [R1, R20]。

---

## 五、应用警示与领域综述

---

### [P11] Soleimani et al. 2022 — JNK-IN-8 在 TNBC 中的应用（off-target 警示）

**文献：** Covalent JNK Inhibitor, JNK-IN-8, Suppresses Tumor Growth in Triple-Negative Breast Cancer by Activating TFEB- and TFE3-Mediated Lysosome Biogenesis and Autophagy. *Mol. Cancer Ther.* 21(10):1847–1858.  
**DOI:** [10.1158/1535-7163.MCT-21-1044](https://doi.org/10.1158/1535-7163.MCT-21-1044)

#### 5.1.1 机制（重要 caveat）

- JNK-IN-8 抑制 TNBC 生长，但 **TFEB/TFE3 激活 / mTOR 抑制 / 溶酶体生物合成** 与 JNK **无关**：
  - JNK1/2 KO **不影响** TFEB 去磷酸化；
  - WT 或 **C116S** JNK 回补 **不能逆转** TFEB 表型。
- JNK 依赖部分仅体现在 **集落形成** 等 readout。

#### 5.1.2 对本项目

- **强警示：** pan-JNK 共价化合物（尤其 JNK-IN-8）的 **细胞表型不可直接归因于 JNK2**；必须配合 C116S / YL5084R / isoform-selective 化合物（YL5084, 56d）做机制拆分。
- 与 Lu 2023 YL5084 MM 增殖 **非 JNK2 完全依赖** [R1] 形成呼应。

---

### [P12] 领域综述与专利 landscape

| 文献 | 角色 |
|------|------|
| Wang 2022 [R9] | JNK 抑制剂设计综述（含共价 + 可逆） |
| Koch 2014 [R34] | JNK 专利 landscape（2010–2014） |
| Feng 2024 [R38] | JNK 专利更新（2015–present；含 CC-90001、共价进展） |
| Dou 2025 [R35] | JNK3 选择性机制综述（**背景**；非 JNK2 主线） |
| Boike 2022 [R13] | 共价药物发现方法学 |

**对本项目：** 背景阅读；**不替代** [P1]–[P4] 共价核心与 [B1]–[B6] 可逆数据。

---

## 六、横向对比表

### 6.1 发现策略对比（共价 + 可逆）

| 文献 | 类型 | 策略 | 起点 | 关键转折 |
|------|------|------|------|---------|
| Zhang 2012 | 共价 | Structure-based | imatinib type-2 | linker 优化 → acrylamide 几何 + flag methyl |
| Lu 2023 | 共价 | Lead optimization | JNK-IN-8/11 | THZ-3-60-1 异构体 hint → (R)-aminopyrrolidine |
| Wydra 2025 | 共价 | **Ligand-first** | aminopyrazole (A-1/21b) | 148× 可逆 JNK2/3 → inverted amide + acrylamide → 56d |
| Park/Zheng | 可逆 | Reversible SAR | aminopyrazole | Leu144/Ile106 → Wydra 前体 |
| Wydra 可逆 | 可逆 | Ligand-first | A-1/16a | 21b 148×; 51d >340× vs JNK1 |
| Tóth 2024 | 共价 | **Warhead-first** | cyclohexenone | 嫁接到 IN-8；C4 调 isoform |
| Nagy 2021 | 可逆 | Clinical SAR | CC-930 系列 | JNK1 纤维化生物学 → CC-90001 **12.9× JNK1** |
| Chen 2026 | 可逆 | VS + SAR | 新 pyridine 骨架 | JNK2 激酶 + MKK7–JNK2 PPI 双靶 → 6l |
| Du 2019 | 工具 | Chemical genetics | JNK-IN-8 + C116S | 功能分离 JNK1 vs JNK2 |

### 6.2 共价确证实验矩阵（仅共价化合物）

| 实验类型 | Zhang/JNK-IN-8 | Lu/YL5084 | Wydra/56d | Tóth/1aR-IN-8 |
|---------|---------------|-----------|-----------|--------------|
| 全蛋白 LC-MS / IPMS | ✓ LC-MS | ✓ LC-MS | ✓ **IPMS 单标记** | — |
| 肽段 MS/MS 定位 Cys | ✓ (C116) | ✓ (C116) | — | — |
| Cys→Ser/Ala 突变 | ✓ C116S ≥100× | ✓ (系列) | — (washout 替代) | ✓ C116S SPR |
| 共晶共价键密度 | ✓ C154 (JNK3 结构) | ✓ C116 (8ELC) | — (叠合设计) | ✓ C116 (8PTA) |
| 非共价对照 (R/-H) | ✓ JNK-IN-6 | ✓ YL2056R | ✓ 56a/56c 丙酰胺 | ✓ IN-8a |
| 预孵育 IC50 右移 | — | 固定 1 h | ✓ PhosphoSens 60 min | — |
| GSH 稳定性 | — | t½ 46 min | ✓ t½ **88.9 h** | ✓ 10 mM 挑战 |
| NanoBRET washout | — | TE | ✓ **不可逆** | ✓ 可逆 |

### 6.3 亚型选择性 readout 对比

| Readout | 用途 | 最佳应用文献 |
|---------|------|-------------|
| **kinact/KI** | 共价效率 + 非共价预定位综合 | YL5084 (21×), 56d (38200/70100 M⁻¹s⁻¹) |
| **预孵育 IC50 (PhosphoSens)** | 共价 TDI；无 ATP 竞争 | 56d (JNK2 25 nM @ 60 min) |
| **固定时间 IC50** | 简单快速；受 incubation time 影响 | YL5084 (31× @ 1 h) |
| **C116S 突变** | 共价必要性 | JNK-IN-8 (≥100×), Du 2019 化学遗传 |
| **NanoBRET washout** | 细胞内不可逆结合 | 56d, 1aR-IN-8 |
| **Biotin 探针竞争** | 细胞内 JNK1 vs JNK2 占据 | YL5084 (500 nM JNK2 only) |
| **NanoBRET TE** | 活细胞靶向 EC50 | YL5084, 56d (883 nM JNK2), 1aR-IN-8 |
| **JNK isoform KO 细胞** | 可逆亚型选择性 | CC-90001 (JNK1 **12.9×** JNK2) |
| **KiNativ / KINOMEscan** | 组学选择性 | JNK-IN-8, YL5084, 56d (97 kin @ 500 nM) |
| **CpHMD pKa** | Cys 反应性预测（计算） | Liu 2021/2022 (Cys116 pKa 7.5–8.0) |
| **同源建模 + MD** | Leu106/Ile106 机制解释 | YL5084 (500 ns) |

### 6.4 活性数据 master 表

#### 6.4.1 共价 JNK2 相关

| 化合物 | JNK1 | JNK2 | JNK3 | 选择性声明 | 细胞 EC50 / TE |
|--------|------|------|------|-----------|---------------|
| JNK-IN-8 | 4.7 nM | 18.7 nM | 1.0 nM | Pan-JNK; JNK vs kinome | p-c-Jun ~300–500 nM |
| YL2056 | 166 nM | 5 nM | — | 33× (1 h IC50) | — |
| **YL5084** | 2173 nM | 70 nM | 84 nM† | **21× kinact/KI; 31× IC50** | JNK2 ~500 nM TE |
| **56d** | >10 μM‡ | 25 nM§ | 40 nM§ | kinact/KI=38200/70100 | NanoBRET JNK2 **883 nM** (>11× vs JNK1) |
| 1aR-IN-8 | ~等 | ~等 | ~等 | Pan-JNK; kinome 更 clean | NanoBRET ~10 nM |
| 1bR-IN-8 | 优 | 10× 弱 | — | JNK1 > JNK2 | — |

#### 6.4.2 可逆 JNK2 相关

| 化合物 | JNK1 | JNK2 | JNK3 | 选择性声明 | 细胞 / 功能 |
|--------|------|------|------|-----------|------------|
| **21b** | 弱 | 活性 | 活性 | **148× vs JNK1** | p-c-Jun ~1.11 μM |
| **51d** | 无 | 活性 | 活性 | **>340× vs JNK1** | — |
| **45a** | — | 4 nM | 6 nM | JNK2/3 双强 | 可逆对照 |
| SP600125 | 40 nM | 40 nM | 90 nM | Pan-JNK 工具 | 广泛引用；off-target 多 |
| CC-930 | Ki 5–61 nM | 近等 | 近等 | Pan-JNK 临床 | Phase II IPF |
| **CC-90001** | Ki 0.3 nM | Ki 0.4 nM | Ki 0.6 nM | 生化近等；**细胞 JNK1 12.9×** | IPF Phase II |
| **6l** | — | IC50 **0.99 μM** | — | JNK2 + MKK7 PPI | IL-6 **0.14 μM** (THP-1) |

† JNK3 Z'-LYTE IC50；最大抑制低于 JNK2。  
‡ 33PanQuinase 无预孵育。  
§ PhosphoSens **60 min 预孵育**（无 ATP 竞争）。

---

## 七、对本项目（JNK2 选择性共价筛选）的启示

### 7.1 阳性对照分级

| Tier | 化合物 | 类型 | 预期 |
|------|--------|------|------|
| **Tier-0** | YL5084 | 共价 | kinact/KI JNK2/JNK1 ≥ 10×；C116S 右移 |
| **Tier-0b** | 56d | 共价 | kinact/KI(JNK2) ≥ 10⁴ M⁻¹s⁻¹；PhosphoSens 预孵育 JNK2 双位数 nM；NanoBRET washout |
| **Tier-1 共价** | YL2056, JNK-IN-8 | 共价 | pan-JNK / 弱 JNK2 偏好参照 |
| **Tier-1 可逆** | 21b, 51d, CC-930, SP600125 | 可逆 | JNK2/3> JNK1 或 pan-JNK 阴性对照 |
| **Tier-1 反向** | CC-90001 | 可逆 | JNK1> JNK2 细胞 12.9×；纤维化参照 |
| **Tier-2** | YL5084R / YL2056R / 56a | 共价阴性 | 非共价对照 |
| **Tier-3** | 1aR-IN-8 | 可逆共价 | GSH 稳定参照 |
| **Tier-4** | 6l | 可逆双靶 | JNK2+PPI；ALI 机制参照（非共价 gate） |

### 7.2 最小湿实验 package（复制文献共识）

1. **生化：** 固定时间 IC50 (JNK1/2/3) + **kinact/KI (Sox sensor)**
2. **共价：** LC-MS 全蛋白 + 肽段 MS/MS；**JNK2 C116S** ≥10× 右移
3. **结构：** 8ELC 对接/共晶（若有条件）
4. **细胞：** biotin-JNK-IN-7 竞争 WB 或 NanoBRET (JNK1 vs JNK2 双 construct)
5. **选择性：** KINOMEscan @ 1 μM 或 KiNativ

### 7.3 计算筛选注意

- **主受体：** 8ELC（DFG-in 共价）；**禁用** 3NPC (DFG-out) [R1, R3]
- **不依赖 Δsel 对接排名**（JNK1 非共价项目已证伪）[R17]
- AF3 共价虚筛需 Phase 0 阳性对照校准 [R19]

---

## 八、参考文献（本综述编号）

| 编号 | 文献 | 角色 |
|------|------|------|
| [R1] | Lu et al., *J. Med. Chem.* 2023 — YL5084 | **JNK2> JNK1 共价核心** |
| [R2] | Zhang et al., *Chem. Biol.* 2012 — JNK-IN-8 | Pan-JNK 共价起点 |
| [R9] | Wang et al., *J. Med. Chem.* 2022 | JNK 抑制剂综述 |
| [R10] | Park et al., *Sci. Rep.* 2015 — aminopyrazole | 可逆 JNK2/3；Leu144 机制 |
| [R11] | Bennett et al., *PNAS* 2001 — SP600125 | Pan-JNK 工具 |
| [R12] | Plantevin Krenitsky et al., *BMCL* 2012 — CC-930 | 可逆 pan-JNK 临床 |
| [R17] | JNK1 项目报告 v2.8 | Δsel 证伪 |
| [R19] | NLRP3 AF3 难度报告 | AF3 gate 警示 |
| [R20] | Wydra et al., *J. Med. Chem.* 2025 — 56d | JNK2/3 ligand-first 共价 + 可逆系列 |
| [R21] | Tóth et al., *Nat. Commun.* 2024 | 可逆共价 + isoform tuning |
| [R25] | Du et al., *ACS Chem. Biol.* 2019 | JNK2C116S 化学遗传 |
| [R26] | Qian et al., *Chem. Commun.* 2019 | JNK-IN-8 live-cell 探针 |
| [R27] | Liu et al., *RSC Med. Chem.* 2021 | MAPK Cys CpHMD |
| [R28] | Liu et al., *J. Med. Chem.* 2022 | FP Ncap Cys 反应性 |
| [R29] | Zheng et al., *J. Med. Chem.* 2014 | 氨基吡唑 SAR；4WHZ |
| [R33] | Soleimani et al., *Mol. Cancer Ther.* 2022 | JNK-IN-8 TNBC；off-target 警示 |
| [R34] | Koch/Gehringer 2014 等 | JNK 专利 landscape |
| [R35] | Dou et al., *J. Med. Chem.* 2025 | JNK3 选择性综述（背景） |
| [R36] | Nagy et al., *J. Med. Chem.* 2021 — CC-90001 | **JNK1> JNK2 可逆** |
| [R37] | Chen et al., *J. Med. Chem.* 2026 — 6l | **JNK2 + MKK7 PPI 可逆** |
| [R38] | Feng et al., *Expert Opin. Ther. Pat.* 2024 | JNK 专利更新 2015–present |

> **已排除（JNK3 共价，不纳入本综述）：** R22 Muth 2016、R23 Wen/JC16I 2024、R30 Wen/25c、R31 Reynders 2021、R32 Hoffelner 2023 — 详见 `REFERENCES.md` 附录。

完整 DOI 列表见 `REFERENCES.md`。

---

## 九、文献空白与未来方向

1. **严格 JNK2> JNK1 且 kinome-clean 的化学 probe** 尚未完全达成：YL5084 有 p38α 非共价 off-target；56d 为 JNK2/3 双轴非单 JNK2；可逆侧 **尚无** 临床级严格 JNK2> JNK1（CC-90001 反向偏 JNK1）[R1, R20, R36]。
2. **JNK2 选择性共晶** 仅 8ELC（YL2056）；YL5084 与 56d **均无独立 PDB**。
3. **56d C116S 突变** 数据未在正文报告（以 washout/IPMS 替代）；建议 follow-up 补充。
4. **可逆 JNK2-selective**（非 JNK2/3 双轴）仍缺系统 SAR；6l 代表 **激酶+PPI** 新方向但 JNK2 激酶 IC50 仅 μM 级 [R37]。
5. **AF3/对接** 对 JNK2 Cys116 共价虚筛的 prospective 验证（COValid 式）仍缺系统 benchmark [R15, R19]。

---

*本文档 v2.0 为 `JNK2选择性共价抑制剂筛选方案.md` 的文献支撑章节，覆盖 **共价 + 可逆 JNK2 抑制剂**（不含 JNK3 共价），与 Phase 0 实验设计直接对应。*
