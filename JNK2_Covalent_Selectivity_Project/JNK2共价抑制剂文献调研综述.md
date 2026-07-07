# JNK2 小分子共价抑制剂 — 已发表文献调研综述

**编制日期：** 2026 年 7 月  
**关联项目：** `JNK2_Covalent_Selectivity_Project/`  
**版本：** v1.0

---

## 执行摘要

截至 2026 年，已发表的小分子 JNK **共价**抑制剂文献可归纳为 **五条主线**：

| 分类 | 代表化合物 | 亚型选择性 | 共价位点 | 文献 |
|------|-----------|-----------|---------|------|
| **A. JNK2 > JNK1 共价（核心）** | YL5084 / YL2056 | kinact/KI ~21×；固定时间 IC50 ~31× | Cys116 | Lu et al. 2023 [R1] |
| **B. JNK2/3 > JNK1 共价** | 56d | 可逆先导 >100× vs JNK1；kinact/KI(JNK2)=38,200 M⁻¹s⁻¹ | Cys116/154 | Wydra et al. 2025 [R20] |
| **C. Pan-JNK 共价起点** | JNK-IN-8 | 三亚型近等效；JNK 家族选择性高 | Cys116/154 | Zhang et al. 2012 [R2] |
| **D. 可逆共价 + 可调亚型** | 1aR-IN-8, 1cR-IN-8 | 弹头侧链可偏 JNK1 或 JNK3 | Cys116 | Tóth et al. 2024 [R21] |
| **E. JNK3 选择性（非 JNK2 主线）** | 化合物 7, JC16I | JNK3 >> JNK1/2 | Cys154 | Muth 2016 [R22]; Wen 2024 [R23] |

**关键结论：**

1. 严格 **JNK2 > JNK1** 的共价小分子，目前以 **Gray/Westover 组 YL5084 系列** 为唯一充分表征的先导；**Wydra 2025** 提供独立 scaffold（氨基吡唑）的 JNK2/3 共价验证。
2. 共价机制确证的标准组合：**LC-MS 全蛋白质量位移 → 胰蛋白酶消化肽段 MS/MS 定位 Cys → C116S/C154A 突变 ≥10–100× IC50 右移 → 共晶 Fo−Fc 连续密度**。
3. 亚型选择性不能仅靠单点 IC50；文献共识 readout 为 **kinact/KI（Sox 荧光肽底物）+ 固定孵育时间 IC50 + 细胞 biotin 探针竞争 / NanoBRET + 同源建模 Leu106/Ile106**。
4. JNK-IN-8 是 pan-JNK 工具化合物，**不是** JNK2 选择性参照；但其 Cys116 共价化学与选择性优化逻辑是 YL5084 的直接前体。

---

## 一、文献全景与 JNK2 相关性分级

### 1.1 为何 JNK2 亚型选择性极难

JNK1 与 JNK2 ATP 口袋 **98% 序列同一**，仅 **Met77/Leu77** 与 **Ile106/Leu106** 两处差异；JNK2 与 JNK3 亦仅 **Leu77/Met115** 一处差异 [R1, R8]。传统可逆抑制剂（如 CC-930、SP600125）几乎无法实现 JNK1 vs JNK2 选择性 [R9, R11]。

共价策略利用 **JNK 家族独有的 Cys116（JNK1/2）/ Cys154（JNK3）**（MAPK 激酶组中罕见），通过 **非共价预定位（KI）× 成键速率（kinact）** 两个维度实现亚型差异 [R1, R2, R13]。

### 1.2 本综述纳入标准

| 纳入 | 排除 |
|------|------|
| 明确共价机制（丙烯酰胺/环己烯酮等 Michael 受体） | 纯可逆 pan-JNK（CC-930 等） |
| 有 Cys116/154 修饰或突变确证 | 仅计算/专利、无湿实验 |
| 报告 JNK1/2/3 至少两亚型活性 | 仅 JNK 通路抑制剂、未证共价 |

---

## 二、分篇深度复盘

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

#### 2.3.3 共价与选择性确证（摘要级 + SI）

| 指标 | 数据 |
|------|------|
| kinact/KI (JNK2) | **38,200 M⁻¹s⁻¹** |
| 可逆先导 vs JNK1 | **>100×** |
| 细胞 | 共价化合物保持 JNK2/3 亚型选择性活性 |
| Kinome | 56d **clean profile**（正文摘要） |

> **注：** 全文 IC50 表格、C116S 突变、共晶 PDB 需查阅 SI（Figshare [10.1021/acs.jmedchem.5c00884.s001](https://doi.org/10.1021/acs.jmedchem.5c00884.s001)）。与 YL5084 系列 **独立验证** 了「氨基吡唑/不同 scaffold + Cys116 共价 → JNK2/3 > JNK1」策略。

#### 2.3.4 湿实验（论文声明）

- 三亚型并行生化/细胞测活
- kinact/KI 动力学
- 细胞 isoform 选择性
- Kinome panel
- PK / 代谢稳定性（可逆系列）

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
| **1cR-IN-8** | **JNK3 > JNK2 >10×** |

机制：C4 取代基投射至 **exon 6 编码的底物结合沟**（JNK1/2/3 序列差异区），环己烯酮 3D 形状感知 subtle 差异。

**对本项目的意义：** 证明 **共价 + 弹头侧链** 可实现亚型 tuning；若目标 JNK2> JNK1，需避免 1bR 型 C4 延伸方向，可参考 1aR/1cR SAR。

#### 2.4.5 湿实验清单

PhALC（自研 MAPK 活性 luciferase 互补）；NanoBRET TE + washout；SPR（WT/C116S）；X-ray；kinome scan ×2 平台交叉验证；AP-1 reporter；MKK7 ACT 工程细胞系；PROTAC (PRT_1) 降解 JNK1；大鼠 hepatocyte / plasma PK。

---

### [P5] Muth et al. 2016 — JNK3 吡啶咪唑共价抑制剂（Tübingen，非 JNK2 主线）

**文献：** Tri- and Tetrasubstituted Pyridinylimidazoles as Covalent Inhibitors of c-Jun N-Terminal Kinase 3. *J. Med. Chem.* 60(2):594–607.  
**DOI:** [10.1021/acs.jmedchem.6b01180](https://doi.org/10.1021/acs.jmedchem.6b01180)

#### 2.5.1 发现与优化

- **起点：** p38α 吡啶咪唑先导 → 取代模式调整 **shift 至 JNK3**（疏水区 I 小甲基 + imidazole N1 甲基化）。
- **共价化：** 骨架嫁接 **丙烯酰胺** → 化合物 **7**（四取代咪唑）。

#### 2.5.2 活性与共价确证

| 指标 | 数据 |
|------|------|
| JNK3 IC50 | **0.3 nM** |
| Kinome (410) | Excellent selectivity |
| 共价确证 | WT vs **JNK3-C154A** 孵育 + MS；突变体活性骤降 |
| 代谢 | 人肝微粒体稳定 |

**亚型：** JNK3 选择性（非 JNK2> JNK1）；与 Wydra 2025 同组，代表 Tübingen **JNK3 共价** 技术路线。

---

### [P6] Wen et al. 2024 — JC16I：JNK3 > JNK1/2 共价抑制剂（欧阳亮组）

**文献：** Selective Covalent Inhibiting JNK3 by Small Molecules for Parkinson's Diseases. *Angew. Chem. Int. Ed.* e202411037.  
**DOI:** [10.1002/anie.202411037](https://doi.org/10.1002/anie.202411037)

#### 2.6.1 发现路径

- 结构引导 + 共价药物设计；靶向 JNK 保守 Cys（JNK3 **Cys154**）。
- **JC16I：** 高效 JNK3 共价抑制剂；**JC-P1** 炔烃探针用于 ABPP。

#### 2.6.2 活性与选择性

| 指标 | 数据 |
|------|------|
| JNK3 IC50 | **5.31 nM** |
| vs JNK1/2 | **>160×** 选择性 |
| 细胞 | 低浓度长效 JNK3 抑制；washout 后持久 |
| 突变 | 关键氨基酸（含 Cys154 及口袋残基）→ 选择性机制 |
| ABPP | JC-P1 特异性标记 SH-SY5Y 中 JNK3；proteome-wide 干净 |
| 体内 | PD 模型神经保护 |

**对本项目：** 反向案例 — 同一 Cys 共价策略可实现 **JNK3 >> JNK1/2**；说明亚型选择性来自 **非共价口袋 fit**，非 Cys 本身。

---

## 三、横向对比表

### 3.1 发现策略对比

| 文献 | 策略 | 起点 | 关键转折 |
|------|------|------|---------|
| Zhang 2012 | Structure-based | imatinib type-2 | linker 优化 → acrylamide 几何 + flag methyl |
| Lu 2023 | Lead optimization | JNK-IN-8/11 | THZ-3-60-1 异构体 hint → (R)-aminopyrrolidine |
| Wydra 2025 | **Ligand-first** | aminopyrazole reversible | >100× JNK2/3 可逆 → +acrylamide |
| Tóth 2024 | **Warhead-first** | cyclohexenone | 嫁接到 IN-8；C4 调 isoform |
| Muth 2016 | Scaffold hopping | p38 pyridinylimidazole | → JNK3 + acrylamide |
| Wen 2024 | Structure-guided covalent | 自有 scaffold | JNK3 口袋 multi-residue fit |

### 3.2 共价确证实验矩阵

| 实验类型 | Zhang/JNK-IN-8 | Lu/YL5084 | Wydra/56d | Tóth/1aR-IN-8 | Muth/7 | Wen/JC16I |
|---------|---------------|-----------|-----------|--------------|--------|-----------|
| 全蛋白 LC-MS | ✓ | ✓ | SI | — | ✓ | — |
| 肽段 MS/MS 定位 Cys | ✓ (C116) | ✓ (C116) | SI | — | ✓ (C154) | — |
| Cys→Ser/Ala 突变 | ✓ C116S ≥100× | ✓ (系列) | SI | ✓ C116S SPR | ✓ C154A | ✓ 多点 |
| 共晶共价键密度 | ✓ C154 (JNK3) | ✓ C116 (8ELC) | SI | ✓ C116 (8PTA) | — | — |
| 非共价对照 (R/-H) | ✓ JNK-IN-6 | ✓ YL2056R | — | ✓ IN-8a | — | — |
| GSH 稳定性 | — | t½ 46 min | — | ✓ 10 mM 挑战 | — | — |

### 3.3 亚型选择性 readout 对比

| Readout | 用途 | 最佳应用文献 |
|---------|------|-------------|
| **kinact/KI** | 共价效率 + 非共价预定位综合 | YL5084 (21×), 56d (38200 M⁻¹s⁻¹) |
| **固定时间 IC50** | 简单快速；受 incubation time 影响 | YL5084 (31× @ 1 h) |
| **C116S/C154A 突变** | 共价必要性 | JNK-IN-8 (≥100×), Tóth SPR |
| **Biotin 探针竞争** | 细胞内 JNK1 vs JNK2 占据 | YL5084 (500 nM JNK2 only) |
| **NanoBRET TE** | 活细胞靶向；可 washout | YL5084, 1aR-IN-8 (~10 nM) |
| **KiNativ / KINOMEscan** | 组学选择性 | JNK-IN-8, YL5084, 1aR-IN-8 |
| **同源建模 + MD** | Leu106/Ile106 机制解释 | YL5084 (500 ns) |

### 3.4 活性数据 master 表

| 化合物 | JNK1 | JNK2 | JNK3 | 选择性声明 | 细胞 EC50 / TE |
|--------|------|------|------|-----------|---------------|
| JNK-IN-8 | 4.7 nM | 18.7 nM | 1.0 nM | Pan-JNK; JNK vs kinome | p-c-Jun ~300–500 nM |
| YL2056 | 166 nM | 5 nM | — | 33× (1 h IC50) | — |
| **YL5084** | 2173 nM | 70 nM | 84 nM† | **21× kinact/KI; 31× IC50** | JNK2 ~500 nM TE |
| 56d | — | kinact/KI=38200 | JNK2/3 | **>100× vs JNK1 (reversible lead)** | 细胞 JNK2/3 选择性 |
| 1aR-IN-8 | ~等 | ~等 | ~等 | Pan-JNK; kinome 更 clean | NanoBRET ~10 nM |
| 1bR-IN-8 | 优 | 10× 弱 | — | JNK1 > JNK2 | — |
| 1cR-IN-8 | — | 10× 弱 | 优 | JNK3 > JNK2 | — |
| Muth 7 | — | — | 0.3 nM | JNK3; 410 kinome | — |
| JC16I | 弱 | 弱 | 5.3 nM | **>160× vs JNK1/2** | 长效 cellular JNK3 |

† JNK3 Z'-LYTE IC50；最大抑制低于 JNK2。

---

## 四、对本项目（JNK2 选择性共价筛选）的启示

### 4.1 阳性对照分级

| Tier | 化合物 | 预期 |
|------|--------|------|
| **Tier-0** | YL5084 | kinact/KI JNK2/JNK1 ≥ 10×；C116S 右移 |
| **Tier-1** | YL2056, JNK-IN-8 | 共价参照；JNK-IN-8 为 pan-JNK 阴性选择性对照 |
| **Tier-2** | YL5084R / YL2056R | 非共价对照 |
| **Tier-3** | 1aR-IN-8 | 可逆共价 + GSH 稳定参照 |

### 4.2 最小湿实验 package（复制文献共识）

1. **生化：** 固定时间 IC50 (JNK1/2/3) + **kinact/KI (Sox sensor)**
2. **共价：** LC-MS 全蛋白 + 肽段 MS/MS；**JNK2 C116S** ≥10× 右移
3. **结构：** 8ELC 对接/共晶（若有条件）
4. **细胞：** biotin-JNK-IN-7 竞争 WB 或 NanoBRET (JNK1 vs JNK2 双 construct)
5. **选择性：** KINOMEscan @ 1 μM 或 KiNativ

### 4.3 计算筛选注意

- **主受体：** 8ELC（DFG-in 共价）；**禁用** 3NPC (DFG-out) [R1, R3]
- **不依赖 Δsel 对接排名**（JNK1 非共价项目已证伪）[R17]
- AF3 共价虚筛需 Phase 0 阳性对照校准 [R19]

---

## 五、参考文献（本综述编号）

| 编号 | 文献 | 角色 |
|------|------|------|
| [R1] | Lu et al., *J. Med. Chem.* 2023 — YL5084 | **JNK2> JNK1 核心** |
| [R2] | Zhang et al., *Chem. Biol.* 2012 — JNK-IN-8 | Pan-JNK 共价起点 |
| [R10] | Park et al., *Sci. Rep.* 2015 — aminopyrazole | Wydra 可逆先导基础 |
| [R17] | JNK1 项目报告 v2.8 | Δsel 证伪 |
| [R19] | NLRP3 AF3 难度报告 | AF3 gate 警示 |
| [R20] | Wydra et al., *J. Med. Chem.* 2025 — 56d | JNK2/3 ligand-first 共价 |
| [R21] | Tóth et al., *Nat. Commun.* 2024 | 可逆共价 + isoform tuning |
| [R22] | Muth et al., *J. Med. Chem.* 2016 | JNK3 共价（同组前作） |
| [R23] | Wen et al., *Angew. Chem.* 2024 — JC16I | JNK3 选择性（反向案例） |

完整 DOI 列表见 `REFERENCES.md`。

---

## 六、文献空白与未来方向

1. **严格 JNK2> JNK1 且 kinome-clean 的化学 probe** 尚未达成（YL5084 仍有 p38α 非共价 off-target；MM 增殖不完全依赖 JNK2）。
2. **JNK2 选择性共晶** 仅 8ELC（YL2056），YL5084 本身无独立 PDB entry。
3. **Wydra 56d** 全文/SI 数据（IC50 表、C116S、PDB）待完整入库。
4. **可逆共价 JNK2- selective warhead**（Tóth 型 C4 优化指向 JNK2> JNK1）尚未见专门报道。

---

*本文档为 `JNK2选择性共价抑制剂筛选方案.md` 的文献支撑章节，与 Phase 0 实验设计直接对应。*
