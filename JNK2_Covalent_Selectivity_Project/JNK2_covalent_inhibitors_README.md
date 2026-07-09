# JNK2 抑制剂主表

本目录含 **两张主表**：

| 文件 | 内容 | 条目数 |
|------|------|--------|
| `JNK2_covalent_inhibitors_master_table.csv` | Cys116 **共价** 抑制剂 | 12 |
| `JNK2_reversible_inhibitors_master_table.csv` | **可逆** ATP 竞争 JNK2 相关抑制剂 | 13 |

格式对齐 `NLRP3_covalent_inhibitors_master_table.csv`（分支 `cursor/nlrp3-covalent-master-table-cd8a`）。

---

# 表一：共价抑制剂

数据文件：`JNK2_covalent_inhibitors_master_table.csv`

## 纳入标准

实验明确 **Cys116（JNK1/2）共价** 机制，且满足以下 **至少两项**：

1. 位点突变挽救（C116S → 结合或抑制显著右移）
2. LC-MS / IPMS 全蛋白或肽段 MS/MS 直接鉴定 Cys116 修饰
3. 不可逆性对照（washout、还原弹头 / 丙酰胺对照、Dialysis 等）

**排除：** 纯 JNK3 Cys154 共价（Muth 7、JC16I 等）— 见 `REFERENCES.md` §九。

## 亚型编号说明（重要）

JNK 家族中，**结构同源半胱氨酸** 在不同亚型 UniProt 编号不同。**这不是录入错误。**

| 亚型 | UniProt | 共价位点 | 代表验证 |
|------|---------|---------|---------|
| JNK1 | MAPK8 | **Cys116** | JNK-IN-2 MS 肽段（Zhang 2012） |
| JNK2 | MAPK9 | **Cys116** | YL5084 CE-MS/MS；8ELC 共晶 |
| JNK3 | MAPK10 | **Cys154** | JNK-IN-2 共晶 PDB **3V6R**（2.60 Å）；JNK-IN-7 共晶 PDB **3V6S**（2.97 Å） |

**统一写法：** 表格保留原文 `covalent_site` + `site_species`；跨亚型比较时 `human_ortholog_site` 统一写 **Cys116**（JNK1/2 编号），JNK3 晶体数据在 `notes` 标注 Cys154。

> JNK3 Cys154 反应性高于 JNK1/2 Cys116（CpHMD pKa 6.3 vs 7.5–8.0），但 **亚型选择性主要来自非共价口袋 fit（Leu106/Ile106 等）**，非 Cys 位点本身 [Liu 2022; Lu 2023]。

## 字段说明

| 字段 | 说明 |
|------|------|
| `compound_type` | 合成药物 / 合成药物（阴性对照） |
| `site_species` | 共价位点验证所用蛋白/实验体系 |
| `human_ortholog_site` | 人源 JNK1/2 统一编号 Cys116 |
| `mechanism_route` | **共价表：** `A_structure_based` / `C_isoform_medchem` / `D_ligand_first` / `B_reversible_covalent`；**可逆表：** `E_reversible_DFG_in` / `F_reversible_DFG_out` / `G_reversible_pan_JNK` / `H_reversible_jnk1_bias` / `I_reversible_dual_target` |
| `irreversible` | `是` / `否` / `可逆共价` |
| `af3_calibration_priority` | 8ELC 共价回顾校准推荐优先级 |
| `pubmed_status` | 引用前建议复核；主文献 PMID 已标注 |

## 机制路线与必要验证

| 路线 | 代表化合物 | 必要 readout |
|------|-----------|-------------|
| A 结构导向 acrylamide | JNK-IN-7/8 | C116S ≥100×；kinome S(10) |
| C 亚型 MedChem | YL5084 | **kinact/KI** + biotin 竞争 / NanoBRET JNK1 vs JNK2 |
| D ligand-first | 56d | PhosphoSens 预孵育 IC50 + washout + IPMS |
| B 可逆共价 | 1aR-IN-8 | GSH 10 mM 挑战 + Dialysis 可逆 |

## 命名警示

| 名称 | 说明 |
|------|------|
| **JNK-IN-12 (Zhang 2012)** | 苯并噻唑乙腈共价 pan-JNK；S(10)=0.025 |
| **JNK-IN-12 (Zhou 2023)** | 线粒体靶向 SP600125 偶联物（**非共价 Cys116 TCI**）— **未纳入本表** |

## 统计（2026-07）

- 总条目：**12**（含 1 条阴性对照 JNK-IN-6）
- 严格 **JNK2 > JNK1** 共价：**YL5084**（Tier-0）
- JNK2/3 > JNK1 共价：**YL2056、56d**
- Pan-JNK 共价：**JNK-IN-7/8/11/12、1aR-IN-8**
- 可逆共价：**1aR-IN-8、1bR-IN-8**

## 关键文献 DOI / PMID（2026-07 复核）

| 文献 | DOI | PMID | 期刊卷期 |
|------|-----|------|---------|
| Zhang 2012 | 10.1016/j.chembiol.2011.11.010 | **22284361** | *Chem Biol* 2012;19(1):140-154 |
| Lu 2023 | 10.1021/acs.jmedchem.2c01834 | 36826833 | *J Med Chem* 2023;66(5):3356-3371 |
| Wydra 2025 | 10.1021/acs.jmedchem.5c00884 | 40404564 | *J Med Chem* 2025;68(11):12004-12028 |
| Tóth 2024 | 10.1038/s41467-024-52573-2 | 39366946 | *Nat Commun* 2024;15:8269 |

## 数据核对状态（2026-07-09）

主表数值已对照 **PMC 全文 Table 1 / 正文** 逐项审计。核对来源如下：

| 化合物组 | 核对来源 | 状态 |
|---------|---------|------|
| JNK-IN-2/6/7/8/11/12 | Zhang 2012 Table 1（[PMC3270411](https://pmc.ncbi.nlm.nih.gov/articles/PMC3270411/)） | ✅ 已修正并填入精确 IC50 |
| YL2056 / YL5084 | Lu 2023 Table 1–2 + 正文（[PMC11190964](https://pmc.ncbi.nlm.nih.gov/articles/PMC11190964/)） | ✅ 已确认 |
| 56d / 56b | Wydra 2025 正文 + 综述 [P3]（PMC12169684） | ✅ 已确认 |
| 1aR / 1bR-IN-8 | Tóth 2024 Table 1 + Fig 6 / Supp Fig 14（[PMC11452492](https://pmc.ncbi.nlm.nih.gov/articles/PMC11452492/)） | ✅ 已确认 |

**本次修正要点：**

1. **PDB 编号：** 3V6R = JNK-IN-2（2.60 Å）；3V6S = JNK-IN-7（2.97 Å）——原 JNK-IN-2 行误写「3V6R 系列」、JNK-IN-7 行误绑 3V6R，已更正。
2. **JNK-IN-11 生化 IC50：** 1.34/0.50/0.50 nM（原 1.3/0.5/0.5 为四舍五入）；A375 EC50 **8.6 nM** 与原文一致。
3. **JNK-IN-12：** 补全 Table 1 数值 13/11.3/11.0 nM 及细胞 EC50。
4. **JNK-IN-2：** 补全 Table 1 精确 IC50（809/1140/709 nM），非仅「~1 μM」。
5. **1aR-IN-8：** 补全 PhALC 22±6 nM、NanoBRET 11.5 nM 等 Tóth Table 1 数据。

## 蛋白结构 / 对接信息（逐分子）

本表说明每个共价分子是否有 **实验共晶**、文献或项目内是否有 **分子对接/叠合结果**，以及对接范式是否与本项目 **3ELJ 式非共价 DFG-in 对接** 相同。

### 对接范式说明（读表前必读）

| 范式 | 代表 PDB | 软件 / 方法 | 用途 | 与本表关系 |
|------|---------|------------|------|-----------|
| **A. 3ELJ 式非共价 DFG-in** | 3ELJ、4L7F（JNK1）；3E7O（JNK2 非共价） | Schrödinger Glide SP/XP，**非共价**，Cys116 **不参与** 共价约束 | JNK1 亚型选择性交叉对接；**不作** 共价 hit 主排序 | 仅部分分子可 **借同源模型** 做 JNK1 参照 pose，**≠** 原文发现路径 |
| **B. 8ELC 式共价 DFG-in** | 8ELC（JNK2–YL2056） | Glide **Covalent** / AF3 共价校准；Cys116 为反应位点 | **本项目共价虚筛主受体** | YL2056/YL5084 及衍生物的标准模板 |
| **C. 3V6S 式 JNK3 共价 type-2** | 3V6R（IN-2）、3V6S（IN-7） | 原文分子建模 + 共晶验证；imatinib **type-2** DFG-in | Zhang 2012 系列；56d 设计叠合参照 | JNK-IN-2/7 有共晶；IN-8/11/12 **无** 独立 PDB |
| **D. 8PTA 式 JNK1 可逆共价** | 8PTA/8PT9/8PT8（JNK1–环己烯酮） | X-ray + SPR；**非** Glide 共价虚筛主线 | Tóth 2024 可逆共价 warhead | 与 8ELC **不同 isoform + 不同 warhead** |
| **E. 叠合设计（无独立共晶）** | 4WHZ + 3V6S + 可逆 21b 系列 | 手工 / PyMOL 叠合 → linker 可达 Cys116 | Wydra 56d ligand-first | **非** 3ELJ 式盲对接 |

> **3ELJ 式** = 标准非共价激酶对接（hinge + ATP 口袋，DFG-in），**不包含** Cys116 共价键约束。本仓库 `config/targets.yaml` 规定：**共价筛选主受体为 8ELC**；3ELJ/4L7F 仅作 JNK1 交叉参照；**禁用 3NPC（DFG-out）** 作共价对接。

### 逐分子结构 / 对接索引

| ID | 化合物 | 实验共晶 | 共晶 PDB（配体） | 文献对接 / 建模 | 是否 3ELJ 式 | 项目推荐对接模板 |
|----|--------|---------|-----------------|----------------|-------------|-----------------|
| 1 | JNK-IN-2 | ✅ JNK3 | **3V6R**（2.60 Å，IN-2） | Zhang 2012：JNK3 **分子建模**（type-2 U 构象）→ 共晶验证；**无** Glide/3ELJ 报告 | ❌ | 3V6R（JNK3）；JNK2 无共晶 → 可叠合至 8ELC/3V6S |
| 2 | JNK-IN-7 | ✅ JNK3 | **3V6S**（2.97 Å，IN-7） | 同上；linker 优化基于 IN-2 共晶/模型 | ❌ | **3V6S**（56d 叠合锚点之一） |
| 3 | JNK-IN-8 | ❌ 无独立 PDB | —（与 IN-7 同系列，未单独沉积） | 继承 IN-7 系列结构假设；**无** 共晶 pose | ❌ | 叠合 **3V6S** 或 **8ELC**（IN-8 骨架） |
| 4 | JNK-IN-11 | ❌ | — | MedChem 末端替换；**无** 共晶/对接论文 | ❌ | 叠合 3V6S/8ELC（YL5084 母核来源） |
| 5 | JNK-IN-12 (Zhang) | ❌ | — | 同系列；**无** 共晶 | ❌ | 叠合 3V6S |
| 6 | JNK-IN-6 | ❌ | — | 丙酰胺阴性对照；**无** 结构 | ❌ | 不适用（机制对照） |
| 7 | YL2056 | ✅ **JNK2** | **8ELC**（2.0 Å，共价 Cys116）；**7N8T**（JNK2–AMP，对照） | Lu 2023：**Prime 同源建模**（JNK1/JNK3←8ELC）+ **Glide SP** 对接 YL5084；**500 ns MD** | ❌（共价 DFG-in） | **8ELC（主模板 ★★★）** |
| 8 | YL5084 | ❌ 无独立 PDB | —（pose 与 YL2056 同系） | Lu 2023：Glide 对接 **JNK2:YL2056** 晶体结构；homology JNK1/JNK3 | ❌（共价 DFG-in） | **8ELC** 共价对接（Tier-0 校准） |
| 9 | 56d | ❌ | — | Wydra 2025：**叠合设计**（16a + **3V6S** + **4WHZ**）→ meta-aminobenzamide + acrylamide；**非** 盲对接 | ❌（叠合 E） | 4WHZ + 3V6S 叠合；共价验证用 **8ELC** Cys116 几何 |
| 10 | 56b | ❌ | — | 同 56 系列；para-linker 几何 SAR | ❌（叠合 E） | 同 56d |
| 11 | 1aR-IN-8 | ✅ **JNK1** | **8PTA**（1aR）；**8PT9**（1aS）；**8PT8**（1a'R） | Tóth 2024：X-ray + SPR；引用 **3V6S** 说明 IN-8 ATP 结合；**无** 3ELJ/Glide 共价虚筛 | ❌（JNK1 可逆共价 D） | **8PTA**（JNK1）；跨 JNK2 须同源至 8ELC |
| 12 | 1bR-IN-8 | ❌ 无独立 PDB | —（C4 延伸基于 **1aR-IN-8** 共晶构象） | Tóth 2024：自 **8PTA** 底物沟延伸设计；PhALC isoform panel | ❌ | 叠合 **8PTA** / 8ELC（IN-8 骨架 + C4 矢量） |

### 快速判读

- **有 JNK2 共价共晶、可直接作共价对接受体：** 仅 **YL2056（8ELC）**。
- **有共晶但为 JNK3（Cys154）或 JNK1：** JNK-IN-2/7（3V6R/3V6S）、1aR 系列（8PTA 等）— 跨亚型使用时需同源建模或叠合。
- **原文有 Glide/MD 对接记录：** YL2056、YL5084（Lu 2023，基于 8ELC，**非 3ELJ 式**）。
- **原文为叠合/建模、无共晶：** JNK-IN-8/11/12、YL5084（无 PDB）、56d/56b、1bR-IN-8。
- **3ELJ 式对接：** 本表 **无** 分子以 3ELJ 式作为原文发现方法；3ELJ 仅在本项目 `JNK2选择性共价抑制剂筛选方案.md` 中作 **JNK1 非共价交叉参照**，**不用于共价 hit 排序**。

主表 CSV 的 `notes` 字段已追加 `[结构]` 标签，与本节一一对应。

### 结构选择性总框架（读表前必读）

JNK1/JNK2 ATP 口袋 **98% 同一**，仅两处差异氨基酸驱动可设计的亚型选择性：

| 位点 | JNK1 | JNK2 | JNK3 | 结构角色 |
|------|------|------|------|---------|
| **106**（HR-I 后口袋） | **Ile106** | **Leu106** | **Leu106** | 芳环/scaffold 后口袋体积；JNK1 Ile 导致 **steric clash** → JNK2/3 优先 |
| **77**（P-loop 附近） | Met77 | Leu77 | Met115 | 与 P-loop 动态、主链柔性相关（Lu MD） |
| **54**（主链） | Ile54 | Val54 | Ile54 | JNK2 Val54 主链可移动 ~0.6 Å 容纳配体（8ELC） |
| **116** | Cys116 | Cys116 | Cys154 | **共价成键位点**；三亚型均保守 → **提供 potency，不自带 JNK2> JNK1** |
| **exon6 底物沟** | 序列差异 | 序列差异 | 序列差异 | Tóth C4 延伸可「编程」亚型偏好（1bR 路径） |

**选择性两层逻辑：**

1. **JNK 家族 vs 其他激酶（kinome 选择性）：** flag methyl（JNK-IN-8/YL5084）、末端芳环形状（JNK-IN-12 苯并噻唑）、环己烯酮 3D 形状（1aR）— 与 imatinib 经验类似，**不涉及 JNK2> JNK1**。
2. **JNK1 vs JNK2 vs JNK3（亚型选择性）：** 共价 Cys116 提供锚定后，**非共价预定位（KI / back pocket fit）** 决定亚型；readout 用 **kinact/KI** 或固定时间 IC50，**非** 3ELJ 式 ΔΔG 排名（JNK1 非共价项目已证伪）[筛选方案 §1.3]。

**为何 3ELJ 式对接不能解释/预测 JNK2 亚型选择性：**

- 3ELJ 为 **非共价** DFG-in，不约束 Cys116 共价键，**不含 kinact/TDI** 维度；
- 多数本表分子原文 **未用** 3ELJ 发现；YL/56d 选择性来自 **8ELC 共价 pose** 或 **4WHZ 可逆 fit + 叠合**；
- JNK1 非共价 Δsel/MM-GBSA 对 benchmark 选择性方向 **无判别力**（通过率 43%）。

---

### 逐分子：对接/结构如何解释选择性

| ID | 化合物 | 亚型选择性（实验） | 结构/对接解释 | 能否用 3ELJ 预测？ |
|----|--------|------------------|--------------|-------------------|
| 1 | JNK-IN-2 | **无**；pan-JNK ~μM | 3V6R：type-2 U 构象下 acrylamide 与 Cys154 **几何次优** → 低效共价（~1 μM）；**未利用** Leu106/Ile106 差异；末端无 flag methyl → kinome 宽 | ❌ 未设计亚型选择性 |
| 2 | JNK-IN-7 | **无**；pan-JNK nM | 3V6S：优化 **linker 酰胺几何** 使 acrylamide 逼近 Cys154 → **potency↑**（~500× vs IN-2），三亚型 IC50 近等；选择性仍靠 **kinome 面板** 而非口袋差异 | ❌ |
| 3 | JNK-IN-8 | **JNK vs kinome ✅**；**JNK2≈JNK1≈JNK3 ❌** | 叠合 3V6S/8ELC：**flag methyl** 占据 imatinib 式后口袋 → 剔除 IRAK1/PIK3C3 等（S(10)=0.031）；IN-8 骨架 **不区分** Ile106/Leu106 → 生化 JNK2 18.7 nM 甚至弱于 JNK1 4.67 nM | ❌ 3ELJ 看不到 flag methyl 对 kinome 的贡献 |
| 4 | JNK-IN-11 | **无亚型**；kinome **差** | 大体积 **2-苯基吡唑并[1,5-a]吡啶** 填充满 ATP+后口袋 → pan-JNK nM + **p38/CK1** 等 off-target；C116S 仅 ~10× → 强抑制 **部分靠非共价** 结合，共价非唯一决定因素 | ❌ THZ-3-60-1 区域异构体曾 hint JNK2 趋势但 kinome 失败 |
| 5 | JNK-IN-12 | **无亚型**；**kinome 最优** | 苯并噻唑乙腈末端 **形状互补** JNK 后口袋 → S(10)=0.025（系列最 clean）；仍 **pan-JNK** 生化（13/11.3/11.0 nM 近等） | ❌ |
| 6 | JNK-IN-6 | 共价 **阴性** | 丙酰胺 **无法** Michael 加成 Cys116 → ~100× 失活；证明 Zhang 系列 potency 依赖共价，但 **不提供** 亚型信息 | — |
| 7 | YL2056 | **JNK2/3 > JNK1 ~33×** | **8ELC 共晶直接证据：** (R)-3-氨基吡咯烷+芳环深入 **Leu106 后口袋**（JNK2）；JNK1 **Ile106 clash**；Val54 主链位移容纳配体；**Cys116 共价** 提 potency，**(R)-构型** 必需（YL2012 对映体选择性消失）。Glide/MD：JNK1 ΔG_binding 升高 ~2.5 kcal/mol | ❌ 须 **8ELC 共价 pose**；3ELJ 无共价+无 (R)-吡咯烷 fit |
| 8 | YL5084 | **JNK2 > JNK1 ~21×（kinact/KI）** | 继承 YL2056 **Leu106/Ile106** 逻辑 + **flag methyl** 降 kinome；Glide→8ELC：JNK2 中 warhead **更优轨迹** 指向 Cys116；MD：**P-loop 更接近** 配体；JNK1 **Arg50–Glu109** 盐桥限制 P-loop → 共价效率低 → kinact/KI(JNK2)/kinact/KI(JNK1)≈21。JNK3 仍活性（84 nM）因 Leu106 共享 | ❌ 3ELJ 缺共价 TDI；Δsel 不可用 |
| 9 | 56d | **JNK2/3 >> JNK1**；strict JNK2> JNK1 弱于 YL5084 | **Ligand-first：** 可逆 **21b/16a** 先优化 **氨基吡唑 + Leu106/Leu144** 后口袋 fit（4WHZ/Park 2015 机制）→ JNK1 无 reversible IC50；再叠合 **3V6S** 引入 meta-**acrylamide** 可达 Cys116 → **预定位 + 共价** → kinact/KI(JNK2)=38200。JNK1：无 fit → PhosphoSens 无活性 | ❌ 叠合设计，非 3ELJ 盲对接 |
| 10 | 56b | **无有效选择性**（>T10 μM） | 同 56d scaffold 但 **para-** 取代 linker → 叠合/model 显示 acrylamide **偏离** Cys116 最优 Michael 角度 → 共价标记成功但 **TDI 效率极低**；说明 **warhead 矢量** 与选择性同等关键 | ❌ meta vs para 几何 lesson |
| 11 | 1aR-IN-8 | **Pan-JNK 三亚型近等**；**kinome > JNK-IN-8** | 8PTA：**IN-8 ATP 口袋**（hinge H-bond，参照 3V6S）+ **环己烯酮** 可逆锁定 Cys116；刚性环己烯酮 **空间位阻** 降低 off-target 半胱氨酸反应；**未延伸** 至 Leu106 后口袋差异区 → 不区分 JNK1/2 | ❌ JNK1 共晶，跨 JNK2 需同源至 8ELC |
| 12 | 1bR-IN-8 | **JNK1 > JNK2 ~10×（反向）** | 自 8PTA：**C4 丙炔酯** 投射入 **exon6 底物结合沟**（GGV 等残基三亚型差异）；环己烯酮 3D 形状 **感知** 沟槽 subtle 差异 → **可编程** 亚型偏好；与 Leu106 轴 **正交** — 设计 JNK2> JNK1 应 **避免** 此 C4 方向 | ❌ 底物沟非 3ELJ 标准 ATP 对接视野 |

### 选择性机制速查（按设计策略）

```
                    ┌─ kinome 选择性（JNK vs 其他激酶）
                    │   → flag methyl / 末端形状 / 环己烯酮 bulk
                    │
共价 Cys116 锚定 ───┼─ 亚型选择性 JNK2 vs JNK1
                    │   → Leu106(Ile106) 后口袋 + P-loop 动态 [YL5084, 56d]
                    │   → 须 kinact/KI 或 预孵育 IC50，非 Δsel
                    │
                    └─ 底物沟编程 [1bR] / 反向 JNK1> JNK2
                        → exon6 差异，与 Leu106 轴独立
```

主表 CSV `notes` 已追加 `[选择性]` 标签，与本节对应。

---

# 表二：可逆 JNK2 抑制剂

数据文件：`JNK2_reversible_inhibitors_master_table.csv`

## 纳入标准

**ATP 竞争性可逆** 小分子，且满足以下 **至少一项**：

1. 报告 **JNK2 生化/细胞** 活性或 **JNK2/3 > JNK1** 亚型选择性数据
2. 作为 **pan-JNK 可逆阴性对照**（CC-930、SP600125）或 **反向 JNK1> JNK2 参照**（CC-90001）
3. 作为 **DFG-out / Type II** 结构-选择性经典案例（BIRB796 + 3NPC）

**排除：** 纯 JNK3 选择性可逆（indazole 25c 等）；JNK3 共价文献。

## 用户常见问题：哪个化合物用了 DFG-out（out 蛋白）对接？

**是的——经典案例是 BIRB796（doramapimod）+ PDB 3NPC（JNK2 DFG-out 共晶）。**

| 维度 | BIRB796 / 3NPC（DFG-out） | 21b / 51d / 26k（DFG-in） |
|------|---------------------------|---------------------------|
| 构象 | **Type II DFG-out**；Phe170 翻出，扩展后口袋 | **Type I DFG-in**；标准 ATP 口袋 |
| 代表 PDB | **3NPC**（JNK2–BIRB796） | **4WHZ**（JNK3–26k 氨基吡唑） |
| 对接文献 | MedChemComm 2016：AutoDock Vina 对接 **3NPC** + JNK1/JNK3 DFG-out **同源模型** | Park 2015 / Zheng 2014 / Wydra 2025：**叠合 + SAR**，基于 4WHZ 共晶 |
| 选择性机制 | DFG-out 扩展口袋 **仅 JNK2 实验解析**；BIRB796 在 JNK2 中 fit 良好；JNK1 **Ile106+Met77** / JNK3 **Leu144+Met115** 使 DFG-out 口袋更小 → Vina 无法合理 pose → JNK2 IC50 ~6 nM vs JNK1 >10 μM | **Leu106(Ile106)** HR-I 后口袋：更大芳环/scaffold 在 JNK2/3 fit，JNK1 **steric clash**；51d **倒置酰胺** 完全消除 JNK1 |
| 与本项目关系 | **共价筛选禁用** 3NPC（与 8ELC DFG-in 几何不一致） | 56d 可逆前体；Tier-1 可逆对照 |

> **注意：** Wydra **21b/51d** 等 JNK2/3 选择性氨基吡唑走的是 **DFG-in（4WHZ）** 路线，**不是** 3NPC DFG-out。用户若记得「out 蛋白对接」，通常指 **BIRB796/3NPC** 或 MedChemComm 2016 同源建模对接研究。

### 可逆对接范式（读表前必读）

| 范式 | 代表 PDB | 代表化合物 | 选择性 readout |
|------|---------|-----------|---------------|
| **F. DFG-out Type II** | **3NPC** | BIRB796 | JNK2 >> JNK1/JNK3；扩展口袋体积 |
| **E. DFG-in 氨基吡唑** | **4WHZ** | 26k, 21b, 51d | Leu106/Ile106 后口袋 fit；148–340× vs JNK1 |
| **G. Pan-JNK 可逆** | 3TTI 等 | CC-930, SP600125 | 无亚型选择性（阴性对照） |
| **H. 反向 JNK1 偏好** | — | CC-90001 | 生化 Ki 近等；**细胞 KO 12.9× JNK1** |
| **I. 激酶 + PPI 双靶** | VS 模型 | 6l | JNK2 功能导向；非经典 Leu106 轴 |

### 逐分子：可逆选择性机制摘要

| ID | 化合物 | 亚型选择性 | 如何实现选择性（结构/对接） |
|----|--------|-----------|---------------------------|
| 1 | BIRB796 | JNK2 >> JNK1/JNK3 | **3NPC DFG-out**：tolyl + morpholino-naphthyl 填满 Type II 口袋；JNK1/JNK3 DFG-out 同源模型口袋更小（Ile106/Met77, Leu144/Met115）→ 对接失败/无结合 |
| 2 | 26k | JNK2/3 > JNK1 | **4WHZ DFG-in**：Leu144 后口袋容受芳环；hinge H-bond 锚定 |
| 3 | 26n | JNK3 导向 | 同 Leu144 轴 + kinome 优化 |
| 4 | A-1 | ~30× vs JNK1 | Park 2015：氨基吡唑 + Leu106 后口袋；Wydra 化学起点 |
| 5 | 16a | 114× | 倒置酰胺 + naphthyl → Leu106 fit |
| 6 | **21b** | **148×** | Ligand-first：并行 JNK1/2/3 SAR → HR-II 苯甲酰胺深度占 Leu106 后口袋 |
| 7 | 21h | 83× | 同 21 系列；PK 优化（t½ 3.33 h） |
| 8 | 45a | JNK2/3 双强 (4/6 nM) | 可逆 potency 上限；选择性 ~12× |
| 9 | **51d** | **>340×** | **倒置酰胺** 重定向 H-bond → JNK1 **完全无结合** |
| 10 | SP600125 | Pan-JNK | 无亚型设计（40/40/90 nM） |
| 11 | CC-930 | Pan-JNK | 临床 pan-JNK；无 Ile106/Leu106 利用 |
| 12 | CC-90001 | **JNK1 > JNK2 12.9×（细胞）** | 纤维化生物学 + SAR；生化 Ki 近等 → 须 **isoform KO** readout |
| 13 | 6l | JNK2 功能导向 | JNK2 激酶 + MKK7–JNK2 PPI 双靶；与共价 Leu106 轴正交 |

### 可逆表统计（2026-07）

- 总条目：**13**
- **JNK2/3 > JNK1 可逆：** 21b（148×）、51d（>340×）、16a、26k 系列
- **DFG-out 结构案例：** BIRB796（3NPC）
- **Pan-JNK 阴性对照：** SP600125、CC-930
- **反向参照：** CC-90001（JNK1 12.9×）
- **双靶新策略：** 6l（2026）

可逆表 CSV `notes` 已追加 `[结构]` / `[选择性]` 标签。

## 关联文档

- `JNK2共价抑制剂文献调研综述.md`（v2.0 逐篇复盘）
- `JNK2选择性共价抑制剂筛选方案.md`（Phase 0 gate；附录 B PDB 索引）
- `REFERENCES.md`（R2–R3 PDB DOI）
- `config/targets.yaml`（8ELC / 3ELJ 受体配置）
