# JNK2 共价抑制剂主表

本目录数据文件：`JNK2_covalent_inhibitors_master_table.csv`

格式对齐 `NLRP3_covalent_inhibitors_master_table.csv`（分支 `cursor/nlrp3-covalent-master-table-cd8a`）。

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
| `mechanism_route` | `A_structure_based`（Zhang 咪唑啉骨架）/ `C_isoform_medchem`（YL 系列）/ `D_ligand_first`（56d）/ `B_reversible_covalent`（Tóth 环己烯酮） |
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

## 关联文档

- `JNK2共价抑制剂文献调研综述.md`（v2.0 逐篇复盘）
- `JNK2选择性共价抑制剂筛选方案.md`（Phase 0 gate；附录 B PDB 索引）
- `REFERENCES.md`（R2–R3 PDB DOI）
- `config/targets.yaml`（8ELC / 3ELJ 受体配置）
