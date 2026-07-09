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

## 关联文档

- `JNK2共价抑制剂文献调研综述.md`（v2.0 逐篇复盘）
- `JNK2选择性共价抑制剂筛选方案.md`（Phase 0 gate）
- `REFERENCES.md`
