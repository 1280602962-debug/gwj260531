# JNK2 共价选择性抑制剂 — 图文报告

> 本报告从项目参考文献与公开结构库整理示意图，供 GitHub 浏览与汇报使用。  
> 文献编号对应根目录 [`REFERENCES.md`](../../REFERENCES.md)。

---

## 1. 项目决策五阶段总览

![五阶段决策流程](figures/fig01_five_stages_flow.png)

**图 1.** 从靶点与结构选择 → 共价种子与诱饵 → AF3 置信度门控 → 共价对接与 MMGBSA → 选择性验证与候选推进。  
**依据：** 项目 [`JNK2项目决策五阶段.md`](../../JNK2项目决策五阶段.md)；AF3 门控参考 COValid [R7b]。

---

## 2. 关键共价/对照化合物（2D 结构）

![关键化合物 2D 结构](figures/fig02_key_compounds_2d.png)

**图 2.** 种子与对照分子的 2D 结构（RDKit 自项目 SMILES 绘制）。

| 分子 | 角色 | 主要文献 |
|------|------|----------|
| YL5084 | JNK2 共价 hit（丙烯酰胺） | Wydra 2014 [R1] |
| YL2056 | 共价 hit / 56d 前体 | Wydra 2014 [R1] |
| JNK-IN-8 | 共价 hit（Chloroacetamide） | Ward 2012 [R3] |
| 56d | 共价 hit（Chloroacetamide） | Wydra 2014 [R1] |
| JNK-IN-6 | 共价 hit（Chloroacetamide） | Ward 2012 [R3] |
| 26k | 可逆 DFG-in 选择性参照 | 4WHZ [R7b] |
| BIRB796 | 可逆 DFG-out 参照（非共价） | 3NPC [R7] |

SMILES 来源：[`data/phase0_af3/phase0_compounds_seed.csv`](../../data/phase0_af3/phase0_compounds_seed.csv)、[`COMPOUNDS.md`](../../data/phase0_af3/COMPOUNDS.md)。

---

## 3. 为何共价初筛不用 DFG-out（3NPC）

![DFG-in vs DFG-out 示意](figures/fig03_dfg_in_vs_out_schematic.png)

**图 3.** 共价初筛以 **DFG-in / 8ELC** 为主：Cys116 在 P-loop 内可及；**3NPC（DFG-out）** 中 Cys116 被 Phe169 等遮挡，且共价 hit 晶体多来自 DFG-in 系列 [R1][R3][R7]。

| 结构 | 构象 | 共价初筛 | 说明 |
|------|------|----------|------|
| **8ELC** | DFG-in | ✅ 主模板 | YL5084–JNK2 共晶 [R1] |
| **4WHZ** | DFG-in | ⚪ 可逆选择性参照 | JNK3–26k 氨基吡唑 [R7b] |
| **3NPC** | DFG-out | ❌ 排除共价初筛 | BIRB796 Type II；与 hit 构象不一致 |

---

## 4. JNK2 vs JNK1 选择性（Leu106）

![Leu106 选择性示意](figures/fig04_leu106_selectivity.png)

**图 4.** JNK2 **Leu106** 对应 JNK1 **Met146**；Ward 2012 共价 hit 在 JNK2 上选择性显著高于 JNK1 [R3]。项目优先 **JNK2 共价 + DFG-in** 路线。

---

## 5. 参考文献 PDB 结构面板

![PDB 结构面板](figures/fig05_pdb_structure_panel.png)

**图 5.** 自 RCSB 下载的组装图（assembly 1）。点击下图可跳转 RCSB 条目页。

| 面板 | PDB | 配体/说明 | 文献 |
|------|-----|-----------|------|
| A | [8ELC](https://www.rcsb.org/structure/8ELC) | YL2056，DFG-in 共价 | [R3] |
| B | [3NPC](https://www.rcsb.org/structure/3NPC) | BIRB796，DFG-out 可逆 | [R7] |
| C | [4WHZ](https://www.rcsb.org/structure/4WHZ) | 26k，DFG-in 选择性 | [R7b] |
| D | [3V6S](https://www.rcsb.org/structure/3V6S) | JNK-IN-7（JNK3 共晶） | [R2] |
| E | [7N8T](https://www.rcsb.org/structure/7N8T) | JNK2–AMP，1.6 Å | [R4] |

原始 JPEG：`figures/*_assembly.jpeg`（来源 [RCSB PDB](https://www.rcsb.org/)）。

---

## 6. Phase 0 AF3 数据包（输入侧）

| 内容 | 路径 |
|------|------|
| 种子/诱饵 SMILES | `data/phase0_af3/phase0_compounds_seed.csv` |
| 分子角色说明 | `data/phase0_af3/COMPOUNDS.md` |
| 结构审计 | `data/phase0_af3/STRUCTURE_AUDIT.md` |
| AF3 FASTA / 模板 | `data/phase0_af3/sequences/`、`templates/` |

AF3 **置信度判断**需在跑完结构后合并 `summary_confidences.json`，计算 mPAE 与 EF@1%（见 [`JNK2项目决策五阶段.md`](../../JNK2项目决策五阶段.md) 阶段 3）。

---

## 7. 核心参考文献（本报告引用）

| ID | 内容 |
|----|------|
| [R1] | Wydra et al., *ACS Chem. Biol.* 2014 — YL5084 / YL2056 / 56d，8ELC |
| [R2] | Zhang et al., *Chem. Biol.* 2012 — JNK-IN-8 系列，3V6S |
| [R3] | Lu et al. / PDB 8ELC — JNK2–YL2056 共价共晶 |
| [R4] | PDB 7N8T — JNK2–AMP |
| [R7] | BIRB796 / 3NPC — DFG-out Type II |
| [R7b] | PDB 4WHZ — 26k DFG-in 选择性 |
| [R7b] | COValid — AF3 共价对接置信度与 EF@1% 门控 |

完整列表：[`REFERENCES.md`](../../REFERENCES.md)。

---

## 图件清单

| 文件 | 说明 |
|------|------|
| `fig01_five_stages_flow.png` | 五阶段流程图 |
| `fig02_key_compounds_2d.png` | 关键分子 2D |
| `fig03_dfg_in_vs_out_schematic.png` | DFG-in/out 与 Cys116 |
| `fig04_leu106_selectivity.png` | JNK1/JNK2 选择性 |
| `fig05_pdb_structure_panel.png` | PDB 组装图面板 |
| `*_assembly.jpeg` | RCSB 原图 |

*生成日期：2026-07-10 · 项目分支 `cursor/jnk2-covalent-master-table-26c6`*
