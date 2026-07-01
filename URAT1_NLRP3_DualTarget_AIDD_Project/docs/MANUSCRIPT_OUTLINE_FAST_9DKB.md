# 快速发表线论文大纲 — URAT1 抑制剂 @ 9DKB（论文 A）

> **双靶叙事 + JMGM/非OA 选刊**：`DUAL_TARGET_AND_FAST_JOURNALS.md`  
> **目标期刊**：*J. Comput.-Aided Mol. Des.*（首选）或 *Molecular Diversity*  
> **类型**：Computational characterization / Case study（**非** hit discovery，**非** 新算法）  
> **与论文 B（JCIM）关系**：仅覆盖 **9DKB 单态**；三态失败与 decoy benchmark 全部留给论文 B

---

## 标题（推荐）

**Computational characterization of clinical URAT1 inhibitors against the inward-open cryo-EM structure 9DKB: Glide docking, redocking validation, and molecular dynamics**

备选：
- *Binding mode analysis of lesinurad-class URAT1 inhibitors using the Suo 2025 inward-open structure*
- *Structure-based assessment of marketed URAT1 inhibitors: a docking and MD study on PDB 9DKB*

---

## Abstract（~200 words）

1. **Background**：URAT1 抑制剂；2025 年 9DKB inward-open 高分辨结构可用。  
2. **Objective**：在 9DKB 上 **表征** 四种临床/临床阶段抑制剂结合模式与动力学稳定性（**不是** 筛选新化合物）。  
3. **Methods**：Schrödinger Glide SP→XP；lesinurad redock；100 ns MD；MM-GBSA。  
4. **Results**（据实填）：redock RMSD = X Å；四药 GlideScore 范围；关键相互作用；MD 稳定性排序。  
5. **Conclusion**：9DKB 适用于 inward 抑制剂结合模式研究；转运体柔性限制单结构外推——**不声称** occluded/outward 对接。

---

## 1. Introduction（~900 words）

- 痛风 / URAT1 机制（1 段）  
- Alternating access 与 inward 抑制剂假说（1 段）  
- **9DKB vs 9B1H / 9JDZ** 选型理由（引用 Suo 2025；说明 9JDZ 非 occluded/outward）  
- 四药临床背景（lesinurad、benzbromarone、verinurad、dotinurad）  
- 本文范围：**已知抑制剂在 9DKB 的计算表征**

---

## 2. Methods（~1500 words）

### 2.1 Structure preparation
- PDB 9DKB；Protein Prep Wizard；去配体；pH 7.4

### 2.2 Grid and Glide
- Grid 中心：lesinurad 共晶；box 22 Å  
- SP → XP；记录 GlideScore

### 2.3 Redocking validation
- lesinurad 回 dock；RMSD 阈值 2.0 Å

### 2.4 Test set
- 四药 SMILES（`teacher_gate_qc_panel_b_direction.csv`）  
- **不写** ChEMBL 822 训练集

### 2.5 Molecular dynamics
- 力场、水模型、膜环境（若简化则明确写真空/implicit）  
- 100 ns；RMSD、RMSF

### 2.6 MM-GBSA
- 轨迹提取；结合自由能

### 2.7 （可选）ADMET
- QikProp / SwissADME

---

## 3. Results（~1500 words）

### 3.1 Redocking
- **Fig 1**：lesinurad 晶体 vs docked（RMSD 数值）

### 3.2 Docking poses of four inhibitors
- **Fig 2**：口袋 overlay  
- **Table 1**：GlideScore、关键 H-bond/疏水接触

### 3.3 Comparison with co-crystal structures（异源 PDB）
- 9DKA（benzbromarone）、9JDY（verinurad）、9JE1（dotinurad）—— **cross-structure 比较，非 redock**

### 3.4 MD stability
- **Fig 3**：复合物 RMSD  
- **Fig 4**：代表帧相互作用

### 3.5 MM-GBSA ranking
- **Table 2**：ΔG_bind 相对排序 vs 实验 IC50 趋势（**相关性弱则诚实写**）

---

## 4. Discussion（~800 words）

- 四药结合模式异同  
- 9DKB 作为 inward 对接结构的适用边界  
- **1 段局限**：未研究 occluded/outward（9B1K/9B1L 刚性对接对药物不可行——详见进行中工作 / 不展开）  
- 与激酶单结构筛选对比  
- **不写** ML、双靶、新 hit

---

## 5. Conclusions（~150 words）

---

## 图表清单

| ID | 内容 |
|----|------|
| Fig 1 | lesinurad redock |
| Fig 2 | 四药 9DKB 口袋 overlay |
| Fig 3 | MD RMSD（四复合物） |
| Fig 4 | MM-GBSA 或相互作用代表帧 |
| Table 1 | 对接分数与相互作用摘要 |
| Table 2 | MM-GBSA + 实验 IC50 对照 |
| Fig S1 | 9DKB vs 9B1H 叠合（可选） |
| Table S1 | ADMET（可选） |

---

## 投稿前 Checklist

- [ ] lesinurad redock RMSD ≤ 2.0 Å  
- [ ] 四药 pose 经人工目检（无明显 clash）  
- [ ] MD 平衡后至少 50 ns 用于分析  
- [ ] MM-GBSA 方法段可复现  
- [ ] 摘要无 discovery / novel algorithm 措辞  
- [ ] 与论文 B 图表 **零重复**（B 才放 B1K/B1L 失败表）  
- [ ] 未引用 URAT1 ML 2/4 作为支持证据

---

## 预计计算剩余量

| 任务 | 估计 |
|------|------|
| Redock 分析整理 | 0.5 天 |
| MD 4 × 100 ns | 1–2 周（取决于硬件） |
| MM-GBSA | 2–3 天 |
| 作图 + 初稿 | 2–3 周 |

**总计**：主要瓶颈是 **MD**，不是对接。
