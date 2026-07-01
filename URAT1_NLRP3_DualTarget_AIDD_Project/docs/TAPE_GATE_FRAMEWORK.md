# TAPE-GATE 框架总览

> **⚠️ 已归档（2026-07）**：旧版双路径融合路线，**请勿再执行**。当前流程：[`WORKFLOW_CURRENT.md`](WORKFLOW_CURRENT.md)

> **T**ransporter-**A**ware **P**aired-path **E**vidence fusion with **G**enerative **A**nd library screening for dual-**T**arget **E**valuation  
> **版本**: 2.0 | **更新**: 2026-06 | **前身**: STAD-AIDD v1.0

---

## 一、框架定位

TAPE-GATE 是面向 **URAT1 + NLRP3** 双靶、**纯计算**、**小数据不对称** 场景的 AI 辅助药物发现框架。相对 STAD-AIDD v1.0 与 PLK1/NLRP3 类工作，v2.0 的核心升级：

1. **双路径候选发现**：库筛（Path A）+ 生成式优化（Path B）并行，统一融合排序
2. **不对称双证据建模**：URAT1 用转运体感知 QSAR + $S_{\text{trap}}$；NLRP3 用 **assay-conditioned** 分类 + 结构证据（非锚点相似性）
3. **可靠性加权融合 + Pareto 排序**：替代固定 0.5/0.5 线性加权
4. **独立双靶模型优先**：ChEMBL 实测 0 重叠 SMILES → MTL 降为可选，分靶训练 + 证据融合为主路径

差异化对照见 [`DIFFERENTIATION_VS_PLK1_NLRP3.md`](DIFFERENTIATION_VS_PLK1_NLRP3.md)。

---

## 二、架构图

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         TAPE-GATE Pipeline                                │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 0  数据层                                                          │
│    URAT1 (822) + NLRP3 (513, IL-1β+B) + SLC22 辅助 + 专利扩充             │
│    Murcko 骨架 GroupKFold + assay 元数据保留（NLRP3）                      │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 1  不对称双证据建模层                                               │
│    URAT1 臂: MiniMol/Chemprop 回归 + Conformal UQ + SLC22 迁移            │
│    NLRP3 臂: Assay-conditioned CLAMP/TwinBooster 分类 + 结构先验          │
│    （可选）MTL 共享编码器 — 仅作消融对照                                   │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 2  双路径候选生成                                                   │
│    Path A 库筛: Enamine/ChEMBL → ML/UQ 过滤 → 系综对接                     │
│    Path B 生成: CLM cross-fine-tune → RL 双靶奖励 → 有效性/新颖性过滤      │
│    合并候选池 C_union                                                       │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 3  结构约束层（转运体核心）                                           │
│    URAT1: 多构象 $S_{\text{trap}}$ + Arg477 关键接触                        │
│    NLRP3: NACHT 变构口袋 + MM-GBSA + 短程 MD                                │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 4  可靠性融合 + Pareto 排序                                         │
│    ω_U, ω_N 动态权重 → 双靶融合分 → 非支配排序 → Top 50–100                 │
├──────────────────────────────────────────────────────────────────────────┤
│  Stage 5  回顾性验证 + 消融                                               │
│    Benchmark 回收 + vs PLK1-style baseline + 双路径贡献分解                 │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## 三、数据现实（驱动设计变更）

基于用户 ChEMBL 导出（2026-06 分析）：

| 数据集 | 清洗后独特 SMILES | 建模策略 |
|--------|------------------|---------|
| URAT1 | **822** | 监督回归 + conformal UQ（数据充足） |
| NLRP3 | **513**（IL-1β + Assay B） | Assay-conditioned 分类（避免全局回归） |
| 双靶重叠 | **0** | 独立模型 + 证据融合（MTL 非主路径） |
| NLRP3 assay 异质性 | 39 assays，7.2% 化合物跨 assay pActivity 极差 >1 log（37/513；见 data_summary.json） | 必须条件化，禁用锚点相似性主策略 |

---

## 四、Stage 1：不对称双证据建模

### 4.1 URAT1 臂（数据丰富 + 转运体）

| 组件 | 方法 |
|------|------|
| 表示 | MiniMol 指纹（冻结）或 Chemprop D-MPNN |
| 任务 | pIC50 回归 |
| 迁移 | **OAT1/OAT3** 摄取抑制 → URAT1 微调；**OCT1/2** 仅脱靶讨论 |
| 不确定性 | **Split conformal prediction**（90% 预测区间） |
| 结构证据 | 构象系综 $S_{\text{trap}}$（见算法文档 §4） |

### 4.2 NLRP3 臂（数据中等 + 高异质性）

**主策略：Assay-conditioned 分类**

对 assay $a$ 嵌入向量 $\mathbf{e}_a$（测定类型、细胞系、终点 one-hot + 可学习 embedding）：

$$
P_{\text{active}}(x \mid a) = \sigma\big( f(\phi(x), \mathbf{e}_a) \big)
$$

实现选项（按算力）：
- **CLAMP**（assay 条件化对比学习 + 分类头）
- **TwinBooster**（assay-aware 梯度提升）
- **Chemprop + assay ID 作为 graph-level 条件**（轻量实现）

**辅策略：结构主导证据**

对 ML 不确定或 assay 外推样本，提升 $S_{\text{NLRP3}}^{\text{struct}}$ 权重：

$$
\omega_N^{\text{struct}} = 1 - c_N, \quad c_N = \max_a P_{\text{active}}(x \mid a)
$$

**明确不做**：5-anchor ECFP max-pooling、ESM-2 口袋为主表征（留给 PLK1-style 消融对照）。

### 4.3 可选 MTL（消融用）

共享 MiniMol 编码器 + 掩码双任务头，用于证明「0 重叠时 MTL 不优于独立模型 + 融合」。

---

## 五、Stage 2：双路径候选生成

### Path A — 库筛（主路径，稳健）

```
Library (~10⁶)
  → ADMET / PAINS / Lipinski
  → URAT1: ŷ^U ≥ t_U 且 conformal 下界 ≥ t_U^lo
  → NLRP3: P_active ≥ t_N（assay 条件化集成）
  → 双靶初筛 (~10⁴)
  → URAT1/NLRP3 系综对接
  → 可靠性融合 + 多样性聚类
  → Path A 候选 C_A (~500)
```

### Path B — 生成式（创新路径）

```
ChEMBL 预训练 CLM
  → URAT1 高活性集 fine-tune → M_U
  → NLRP3 高活性集 fine-tune → M_N
  → Cross-sampling: logits = α·logits_U + (1-α)·logits_N
  → RL 优化（REINFORCE/PPO, 3000–5000 steps）
  → 奖励 R(x): 双靶 ML + S_trap + NLRP3 struct + QED - SA + Nov
  → 后处理（有效/唯一/新颖/可合成）
  → Path B 候选 C_B (~500–2000)
```

### 合并策略

$$
\mathcal{C}_{\text{union}} = \mathcal{C}_A \cup \mathcal{C}_B \setminus \text{duplicates}
$$

对合并池重新计算融合分，标注来源标签 `source ∈ {library, generative}` 用于消融。

---

## 六、Stage 4：可靠性加权融合 + Pareto

### 6.1 动态权重

| 证据臂 | 可靠性指标 | 权重 |
|--------|-----------|------|
| URAT1 ML | Conformal 区间宽度 $w_U$ | $\omega_U^{\text{ml}} \propto 1/w_U$ |
| NLRP3 ML | Assay-conditioned 置信度 $c_N$ | $\omega_N^{\text{ml}} \propto c_N$ |
| URAT1 结构 | 系综 pose 一致性 | $\omega_U^{\text{struct}}$ |
| NLRP3 结构 | MD 稳定性 + MM-GBSA | $\omega_N^{\text{struct}}$ |

### 6.2 融合分

$$
S_{\text{dual}} = \sum_{i \in \{U,N\}} \omega_i^{\text{ml}} \cdot s_i^{\text{ml}} + \gamma \cdot \sqrt{S_U^{\text{struct}} \cdot S_N^{\text{struct}}}
$$

### 6.3 Pareto 多目标排序

目标向量：$(S_{\text{dual}}, \text{QED}, -\text{SA}, \text{Novelty})$

取 **非支配前沿** 第 1–2 层作为最终候选，避免单标量加权掩盖双靶均衡性。

---

## 七、必做消融实验（含 PLK1-style 对照）

| ID | 变体 | 目的 |
|----|------|------|
| Abl-1 | 无 $S_{\text{trap}}$（单 PDB URAT1） | 转运体构象必要性 |
| Abl-2 | NLRP3 锚点相似性（PLK1-style） | 证明 assay-conditioned 更优 |
| Abl-3 | 0.5/0.5 固定融合 | 证明可靠性加权更优 |
| Abl-4 | 仅 Path A（无生成） | 量化生成路径贡献 |
| Abl-5 | 仅 Path B（无库筛） | 量化库筛路径贡献 |
| Abl-6 | MTL vs 独立双模型 | 0 重叠数据下架构选择 |
| Abl-7 | 无 SLC22 迁移 | URAT1 小样本策略 |

---

## 八、与 STAD-AIDD v1.0 的关系

| STAD-AIDD v1.0 | TAPE-GATE v2.0 |
|----------------|----------------|
| MTL 为核心 | 独立双模型 + 融合为核心，MTL 可选 |
| 单路径库筛 + 可选生成 | **双路径并行**为正式设计 |
| NLRP3 全局回归 MTL 头 | **Assay-conditioned 分类** |
| 简单几何平均融合 | **可靠性加权 + Pareto** |
| 未对标 PLK1/NLRP3 | 显式差异化 + PLK1-style 消融 |

技术细节见 [`ALGORITHM_FRAMEWORK.md`](ALGORITHM_FRAMEWORK.md)。

---

## 九、推荐论文题目

1. *TAPE-GATE: Transporter-aware paired-path evidence fusion for URAT1/NLRP3 dual-target discovery under assay-heterogeneous conditions*
2. *Assay-conditioned and conformation-ensemble dual evidence for hyperuricemia dual-target virtual screening with generative augmentation*
3. *Beyond asymmetric similarity: library and generative paired screening for URAT1 and NLRP3 co-inhibition*
