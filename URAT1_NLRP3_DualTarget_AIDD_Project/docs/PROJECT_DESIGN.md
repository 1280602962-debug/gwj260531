# TAPE-GATE 项目总体设计

> **版本**: 2.0（TAPE-GATE）  
> **更新**: 2026-06  
> **类型**: 纯计算方法学（无湿实验）  
> **前身**: STAD-AIDD v1.0

---

## 一、科学背景与论文叙事逻辑

### 1.1 疾病分层：代谢层 + 炎症层

高尿酸血症（HUA）的危害来自两条并行通路：

```
高血尿酸 → MSU 晶体沉积 → 巨噬细胞吞噬
                              ↓
                    NLRP3 炎症小体激活
                              ↓
              Caspase-1 → IL-1β / IL-18 成熟释放
                              ↓
           急性痛风发作 + 慢性代谢炎症（CKD、AS、MS）
```

- **代谢调控层**：肾脏近曲小管 **URAT1（SLC22A12）** 重吸收约 90% 滤过尿酸（Dai et al., *Cell Res* 2024; Fedor et al., *Nat Commun* 2025）。
- **炎症效应层**：**MSU 晶体** 激活 **NLRP3 炎症小体** 是痛风急性发作的核心机制（Martinon et al., *Nature* 2006; Chen et al., *J Inflamm Res* 2023）。

### 1.2 双靶协同的合理性

| 靶点 | 层级 | 作用 | 代表药物/工具 |
|------|------|------|--------------|
| URAT1 | 代谢 | 减少尿酸重吸收，降低 MSU 形成风险 | lesinurad, verinurad, dotinurad |
| NLRP3 | 炎症 | 阻断 IL-1β 释放，抑制急性/慢性炎症 | MCC950, GDC-2394, NT-0796 |

**论文核心论点**：URAT1 与 NLRP3 是疾病网络中 **上游代谢压力** 与 **下游炎症放大** 的耦合节点；双靶小分子可实现协同治疗。

### 1.3 研究空白（Gap）

1. URAT1/NLRP3 双靶小分子系统计算研究极少。
2. 公开数据 **无重叠化合物**（ChEMBL 实测 0 shared SMILES），标准双靶 QSAR/MTL 失效。
3. URAT1 是 **膜转运蛋白**，须构象系综评分，不能用激酶式单结构对接。
4. NLRP3 存在 **严重 assay 异质性**（47% 多 assay 化合物 >1 log 差），全局回归不可靠。
5. 现有 PLK1/NLRP3 类不对称框架采用锚点相似性 + 固定融合，**不适用于 URAT1 转运体场景**，且方法学创新性不足。

---

## 二、TAPE-GATE 框架总览

```
┌─────────────────────────────────────────────────────────────────┐
│                    TAPE-GATE 五阶段 + 双路径                        │
├─────────────────────────────────────────────────────────────────┤
│  Stage 0: 数据层                                                  │
│    URAT1(822) + NLRP3(503) + assay 元数据 + SLC22 辅助            │
├─────────────────────────────────────────────────────────────────┤
│  Stage 1: 不对称双证据建模                                        │
│    URAT1: 回归 + Conformal UQ + SLC22 迁移                        │
│    NLRP3: Assay-conditioned 分类（非锚点相似性）                   │
├─────────────────────────────────────────────────────────────────┤
│  Stage 2: 双路径候选生成 ★                                        │
│    Path A 库筛: Enamine ~10⁶ → ML/UQ → 对接                       │
│    Path B 生成: CLM cross-fine-tune + RL 双靶奖励                  │
│    → 合并候选池 C_union                                            │
├─────────────────────────────────────────────────────────────────┤
│  Stage 3: 结构约束层                                              │
│    URAT1 $S_{\text{trap}}$ 构象系综 + NLRP3 NACHT 变构对接         │
├─────────────────────────────────────────────────────────────────┤
│  Stage 4: 可靠性加权融合 + Pareto 排序                             │
├─────────────────────────────────────────────────────────────────┤
│  Stage 5: 回顾性验证 + 消融（含 PLK1-style 阴性对照）              │
└─────────────────────────────────────────────────────────────────┘
```

详细算法见 [`TAPE_GATE_FRAMEWORK.md`](TAPE_GATE_FRAMEWORK.md)、[`ALGORITHM_FRAMEWORK.md`](ALGORITHM_FRAMEWORK.md)。  
与 PLK1/NLRP3 差异化见 [`DIFFERENTIATION_VS_PLK1_NLRP3.md`](DIFFERENTIATION_VS_PLK1_NLRP3.md)。

---

## 三、为什么这是「稳健 + 创新」平衡型

| 设计 | 做法 |
|------|------|
| 不夸大 ML | 骨架 CV + conformal 区间；vs XGBoost / PLK1-style baseline |
| 转运体专属 | $S_{\text{trap}}$ 构象系综，单独章节论证 |
| NLRP3 异质性 | Assay-conditioned 分类，保留 assay 元数据 |
| 双路径 | 库筛（稳健）+ 生成式（创新），分别报告贡献 |
| 融合可辩护 | 可靠性动态权重 + Pareto，非黑箱 0.5/0.5 |
| 可复现 | 固定种子、公开配置、PLK1-style 消融对照 |
| 无湿实验 | lesinurad、MCC950 等 benchmark 强制回收 |

---

## 四、分阶段实施计划

### Phase 0：数据准备

| 数据源 | URAT1 | NLRP3 |
|--------|-------|-------|
| ChEMBL（用户 CSV） | **822** 独特 SMILES | **503**（IL-1β + Assay B） |
| 专利/文献 | verinurad 系列 | WO2021214284A1（外部验证） |
| 辅助 | SLC22A1/A2 摄取数据 | THP-1 子集 ~359 |

**关键**：NLRP3 导出时保留 `assay_id`, `assay_type`, `cell_line`

### Phase 1：数据集表征

- UMAP 化学空间（URAT1 vs NLRP3，证明正交）
- NLRP3 assay 冲突热图（支撑条件化建模）
- **0 重叠** → 论证独立模型 + 双路径必要性

### Phase 2：不对称双证据模型

**URAT1**：MiniMol/Chemprop + conformal + SLC22 迁移  
**NLRP3**：Assay-conditioned Chemprop/CLAMP  
**对照**：PLK1-style（SVR + 锚点相似性 + 0.5 融合）

### Phase 3：双路径候选生成

**Path A**：Enamine REAL 库筛漏斗  
**Path B**：CLM + RL 生成 500–2000 分子  
**合并**：去重 + 来源标注

### Phase 4：结构约束与融合排序

- URAT1 9B1H/9DKB/9JDZ 系综
- NLRP3 7ALV/8ETR + MM-GBSA/MD
- 可靠性加权 + Pareto → Top 50–100

### Phase 5：回顾性验证与论文

- Benchmark 回收（分 Path A/B/union）
- 7 组消融（见 TAPE_GATE_FRAMEWORK.md）
- 撰写与开源

---

## 五、目标期刊策略

### 5.1 稳健型首选

| 期刊 | 适配理由 |
|------|---------|
| **Journal of Cheminformatics** | 开源 pipeline + 双路径方法学 |
| **JCIM** | 小数据 ML + 转运体对接经典阵地 |
| **Briefings in Bioinformatics** | 强调 assay-conditioned + 消融协议 |

### 5.2 冲高创新

| 期刊 | 条件 |
|------|------|
| **Artificial Intelligence in Chemistry** | Path B 生成模块完整 + PLK1-style 对照显著优 |
| **Computers in Biology and Medicine** | 疾病网络叙事 + 双路径结果 |

### 5.3 标题角度

强调 **transporter-aware conformation-ensemble** + **assay-conditioned** + **paired-path generative**，避免仅用 "reliability-driven asymmetric"（与 PLK1/NLRP3 撞车）。

---

## 六、无湿实验时的结果表述

### 可以写

- TAPE-GATE 在骨架 CV 与 benchmark 回收上优于 PLK1-style baseline
- Path B 生成候选的化学空间覆盖度与新颖性
- $S_{\text{trap}}$ 消融提升 URAT1 药物回收率
- Assay-conditioned NLRP3 优于锚点相似性（Abl-2）

### 不能写

- ❌ 「发现了新型双靶先导化合物」并暗示已实验验证
- ✅ 「computational prioritization via paired-path evidence fusion」

---

## 七、与 PLK1/NLRP3 及 JNK1 项目的差异

| 维度 | PLK1/NLRP3 | JNK1 项目 | **TAPE-GATE** |
|------|-----------|----------|---------------|
| 主靶类型 | 激酶 | 激酶亚型 | **转运体** |
| NLRP3 策略 | 锚点相似性 | N/A | **Assay-conditioned** |
| 融合 | 0.5/0.5 | 选择性评分 | **可靠性 + Pareto** |
| 候选来源 | 仅库筛 | 库筛 | **库筛 + 生成式** |
| 结构核心 | 单结构对接 | 激酶口袋 | **$S_{\text{trap}}$ 系综** |

---

## 八、关键风险与缓解

| 风险 | 缓解 |
|------|------|
| 与 PLK1/NLRP3 审稿撞车 | 差异化文档 + PLK1-style 消融对照 |
| NLRP3 assay 噪声 | 条件化建模 + THP-1 子集 + 结构加权 |
| 0 双靶重叠 | 独立模型 + 双路径 + 融合评分 |
| 生成不可合成 | SA/QED + 逆合成 + 与库筛候选对比 |
| 算力 | Path B 可降采样；Path A 单独可发 JCIM |

---

## 九、下一步行动

详见 [`PREPARATION_CHECKLIST.md`](PREPARATION_CHECKLIST.md)。

**MVP 路径**：
1. Phase 0–2 + Path A + 验证 → JCIM 级别
2. 加入 Path B + PLK1-style 消融 → Briefings in Bioinformatics / AI in Chemistry
