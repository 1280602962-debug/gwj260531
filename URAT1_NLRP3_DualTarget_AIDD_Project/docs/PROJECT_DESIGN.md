# STAD-AIDD 项目总体设计

> **版本**: 1.0  
> **更新**: 2026-06  
> **类型**: 纯计算方法学（无湿实验）

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

临床流行病学显示约 **80%–90%** HUA 患者为尿酸排泄障碍型，促尿酸排泄（URAT1 抑制）符合多数患者病理生理；但 **单纯降尿酸不能迅速终止已启动的炎症级联**，而单独抗炎不能消除 MSU 持续刺激（Frontiers in Immunology 2023 综述）。

### 1.2 双靶协同的合理性

| 靶点 | 层级 | 作用 | 代表药物/工具 |
|------|------|------|--------------|
| URAT1 | 代谢 | 减少尿酸重吸收，降低 MSU 形成风险 | lesinurad, verinurad, dotinurad |
| NLRP3 | 炎症 | 阻断 IL-1β 释放，抑制急性/慢性炎症 | MCC950, GDC-2394, NT-0796 |

**论文核心论点**：URAT1 与 NLRP3 不是孤立靶点，而是疾病网络中 **上游代谢压力** 与 **下游炎症放大** 的耦合节点；双靶小分子可实现协同治疗，优于单靶或简单联合用药（依从性、DDI 风险）。

### 1.3 研究空白（Gap）

1. URAT1 与 NLRP3 抑制剂均有临床/临床前进展，但 **URAT1/NLRP3 双靶小分子系统研究极少**。
2. 公开活性数据均有限（URAT1 ChEMBL 约 **数十至百余条**；NLRP3 合并专利后约 **400–1200 条**），难以直接训练大模型。
3. **URAT1 是膜转运蛋白**，抑制依赖构象捕获而非催化位点阻断，传统酶导向对接易失效。
4. 缺乏面向 **小数据 + 多靶点 + 转运体结构约束** 的统一 AI 框架。

---

## 二、STAD-AIDD 框架总览

```
┌─────────────────────────────────────────────────────────────────┐
│                    STAD-AIDD 四阶段流水线                          │
├─────────────────────────────────────────────────────────────────┤
│  Stage 0: 数据层                                                  │
│    ChEMBL + 专利 + 文献 → 清洗 → 骨架分组划分                      │
│    + SLC22 家族辅助数据（迁移学习）                                 │
├─────────────────────────────────────────────────────────────────┤
│  Stage 1: 表示学习层（AI 核心）                                    │
│    MiniMol/Chemprop 预训练指纹 → 多任务头（URAT1 + NLRP3 + dual）  │
│    小样本鲁棒微调 + 不确定性估计                                    │
├─────────────────────────────────────────────────────────────────┤
│  Stage 2: 结构约束层（转运体关键）                                  │
│    URAT1 构象系综对接 + 转运循环阻断评分                            │
│    NLRP3 NACHT 变构口袋对接 + MM-GBSA/MD 稳定性                    │
├─────────────────────────────────────────────────────────────────┤
│  Stage 3: 生成优化层（创新亮点）                                    │
│    RL 微调化学语言模型 → 双靶奖励函数 → 候选 de novo 设计           │
├─────────────────────────────────────────────────────────────────┤
│  Stage 4: 验证层（无湿实验的替代）                                  │
│    文献 benchmark 回顾 / 骨架分组 CV / 外部专利集 / 消融实验        │
└─────────────────────────────────────────────────────────────────┘
```

---

## 三、为什么这是「稳健型」而非「空想型」

| 稳健性设计 | 具体做法 |
|-----------|---------|
| 不夸大 ML 精度 | 报告骨架分组 CV + 置信区间；与 XGBoost/单任务 Chemprop 严格对比 |
| 转运体不被当酶处理 | 单独章节 + 构象系综 + 转运抑制评分（见 `URAT1_TRANSPORTER_VALIDATION.md`） |
| 无湿实验仍有说服力 | 已知药物（lesinurad, MCC950 等）必须被漏斗 **回顾性回收** |
| 可复现 | 固定随机种子、公开配置、WelQrate 式数据清洗协议 |
| 创新可辩护 | 框架是 **组合创新**（foundation model + transporter ensemble + dual-target RL），非黑箱 |

---

## 四、分阶段实施计划

### Phase 0：数据准备（1–2 周）

**目标**：构建可发表质量的训练/测试集。

| 数据源 | URAT1 | NLRP3 |
|--------|-------|-------|
| ChEMBL | CHEMBL6120，约 20–100 化合物级记录 | CHEMBL1741208，约 400–530 条（清洗后） |
| 专利 WO | lesinurad/verinurad 系列 | NLRP3 抑制剂专利（参考 Zhao et al. 2024 九项专利） |
| 文献补充 | Dai 2024 共晶配体、Verinurad 系列 SAR | MCC950、GDC-2394、NP3-562 等 |

**清洗规则**（与 JNK1 项目一致，见 `config/targets.yaml`）：
- 仅保留 `Standard Relation = '='`
- pActivity 4.0–10.0
- 同 SMILES 冲突丢弃
- **Murcko 骨架分组** 划分 train/val/test

**预期样本量**（清洗后）：
- URAT1：**80–150** 独特分子（需专利扩充）
- NLRP3：**350–800** 独特分子
- 双靶重叠：**< 10**（预期极少 → 这正是方法学动机）

### Phase 1：数据集表征（3–5 天）

- 化学空间 PCA/UMAP（URAT1 vs NLRP3 vs 双活性交集）
- 骨架多样性、活性分布、Murcko 唯一性
- **关键图表**：证明两靶点化学空间部分正交 → 需要生成式/结构约束而非简单拼药效团

### Phase 2：多任务活性模型（1–2 周）

**主模型**：MiniMol 指纹 + MLP 多任务头  
**Baseline**：XGBoost (ECFP4)、Chemprop 单任务、Random Forest

**评估协议**（参考 ChemRxiv 2024 方法比较指南）：
- 5-fold 骨架 GroupKFold
- 指标：RMSE, MAE, R², Spearman, EF@1%
- Wilcoxon 符号秩检验 vs baseline

**小数据策略**：
1. 冻结 MiniMol，仅训 MLP head
2. 辅助任务：SLC22A1/A2 摄取抑制数据预训练
3. 不确定性：ensemble 或 conformal prediction 区间

### Phase 3：结构约束虚拟筛选（2–3 周）

**URAT1**（重点）：
- 下载 PDB 9B1H / 9DKB / 9JDZ
- 准备 outward-open、inward-open、occluded 代表构象
- 对接 + **构象捕获评分**（见算法文档）

**NLRP3**：
- PDB 7ALV、8ETR，NACHT 域口袋
- Glide SP/XP 或 AutoDock Vina + 50 ns MD + MM-GBSA

**漏斗**：ML 预测 → 双靶几何平均 → 结构评分 → 多样性聚类 → Top 50–100 候选

### Phase 4：生成式双靶优化（2–4 周，可选但建议做）

- 基于已知 URAT1 与 NLRP3 活性分子 fine-tune 化学语言模型
- RL 奖励：双靶预测活性 + 双靶对接 + QED + SA
- 参考：POLYGON (Nat Commun 2024)、CLM dual-target (Nat Commun 2024)

### Phase 5：回顾性验证与论文撰写（2–3 周）

- Benchmark 回收率：lesinurad, benzbromarone, MCC950, GDC-2394 等
- 消融：去掉结构约束 / 去掉迁移学习 / 单靶 vs 双靶
- 输出最终候选 + 可合成性评估

---

## 五、目标期刊策略

### 5.1 稳健型首选（方法学 + 应用平衡）

| 期刊 | IF 区间 | 适配理由 |
|------|---------|---------|
| **Journal of Cheminformatics** | ~7 | 开放获取，接受计算流程 + 开源代码 |
| **Journal of Chemical Information and Modeling** | ~5 | 小数据 ML + 对接方法经典阵地 |
| **Briefings in Bioinformatics** | ~7 | 若强调 AI 框架与 benchmark 协议 |
| **Pharmaceutics** (MDPI) | ~5 | 痛风治疗背景 + 计算药学，审稿相对快 |

### 5.2 冲高创新（AI 偏重）

| 期刊 | 风险 | 条件 |
|------|------|------|
| **Nature Communications** (子刊级方法) | 高 | 需强消融 + 可能需合作方验证 |
| **Artificial Intelligence in Chemistry** | 中 | 生成式模块完整、对比充分 |
| **Computers in Biology and Medicine** | 中 | 疾病网络叙事 + 完整 pipeline |

### 5.3 建议策略

**主投**：Journal of Cheminformatics 或 JCIM  
**备投**：Pharmaceutics / Briefings in Bioinformatics  
**标题角度**：强调 **transporter-aware** + **small-data dual-target** + **reproducible benchmark**，而非「发现了新药」。

---

## 六、无湿实验时的「结果」应如何表述

### 可以写的（Computational findings）

- 框架在骨架 CV 上优于 baseline（统计显著）
- 回顾性回收已知 URAT1/NLRP3 抑制剂
- 生成候选具有合理 ADMET、合成可及性、双靶对接模式
- 构象系综对接比单结构对接提升 benchmark 回收率（消融证明）

### 不能写的（需避免过度声称）

- ❌ 「发现了新型双靶先导化合物」并暗示已实验验证
- ❌ 「体内有效」
- ✅ 「computational prioritization of dual-target candidates warranting experimental validation」

### 增强可信度的补充（仍无需自建 lab）

1. **SwissADME / pkCSM** 预测理化性质
2. **SwissTargetPrediction** 脱靶风险
3. **合成可及性 SA score** + 逆合成分析（AiZynthFinder 等）
4. 在 Discussion 中明确列出 **建议实验验证方案**（见 `URAT1_TRANSPORTER_VALIDATION.md`）

---

## 七、与现有 JNK1 项目的差异

| 维度 | JNK1 选择性项目 | URAT1/NLRP3 双靶项目 |
|------|----------------|---------------------|
| 靶点关系 | 同源激酶亚型 | 跨通路（代谢 + 炎症） |
| 数据量 | 444–1147/靶点 | 80–800/靶点（更少） |
| 结构挑战 | 激酶保守 ATP 口袋 | 转运体构象动态 + NLRP3 变构 |
| 核心创新 | 选择性 ML + 对接 | 转运体感知 + 双靶 MTL + 生成式 |
| 验证重点 | 亚型选择性 | 构象捕获 + 双靶协同评分 |

---

## 八、关键风险与缓解

| 风险 | 缓解 |
|------|------|
| URAT1 数据过少 | 专利数据 + SLC22 迁移 + 结构对接主导后期漏斗 |
| 双靶无重叠训练样本 | 多任务学习 + 生成式 RL + 分别验证后融合评分 |
| 审稿人质疑无实验 | 强回顾性 benchmark + 消融 + 透明局限性讨论 |
| 对接假阳性 | 系综 + MD + 与共晶配体 RMSD 对照 |
| 生成分子不可合成 | SA/QED 约束 + 文献相似性过滤 |

---

## 九、下一步行动

详见 [`PREPARATION_CHECKLIST.md`](PREPARATION_CHECKLIST.md) 与 [`ALGORITHM_FRAMEWORK.md`](ALGORITHM_FRAMEWORK.md)。

**最小可行发表（MVP）路径**：
1. 完成 Phase 0–3 + 回顾验证（不做生成式也可发 JCIM 级别）
2. 若加入 Phase 4 生成模块 + 完整消融 → 可冲 Briefings in Bioinformatics / 更高档期刊
