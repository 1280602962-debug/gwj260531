# SCI 论文大纲（稳健型 + AI 方法偏置）

> 建议目标：**Journal of Cheminformatics** 或 **Journal of Chemical Information and Modeling**  
> 文章类型：Original Research / Methodology  
> 预计字数：6000–8000 words + SI

---

## 标题（推荐）

**STAD-AIDD: A structure-constrained transporter-aware framework for dual-target URAT1/NLRP3 inhibitor discovery in hyperuricemia**

---

## Abstract（250 words 结构）

1. **Background**（2–3 句）：HUA/痛风；URAT1 代谢 + NLRP3 炎症双通路；小数据 + 转运体机制挑战
2. **Methods**（3–4 句）：STAD-AIDD = foundation model MTL + conformational ensemble + RL generation
3. **Results**（3–4 句）：CV 性能；benchmark 回收率；Top 候选性质
4. **Conclusions**（1–2 句）：计算优先双靶策略；待实验验证

---

## 1. Introduction（~1200 words）

### 1.1 高尿酸血症与痛风流行病学
- 患病率上升，代谢综合征共病
- 引用：G-CAN 2023 review; Frontiers Immunology 2023 NLRP3-gout

### 1.2 URAT1：尿酸重吸收与促排泄治疗
- SLC22A12，~90% 重吸收
- 80–90% 排泄障碍型 HUA
- 引用：Dai 2024; Fedor 2025

### 1.3 NLRP3 炎症小体与 MSU 晶体
- Martinon 2006 经典；痛风急性发作机制
- 降尿酸 vs 抗炎的时间差问题
- 引用：Chen 2023; Li 2023 NLRP3 review (J Med Chem)

### 1.4 双靶协同治疗 rationale
- 疾病网络图（Fig 1A）
- 单靶/联合用药局限

### 1.5 AI 药物发现现状与 Gap
- 小数据、转运体、双靶挑战
- 本文贡献（5 bullet points，对应 INNOVATION_POINTS.md）

---

## 2. Materials and Methods（~2500 words）

### 2.1 Data collection and curation
- ChEMBL, patents, assay filtering
- 清洗规则、样本量统计表

### 2.2 Molecular representation and multi-task learning
- MiniMol 预训练
- MTL 架构、损失函数、SLC22 迁移
- Baselines、CV 协议

### 2.3 Conformational ensemble docking
- **URAT1**：构象态、$S_{\text{trap}}$ 定义（核心公式）
- **NLRP3**：NACHT 口袋、MM-GBSA
- 双靶融合评分

### 2.4 Generative dual-target optimization
- CLM + RL
- 奖励函数各组分

### 2.5 Retrospective validation and ablation
- Benchmark 化合物列表
- 消融 5 变体

### 2.6 ADMET and synthetic accessibility
- SwissADME, SA score, PAINS

---

## 3. Results（~2000 words）

### 3.1 Dataset characterization
- **Fig 2**：UMAP 化学空间（URAT1 vs NLRP3）
- 骨架多样性、活性分布
- 双靶重叠极少 → 方法动机

### 3.2 Multi-task model performance
- **Fig 3A**：5-fold CV parity plot
- **Fig 3B**：消融柱状图（MTL vs single-task vs no-transfer）
- **Table S1**：完整指标 RMSE/R²/EF@1%

### 3.3 Structure-based screening
- **Fig 4**：URAT1 构象系综对接 — lesinurad 重现
- **Fig 5**：NLRP3 对接模式 vs MCC950/GDC-2394
- 漏斗统计

### 3.4 Retrospective benchmark recovery
- **Fig 6**：已知药物在排名分布中的位置
- Recall@100, @500

### 3.5 Generated dual-target candidates
- **Fig 7**：生成候选评分分布、代表结构
- **Table 2**：Top 20 候选（SMILES, 预测 pActivity, 结构分, QED, SA）

### 3.6 Case studies
- 选 2–3 个候选详细分析相互作用（URAT1 Arg477, NLRP3 NACHT glue）

---

## 4. Discussion（~1200 words）

### 4.1 转运体感知对接的必要性
- 对比单结构对接消融结果

### 4.2 小数据 MTL 的有效性边界
- 适用域分析；何时预测不可靠

### 4.3 双靶计算设计的局限
- 无实验验证；双靶分子成药性挑战
- 与联合用药的比较

### 4.4 建议实验验证方案
- 引用 URAT1_TRANSPORTER_VALIDATION.md §6

### 4.5 未来方向
- 主动学习闭环；FEP；湿实验合作

---

## 5. Conclusions（~200 words）

---

## 图表清单

| ID | 类型 | 内容 |
|----|------|------|
| Fig 1 | 示意图 | STAD-AIDD 全流程 + 疾病网络 |
| Fig 2 | UMAP | 数据集化学空间 |
| Fig 3 | 性能图 | MTL CV + 消融 |
| Fig 4 | 3D | URAT1 系综对接 |
| Fig 5 | 3D | NLRP3 对接 |
| Fig 6 | 柱状/散点 | Benchmark 回收 |
| Fig 7 | 化学结构 | 生成候选 |
| Table 1 | 对比 | 相关工作 |
| Table 2 | 数据 | Top 候选 |

### Supplementary

- SI Table S1–S5：完整候选、对接参数、超参、外部验证
- SI Fig S1–S4：MD 轨迹、更多对接 pose、SHAP 解释

---

## 写作注意事项

### 语言

- 全文英文；可另附中文导读（非投稿用）
- 避免 "we discovered a new drug" → 用 "computational prioritization"

### 审稿应对

- 预备 Response Letter 模板：针对无实验、数据少、转运体验证

### 开源要求

- J Cheminformatics 要求代码可用；GitHub + Zenodo DOI
- 提供 `requirements.txt` 与复现脚本

---

## 投稿前 Checklist

- [ ] 所有 PDB 引用含 DOI
- [ ] 骨架 CV 非随机划分
- [ ] Benchmark 化合物 SMILES 准确
- [ ] 消融完整
- [ ] Limitations 段落诚实
- [ ] 作者贡献、利益冲突、数据可用性声明
- [ ] Graphical Abstract（J Cheminformatics 需要）
