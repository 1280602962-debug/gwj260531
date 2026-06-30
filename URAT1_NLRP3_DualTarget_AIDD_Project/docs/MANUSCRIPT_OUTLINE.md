# SCI 论文大纲（TAPE-GATE v2.0）

> 建议目标：**Journal of Cheminformatics** 或 **JCIM**  
> 文章类型：Original Research / Methodology  
> 预计字数：6500–8500 words + SI

---

## 标题（推荐）

**TAPE-GATE: Transporter-aware paired-path evidence fusion for URAT1/NLRP3 dual-target discovery under assay-heterogeneous conditions**

备选：
- *Assay-conditioned and conformation-ensemble dual evidence for hyperuricemia dual-target screening with generative augmentation*
- *Beyond asymmetric similarity: library and generative paired screening for URAT1 and NLRP3 co-inhibition*

---

## Abstract（250 words 结构）

1. **Background**：HUA/痛风；URAT1 代谢 + NLRP3 炎症；数据无重叠 + assay 异质性 + 转运体机制
2. **Methods**：TAPE-GATE = 不对称双证据（conformal URAT1 + assay-conditioned NLRP3）+ $S_{\text{trap}}$ 系综 + **库筛/生成双路径** + 可靠性 Pareto 融合
3. **Results**：vs PLK1-style baseline 的 benchmark 回收；Path A/B 贡献；Top 候选性质
4. **Conclusions**：paired-path 计算优先策略；待实验验证

---

## 1. Introduction（~1300 words）

### 1.1–1.4 疾病与双靶 rationale
（同 v1.0，强调痛风场景）

### 1.5 计算药物发现的挑战
- URAT1 转运体 vs 激酶对接范式
- NLRP3 assay 异质性（引用用户数据统计：513 化合物、7.2% 跨 assay 离散）
- **0 重叠 SMILES** → 标准 MTL 失效
- PLK1/NLRP3 类 asymmetric 框架的局限（**不点名攻击，客观对比**）

### 1.6 本文贡献（7 bullet，对应 INNOVATION_POINTS.md C1–C7）

---

## 2. Materials and Methods（~2800 words）

### 2.1 Data collection and curation
- ChEMBL 822 URAT1 / 513 NLRP3（IL-1β + Assay B）
- Assay 元数据保留策略
- 0 重叠统计

### 2.2 Asymmetric dual-evidence modeling
- **URAT1**：MiniMol/Chemprop + conformal prediction + SLC22 transfer
- **NLRP3**：Assay-conditioned classification（公式 + 实现）
- Independent vs MTL（消融）

### 2.3 Paired-path candidate generation ★
- **Path A**：Library screening funnel
- **Path B**：CLM cross-fine-tune + RL reward（含 $S_{\text{trap}}$）
- Merge protocol

### 2.4 Conformational ensemble docking
- URAT1 $S_{\text{trap}}$（核心公式）
- NLRP3 NACHT + MM-GBSA/MD

### 2.5 Reliability-weighted Pareto fusion
- 动态权重公式
- vs 0.5/0.5 fixed fusion

### 2.6 PLK1-style baseline（Methods 中明示）
- SVR(URAT1) + 5-anchor similarity(NLRP3) + equal-weight fusion
- 用于公平对比，非本文方法

### 2.7 Retrospective validation and ablation
- 7 组消融 + benchmark 列表

---

## 3. Results（~2200 words）

### 3.1 Dataset characterization
- **Fig 2A**：UMAP 化学空间
- **Fig 2B**：NLRP3 assay 冲突热图
- **Fig 2C**：0 重叠示意

### 3.2 Asymmetric model performance
- **Fig 3A**：URAT1 conformal parity plot
- **Fig 3B**：NLRP3 assay-conditioned ROC/PR（分 assay）
- **Table S1**：vs XGBoost / PLK1-style baseline

### 3.3 PLK1-style vs TAPE-GATE comparison ★
- **Fig 4**：Benchmark 回收率柱状图（双方法）
- **Table 1**：相关工作 + 本框架对比

### 3.4 Structure-based screening
- **Fig 5**：URAT1 $S_{\text{trap}}$ — lesinurad 重现
- **Fig 6**：NLRP3 对接模式

### 3.5 Dual-path contribution
- **Fig 7A**：Path A vs Path B vs Union 回收率
- **Fig 7B**：生成候选新颖性（Tc 分布）
- **Table 2**：Top 20 候选（含 source 标签）

### 3.6 Ablation summary
- **Fig 8**：7 组消融雷达图或热力图

---

## 4. Discussion（~1300 words）

### 4.1 为何转运体需要 $S_{\text{trap}}$（非 PLK1 式对接）
### 4.2 Assay-conditioned 对 NLRP3 的必要性
### 4.3 双路径的互补性（库筛稳健 vs 生成新颖）
### 4.4 与 PLK1/NLRP3 框架的方法学边界
### 4.5 局限性与建议实验验证
### 4.6 未来：主动学习、FEP、湿实验合作

---

## 5. Conclusions（~200 words）

---

## 图表清单

| ID | 内容 |
|----|------|
| Fig 1 | TAPE-GATE 全流程 + 双路径 + 疾病网络 |
| Fig 2 | 数据集表征（UMAP + assay 冲突 + 0 重叠） |
| Fig 3 | 不对称双证据模型性能 |
| Fig 4 | **TAPE-GATE vs PLK1-style benchmark 回收** |
| Fig 5 | URAT1 构象系综 |
| Fig 6 | NLRP3 对接 |
| Fig 7 | 双路径贡献 + 生成新颖性 |
| Fig 8 | 消融实验 |
| Table 1 | 相关工作对比（含 PLK1/NLRP3） |
| Table 2 | Top 候选 |

---

## 写作注意事项

### 差异化表述

- 用 "assay-conditioned evidence" 而非仅 "reliability-driven"
- 用 "paired-path" 而非仅 "virtual screening"
- Methods 中 **单独小节** 描述 PLK1-style baseline，避免审稿人认为抄袭

### 开源

- 脚本：`run_tape_gate_pipeline.py`
- 配置：`config/dual_path.yaml`, `config/model_hierarchy.yaml`

---

## 投稿前 Checklist

- [ ] PLK1-style baseline 已运行并写入 Results
- [ ] Path A / Path B / Union 分路径回收率
- [ ] NLRP3 assay 元数据在 SI 中完整说明
- [ ] 7 组消融完整
- [ ] DIFFERENTIATION 文档中的模块均未作为主方法误用
