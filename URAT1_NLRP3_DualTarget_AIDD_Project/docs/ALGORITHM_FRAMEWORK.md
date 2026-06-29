# STAD-AIDD 算法框架详解

> 本文档是论文 **Methods** 部分的技术核心，逐模块说明算法原理、输入输出、超参与公式。

---

## 一、符号与问题定义

设分子集合 $\mathcal{D} = \{(x_i, y_i^{U}, y_i^{N})\}$，其中：

- $x_i$：分子图或 SMILES
- $y_i^{U}$：URAT1 活性（pIC50 或 pKi），可能缺失
- $y_i^{N}$：NLRP3 活性（pIC50，IL-1β 终点），可能缺失

**目标**：学习映射 $f: x \mapsto (\hat{y}^{U}, \hat{y}^{N})$，并生成/筛选 **双靶优先** 候选集 $\mathcal{C}^*$：

$$
\mathcal{C}^* = \arg\max_{x \in \mathcal{X}} S_{\text{dual}}(x)
$$

其中 $S_{\text{dual}}$ 融合 ML 预测、结构对接、药物相似性与合成可及性。

**小数据约束**：$|\mathcal{D}^{U}| \sim 10^2$，$|\mathcal{D}^{N}| \sim 10^3$，$|\mathcal{D}^{U \cap N}| \ll 10$。

---

## 二、Stage 0：数据清洗算法

### 2.1 活性标准化

对每条记录 $r$：

1. 过滤 `standard_relation != '='`
2. 单位统一为 nM，计算 $\text{pActivity} = -\log_{10}(\text{value}_{\text{nM}} \times 10^{-9})$
3. 保留 $\text{pActivity} \in [4, 10]$

### 2.2 冲突消解（per SMILES）

对同一 $(\text{SMILES}, \text{target})$ 的多个测定值 $\{v_j\}$：

```
if std(v) > 0.5 OR (max(v) - min(v)) > 1.0:
    discard SMILES
else:
    y = median(v)
```

### 2.3 骨架分组划分

使用 **Murcko 骨架** $s(x)$ 作为分组标签（Bemis & Murcko, 1996）：

- GroupKFold($k=5$) 在骨架级别划分
- 保证测试集骨架与训练集 **零重叠** → 模拟真实新药发现外推能力

### 2.4 NLRP3 专用过滤

仅保留 assay description 含 `IL-1β` / `IL-1beta` / `interleukin-1` 的条目（参考 Zhao et al., *BMC Chemistry* 2024）。

### 2.5 URAT1 数据扩充

ChEMBL 数据不足时，从以下来源合并：
- 专利 WO2018015445（verinurad 系列）
- Dai et al. 2024 共晶配体及其 SAR 类似物
- 文献报道的 benzbromarone 系列

---

## 三、Stage 1：分子表示学习

### 3.1 分子基础模型：MiniMol

**原理**（Ahmadi et al., 2024; arXiv:2404.14986）：

- 骨干：GNN（MPNN++ / GINE）
- 预训练：Graphium LargeMix，**3300+ 稀疏多任务**（量子化学 + 生物活性）
- 输出：分子指纹 $\phi(x) \in \mathbb{R}^{d}$，$d \approx 512$–768
- 参数量约 **10M**（比 MolE 小 10×，适合小数据微调）

**为何选 MiniMol 而非从头训 GNN**：
- URAT1 样本 < 200 时，端到端 GNN 严重过拟合
- 预训练指纹已编码子结构–性质关系，下游只需浅层 MLP（Beaini et al., ICLR 2024 证明量子任务辅助生物任务）

### 3.2 多任务学习（MTL）架构

```
SMILES x
   ↓
[MiniMol encoder]  (frozen or partial unfreeze)
   ↓
φ(x) ∈ R^d
   ↓
┌──────────────┬──────────────┬──────────────┐
│  Head_U      │  Head_N      │  Head_dual   │
│  MLP → ŷ^U   │  MLP → ŷ^N   │  MLP → p̂    │
└──────────────┴──────────────┴──────────────┘
```

**损失函数**（缺失标签可掩码）：

$$
\mathcal{L} = \frac{1}{|\mathcal{M}_U|}\sum_{i \in \mathcal{M}_U} (y_i^U - \hat{y}_i^U)^2 + \frac{1}{|\mathcal{M}_N|}\sum_{i \in \mathcal{M}_N} (y_i^N - \hat{y}_i^N)^2 + \lambda \mathcal{L}_{\text{dual}}
$$

其中 $\mathcal{M}_U, \mathcal{M}_N$ 为有标签索引集；$\mathcal{L}_{\text{dual}}$ 为双活性分类 BCE（标签：$y^U \geq 6$ AND $y^N \geq 6$）。

**任务权重**：$\lambda = 0.5$（双活性样本极少，避免主导梯度）。

### 3.3 迁移学习：SLC22 辅助预训练

**动机**：URAT1 属于 SLC22 家族，OCT1/OCT2 有更多 ChEMBL 摄取数据。

**两阶段微调**：

1. **Stage A**：在 SLC22A1/A2 摄取抑制数据上微调 MLP head（冻结 MiniMol）
2. **Stage B**：在 URAT1 数据上继续微调 head
3. **Stage C**（可选）：解冻 MiniMol 最后一层 GNN，lr = 1e-5

参考 OCT1/OCT2 交替开放机制（Pan et al., *Nature* 2023; Feng et al., *Cell Res* 2023）。

### 3.4 Baseline 模型

| 模型 | 特征 | 用途 |
|------|------|------|
| XGBoost | ECFP4 (2048-bit) | 小数据强 baseline |
| Chemprop | D-MPNN, 单任务 | 图学习对照 |
| Random Forest | Morgan 2048 | 传统 QSAR |
| kNN + Tanimoto | 相似性外推 | 适用域分析 |

### 3.5 评估指标

**回归**：RMSE, MAE, $R^2$, Spearman $\rho$

**虚拟筛选**：EF@1%, EF@5%, BEDROC $\alpha=20$

**统计检验**：5-fold CV 上 Wilcoxon signed-rank test vs 最佳 baseline（ChemRxiv 2024 协议）

### 3.6 不确定性估计

**Ensemble**：5 个不同随机种子的 MTL 模型，预测区间 $[\hat{y} - 1.96\sigma, \hat{y} + 1.96\sigma]$

**应用**：筛选时丢弃高不确定性分子（$\sigma > 0.5$ log unit）

---

## 四、Stage 2：结构约束虚拟筛选（核心差异化）

### 4.1 为什么 URAT1 不能用「酶式对接」

酶抑制剂：阻断活性位点催化 → 单构象对接通常有效。

**URAT1 抑制机制**（Dai et al., *Cell Res* 2024）：
1. 配体占据 **Phe-rich 底物口袋**，与尿酸竞争
2. 稳定 **inward-facing** 或 **occluded** 构象
3. **阻断 alternating access 转运循环** → 转运抑制

因此评分必须回答：**该分子能否在不同构象态中选择性稳定抑制态？**

### 4.2 URAT1 构象系综对接

**输入结构**（见 `data/structures/docking_ensemble_pdb.csv`）：

| PDB | 构象态 | 意义 |
|-----|--------|------|
| 9B1H | inward-open + lesinurad | 药物结合参考 |
| 9JDZ | outward-open + urate | 底物转运入口态 |
| 9JDZ | occluded + urate | 中间态 |

**对接流程**：

```
for each conformer C_k in ensemble:
    pose_k = dock(ligand, C_k, binding_site_k)
    score_k = VinaScore(pose_k)

S_trap = w_in * score_inward + w_occ * score_occluded - w_out * score_outward
```

**构象捕获分** $S_{\text{trap}}$：inward/occluded 结合强、outward 结合弱 → 高分配体。

**关键残基接触分** $S_{\text{key}}$：

检查与文献报道关键位点距离 < 4 Å 的比例：
- Phe-rich cage 残基（底物识别）
- **Arg477**（verinurad/dotinurad 高亲和力关键，Dai 2024 Fig.6）

**综合 URAT1 结构分**：

$$
S_{\text{URAT1}}^{\text{struct}} = 0.30 S_{\text{trap}} + 0.25 S_{\text{vina}} + 0.25 S_{\text{key}} + 0.20 S_{\text{compete}}
$$

$S_{\text{compete}}$：与尿酸共晶位点空间重叠度（竞争性抑制）。

### 4.3 NLRP3 结构对接

**结合位点**：NACHT 域 Walker B 变构口袋（Coll et al., *Nat Commun* 2019; PDB 7ALV）

**机制**：小分子作为 **分子内胶（intramolecular glue）** 锁定 NLRP3 非活性构象。

**评分**：

$$
S_{\text{NLRP3}}^{\text{struct}} = 0.30 S_{\text{vina}} + 0.35 \Delta G_{\text{MMGBSA}} + 0.35 S_{\text{stability}}
$$

$S_{\text{stability}}$：50–100 ns MD 后：
- 配体 RMSD < 2.5 Å
- NACHT 亚域间距离波动 < 阈值（变构锁定）

### 4.4 双靶结构协同分

$$
S_{\text{dual}}^{\text{struct}} = \sqrt{S_{\text{URAT1}}^{\text{struct}} \cdot S_{\text{NLRP3}}^{\text{struct}}} + \beta \cdot \min(S_{\text{URAT1}}^{\text{struct}}, S_{\text{NLRP3}}^{\text{struct}})
$$

几何平均惩罚「单靶极强、另一靶极弱」的分子；$\min$ 项奖励均衡双靶（$\beta = 0.2$）。

### 4.5 ML + 结构融合漏斗

```
Library (10^6 SMILES)
    │  ADMET / PAINS / Lipinski
    ▼
ML filter: ŷ^U ≥ t_U AND ŷ^N ≥ t_N  (约 10^4)
    ▼
Docking ensemble: S_dual^struct top 5%  (约 500)
    ▼
Diversity clustering (Butina, Tc=0.4)
    ▼
Top 50–100 candidates
```

阈值 $t_U, t_N$ 在文献 benchmark 上校准（参考 JNK1 项目 `calibrate_threshold.py` 思路）。

---

## 五、Stage 3：生成式双靶优化

### 5.1 方法选择

| 方法 | 代表文献 | 适用性 |
|------|---------|--------|
| CLM + pooled fine-tuning | Schneider et al., *Nat Commun* 2024 (52060) | 小样本双靶偏置生成 |
| POLYGON (RL + VAE) | Ferreira et al., *Nat Commun* 2024 (47120) | 可定制双靶奖励 |
| MTMol-GPT (IRL) | *PLOS Comput Biol* 2024 | 多靶对抗生成 |

**推荐**：**CLM + RL 混合**（稳健且算力适中）

### 5.2 化学语言模型（CLM）

**预训练**：ChEMBL 全库 SMILES（或 MolGPT/ChemGPT 公开权重）

**Fine-tuning 集**：
- 正例：URAT1 高活性分子（p ≥ 6）+ NLRP3 高活性分子（p ≥ 6）的 **并集**
- 可选：已知双药理分子（如某些 NSAID，仅作弱监督）

**Cross fine-tuning**（Schneider 2024 技术）：
- 对 URAT1 活性分子集 fine-tune → 模型 $M_U$
- 对 NLRP3 活性分子集 fine-tune → 模型 $M_N$
- 生成时交替采样或 logits 插值：$\text{logits} = \alpha \text{logits}_U + (1-\alpha) \text{logits}_N$

### 5.3 强化学习奖励函数

参考 POLYGON，定义生成分子 $x$ 的奖励：

$$
R(x) = w_1 \hat{y}^U(x) + w_2 \hat{y}^N(x) + w_3 S_{\text{URAT1}}^{\text{struct}}(x) + w_4 S_{\text{NLRP3}}^{\text{struct}}(x) + w_5 \text{QED}(x) - w_6 \text{SA}(x) + w_7 \text{Nov}(x)
$$

| 项 | 权重建议 | 说明 |
|----|---------|------|
| $\hat{y}^U, \hat{y}^N$ | 0.20 each | MTL 模型预测 |
| 结构分 | 0.15 each | 系综对接（可每 100 步算一次以省算力） |
| QED | 0.10 | 药物相似性 |
| SA | -0.10 | 合成惩罚（Ertl 2009） |
| Nov | 0.05 | 与训练集最大 Tanimoto < 0.85 奖励 |

**RL 算法**：REINFORCE 或 PPO，训练 3000–5000 步。

### 5.4 生成候选质量控制

1. **有效性**：RDKit 解析成功率 > 95%
2. **唯一性**：去重 SMILES
3. **新颖性**：与 ChEMBL 最大 Tc < 0.85
4. **可合成性**：SA < 6，AiZynthFinder 至少一条合成路线
5. **PAINS 过滤**

---

## 六、Stage 4：回顾性验证算法

### 6.1 Benchmark 回收测试

对 `data/benchmarks/literature_benchmarks.csv` 中 **must_recover** 化合物：

```
Recall@K = |benchmark ∩ top_K_candidates| / |benchmark|
```

**通过标准**（无湿实验时的最低可信度）：
- URAT1 药物（lesinurad, benzbromarone）在 Top 500 回收
- NLRP3 工具化合物（MCC950）在 Top 500 回收
- 阴性对照（allopurinol 作为 URAT1 非抑制剂）排名靠后

### 6.2 消融实验（论文必做）

| 变体 | 去掉的模块 | 预期 |
|------|-----------|------|
| Abl-1 | 无构象系综（单 PDB） | URAT1 benchmark 回收下降 |
| Abl-2 | 无 MiniMol 预训练 | 小数据 CV $R^2$ 下降 |
| Abl-3 | 无 SLC22 迁移 | URAT1 任务 RMSE 上升 |
| Abl-4 | 单靶独立筛选再交集 | 双靶均衡候选减少 |
| Abl-5 | 无 RL 生成 | 化学空间覆盖度下降 |

### 6.3 外部验证集

- 从 NLRP3 专利 WO2021214284A1 留出一部分 **不参与训练** 的化合物作外部测试
- URAT1 用 Dai 2024 论文 SI 中部分类似物作外部测试

---

## 七、伪代码：端到端 Pipeline

```python
# run_stad_pipeline.py 逻辑摘要

# Stage 0
D_urat1, D_nlrp3 = prepare_data(config)
splits = murcko_group_split(D_urat1, D_nlrp3, k=5)

# Stage 1
encoder = MiniMolEncoder(frozen=True)
mtl_model = MTLHead(encoder, tasks=["urat1", "nlrp3", "dual"])
mtl_model = finetune(mtl_model, splits, auxiliary="SLC22")

# Stage 2
library = load_library("enamine_real.smi")  # ~10^6
hits_ml = filter_ml(library, mtl_model, thresholds)
hits_struct = ensemble_dock(hits_ml, urat1_ensemble, nlrp3_ensemble)
candidates = diversity_select(hits_struct, top_n=100)

# Stage 3 (optional)
clm = load_pretrained_clm()
clm = rl_finetune(clm, reward_fn=build_reward(mtl_model, ensembles))
generated = sample(clm, n=5000)
candidates_gen = postprocess(generated)

# Stage 4
report = retrospective_validate(candidates + candidates_gen, benchmarks)
ablation = run_ablations(config)
write_manuscript_figures(report, ablation)
```

---

## 八、算法创新点小结（对应论文 Contribution）

1. **Transporter-aware ensemble scoring** $S_{\text{trap}}$：首次在 URAT1 双靶框架中系统使用构象捕获评分（非传统酶口袋对接）。
2. **Hierarchical transfer for electrochemical transporter QSAR**：SLC22 家族 → URAT1 的分层迁移。
3. **Dual-target MTL with missing labels**：处理双靶数据极少重叠的实际场景。
4. **Structure-constrained generative RL**：奖励函数同时嵌入转运体系综与 NLRP3 变构对接。
5. **Reproducible small-data benchmark protocol**：骨架 CV + 文献回收 + 消融，符合 ChemRxiv 2024 / WelQrate 精神。

---

## 九、推荐阅读的实现参考

| 组件 | 开源实现 |
|------|---------|
| MiniMol / Graphium | https://github.com/datamol-io/graphium |
| Chemprop MTL | https://github.com/chemprop/chemprop |
| AutoDock Vina | https://github.com/ccsb-scripps/AutoDock-Vina |
| RDKit 指纹/骨架 | https://www.rdkit.org |
| POLYGON 思路 | Nat Commun 2024, 代码需自实现或参考 MolGen 库 |
| 对接系综 | 参考 JNK1 项目 `config/docking_ensemble.yaml` 模式 |

详细验证要求见 [`URAT1_TRANSPORTER_VALIDATION.md`](URAT1_TRANSPORTER_VALIDATION.md)。
