# TAPE-GATE 算法框架详解

> 本文档是论文 **Methods** 部分的技术核心。框架名：**TAPE-GATE** v2.0  
> 总览见 [`TAPE_GATE_FRAMEWORK.md`](TAPE_GATE_FRAMEWORK.md)；与 PLK1/NLRP3 差异化见 [`DIFFERENTIATION_VS_PLK1_NLRP3.md`](DIFFERENTIATION_VS_PLK1_NLRP3.md)。

---

## 一、符号与问题定义

设两个 **无标签重叠** 的单靶数据集：

- $\mathcal{D}_U = \{(x_i, y_i^U)\}$，URAT1 pIC50，$|\mathcal{D}_U| \approx 822$
- $\mathcal{D}_N = \{(x_j, y_j^{N}, a_j)\}$，NLRP3 活性 + **assay 条件** $a_j$，$|\mathcal{D}_N| \approx 513$（609 条记录）

**目标**：从库筛与生成式两条路径得到候选池 $\mathcal{C}$，最大化双靶优先排序：

$$
\mathcal{C}^* = \text{ParetoTop}\big(\{x : S_{\text{dual}}(x)\}\big)
$$

其中 $S_{\text{dual}}$ 融合 ML 证据、结构证据、药物相似性与各臂 **可靠性权重**。

---

## 二、Stage 0：数据清洗算法

### 2.1 URAT1 清洗

1. 保留 `standard_relation = '='`，单位统一 nM
2. $\text{pActivity} = -\log_{10}(\text{value}_{\text{nM}} \times 10^{-9}) \in [4, 10]$
3. 同 SMILES 冲突：std > 0.5 或 range > 1.0 → 丢弃，否则 median
4. 终点：urate uptake / urate transport 相关 assay
5. Murcko 骨架 GroupKFold($k=5$)

### 2.2 NLRP3 清洗（关键：assay 元数据保留）

1. 仅保留 IL-1β / IL-1beta / interleukin-1 终点
2. 推荐子集：**Assay Type B**（用户数据验证，同质性更高）
3. 可选 THP-1 子集（curated 约 **302** 独特 SMILES / 313 条记录）作高置信训练域
4. **保留 assay_id、cell_line、assay_type** 用于条件化建模
5. 活性二值化：$y^N_{\text{bin}} = \mathbb{1}[y^N \geq 6]$（分类主任务）

### 2.3 专利与文献扩充

- URAT1：WO2018015445、Dai 2024 共晶 SAR
- NLRP3：WO2021214284A1 等（外部验证集，不参与训练）

---

## 三、Stage 1：不对称双证据建模

### 3.1 URAT1 臂：监督回归 + Conformal UQ

**编码器**：MiniMol $\phi(x) \in \mathbb{R}^d$（冻结或末层微调）

**回归头**：

$$
\hat{y}^U = f_U(\phi(x)), \quad \mathcal{L}_U = \frac{1}{|\mathcal{D}_U|}\sum_i (y_i^U - \hat{y}_i^U)^2
$$

**Conformal prediction**（非 PLK1 式 bootstrap/kNN）：

1. 在 calibration set 上计算残差 $r_i = |y_i - \hat{y}_i|$
2. 取 $(1-\alpha)$ 分位数 $q_{1-\alpha}$
3. 预测区间 $[\hat{y}^U - q, \hat{y}^U + q]$，区间宽度 $w_U = 2q$

**筛选规则**：$\hat{y}^U \geq t_U$ 且 **下界** $\geq t_U^{\text{lo}}$（保守策略）

**SLC22 迁移**（三阶段；**OAT 优先**）：

1. **OAT1/OAT3**（SLC22A6/A8）摄取/抑制数据上微调 MLP head — 阴离子底物化学空间，与 URAT1 亚家族一致
2. URAT1 数据上继续微调 head
3. 可选解冻 MiniMol 末层，lr = 1e-5

**OCT1/OCT2** 数据 **不** 作为主迁移源；仅用于脱靶对接与 $R_{\text{sel}}$（Tier 3）。须做消融：无 OAT 迁移 vs 无 OCT 脱靶特征。

### 3.2 NLRP3 臂：Assay-conditioned 分类

**禁止作为主路径**：多锚点 ECFP max-pooling（见差异化文档）

**Assay 嵌入** $\mathbf{e}_a$：

```
e_a = [onehot(assay_type); onehot(cell_line); Embedding(assay_id)]
```

**条件化分类器**：

$$
P_{\text{active}}(x \mid a) = \sigma\big( \text{MLP}([\phi(x); \mathbf{e}_a]) \big)
$$

**训练**：加权 BCE，按 assay 样本量逆频率加权，缓解 39 assays 长尾

**推理（库筛）**：对关键 assay 集合 $\mathcal{A}^*$（如 THP-1 + IL-1β）取：

$$
P_{\text{active}}^{\text{ens}}(x) = \frac{1}{|\mathcal{A}^*|}\sum_{a \in \mathcal{A}^*} P_{\text{active}}(x \mid a)
$$

**置信度**：$c_N = \max\big(P_{\text{ens}}, 1 - \text{entropy}(\{P(x|a)\})\big)$

**实现备选**：

| 方法 | 复杂度 | 适用 |
|------|--------|------|
| Chemprop + assay_id 条件 | 低 | MVP 首选 |
| TwinBooster assay-aware | 中 | 表格特征强时 |
| CLAMP 对比学习 | 高 | 冲高期刊 |

### 3.3 独立模型 vs MTL

**默认**：URAT1、NLRP3 **独立训练**，Stage 4 融合

**可选 MTL**（消融 Abl-6）：

$$
\mathcal{L}_{\text{MTL}} = \mathcal{L}_U + \mathcal{L}_N + \lambda \mathcal{L}_{\text{dual-bce}}
$$

$\mathcal{L}_{\text{dual-bce}}$ 仅对极少数双标签样本有效；预期在 0 重叠下 **不优于独立模型**。

### 3.4 Baselines

| 模型 | 用途 |
|------|------|
| XGBoost + ECFP4 | URAT1 强 baseline |
| Chemprop 单任务 | 图学习对照 |
| **PLK1-style pipeline** | SVR(URAT1) + 锚点相似性(NLRP3) + 0.5 融合 — **阴性对照** |
| Docking-only | 结构单独特证 |

---

## 四、Stage 2：双路径候选生成

### 4.1 Path A — 库筛漏斗

```
Input: Library L (~10⁶ SMILES, e.g. Enamine REAL)
  Step 1: RDKit valid + PAINS + Lipinski + reactive filters
  Step 2: URAT1 conformal filter (ŷ^U, lower bound)
  Step 3: NLRP3 P_active^ens ≥ t_N
  Step 4: Dual-target pre-score → top 10⁴
  Step 5: Ensemble docking (§5)
  Step 6: Reliability fusion (§6) + Butina diversity (Tc=0.4)
Output: C_A, |C_A| ~ 300–500
```

阈值 $t_U, t_N$ 在文献 benchmark 上校准（lesinurad, MCC950 分位数法）。

### 4.2 Path B — 生成式优化

**CLM 双靶偏置**（Schneider et al., Nat Commun 2024）：

1. 预训练：ChEMBL SMILES
2. $M_U$：URAT1 高活性 fine-tune；$M_N$：NLRP3 高活性 fine-tune
3. 采样：$\text{logits} = \alpha \text{logits}_U + (1-\alpha)\text{logits}_N$，$\alpha \in [0.4, 0.6]$

**RL 奖励**（POLYGON 思路，嵌入转运体项）：

$$
R(x) = \sum_i w_i r_i(x)
$$

| 组分 $r_i$ | 权重 $w_i$ | 说明 |
|-----------|-----------|------|
| $\hat{y}^U(x)$ | 0.18 | URAT1 回归模型 |
| $P_{\text{active}}^{\text{ens}}(x)$ | 0.18 | NLRP3 条件化分类 |
| $S_{\text{trap}}(x)$ | 0.14 | URAT1 系综（每 100 RL 步算一次） |
| $S_{\text{NLRP3}}^{\text{struct}}(x)$ | 0.14 | NLRP3 对接 |
| QED | 0.12 | 类药性 |
| $-\text{SA}$ | 0.10 | 合成惩罚 |
| Novelty | 0.08 | Tc < 0.85 vs ChEMBL |
| Conformal penalty | 0.06 | $w_U$ 过大则惩罚 |

**后处理**：有效 > 95%，去重，SA < 6，AiZynthFinder 至少 1 条路线

**输出**：$\mathcal{C}_B$，$|\mathcal{C}_B| \sim 500$–2000

### 4.3 候选池合并

```python
C_union = dedupe_by_smiles(C_A ∪ C_B)
for x in C_union:
    x.source = "library" | "generative" | "both"
```

---

## 五、Stage 3：结构约束评分

### 5.1 URAT1 构象系综与 $S_{\text{trap}}$

**机制**（Dai et al., Cell Res 2024）：配体稳定 inward/occluded 态，阻断 alternating access

**构象系综**：

| PDB | 态 | 权重 |
|-----|-----|------|
| 9B1H / 9DKB | inward-open + inhibitor | $w_{\text{in}}$ |
| 9JDZ | occluded | $w_{\text{occ}}$ |
| 9JDZ | outward-open | $w_{\text{out}}$ |

$$
S_{\text{trap}} = w_{\text{in}} s_{\text{in}} + w_{\text{occ}} s_{\text{occ}} - w_{\text{out}} s_{\text{out}}
$$

**关键残基** $S_{\text{key}}$：Phe-rich cage、**Arg477**（< 4 Å 接触比例）

$$
S_{\text{URAT1}}^{\text{struct}} = 0.30 S_{\text{trap}} + 0.25 S_{\text{vina}} + 0.25 S_{\text{key}} + 0.20 S_{\text{compete}}
$$

### 5.2 NLRP3 结构评分

NACHT Walker B 变构口袋（7ALV, 8ETR）；分子内胶机制

$$
S_{\text{NLRP3}}^{\text{struct}} = 0.30 S_{\text{vina}} + 0.35 \Delta G_{\text{MMGBSA}} + 0.35 S_{\text{stability}}^{\text{MD}}
$$

MD：50–100 ns，配体 RMSD < 2.5 Å，NACHT 亚域距离波动阈值

### 5.3 结构协同（非简单几何平均）

$$
S_{\text{struct}}^{\text{synergy}} = \sqrt{S_{\text{URAT1}}^{\text{struct}} \cdot S_{\text{NLRP3}}^{\text{struct}}} + \beta \cdot \min(\cdot)
$$

$\beta = 0.2$；当 $c_N < 0.5$ 时，提升结构项总权重（ML 不可靠）

---

## 六、Stage 4：可靠性加权融合 + Pareto

### 6.1 ML 可靠性权重

$$
\omega_U = \frac{1/w_U}{(1/w_U) + c_N}, \quad \omega_N = \frac{c_N}{(1/w_U) + c_N}
$$

### 6.2 综合融合分

$$
S_{\text{dual}}(x) = \omega_U \cdot \tilde{y}^U + \omega_N \cdot P_{\text{active}}^{\text{ens}} + \gamma \cdot S_{\text{struct}}^{\text{synergy}}
$$

$\tilde{y}^U$：归一化至 $[0,1]$；$\gamma = 0.35$

### 6.3 Pareto 排序

目标向量 $\mathbf{o}(x) = (S_{\text{dual}}, \text{QED}, -\text{SA}, \text{Nov})$  
取 rank 1–2 非支配解；同 rank 内按 $S_{\text{dual}}$ 降序

**对比 PLK1-style**（消融 Abl-3）：

$$
S_{\text{PLK1-style}} = 0.5 \tilde{s}_U + 0.5 \tilde{s}_N \quad \text{（固定等权，本框架明确弃用为主方法）}
$$

---

## 七、Stage 5：回顾性验证

### 7.1 Benchmark 回收

`data/benchmarks/literature_benchmarks.csv`：

- URAT1：lesinurad, benzbromarone, verinurad（must_recover）
- NLRP3：MCC950, GDC-2394, NT-0796
- 阴性：allopurinol（URAT1 非抑制剂，应排名靠后）

分别报告 **Path A only**、**Path B only**、**C_union** 的 Recall@100/@500

### 7.2 必做消融

见 `TAPE_GATE_FRAMEWORK.md` §七；额外报告：

- PLK1-style baseline vs TAPE-GATE：双靶均衡性（$\min$ 分位数）、benchmark 回收率
- 生成路径新颖性：与 ChEMBL 最大 Tc 分布

---

## 八、端到端伪代码

```python
# run_tape_gate_pipeline.py

# Stage 0
D_U, D_N = prepare_data(config)  # assay metadata preserved for NLRP3
splits = murcko_group_split(D_U, D_N, k=5)

# Stage 1 — asymmetric dual evidence
model_U = train_urat1_regressor(D_U, conformal=True, transfer="SLC22")
model_N = train_nlrp3_assay_conditioned(D_N, method="chemprop_cond")

# Stage 2 — dual paths
# Path A
library = load_library("enamine_real.smi")
C_A = library_screen(library, model_U, model_N, ensembles)

# Path B
clm = cross_finetune_clm(D_U, D_N)
C_B = rl_generate(clm, reward=build_dual_reward(model_U, model_N, ensembles))
C_union = merge_candidates(C_A, C_B)

# Stage 3–4
scored = reliability_pareto_rank(C_union, model_U, model_N, ensembles)

# Stage 5
report = retrospective_validate(scored, benchmarks)
ablation = run_ablations(include_plk1_style_baseline=True)
```

---

## 九、算法创新点（论文 Contribution）

1. **Transporter-aware $S_{\text{trap}}$** 与 PLK1 式激酶对接本质区分
2. **Assay-conditioned NLRP3 证据** 替代锚点相似性，应对 7.2% 跨 assay >1 log 离散（37/513）
3. **Paired-path discovery**：库筛 + 生成式并行，扩展双靶化学空间
4. **Reliability-weighted Pareto fusion** 替代固定 0.5/0.5 线性加权
5. **PLK1-style 阴性对照消融**：可证伪的差异化验证协议

---

## 十、实现参考

| 组件 | 资源 |
|------|------|
| MiniMol / Graphium | https://github.com/datamol-io/graphium |
| Conformal prediction | `mapie` 或自实现 split conformal |
| Chemprop v2 | https://github.com/chemprop/chemprop |
| AutoDock Vina | https://github.com/ccsb-scripps/AutoDock-Vina |
| CLAMP 思路 | 参考 assay-conditioned QSAR 文献 |

转运体验证要求见 [`URAT1_TRANSPORTER_VALIDATION.md`](URAT1_TRANSPORTER_VALIDATION.md)。
