# MASFL v3.1 完整流程（可证伪修订版）

> **MASFL** = **M**echanism-**A**ligned **S**tate and **F**unction **L**earning  
> **版本**: 3.1-Verifiable | **前身**: v3.0-Practical → TAPE-GATE v2.0  
> **定位**: 蛋白条件化多态蒸馏 + 去混杂功能学习 + 可行域 Pareto + **模块门槛降级**

---

## 零、问题定义与设计哲学

### 0.1 双靶目标（模式 A + 可行域）

ChEMBL 实测 **URAT1 / NLRP3 零重叠 SMILES**。主目标不是宣称单一分子双靶，而是：

1. 分别在 URAT1、NLRP3 上得到可信排序 $s_U, s_N$
2. 在 **可行域** $\mathcal{F}$ 内做 Pareto 推荐
3. 输出 **三级推荐**（双靶兼顾 / 单靶优选 / 仅讨论）

$$
\mathcal{F} = \{x:\; s_U(x) \geq t_U \;\wedge\; s_N(x) \geq t_N\}
$$

$$
\mathcal{C}^* = \arg\max_{x \in \mathcal{F}} \min(s_U(x), s_N(x)) \quad \text{（Pareto knee on Chebyshev）}
$$

阈值 $t_U, t_N$ 由 benchmark 分位数校准，**非**人工拍脑袋。

### 0.2 六条原则（回应 v3.0 第二轮质疑）

| # | 原则 | 回应的问题 |
|---|------|-----------|
| P1 | **算力分层** | Teacher 仅蒸馏集 + L2；Student 用于 L1 triage |
| P2 | **监督诚实** | π 为 surrogate index；MM-GBSA 不进 loss |
| P3 | **蛋白条件化** | Student 输入受体三态编码，破解信息瓶颈 |
| P4 | **蒸馏集扩充** | 5×10³–2×10⁴ 多源标注，含无活性负样本 |
| P5 | **不确定性拒绝** | 高 $u_{\text{epi}}$ 分子不凭 Student 淘汰 |
| P6 | **模块门槛降级** | 任一模块未过线 → 自动旁路，不硬撑主张 |

### 0.3 架构总览

```
Stage 0    数据、结构、蒸馏集分层构建
Stage 1    TC-Encoder（L0 粗筛 only）
Stage 2    Teacher M-CPDL（多任务标签，5k–20k 分子）
Stage 2′   PC-Student（蛋白条件化 + 排序蒸馏 + 集成不确定性）
Stage 2″   主动学习（可选，1–3 轮）
Stage 3    DFIM（去混杂功能抑制 + NLRP3 药效团 L1）
Stage 4    Path A 四层漏斗（L0→L1→L2→L3）
Stage 5    Path B 离线 RL（可选，默认降级）
Stage 6    LEF-Cal + 可行域 Pareto + 三级推荐
Stage 7    L2 Teacher 精修 + 结构 sanity
Stage 8    模块门槛评估 + 多 split 验证 + 算力-精度曲线
```

### 0.4 模块门槛总表（Gate）

| 模块 | 必过指标 | 未过线处置 |
|------|---------|-----------|
| Teacher M-CPDL | redock ≤2.0 Å；benchmark 方向 4/4 | 停用 CPDL，回 v2 $S_{\text{trap}}$ |
| PC-Student | scaffold-novel Spearman≥0.45；benchmark 方向 ≥3/4 | Student 仅弱筛；L1 阈值放宽；L2 全量 Teacher |
| DFIM | assay hold-out AUROC ≥0.75 | 降 FIM 权重；L1 改药效团主导 |
| LEF-Cal | OOF 等渗 Brier 改善 | 固定解析权重 |
| 整管线 | 回溯 EF@1% > 3×随机；URAT1 benchmark ≥3/4 | 收缩为「Teacher 晚阶段筛选」论文 |

---

## Stage 0：数据、结构与蒸馏集构建

### 0.1 生物活性数据

| 数据集 | 规模 | 操作 |
|--------|------|------|
| URAT1 | 822 | pIC50；Murcko GroupKFold($k=5$) |
| NLRP3 | 513 | IL-1β；粗粒度 $a_{\text{coarse}}$ |
| NLRP3 主训子集 | ~359 | THP-1 + IL-1β（DFIM 主域） |
| 重叠 | 0 | 独立双模型 |

```text
a_coarse = [assay_type, cell_line, readout_type]   # 5–15 类
n(assay_id) < 10 → "other"
```

### 0.2 结构数据与受体编码（PC-Student 用）

| 靶点 | PDB | 用途 |
|------|-----|------|
| URAT1 inward | 9B1H, 9DKB | Grid + 口袋图 $E_{\text{in}}$ |
| URAT1 occluded | 9JDZ | $E_{\text{occ}}$ |
| URAT1 outward | 9JDZ | $E_{\text{out}}$ |
| NLRP3 | 7ALV, 8ETR | 药效团模板 + L2 Grid |

**离线一次**：从三态结合位点构建口袋图（残基类型、坐标、部分电荷）→ `structures/pocket_graph_{in,occ,out}.pt`

### 0.3 蒸馏集 $\mathcal{D}_{\text{distill}}$（v3.1 核心扩充）

| 子集 | 规模 | 来源 | 活性标签 |
|------|------|------|---------|
| A 训练活性 | ~650/fold | URAT1 训练折 | pIC50 |
| B 骨架代表 | 500–1000 | Murcko 聚类中心 | 有 |
| C 化学空间采样 | 2000–5000 | ChEMBL SLC22 邻域 FPS | 有 |
| D 无活性多样性 | 3000–10000 | Enamine/ChEMBL 随机 | **无活性** |
| E 边界集 | ~200 | benchmark 类似物 | 有/参考 |

**总量目标**：5×10³ – 2×10⁴（一次性 Teacher 标注投资）

负样本（D）蒸馏目标：$\pi_{\text{in}}+\pi_{\text{occ}}$ 应低于活性中位数或方向随机——用于抑制 Student 假阳性。

### 0.4 NLRP3 药效团模板（L1 用）

从 7ALV + MCC950 共晶提取 Catalyst/Phase 药效团 → `structures/NLRP3_pharmacophore.hyp`

### 0.5 输出

```text
data/processed/urat1_curated.csv
data/processed/nlrbp3_curated.csv
data/distill/distill_manifest.csv          # 分层来源标记
data/splits/scaffold_fold_{0-4}.csv
structures/pocket_graph_{in,occ,out}.pt
structures/NLRP3_pharmacophore.hyp
```

---

## Stage 1：TC-Encoder（仅 L0）

### 1.1 模型

- Backbone：Chemprop（ChEMBL 预训练）或 MiniMol（冻结）
- 双任务头：URAT1 回归 + NLRP3 分类（$a_{\text{coarse}}$）
- **对比学习：关闭**（ablation only）

### 1.2 Conformal 不确定性 $w_U$

- 在 **全 URAT1 OOF 残差**（~822 点）上算 90% 区间
- **评价**：覆盖率是否 ≈90%；若平均宽度 >1.5 log → **$w_U$ 不用于路由**

### 1.3 角色边界

```text
TC-Encoder → L0 物理化学 + 活性粗筛
≠ URAT1 主排序
≠ NLRP3 主排序（DFIM 取代）
```

---

## Stage 2：Teacher M-CPDL（多任务构象教师）

> 对 $\mathcal{D}_{\text{distill}}$ 全量运行（5k–20k），不对百万库运行。

### 2.1 三态对接

```text
for x in D_distill:
  for s in {in, occ, out}:
    Glide SP → XP → best pose
    记录：GlideScore_s, pose_xyz
```

### 2.2 多任务 Teacher 标签（不止 π）

| 标签 | 符号 | 定义 |
|------|------|------|
| 构象偏好分布 | $\pi^T(s\|x)$ | softmax(Glide 归一化；out 取负) |
| Arg477 接触 | $c_{\text{Arg}}(x)$ | 任意态 inward pose 距离 <4 Å |
| Phe cage 接触 | $c_{\text{Phe}}(x)$ | 同上 |
| 方向符号 | $d(x)$ | sign($u_{\text{in}} - u_{\text{out}}$) |
| IFP 相似 | $f_{\text{IFP}}(x)$ | 与 lesinurad consensus IFP Tanimoto |

MM-GBSA **仅存档**，不进任何 loss。

### 2.3 Teacher 可信 Gate（必须先过）

| 检验 | 标准 |
|------|------|
| lesinurad redocking | RMSD ≤ 2.0 Å |
| 四药 inward>outward | **4/4** |
| 负样本集 D 方向 | $\pi_{\text{in}}+\pi_{\text{occ}}$ 中位数 < 活性集 |

**任一失败 → 禁止进入 Stage 2′，回退 v2 手工 $S_{\text{trap}}$。**

### 2.4 输出

```text
teacher/labels_fold{k}.csv    # π, c_Arg, c_Phe, d, f_IFP
teacher/poses/                # 可选
```

---

## Stage 2′：PC-Student（蛋白条件化排序蒸馏）

### 2′.1 模型架构

```
输入：
  配体图 L(x)  — Chemprop / GNN
  口袋图 E_s   — 预计算，s ∈ {in, occ, out}

交互：
  h_s = CrossAttention(L, E_s)   或  Ligand-Pocket GNN pair

输出：
  π_student(s|x) = softmax_s( w · h_s )
  ĉ_Arg, ĉ_Phe  — 辅助头
  S_π = π_in + π_occ - π_out
```

**禁止输入**：Glide score、MM-GBSA、对接 pose 坐标。

### 2.2 损失函数（排序为主，非绝对 KL）

$$
\mathcal{L} = \mathcal{L}_{\text{rank}} + \lambda_1 \mathcal{L}_{\text{contact}} + \lambda_2 \mathcal{L}_{\text{MIL}} + \lambda_3 \mathcal{L}_{\text{dir}}
$$

| 项 | 定义 |
|----|------|
| $\mathcal{L}_{\text{rank}}$ | pairwise：$\text{sign}(S^T_i - S^T_j)$ 与 $\hat{S}_i - \hat{S}_j$ 一致 |
| $\mathcal{L}_{\text{contact}}$ | BCE($\hat{c}_{\text{Arg}}, c^T_{\text{Arg}}$) + BCE($\hat{c}_{\text{Phe}}, c^T_{\text{Phe}}$) |
| $\mathcal{L}_{\text{MIL}}$ | pIC50 ~ weighted sum of $\hat{S}_s$（仅活性子集 A–C） |
| $\mathcal{L}_{\text{dir}}$ | hinge：benchmark 上 $\pi_{\text{in}} > \pi_{\text{out}}$ |

可选轻量 KL：仅作正则，权重 ≤0.1。

### 2.3 深度集成不确定性

- **5 个独立初始化** PC-Student → $\mu(S_\pi)$, $\sigma(S_\pi)$
- 认知不确定性：$u_{\text{epi}}(x) = \sigma(S_\pi)$

### 2.4 部署规则（必备）

```text
if u_epi(x) > τ_epi（验证集 90th 百分位）:
    不凭 Student 淘汰 x
    标记为 "uncertain" → 保留至 L2 或宽松通过 L1
else:
    可用 Student 强力筛除
```

### 2.5 PC-Student Gate

| 指标 | 通过线 | 失败处置 |
|------|--------|---------|
| Scaffold-novel Spearman($\hat{S}$, $S^T$) | ≥ 0.45 | Student 弱筛模式 |
| Benchmark 方向 | ≥ 3/4 | 禁止 Student 单独筛 URAT1 |
| EF@5% vs decoy | > 2.0 | L1 回退 Teacher 子采样 |
| PC vs 无蛋白 Student ablation | PC 显著优 | 否则改用药效团匹配 Student |

### 2.6 输出

```text
models/pc_student_ensemble/   # 5 checkpoints
student/u_epi_threshold.json
```

---

## Stage 2″：主动学习（可选，强烈建议）

```text
Round r = 1..R（R≤3）:
  1. PC-Student 对 Enamine 10^5 子集推理
  2. 选取：Top u_epi 的 500 + Top S_π 的 500 + 随机 500
  3. Teacher 标注 → 加入 D_distill
  4. 重训 PC-Student
  5. 检查 Gate 是否改善
```

**评价**：每轮 scaffold-novel Spearman 是否单调升；若不升 → 停止，避免过拟合蒸馏循环。

---

## Stage 3：DFIM（去混杂功能抑制 + 药效团 L1）

### 3.1 适用域声明

> DFIM 主训练域 = THP-1 + IL-1β + 细胞实验。  
> 其他 assay 仅 hold-out 外推报告，**不**用于主漏斗阈值。

### 3.2 去混杂模型

$$
P_\phi(\text{block} \mid x, a) = \sigma\Big( f(\phi(x), e_a) - \gamma \cdot \hat{T}(x) \Big)
$$

- $\hat{T}(x)$：cLogP、TPSA、预测 cytotox（pkCSM 或简易 RF）
- $\gamma$：交叉验证选取

### 3.3 NLRP3 L1 综合分（非纯 QSAR）

$$
s_N^{\text{L1}} = \beta_1 P_\phi + \beta_2 \text{Pharm}_{\text{MCC950}}(x) + \beta_3 \text{Shape}_{7\text{ALV}}(x)
$$

- $\text{Pharm}$：毫秒级药效团 FitValue
- $\beta$：在 assay hold-out CV 上选

### 3.4 DFIM Gate

| 指标 | 通过线 | 失败处置 |
|------|--------|---------|
| Scaffold CV AUROC | ≥ 0.85 | — |
| Assay hold-out AUROC | ≥ 0.75 | $\beta_1$ 下调 0.3；药效团权重升至 0.5 |
| MCC950 类似物 EF@1% | > 5 | — |
| colchicine 排名百分位 | < 20% | — |
| 去混杂 vs 未去混杂 hold-out | 提升或持平 | 否则去掉 $\hat{T}$ 项 |

### 3.5 输出

```text
models/dfim.pt
fim/oof_predictions.csv
fim/assay_holdout_report.json
fim/deconfound_ablation.json
```

---

## Stage 4：Path A 四层漏斗

### 4.1 漏斗（含不确定性路由）

```
L0  (~10⁶ → ~10⁵)
    RDKit + PAINS + Lipinski
    TC-Encoder：ŷ^U 下界 > 4.5（宽松）；非强淘汰

L1  (~10⁵ → ~10⁴)
    URAT1：PC-Student S_π 排序
           u_epi > τ → 保留不淘汰（uncertain 池）
           负样本方向约束：π_in+π_occ 过低可淘汰
    NLRP3：s_N^L1（DFIM + 药效团 + 形状）
    双靶各 Top 2×10⁴ 并集

L2  (~10⁴ → ~10²)   ← 权威决策层
    URAT1：Teacher M-CPDL 三态 XP（并集 + uncertain 池）
    NLRP3：Glide SP/XP 7ALV
    重算 s_U^L2, s_N^L2

L3  (~10² → 50–100)
    关键残基接触 + 可选 MM-GBSA
    Butina 多样性 Tc=0.4
```

### 4.2 L2 分数

$$
s_U^{\text{L2}} = 0.35 S_\pi^T + 0.30 c_{\text{Arg}} + 0.20 c_{\text{Phe}} + 0.15 S_{\text{Glide}}^{\text{in}}
$$

$$
s_N^{\text{L2}} = 0.40 P_\phi + 0.35 S_{\text{Glide}}^N + 0.15 \text{Pharm} + 0.10 S_{\text{key}}^{\text{NACHT}}
$$

### 4.3 算力-精度曲线（必报）

横轴：L2 对接分子数；纵轴：benchmark 回收率、EF@1%。  
证明 Student 的价值 = **在相同回收率下减少 L2 对接量**。

### 4.4 输出

```text
screening/L1_union.csv          # 含 u_epi, uncertain_flag
screening/L2_refined.csv
screening/cost_precision_curve.json
```

---

## Stage 5：Path B 离线 RL（可选，默认降级）

### 5.1 定位

- **默认**：Path B 为探索性模块，**不作**主 Contribution
- 仅当 PC-Student Gate 全过 + 主动学习 ≥1 轮后启用

### 5.2 离线 RL（非在线）

```text
数据集：D_distill 上已有 Teacher 标签的 10^4 分子
算法：保守 Q-learning / CQL 或行为克隆 + 奖励加权
奖励（缓存，无需实时对接）：
  R = 0.35·rank(S^T) + 0.25·S_π^student + 0.25·P_φ + 0.15·QED - 0.10·SA
```

### 5.3 硬约束

```text
必须通过：NLRP3 药效团 + SA<6 + Murcko Tc(已知活性) ∈ [0.2, 0.85]
禁止：在线 RL 纯 Student 奖励
```

### 5.4 RL 评价

| 指标 | 标准 |
|------|------|
| Teacher 通过率（L2 后 Top 10%） | > 30% |
| Student–Teacher Spearman（生成集） | > 0.4 |
| 奖励黑客检测 | Student Top100 与 Teacher Top100 重叠 > 20% |

---

## Stage 6：LEF-Cal + 可行域 Pareto + 三级推荐

### 6.1 等渗校准（500+ OOF 点，非 6 benchmark）

$$
\tilde{s}_U = \text{Isotonic}(\hat{y}^U_{\text{OOF}}), \quad
\tilde{s}_N = \text{Isotonic}(P_{\phi,\text{OOF}})
$$

**评价**：Brier score、reliability diagram 单调性。

### 6.2 可行域阈值

```text
t_U = percentile({s_U^L2(benchmark_active)}, 25)
t_N = percentile({P_φ(MCC950, GDC)}, 50)
```

### 6.3 可行域 Pareto + Chebyshev knee

在 $\mathcal{F}$ 内：

1. Pareto 非支配排序
2. 最大化 $\min(s_U^{\text{L2}}, s_N^{\text{L2}})$ 选 knee
3. 辅助 tie-break：$s_{\text{aux}} = \omega_U \tilde{s}_U + \omega_N \tilde{s}_N$

### 6.4 三级推荐

| Tier | 条件 | 论文角色 |
|------|------|---------|
| **Tier 1** | $x \in \mathcal{F}$ + Pareto knee | 主表 Top 20 |
| **Tier 2** | 单靶 $s > t$ 且另一靶 $\geq t \times 0.7$ | 备选 |
| **Tier 3** | Pareto 但一侧极弱 | 讨论 only，不推荐合成 |

### 6.5 Benchmark（回收检验，非训练）

| 化合物 | 要求 |
|--------|------|
| URAT1 四药 | L2 后 Top-500；≥3/4 |
| NLRP3 两药 | L2 后 Top-500；2/2 |
| allopurinol, colchicine | 百分位 < 20% |

---

## Stage 7：L2 结构精修

对 Tier 1 + Tier 2（~50–100）：

**URAT1**：Teacher XP、Arg477/Phe、可选 50 ns 膜 MD（RMSD <2.5 Å）  
**NLRP3**：7ALV XP、MM-GBSA、可选 MD  
**OCT1/OCT2**：对接比值 $R_{\text{sel}}$ — **Tier 3 参考**，表述为 computational hypothesis

---

## Stage 8：验证、消融与降级决策

### 8.1 必做验证

| 验证 | 报告 |
|------|------|
| Scaffold GroupKFold | 全模块 OOF |
| Temporal split | ChEMBL 年份 |
| Assay hold-out | DFIM 外推 |
| Scaffold-novel benchmark | lesinurad, dotinurad 单独 |
| Retrospective screen | EF@1%, EF@5%, BEDROC |
| Decoy | 性质匹配 |
| Cost-precision curve | Student 价值核心图 |

### 8.2 等预算消融

| ID | 对照 |
|----|------|
| Abl-1 | 无蛋白 Student vs PC-Student |
| Abl-2 | 排序蒸馏 vs KL 蒸馏 |
| Abl-3 | 蒸馏集 1k vs 10k |
| Abl-4 | 无负样本蒸馏 vs 有 |
| Abl-5 | DFIM vs 未去混杂 FIM |
| Abl-6 | 无约束 Pareto vs 可行域 Pareto |
| Abl-7 | v2 vs v3.1 |
| Abl-8 | PLK1-style baseline |

### 8.3 整管线降级树

```text
PC-Student Gate 失败
  → L1 全送 L2 Teacher（算力↑，主张收缩）
DFIM hold-out < 0.65
  → NLRP3 L1 药效团主导
Teacher Gate 失败
  → 回 v2 手工 S_trap
全部 Gate 失败
  → 发表定位为「双靶结构晚阶段筛选案例研究」，非学习框架
```

### 8.4 禁止表述

- ❌ Student 学到构象占据率 / 转移概率
- ❌ 百万库零对接可信排序（应为 triage + L2 权威）
- ❌ 6 分子校准融合网络
- ❌ OCT 对接 = 实验选择性
- ✅ protein-conditioned ranking distillation
- ✅ uncertainty-aware triage with late-stage teacher arbitration
- ✅ feasible-domain Pareto for dual-pathway prioritization

---

## 论文 Contribution（v3.1）

1. **M-CPDL + PC-Student**：多任务三态 Teacher + 蛋白条件化排序蒸馏 + 不确定性拒绝
2. **DFIM**：去混杂 assay 功能抑制 + 药效团 L1 结构锚定
3. **可行域 Pareto + 三级推荐**：药理学约束的双靶多目标决策
4. **模块门槛 + 算力-精度曲线**：可证伪、可降级的验证协议

---

## 与 v3.0 / v2 对照

| 维度 | v2.0 | v3.0 | **v3.1** |
|------|------|------|----------|
| 蒸馏集 | — | ~10³ | **5k–20k + 负样本** |
| Student | — | 盲蒸馏 2D | **PC + 排序蒸馏 + 集成** |
| Student 失败 | — | 无处置 | **Gate + 降级** |
| NLRP3 L1 | docking | 纯 FIM | **DFIM + 药效团** |
| Pareto | 无约束 | 无约束 | **可行域 + Tier** |
| RL | 在线 Student | 在线 | **离线 / 默认关闭** |
| 评价 | benchmark | +多 split | **+门槛 + cost-precision** |

---

## 执行脚本（规划）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
python scripts/00_prepare_data.py
python scripts/00b_build_distill_set.py           # Stage 0.3
python scripts/00c_encode_pocket_graphs.py        # Stage 0.2
python scripts/01_train_tc_encoder.py             # Stage 1
python scripts/02_teacher_mcpdl.py              # Stage 2（Glide）
python scripts/03_train_pc_student.py             # Stage 2′
python scripts/03b_active_learning.py           # Stage 2″（可选）
python scripts/04_train_dfim.py                   # Stage 3
python scripts/05_screen_path_a.py                # Stage 4
python scripts/06_offline_rl.py                   # Stage 5（可选）
python scripts/07_fusion_feasible_pareto.py       # Stage 6
python scripts/08_l2_validate.py                  # Stage 7
python scripts/09_module_gates_and_ablation.py    # Stage 8
```

---

## 相关文档

- [`TAPE_GATE_FRAMEWORK.md`](TAPE_GATE_FRAMEWORK.md) — v2.0
- [`URAT1_TRANSPORTER_VALIDATION.md`](URAT1_TRANSPORTER_VALIDATION.md)
- [`MODEL_QUALITY_REPORT.md`](MODEL_QUALITY_REPORT.md)
