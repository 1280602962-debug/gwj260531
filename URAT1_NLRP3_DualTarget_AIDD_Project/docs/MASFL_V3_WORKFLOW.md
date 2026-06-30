# MASFL v3.0 完整流程（修订版）

> **MASFL** = **M**echanism-**A**ligned **S**tate and **F**unction **L**earning  
> **版本**: 3.0-Practical | **前身**: TAPE-GATE v2.0  
> **定位**: 昂贵结构证据蒸馏 + 表型功能条件化 + 可校准双靶融合（非端到端对接学习）

---

## 零、问题定义与核心原则

### 0.1 双靶目标（模式 A，与当前数据一致）

ChEMBL 实测 **URAT1 / NLRP3 零重叠 SMILES**，故不宣称发现「单一分子同时双靶」为主目标，而采用：

$$
\max_{x \in \mathcal{C}} \;\big(s_U(x),\; s_N(x)\big) \quad \text{via Pareto 非支配排序}
$$

- $s_U$：URAT1 构象偏好 + 摄取抑制证据
- $s_N$：NLRP3 功能抑制 + 结构锚定证据
- 最终推荐：Pareto 前沿上的 knee point 或高 hypervolume 分子

### 0.2 四条工程原则（回应 v3.0 初稿质疑）

| 原则 | 含义 |
|------|------|
| **算力分层** | 对接仅 Teacher / L2；库筛与 RL 用 Student 代理 |
| **监督诚实** | MM-GBSA 不作 regression label；π 为 surrogate preference index |
| **Assay 粗粒度** | FIM 用 assay_type + cell_line，不用 92 维 assay_id 主嵌入 |
| **校准非过拟合** | 6 个 benchmark 仅回收检验；融合参数由 OOF CV 选 |

### 0.3 架构总览

```
Stage 0   数据与结构准备
Stage 1   TC-Encoder（共享分子表征 + 双任务头）
Stage 2   Teacher CPDL（离线，三态对接 → π_teacher）  ─┐
Stage 2′  Student CPDL（蒸馏，2D/3D → π_student）     ─┤ URAT1 臂
Stage 3   FIM-Lite（粗粒度 assay 功能抑制）            ─┤ NLRP3 臂
Stage 4   Path A 分层库筛（L0→L1→L2）                  │
Stage 5   Path B 生成式 RL（Student 奖励 + Teacher 校正）│
Stage 6   LEF-Cal 融合 + Pareto 排序                   │
Stage 7   L2 Teacher 精修（Top 1%）+ 结构 sanity       │
Stage 8   多 split 回溯验证 + 消融                        ┘
```

---

## Stage 0：数据与结构准备

### 0.1 生物活性数据

| 数据集 | 规模 | 操作 |
|--------|------|------|
| URAT1 | 822 SMILES | pIC50 回归；Murcko GroupKFold($k=5$) |
| NLRP3 | 503 SMILES（IL-1β + Assay B） | 二值/概率分类；**粗粒度 assay 编码** |
| 重叠 | 0 | 独立双模型，MTL 仅消融 |
| SLC22 辅助 | OCT1/OCT2 | TC-Encoder 预训练迁移（可选） |

**NLRP3 assay 编码（FIM 用）**：

```text
a_coarse = [assay_type, cell_line, readout_type]   # 5–15 类，非 92 维 assay_id
assay_id 样本 <10 → 合并为 "other"
主训练子集：THP-1 + IL-1β（~359）可选
```

### 0.2 结构数据

| 靶点 | PDB | 角色 |
|------|-----|------|
| URAT1 inward | 9B1H, 9DKB | Teacher Grid $s_{\text{in}}$ |
| URAT1 occluded | 9JDZ | Teacher Grid $s_{\text{occ}}$ |
| URAT1 outward | 9JDZ | Teacher Grid $s_{\text{out}}$（减权） |
| NLRP3 NACHT | 7ALV（主）, 8ETR（辅） | L2 结构锚定 |

**预处理**：Schrödinger Protein Prep Wizard 或等效流程；lesinurad redocking RMSD ≤ 2.0 Å 为 Teacher 可信前提。

### 0.3 输出

```text
data/processed/urat1_curated.csv
data/processed/nlrbp3_curated.csv          # 含 a_coarse 列
data/splits/scaffold_fold_{0-4}.csv
structures/URAT1_{in,occ,out}_grid.zip
structures/NLRP3_7ALV_grid.zip
```

---

## Stage 1：TC-Encoder（共享分子表征）

### 1.1 模型

- **Backbone**：MiniMol（冻结）或 Chemprop D-MPNN（ChEMBL 预训练）
- **双任务头**：URAT1 回归 $\hat{y}^U$；NLRP3 分类 $\hat{p}^N_{\text{raw}}$
- **靶点条件化**：target embedding 拼接至 head 输入
- **对比学习**：默认 **关闭**；仅作 ablation（小样本下不稳）

### 1.2 训练

```text
L_TC = L_reg(URAT1) + L_bce(NLRP3, a_coarse)
CV：scaffold GroupKFold
URAT1 附加：conformal 残差 → 区间宽度 w_U
```

### 1.3 输出

- `models/tc_encoder.pt`
- OOF 预测 → 供 LEF-Cal 校准
- **用途**：L0 快速过滤；**不**作为 URAT1 主排序依据

---

## Stage 2：Teacher CPDL（离线，昂贵）

> 仅对 **训练集 + benchmark + 骨架代表集**（约 10³ 分子量级）运行，不对全库运行。

### 2.1 三态对接（Teacher）

对每个配体 $x$：

```text
for s in {inward, occluded, outward}:
    Glide SP → XP（Top 10%）→ 取最佳 pose
    记录：GlideScore_s, 关键接触(Arg477, Phe cage)
```

**不**将 MM-GBSA 写入训练 label；MM-GBSA 仅存档供分析。

### 2.2 Teacher 偏好分布

手工初始化（可学习微调仅在有方向约束的子集上）：

$$
\pi_{\text{teacher}}(s \mid x) = \text{softmax}_s\big(\beta \cdot u_s(x)\big)
$$

其中 $u_{\text{in}}, u_{\text{occ}}$ 为 Glide 归一化分，$u_{\text{out}}$ 取负（抑制 outward 偏好）。

**方向监督**（仅 benchmark + 共晶，≤10 分子）：

```text
L_dir = Σ max(0, π_teacher(out|x) - π_teacher(in|x) + margin)   # hinge
```

### 2.3 输出

```text
teacher/π_teacher_{fold}.csv      # SMILES, π_in, π_occ, π_out
teacher/poses_{fold}/               # 可选，供 PLIG 分析
```

---

## Stage 2′：Student CPDL（在线，快速代理）

### 2′.1 模型

- **输入**：2D 分子图 / RDKit 3D 构象 / MiniMol 嵌入（**禁止**输入 Glide score）
- **输出**：$\pi_{\text{student}}(s \mid x)$，及标量 $S_{\pi} = \pi_{\text{in}} + \pi_{\text{occ}} - \pi_{\text{out}}$

### 2′.2 训练损失

$$
\mathcal{L}_{\text{student}} = \underbrace{\text{KL}(\pi_{\text{student}} \| \pi_{\text{teacher}})}_{\text{蒸馏}} + \underbrace{\mathcal{L}_{\text{MIL}}(\sum_s \pi(s) f_s,\; y^U)}_{\text{弱监督}} + \lambda_{\text{dir}} \mathcal{L}_{\text{dir}}
$$

- $\mathcal{L}_{\text{MIL}}$：$pIC50$ 与偏好加权对接分的多实例回归
- Teacher 标签按 fold **严格隔离**（当前 fold 的 test 集 Teacher 标签不参与 Student 训练）

### 2′.3 推理成本

| 操作 | 耗时量级 |
|------|---------|
| Student 前向 | ~1 ms/分子 |
| Teacher 三态对接 | ~1–5 min/分子 |

**库筛 / RL 全程使用 Student**；仅 L2 对 Top 1% 调用 Teacher。

### 2′.4 输出

```text
models/cpdl_student.pt
student/π_student_scores.csv      # 全库可批推理
```

---

## Stage 3：FIM-Lite（NLRP3 功能抑制）

### 3.1 模型

$$
P_{\phi}(\text{block activation} \mid x, a_{\text{coarse}}) = \sigma\big(\text{MLP}([\phi(x); \mathbf{e}_{a_{\text{coarse}}}])\big)
$$

- $\phi(x)$：TC-Encoder 共享 backbone
- **不用** fine-grained assay_id 作主特征

### 3.2 结构锚定（非主排序）

NLRP3 综合分（筛库用）：

$$
s_N(x) = \alpha \cdot P_{\phi}(\text{block} \mid x, a^*) + (1-\alpha) \cdot \hat{p}^N_{\text{raw}}
$$

- $a^*$：推理时固定为 THP-1 + IL-1β 条件
- $\alpha$：在 **assay hold-out CV** 上选，默认 0.7–0.9（功能为主）
- **Docking 不参与 L0–L1**；L2 对 Top 候选算 $S_{\text{dock}}^N$

### 3.3 去偏

- PAINS / 反应性过滤
- 可选细胞毒性先验降权
- 泛抑制剂（colchicine 等）纳入阴性 benchmark

### 3.4 输出

```text
models/fim_lite.pt
fim/oof_predictions.csv
fim/assay_holdout_report.json
```

---

## Stage 4：Path A 分层库筛

### 4.1 漏斗定义

```
L0  (~10⁶)  RDKit valid + PAINS + Lipinski
              + TC-Encoder 粗筛（ŷ^U 下界、P_φ 阈值）
              → ~10⁵

L1  (~10⁵→10⁴)  Student CPDL：S_π 排序 + 不确定性边界过采样
              + FIM：P_φ 排序
              双靶预分：各自 Top 2×10⁴ 取并集 → ~10⁴

L2  (~10⁴→10²)  Teacher CPDL：三态 Glide XP（仅并集）
              + NLRP3 Glide SP（7ALV Grid）
              重算 s_U^L2, s_N^L2

L3  (~10²→50)   MM-GBSA（可选）+ 关键残基接触检查
              + Butina 多样性 (Tc=0.4)
```

### 4.2 URAT1 L2 分数

$$
s_U^{\text{L2}} = 0.4\, S_{\pi}^{\text{teacher}} + 0.3\, S_{\text{key}} + 0.3\, S_{\text{Glide}}^{\text{in}}
$$

### 4.3 NLRP3 L2 分数

$$
s_N^{\text{L2}} = 0.5\, P_{\phi} + 0.3\, S_{\text{Glide}}^N + 0.2\, S_{\text{key}}^{\text{NACHT}}
$$

### 4.4 输出

```text
screening/L0_pass.smi
screening/L1_union_top10k.csv
screening/L2_refined_top500.csv
```

---

## Stage 5：Path B 生成式 RL（可选）

### 5.1 CLM 双靶偏置

```text
预训练 ChEMBL SMILES
→ M_U：URAT1 高活性 fine-tune
→ M_N：NLRP3 高活性 fine-tune
采样：logits = α·logits_U + (1-α)·logits_N
```

### 5.2 RL 奖励（全 Student，无对接）

每步生成分子 $x$：

$$
R(x) = w_1 \hat{y}^U(x) + w_2 P_{\phi}(x) + w_3 \Delta\pi_{\text{student}}(x) + w_4 \text{QED} - w_5 \text{SA}
$$

| 组分 | 来源 | 耗时 |
|------|------|------|
| $\Delta\pi_{\text{student}}$ | CPDL Student | ms |
| $\hat{y}^U, P_{\phi}$ | TC-Encoder / FIM | ms |
| Teacher 校正 | 每 500 步对 Top-50 | 小时级 |

### 5.3 后处理

```text
有效 SMILES >95%
SA < 6
去重 + 与 Path A 合并 → C_union
```

---

## Stage 6：LEF-Cal 融合 + Pareto 排序

### 6.1 不用大 MLP 在 6 个 benchmark 上训练

**LEF-Cal** = 解析式可靠性路由 + OOF 等渗校准：

$$
\tilde{s}_U = \text{Isotonic}(\hat{y}^U_{\text{OOF}}), \quad
\tilde{s}_N = \text{Isotonic}(P_{\phi,\text{OOF}})
$$

$$
\omega_U = \frac{1/w_U}{(1/w_U) + c_N}, \quad
\omega_N = \frac{c_N}{(1/w_U) + c_N}
$$

### 6.2 双靶融合（模式 A）

**不**强行标量融合为双靶分数，而：

1. 各臂输出 $(s_U^{\text{L2}}, s_N^{\text{L2}})$
2. **Pareto 非支配排序**
3. Knee point 选取（或 hypervolume 最大）

可选标量辅助（仅排序 tie-break）：

$$
s_{\text{aux}} = \omega_U \tilde{s}_U + \omega_N \tilde{s}_N
$$

### 6.3 Benchmark 角色（非训练集）

| 化合物 | 检验 |
|--------|------|
| lesinurad, benzbromarone, verinurad, dotinurad | URAT1 Top-500 回收 |
| MCC950, GDC-2394 | NLRP3 Top-500 回收 |
| allopurinol, colchicine | 阴性：排名应低 |

### 6.4 输出

```text
fusion/pareto_front.csv
fusion/top50_recommendations.csv
fusion/benchmark_recovery.json
```

---

## Stage 7：L2 Teacher 精修与结构 Sanity

对 Pareto 前沿 Top 50–100：

### 7.1 URAT1

- [ ] Teacher 三态 XP 重对接
- [ ] Arg477 + Phe cage 接触 < 4 Å
- [ ] lesinurad redocking 对照 RMSD
- [ ] 可选：50 ns Desmond MD（膜嵌入），配体 RMSD < 2.5 Å

### 7.2 NLRP3

- [ ] 7ALV Glide XP + MCC950 pose 对照
- [ ] MM-GBSA Top pose
- [ ] 可选：50 ns MD，NACHT 亚域距离稳定

### 7.3 选择性

- [ ] OCT1/OCT2 平行对接（SLC22 脱靶）

### 7.4 输出

```text
final/top_candidates_structurally_validated.csv
final/interaction_diagrams/
```

---

## Stage 8：验证与消融协议

### 8.1 必做验证

| 验证 | 方法 | 通过标准 |
|------|------|---------|
| Scaffold split | GroupKFold | OOF R²/AUROC 报告 |
| Temporal split | ChEMBL 年份切分 | 新化合物泛化 |
| Assay hold-out | 留出一类 cell_line | FIM AUROC 不降 >10% |
| Retrospective screen | ChEMBL 全库或 Enamine 10⁵ | EF@1%, EF@5%, BEDROC |
| Decoy | 性质匹配 decoys | 活性/诱饵 AUC |
| Benchmark 回收 | 6+2 药 | URAT1 ≥3/4；NLRP3 2/2 |
| Scaffold-novel 子集 | lesinurad, dotinurad 排除训练 | 单独报告 |

### 8.2 等计算预算消融

| ID | 对照 | 验证点 |
|----|------|--------|
| Abl-1 | Hand $S_{\text{trap}}$ vs Student $\pi$ | 零对接推理回收率 |
| Abl-2 | Student only vs Student+Teacher L2 | 晚阶段增益 |
| Abl-3 | FIM coarse vs assay_id fine | 过拟合检验 |
| Abl-4 | FIM-only vs FIM+docking L2 | 机制锚定价值 |
| Abl-5 | 固定权重 vs LEF-Cal | 融合校准 |
| Abl-6 | TAPE-GATE v2 vs MASFL v3 | 整体 EF / benchmark |
| Abl-7 | PLK1-style baseline | 阴性对照 |

### 8.3 禁止的表述

- ❌ 「学到真实构象占据率 / 转移概率」
- ❌ 「MM-GBSA 监督的可学习能量函数」
- ❌ 「92 assay 全空间插值」
- ❌ 「6 个药物训练的融合神经网络」
- ✅ 「surrogate conformational preference index」
- ✅ 「assay-coarse functional inhibition proxy」
- ✅ 「distilled multi-state docking evidence」

---

## 软件与算力分工

| 阶段 | 工具 | 规模 |
|------|------|------|
| 结构准备 / Teacher 对接 | Schrödinger Glide XP | 10³ |
| Student 训练 / 融合 | Python, Chemprop, XGBoost | CPU/GPU |
| L0–L1 库筛 | RDKit + PyTorch 批推理 | GPU |
| NLRP3 L2 对接 | Glide SP/XP | 10⁴ |
| MD（可选） | Desmond / GROMACS | Top 50 |
| 可视化 | Maestro, DS | 论文图 |

详见 [`COMPUTE_REQUIREMENTS.md`](COMPUTE_REQUIREMENTS.md)。

---

## 论文 Contribution 表述（修订后）

1. **T-S CPDL**：将三态对接证据蒸馏为可扩展构象偏好代理，实现零对接高通量筛库与 RL
2. **FIM-Lite**：粗粒度 assay 条件化功能抑制学习，结构对接作晚阶段锚定
3. **LEF-Cal + Pareto**：可校准双靶证据路由与多目标推荐，避免小样本融合过拟合
4. **系统验证**：scaffold / temporal / assay-holdout / retrospective EF 协议

---

## 与 TAPE-GATE v2 对照

| 维度 | v2.0 | v3.0-Practical |
|------|------|----------------|
| URAT1 主证据 | 手工 $S_{\text{trap}}$，每分子三态对接 | Student $\pi$ 筛库 + Teacher L2 精修 |
| NLRP3 | Assay-conditioned + docking 并行 | FIM-Lite 主 + docking L2 锚定 |
| 融合 | 固定经验权重 | LEF-Cal + Pareto |
| RL | $S_{\text{trap}}$ 每步对接 | Student $\Delta\pi$ + 周期 Teacher |
| 验证 | Benchmark 2/4 URAT1 | 多 split + 回溯 EF + 等预算消融 |
| 主张 | 转运体感知流水线 | 机制对齐蒸馏学习框架 |

---

## 一键执行（规划）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
python scripts/00_prepare_data.py
python scripts/01_dataset_analysis.py
python scripts/02_train_tc_encoder.py          # Stage 1
python scripts/03_teacher_cpdl.py              # Stage 2（需 Glide）
python scripts/04_train_cpdl_student.py        # Stage 2′
python scripts/05_train_fim_lite.py            # Stage 3
python scripts/06_screen_path_a.py             # Stage 4
python scripts/07_generate_path_b.py           # Stage 5（可选）
python scripts/08_fusion_pareto.py             # Stage 6
python scripts/09_l2_validate.py               # Stage 7
python scripts/10_ablation_benchmark.py        # Stage 8
```

> 注：Stage 2–7 部分脚本待实现；当前仓库已实现 v2 子集，见 `run_tape_gate_pipeline.py`。

---

## 相关文档

- [`TAPE_GATE_FRAMEWORK.md`](TAPE_GATE_FRAMEWORK.md) — v2.0 总览
- [`ALGORITHM_FRAMEWORK.md`](ALGORITHM_FRAMEWORK.md) — v2 算法细节
- [`URAT1_TRANSPORTER_VALIDATION.md`](URAT1_TRANSPORTER_VALIDATION.md) — 结构验证标准
- [`MODEL_QUALITY_REPORT.md`](MODEL_QUALITY_REPORT.md) — 当前 ML 评估（URAT1_NO_GO）
