# URAT1 三态对接 Benchmark 执行计划

> 配套战略文档：`PAPER_PIVOT_BENCHMARK.md`  
> 结构映射权威：`URAT1_THREE_STATE_DOCKING.md`  
> Gate 数据集：`TEACHER_GATE_QC_DATASETS.md`

---

## 1. 目标

在 **可承受的计算预算** 内，产出可支撑 JCIM / J. Cheminformatics benchmark 论文的：

1. 多协议对接结果表（CSV）
2. Gate 2/3 通过或失败的 **定量** 记录
3. 可开源的 SMILES + 参数 + 分数子集

---

## 2. 化合物面板

### 2.1 核心面板（必做）

| 组别 | 化合物 | 数量 | 文件 |
|------|--------|------|------|
| Gate redock | lesinurad | 1 | `teacher_gate_qc_panel_b_direction.csv` |
| Gate 方向 | lesinurad, benzbromarone, verinurad, dotinurad | 4 | 同上 |
| Benchmark 扩展 | `literature_benchmarks.csv` 中 URAT1 Tier1a | ~4–8 | `data/benchmarks/literature_benchmarks.csv` |

### 2.2 Decoy 集（Gate 3 + enrichment）

| 来源 | 数量 | 说明 |
|------|------|------|
| `distill_subset_d.csv` 随机子集 | 50–200 | seed=42，与活性理化性质匹配 |
| 阴性对照 | allopurinol 等 | `literature_benchmarks.csv` Tier_neg |

**不要** 对 8000 条 D 全量跑 B1K/B1L 刚性对接。

### 2.3 活性参照集（Gate 3 统计）

| 集合 | 文件 | 子集建议 |
|------|------|----------|
| A | `distill_subset_a.csv` | 随机 50–200 与 D 配对 |

---

## 3. 协议定义

### Protocol A — 单态 inward（基线）

```
配体 → Glide SP → XP @ 9DKB grid
输出: GlideScore_in, pose_in
```

- **用途**：传统激酶式单结构对接基线
- **预期**：4 药均可产生 pose（已验证）

### Protocol B — 刚性三态 Glide（失败文档化）

```
对每个 grid (9DKB, 9B1K, 9B1L):
  配体 → Glide SP → XP
输出: 每态 pose 数、GlideScore、失败原因
```

- **关键指标**：`pose_viability_rate = n_success / n_total`
- **已知结果**：药物在 9B1K/9B1L **0% pose**；urate 对照有 pose
- **论文价值**：证明 rigid transporter ensemble 对药物不可直接推广

### Protocol C — inward dock + pose 转移 + rescoring（Scheme 2，主推荐）

```
Step 0: 蛋白叠合（一次性）
  9B1K, 9B1L → align to 9DKB (backbone TM 区)

Step 1: 仅 9DKB
  Glide SP → XP → 取 top pose

Step 2: 复制 pose 至 B1K/B1L（在叠合坐标系下，不重叠合蛋白）

Step 3a (9B1K / occluded):
  Prime Minimize (ligand + binding site)
  MM-GBSA → ΔG_occ

Step 3b (9B1L / outward):
  Clash 检测 (vdW overlap)
  若严重 clash: ΔG_out = ΔG_occ + penalty (+5 kcal/mol 初值，敏感性 ±2)
  若可最小化: Prime → MM-GBSA → ΔG_out

Step 4: Boltzmann 权重
  π_s = exp(-ΔG_s / RT) / Σ_s exp(-ΔG_s / RT)
  S_π = π_in + π_occ - π_out
```

**注意**：

- π 基于 **rescored ΔG**，不是 raw GlideScore（B1K/B1L 无 Glide pose 时）
- 叠合 RMSD 大（B1L ~19 Å）必须在 Discussion 承认

### Protocol D — IFD 敏感性（可选，仅 4 药）

```
Glide IFD @ 9B1K（柔性口袋）
```

- **用途**：给出「若允许诱导契合，occluded 能否结合」的上界
- **规模**：4 药 × 1 grid，计算可承受

---

## 4. 评分与 Gate 判据

### 4.1 Boltzmann $\pi$（Protocol C）

```python
import numpy as np

RT = 0.593  # kcal/mol at 298 K

def boltzmann_pi(delta_g: dict[str, float]) -> dict[str, float]:
    """delta_g keys: inward, occluded, outward"""
    states = list(delta_g.keys())
    dg = np.array([delta_g[s] for s in states])
    dg_min = dg.min()
    w = np.exp(-(dg - dg_min) / RT)
    pi = w / w.sum()
    return {s: float(p) for s, p in zip(states, pi)}

def s_pi(pi: dict[str, float]) -> float:
    return pi["inward"] + pi["occluded"] - pi["outward"]
```

实现文件：`scripts/utils_three_state_scoring.py`（待建）

### 4.2 Gate 判据

| Gate | 条件 | 失败时行动 |
|------|------|------------|
| G1 | lesinurad redock RMSD ≤ 2.0 Å @ 9DKB | 检查 grid / 质子化 |
| G2 | 4/4 药 $S_\pi > 0$ under Protocol C | 调 penalty、结合位点或改 rescoring |
| G3 | median($\pi_{in}+\pi_{oc}$)_A > median(...)_D | 扩大 decoy 或改特征；仍失败则只报 A vs C |

---

## 5. 评估指标

| 指标 | 定义 | 用途 |
|------|------|------|
| Pose viability | 有有效 pose 的分子比例 | Protocol B 主图 |
| Redock RMSD | 共晶配体 vs docked pose | Gate 1 |
| $S_\pi$ sign accuracy | 已知抑制剂 $S_\pi>0$ 比例 | Gate 2 |
| Enrichment@k | 活性在 top-k% 占比 vs random | A vs C |
| ROC-AUC | 活性 vs decoy 分类 | decoy 集 |
| Rank of actives | 四药在 decoy 池中的排名 | 主表 |

统计要求（jcim.5c01609）：

- 报告效应量（如 Cliff's delta、Cohen's d）与 95% CI
- Scaffold split 若做 ML 对比（本篇可只做对接 rank）

---

## 6. 计算任务量估算

| 任务 | 分子数 | Grid | 估计 job |
|------|--------|------|----------|
| Gate A+B+C（4 药） | 4 | 3 + rescoring | ~12 Glide + 8 Prime/MM-GBSA |
| Benchmark 扩展 A+C | 8 | 同上 | ~24 Glide + rescoring |
| Protocol B 失败记录 | 4 | 3 | 12（已知大部分失败） |
| Decoy C only | 100 | 1 dock + 2 rescoring | ~100 Glide + 200 rescoring |
| IFD (D) | 4 | 1 | 4 IFD |

**总计**：远低于 26,466 全量三态 Glide。

---

## 7. 输出文件规范

建议目录：`results/three_state_benchmark/`

```
results/three_state_benchmark/
├── protocol_a_scores.csv
├── protocol_b_viability.csv      # 含 failure_reason 列
├── protocol_c_rescored.csv       # ΔG_*, pi_*, S_pi
├── gate_summary.json
├── enrichment_curves.png         # 作图脚本输出
└── README.md                     # 复现步骤
```

### CSV 列（protocol_c_rescored.csv）

| 列 | 说明 |
|----|------|
| `compound_id` | 内部 ID |
| `smiles` | |
| `protocol` | A / B / C / D |
| `dg_inward` | kcal/mol |
| `dg_occluded` | |
| `dg_outward` | |
| `pi_inward` | |
| `pi_occluded` | |
| `pi_outward` | |
| `s_pi` | |
| `clash_outward` | bool |
| `penalty_applied` | kcal/mol or NA |
| `notes` | log 摘要 |

---

## 8. Maestro 操作要点（你已遇到的问题）

### 8.1 叠合参考配体

- Superposition 的 reference 必须在 **Project Table 中选中**
- 先 align 蛋白 backbone，再复制配体 pose

### 8.2 B1K/B1L 刚性 Glide 失败

- Log 关键词：`GRID-ENERGY MIN FAILED`, `GlideScore=10000`
- **不是** grid 建错（urate 能 dock）— 是 **药物构象与刚性口袋不兼容**
- 此结果应写入 `protocol_b_viability.csv`，是论文 **核心发现之一**

### 8.3 SP → XP 流程

- 最终 $\pi$ 若用 GlideScore（仅 Protocol A/B 成功态）：取 **XP GlideScore**
- Protocol C：取 MM-GBSA $\Delta G$

---

## 9. 时间线（技术阶段，非日历）

```
[Done] PDB 映射校正 (9DKB/9B1K/9B1L)
[Done] Gate 数据集与 manifest
[Done] 9DKB 四药对接
[Done] B1K/B1L 刚性失败记录
[Now]  Protocol C pose 转移 + MM-GBSA
[Next] Gate 2 四药 S_π
[Next] Decoy 子集 + Gate 3
[Next] Protocol A vs C enrichment 表
[Next] 初稿图表
```

---

## 10. 参考文献（Methods 必引）

1. Dai Y, Lee CH. *Cell Res* 2024 — 9B1K, 9B1L  
2. Suo Y et al. *Nat Commun* 2025 — 9DKB  
3. Wu C et al. *Cell Discov* 2025 — 9JDZ 局限  
4. Sindt et al. JCIM 2025 (jcim.5c00730) — rescoring 局限  
5. Practically Significant Method Comparison, JCIM 2025 (jcim.5c01609) — 统计报告  
6. Burns et al. 2016 — 细胞 IC50 benchmark  

---

## 11. 与全量 MASFL 的关系

| 条件 | 行动 |
|------|------|
| Gate 2 < 4/4 | 不启动 Teacher；论文止于 benchmark |
| Gate 2 = 4/4, Gate 3 通过 | 可对 **子集**（如 200+200）试 Teacher 标签，写入 **future work** |
| Gate 2 = 4/4, Gate 3 失败 | 发表 benchmark；讨论 $S_\pi$ 对 URAT1 判别力不足 |

全量 8973 × 3 仅在 **Protocol C 被验证** 且 **有集群自动化** 后考虑。
