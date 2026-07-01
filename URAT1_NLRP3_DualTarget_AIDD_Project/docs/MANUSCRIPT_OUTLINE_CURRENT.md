# 论文定稿思路（2026-07 修订）

> **首投**：*Journal of Molecular Modeling*  
> **定位**：痛风 URAT1–NLRP3 双节点下的 **临床药物重定位 + 不对称双证据漏斗**  
> **不是**：双靶新药发现、8973 百万虚筛、Teacher/OAT 创新

---

## 0. 一句话

在全临床药物库上用 **NLRP3 ML 快筛** → **P(active)≥0.5（n≈1588）双靶对接** → **Pareto 短名单**；用 **8973 仅证明 URAT1 应对接**；用 **代表药 MD（2+2）** 解释机制。

---

## 1. 三套数据、三个用途（禁止混用）

| 数据集 | n | 用途 | 是否双靶筛选 |
|--------|---|------|-------------|
| **ChEMBL 临床药物 manifest** | 8319（III期+上市 1283） | **主筛选库**：NLRP3 ML → 对接漏斗 → Pareto | ✅ |
| **8973 distill** | 8928 docked | **仅 URAT1 回顾**：A vs D 富集、四药百分位 | ❌ 不做 NLRP3 ML |
| **Benchmark 六药 + MD** | 4+2 | 方法验证与机制图 | 平行表征，非筛库 |

---

## 2. 计算流程（你采纳的新流程）

```
ChEMBL 临床药物库 (n=8319)
        │
        ▼  【Step 1】NLRP3 assay-conditioned ML（~20 s）
        │   → nlrp3_ml_scores_clinical_all.csv
        ▼
        │  【Step 2】P(active) ≥ 0.5 → n≈1588（对接池）
        │   → URAT1 @ 9DKB XP + NLRP3 @ 8ETR XP
        ▼
        │  【Step 3】整合
        │   • 主图：P(S_U^dock) vs P(S_N^ML) 或双对接 Pareto
        │   • 标注：lesinurad、colchicine、MCC950 等对照
        ▼
        │  【Step 4】MD（主文 2+2，非全库）
        │   URAT1：benzbromarone + dotinurad（或 +lesinurad SI）
        │   NLRP3：MCC950 + GDC-2394
        ▼
   计算重定位短名单 + 机制讨论（无湿实验）
```

**与经典「先 ML 再对接」**：NLRP3 ML 预筛 → 命中分子双靶结构对接（对标 FDA 2697 全库对接，但先用 ML 缩至 ~1588）。

---

## 3. NLRP3 ML 筛选（已跑通）

```bash
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto
```

| 指标 | 数值 |
|------|------|
| 全库打分 | **8319** |
| **P(active)≥0.5（对接池）** | **1588** |
| max_phase≥3 子集（SI） | 247 |
| Top 5%（可选对照） | 416 |

输出：`results/repurposing/docking_pool_p05.csv`

**已知对照药排名（全库）**：

| 药物 | NLRP3 rank | P(active) | 备注 |
|------|------------|-----------|------|
| verinurad | 595 | 0.92 | 高（训练集内药） |
| colchicine | 1038 | 0.92 | 炎症相关，非直接 NLRP3 |
| lesinurad | 1513 | 0.60 | 中等 |
| benzbromarone / dotinurad / allopurinol / febuxostat | >4400 | ~0 | 预期偏低（非 NLRP3 药） |

**写作注意**：Top10 多为 **Phase 1–2 早期候选**，对接前建议加一层：

```bash
# 可选：仅 III 期+上市子集再筛
python3 scripts/screen_repurposing_library.py --panel phase_ge3 ...
```

或在 **P≥0.5 对接池** 中 **过滤 max_phase≥3**（247 条）作 SI 敏感性分析。

---

## 4. 文章结构（五节 Results）

### R1. 数据不对称与方法设计
- ChEMBL 0 重叠；URAT1 ML 2/4 vs NLRP3 AUROC 0.89  
- 流程图：**NLRP3 ML 筛临床库 → 双靶对接 → Pareto**

### R2. NLRP3 ML 筛选临床药物库
- 8319 全库；**1588** 条 P≥0.5 进入对接池  
- 对照药排名表  

### R4. 对接池双靶对接 + Pareto（主图）
- **1588** @ 9DKB + 8ETR（SI：phase≥3 子集 247）  
- $S_U$–$S_N$ 平面；已知药标注  
- SI：NLRP3 对接 vs ML 在对接池上的 Spearman

### R5. 代表药 MD（2+2）
- 结合稳定性；**不声称**全库 hit 验证

---

## 5. 主图方案（6 张）

| 图 | 内容 |
|----|------|
| Fig 1 | 双节点 + **ML→双对接→Pareto** 流程 |
| Fig 2 | NLRP3 ML 全库分数分布 + 对照药 |
| Fig 3 | 8973 URAT1 回顾富集（与 R3 绑定） |
| Fig 4 | 短名单双靶对接 Pareto（**主贡献**） |
| Fig 5 | URAT1 代表药 pose + MD |
| Fig 6 | NLRP3 代表药 pose + MD |

---

## 6. 本地 Agent 下一步命令

```bash
# 1. NLRP3 ML 全库 + 导出 P≥0.5 对接池（若未跑）
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto
# → results/repurposing/docking_pool_p05.csv（n≈1588）

# 2. Maestro：对 docking_pool_p05 跑 URAT1 9DKB XP + NLRP3 8ETR XP
# 导出至 results/repurposing/docking_raw/

# 3. Pareto 整合
python3 scripts/merge_docking_pareto.py \
  --ml-scores results/repurposing/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_8etr_p05.csv \
  --pool results/repurposing/docking_pool_p05.csv

# 4. SI：仅 III 期+上市子集
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel phase_ge3 --export-p05-pool --skip-tanimoto
```

---

## 7. 标题方向

*Clinical drug repurposing for gout-related URAT1 and NLRP3 targets: NLRP3 machine-learning prescreening, dual-target docking, and molecular dynamics of benchmark inhibitors*

---

## 8. 明确不写的

- 8973 上 NLRP3 ML / Pareto  
- OAT 迁移、Teacher 8973×三态  
- 「发现双靶抑制剂」  
- 全库 MD  
