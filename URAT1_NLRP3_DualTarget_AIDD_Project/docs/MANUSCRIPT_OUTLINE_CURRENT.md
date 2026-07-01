# 论文定稿思路（2026-07 修订）

> **首投**：*Journal of Molecular Modeling*  
> **定位**：痛风 URAT1–NLRP3 双节点下的 **临床药物重定位 + 不对称双证据漏斗**  
> **不是**：双靶新药发现、8973 百万虚筛、Teacher/OAT 创新

---

## 0. 一句话

在全临床药物库上用 **NLRP3 ML 快筛** → **Top 命中双靶对接** → **Pareto 短名单**；用 **8973 仅证明 URAT1 应对接**；用 **代表药 MD** 解释机制。

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
        ▼  【Step 1】NLRP3 assay-conditioned ML（~20 s，已完成）
        │   → nlrp3_ml_scores_clinical_all.csv
        │   → Top 5% / Top150：nlrp3_top_for_dual_docking_clinical_all.csv (n=416)
        ▼
        │  【Step 2】对短名单双靶对接（本地 Maestro）
        │   • URAT1 @ 9DKB Glide XP  → S_U^dock
        │   • NLRP3 @ 7ALV 或 8ETR XP → S_N^dock（交叉验证 ML）
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

**与经典「先 ML 再对接」的关系**：  
这是 **单靶 ML 预筛（NLRP3）→ 双靶结构确认（URAT1+NLRP3 对接）**，符合重定位文献（如 2697 FDA 药全库对接；你先用 ML 缩小到 ~400 再对接，**更省算力、更好写**）。

---

## 3. NLRP3 ML 筛选结果（已跑通）

**命令**：
```bash
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --top-n 150 --top-pct 5 --skip-tanimoto
```

| 指标 | 数值 |
|------|------|
| 全库打分 | **8319** |
| P(active)≥0.5 | **1588** |
| 对接短名单 | **416**（Top 5%） |
| 耗时 | ~20 s |

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

或在 416 短名单中 **过滤 max_phase≥3** 再对接（主文更稳）。

---

## 4. 文章结构（五节 Results）

### R1. 数据不对称与方法设计
- ChEMBL 0 重叠；URAT1 ML 2/4 vs NLRP3 AUROC 0.89  
- 流程图：**NLRP3 ML 筛临床库 → 双靶对接 → Pareto**

### R2. NLRP3 ML 筛选临床药物库
- 8319 全库分布、1588 阳性预测  
- 对照药排名表  
- 导出 416 对接短名单的理由（Top 5%）

### R3. URAT1 对接回顾验证（8973，独立一节）
- A vs D EF@5%≈4.2；四药百分位  
- **说明**：为何 URAT1 在 Step 2 用对接、不用 ML  
- **不参与** NLRP3 筛选叙事

### R4. 短名单双靶对接 + Pareto（主图）
- 416（或 phase≥3 子集）@ 9DKB + 7ALV/8ETR  
- $S_U$–$S_N$ 平面；已知药标注  
- SI：NLRP3 对接 vs ML 在短名单上的 Spearman

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

```
# 1. 已完成：NLRP3 ML 全库
# 输出：results/repurposing/nlrp3_top_for_dual_docking_clinical_all.csv

# 2. 可选：只要 III 期+上市进入对接
python3 -c "
import pandas as pd
df=pd.read_csv('results/repurposing/nlrp3_top_for_dual_docking_clinical_all.csv')
sub=df[df.max_phase>=3].copy()
sub.to_csv('results/repurposing/docking_shortlist_phase_ge3.csv', index=False)
print(len(sub), 'compounds')
"

# 3. 对 docking_shortlist 跑 URAT1 9DKB XP + NLRP3 7ALV/8ETR XP（Maestro）

# 4. 合并对接分 → Pareto 脚本（待建 merge_docking_pareto.py）
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
