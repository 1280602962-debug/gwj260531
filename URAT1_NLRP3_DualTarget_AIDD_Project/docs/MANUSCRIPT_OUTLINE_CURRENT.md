# 论文定稿思路（2026-07 修订 · 9DKB + 7ALV）

> **首投**：*Journal of Molecular Modeling*  
> **定位**：痛风 URAT1–NLRP3 双节点下的 **临床药物重定位 + 不对称双证据漏斗**  
> **对接结构**：URAT1 **9DKB** + NLRP3 **7ALV**（主文）；8ETR 不作硬性要求  
> **不是**：双靶新药发现、8973 双靶 Pareto、Teacher/OAT 创新

---

## 0. 一句话

在全临床药物库上用 **NLRP3 ML 快筛** → **P(active)≥0.5（n≈1588）双靶对接（9DKB+7ALV）** → **Pareto 短名单（n=6）**；用 **8973 仅证明 URAT1 应对接**；用 **代表药 MD（2+2）** 解释机制。

---

## 1. 三套数据、三个用途（禁止混用）

| 数据集 | n | 用途 | 是否双靶筛选 |
|--------|---|------|-------------|
| **ChEMBL 临床药物 manifest** | 8319 | **主筛选库**：NLRP3 ML → 9DKB+7ALV → Pareto | ✅ |
| **8973 distill** | 8928 docked | **仅 URAT1 回顾**：A vs D 富集、四药百分位 | ❌ |
| **Benchmark + MD** | 4+2 | 方法验证与机制图 | 平行表征 |

---

## 2. 计算流程（已定稿）

```
ChEMBL 临床药物库 (n=8319)
        │
        ▼  NLRP3 ML（~20 s）
        ▼  P(active) ≥ 0.5 → n=1588
        ▼  URAT1 @ 9DKB XP  +  NLRP3 @ 7ALV XP
        ▼  双靶合并 n=1451 → Pareto 前沿 n=6
        ▼  MD：benzbromarone + dotinurad (URAT1)；MCC950 + 1 Pareto lead (NLRP3)
```

**已完成对接结果**：见 [`RESULTS_DOCKING_9DKB_7ALV.md`](RESULTS_DOCKING_9DKB_7ALV.md)

---

## 3. Results 五节（含 R3）

### R1. 数据不对称与方法设计
- ChEMBL 0 重叠；URAT1 ML 2/4 benchmark fail vs NLRP3 AUROC ~0.89  
- 流程图：ML → 9DKB+7ALV → Pareto；8973 独立回顾轨

### R2. NLRP3 ML 预筛（Fig 2）
- 8319 → 1588；对照药：colchicine/verinurad 高 P（局限）；痛风药低 P  
- Phase 组成偏倚（Phase I 偏高）— Discussion 去混淆

### R3. URAT1 8973 回顾（Fig 3）
- AUC 0.705；EF@5% 4.23；dotinurad 对接高 / ML 极低 → **URAT1 用对接**

### R4. 双靶对接 Pareto（Fig 4）— **主贡献**
- **1451** 双靶合并；**6** Pareto 分子（表 1）  
- 轴：$S_U$（9DKB 百分位）vs $S_N=\max(P_{ML}, P_{7ALV})$  
- 标注：lesinurad、verinurad、colchicine  
- **不写** benzbromarone/dotinurad 在 Pareto 内（未进 P≥0.5 池）

### R5. 代表药 MD（Fig 5–6）
- URAT1：benzbromarone、dotinurad @ 9DKB（单独 pose，因未进 ML 池）  
- NLRP3：MCC950 @ 7ALV（analog 模板说明）+ EGCG 或 FOSIGOTIFATOR（Pareto 代表）  
- 50–100 ns；RMSD、关键相互作用、MM-GBSA 定性

---

## 4. 主图（6 张）

| 图 | 内容 | 状态 |
|----|------|------|
| Fig 1 | 双节点 + ML→9DKB/7ALV→Pareto 流程 | 待画 |
| Fig 2 | NLRP3 ML 分布 + 对照 + 漏斗（1451 双靶完成） | ✅ |
| Fig 3 | 8973 URAT1 回顾 | ✅ |
| Fig 4 | **9DKB–7ALV Pareto** | ✅ |
| Fig 5 | URAT1 MD（benz + dotinurad） | 待本地 MD |
| Fig 6 | NLRP3 MD（MCC950 + lead） | 待本地 MD |

---

## 5. Discussion 要点

1. **不对称双证据**合理：NLRP3 ML 缩库 + URAT1 结构过滤。  
2. **7ALV** 为 MCC950-class analog 口袋；与 8ETR 差异可放 SI（若日后补算）。  
3. **colchicine** 与 **phase 偏倚** 为方法局限，非生物学“Phase I 富集”。  
4. **Pareto 六分子** 为计算假设，需临床阶段与安全性文献支撑。  
5. **不声称**发现首个双靶痛风抑制剂。

---

## 6. 本地下一步（MD + 投稿）

```bash
# 已完成
python3 scripts/merge_docking_pareto.py ...  # 9DKB + 7ALV
python3 scripts/analyze_pareto_benchmarks.py
python3 scripts/plot_available_figures.py

# 本地 Maestro：单独对接 benchmark（未在 P05 池）
# benzbromarone, dotinurad @ 9DKB；MCC950 @ 7ALV

# Desmond/GROMACS：2+2 MD → 导出 Fig 5–6
```

---

## 7. 标题

*Clinical drug repurposing for gout-related URAT1 and NLRP3 targets: NLRP3 machine-learning prescreening, dual-target docking at 9DKB and 7ALV, and molecular dynamics of benchmark inhibitors*

**完整英文稿：** [`MANUSCRIPT_DRAFT_CURRENT.md`](MANUSCRIPT_DRAFT_CURRENT.md)

---

## 8. 明确不写

- 8973 上 NLRP3 ML / Pareto  
- OAT 迁移、Teacher 8973×三态  
- 「发现双靶抑制剂」  
- 全库 MD  

---

*数据与数字：[`RESULTS_DOCKING_9DKB_7ALV.md`](RESULTS_DOCKING_9DKB_7ALV.md)*
