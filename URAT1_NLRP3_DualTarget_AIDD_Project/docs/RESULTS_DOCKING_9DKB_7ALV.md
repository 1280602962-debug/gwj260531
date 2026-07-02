# 双靶对接结果：9DKB + 7ALV（2026-07 定稿）

> **⚠️ 引擎迁移**：下文数值来自 **Glide XP 开发跑**（历史）。投稿前须用 **AutoDock Vina 1.2.5** 全量重对接并更新本文件（见 [`OPEN_SOURCE_DOCKING.md`](OPEN_SOURCE_DOCKING.md)）。Pareto **逻辑与 benchmark 行为**（colchicine 失败、lesinurad 非前沿）在池内百分位框架下应保持稳定。

> **结构选择**：URAT1 **9DKB**（inward-open cryo-EM）+ NLRP3 **7ALV**（NACHT + MCC950-class analog 共晶模板）。  
> **不再要求 8ETR** 作为主对接结构；8ETR 可作为未来 SI。

---

## 1. 数据覆盖

| 阶段 | n | 说明 |
|------|---|------|
| 临床库 NLRP3 ML 打分 | 8319 | Phase 0–1 |
| P(active) ≥ 0.5 对接池 | 1588 | 进入 Vina 批量对接 |
| 9DKB 有分 | 1455 | Glide 开发跑；Vina 待重跑 |
| 7ALV 有分 | 1517 | Glide 开发跑；Vina 待重跑 |
| **双靶均有分（合并）** | **1451** | Pareto 输入 |
| **Pareto 前沿** | **6** | 短名单 |

---

## 2. Pareto 短名单（6）

| 化合物 | S_U (9DKB %) | S_N (max ML, 7ALV %) | max_phase |
|--------|--------------|----------------------|-----------|
| SLV-334 | 99.9 | 92.1 | 2 |
| LANPROSTON | 99.9 | 96.8 | 2 |
| LASALOCID | 99.7 | 98.3 | 2 |
| EPIGALOCATECHIN GALLATE | 99.2 | 99.7 | 3 |
| FOSIGOTIFATOR | 98.7 | 99.8 | 2 |
| FOSRAVUCONAZOLE | 96.9 | 99.9 | 2 |

**解读**：六者均在双百分位上极高，但多为 **Phase 1–2 早期候选**；仅 EGCG 为 Phase 3。主文应强调这是 **计算重定位假说清单**，需文献与实验复核，不宜直接宣称“发现双靶新药”。

---

## 3. 与 Benchmark / 对照药比较

### 3.1 在 P≥0.5 池内、且完成双靶对接者

| 药物 | 角色 | S_U | S_N | Pareto | 结论 |
|------|------|-----|-----|--------|------|
| **lesinurad** | URAT1 阳性 | 91.6 | 95.0 | 否 | URAT1 对接百分位高；NLRP3 轴被 ML 拉高（非 URAT1 药） |
| **verinurad** | URAT1 阳性 | 77.7 | 97.9 | 否 | 同上；URAT1 对接仍较好 |
| **colchicine** | NLRP3 间接 | 30.7 | 50.1 | 否 | 在 P05 池内（benchmark 强制纳入）；ML 高分为已知局限；7ALV 对接一般 → **不应作 NLRP3 直接 hit** |

### 3.2 不在 P≥0.5 池内（被 NLRP3 ML 预筛排除）

| 药物 | P(active) | 说明 |
|------|-----------|------|
| benzbromarone | ~0 | 强 URAT1 药，**不经过 NLRP3 漏斗** → 用 **8973 回顾 + 单独 9DKB pose** 验证 |
| dotinurad | ~0 | 同上；8973 上对接百分位 ~89%，ML 失败 |
| allopurinol / febuxostat | ~0 | 预期阴性 |
| MCC950 / GDC-2394 | 不在临床 manifest | **MD 用单独 redock**，不进 1588 池 |

### 3.3 ML vs 7ALV 对接（对接池内）

- Spearman ρ( P(active), 7ALV XP ) ≈ **−0.04**（p≈0.17）  
- **结论**：NLRP3 纵轴宜采用 **max(ML 百分位, 7ALV 对接百分位)**（`--sn-mode both`），与“ML 筛库 + 结构验证”叙事一致。

---

## 4. 与 8973 URAT1 回顾的一致性

| 证据 | 结果 | 含义 |
|------|------|------|
| 8973 AUC | 0.705 | 对接对 URAT1 活性有中等区分 |
| EF@5% | 4.23 | 前 5% 富集活性 |
| dotinurad 对接 vs ML | ~89% vs ~5% | **URAT1 主证据必须是对接，不是 ML** |
| lesinurad 在 P05 池 S_U | 91.6% | 与“弱 URAT1 抑制剂”一致但对接仍可识别 |

双轨设计 **自洽**：8973 证明 URAT1 对接管线；1588 池 Pareto 在 NLRP3 ML 预筛后叠加 URAT1 结构约束。

---

## 5. 能否进入 MD？——可以，推荐 2+2

### URAT1 @ 9DKB（主文 Fig 5）

| 化合物 | 理由 |
|--------|------|
| **benzbromarone** | 文献强 URAT1；8973 回顾优秀；不在 P05 池 → **单独取 pose** |
| **dotinurad** | 日本上市 SURI；8973 对接高、ML 低 → 支撑“对接主导”论点 |

备选 SI：**lesinurad**（池内 S_U 高，但临床 URAT1 活性弱于 benz/dotinurad）。

### NLRP3 @ 7ALV（主文 Fig 6）

| 化合物 | 理由 |
|--------|------|
| **MCC950** | 金标准工具药；7ALV 为 analog 模板 → **单独 redock + Methods 说明** |
| **EPIGALOCATECHIN GALLATE** 或 **FOSIGOTIFATOR** | Pareto 前沿代表；前者 Phase 3 更易写讨论 |

**不建议**对 6 个 Pareto 分子全部 MD；选 1 个重定位 lead + MCC950 即可。

### MD 报告指标（与提纲一致）

- 50–100 ns；RMSD、关键残基距离、MM-GBSA（定性）
- URAT1：Phe 笼、substrate pocket；NLRP3：Walker B / sulfonylurea 口袋

---

## 6. 论文可写结论（Results / Discussion 要点）

1. **漏斗可行**：8319 → 1588 → 1451 双靶对接 → 6 Pareto，计算成本可控。  
2. **URAT1 轴可信**：8973 富集 + P05 池内 lesinurad/verinurad 高 S_U；dotinurad/benz 走 benchmark 轨。  
3. **NLRP3 轴需去混淆**：colchicine 高 ML 低结构一致性；phase 组成偏倚见 Fig 2c。  
4. **重定位输出**：6 个双高百分位分子为 **假设清单**；EGCG 等需结合适应症与毒性文献讨论。  
5. **局限**：7ALV 非 MCC950 共晶；133+71 缺分；Pareto 前沿少 → 主文强调 **方法 + 代表药 MD**，非全库命中声明。

---

## 7. 复现命令

```bash
python3 scripts/normalize_canvas_docking_export.py --input results/repurposing/docking_raw/9dkb_xp_raw.csv --pdb 9DKB --output results/repurposing/docking_raw/urat1_9dkb_p05.csv
python3 scripts/normalize_canvas_docking_export.py --input results/repurposing/docking_raw/7alv_xp_raw.csv --pdb 7ALV --output results/repurposing/docking_raw/nlrp3_7alv_p05.csv
python3 scripts/merge_docking_pareto.py --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv --nlrp3-dock results/repurposing/docking_raw/nlrp3_7alv_p05.csv --pool data/repurposing/screening/docking_pool_p05.csv --sn-mode both
python3 scripts/analyze_pareto_benchmarks.py
python3 scripts/plot_available_figures.py
```

---

*关联：`docs/MANUSCRIPT_OUTLINE_CURRENT.md` · `docs/MANUSCRIPT_DRAFT_CURRENT.md` · `docs/WORKFLOW_CURRENT.md`*
