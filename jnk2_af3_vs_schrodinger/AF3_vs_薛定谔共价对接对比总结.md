# JNK2 Phase 0：AF3 共价对接 vs 薛定谔 Covalent Docking 对比总结

**体系**：PDB 8ELC / Cys116 共价位点  
**化合物集**：4 个 active（YL5084、YL2056、JNK-IN-8、56d）+ property-matched decoy（约 141–143）  
**本地原始数据**：`D:\CADD paper exercise\JNK2\`（`af_output/`、`docking/JNK2covalent_docking_1/`）  
**本目录**：汇总表与对比结论（不含完整 AF3/Glide 原始轨迹体积）

---

## 一句话结论

薛定谔 **Covalent Docking + Glide GScore** 能部分区分 active/decoy，但**明显弱于 AF3 + mPAE**；两者排序只有**弱–中等相关**，**不能互换**。Phase 0 gate 应以 **AF3 mPAE** 为主，Glide 作辅助验证；**56d anchor** 是 Glide 假阳性高发区。

---

## 1. AF3 侧（共价约束 cofolding）

| 项目 | 结果 |
|------|------|
| 任务数 | 147（4 active + 143 decoy） |
| 最佳排序指标 | **mPAE**（`protein_ligand_pae_min`） |
| mPAE AUC | **1.00**；Cohen's d ≈ **5.3** |
| decoy 压过最差 active | **0%** |
| EF@1% / EF@5% | **36.8** / **15.8–21** |

Active mPAE 约 **0.80–0.90 Å**；最好 decoy 约 **1.25–1.30 Å**，中间有约 **0.4 Å** 安全间隔。

ipTM / ranking_score AUC 也约 0.98，但分辨率不足：约 **18/143** decoy 的 ipTM ≥ 0.85（典型假阳性如 `DEC_56d_032`：ipTM 0.94 但 mPAE 1.25 Å）。与 COValid（Shamir et al., JACS 2025）结论一致：**mPAE 为主指标**。

---

## 2. 薛定谔侧（Covalent Docking + Glide）

| 项目 | 结果 |
|------|------|
| 对接任务 | 148（job 000000–000147） |
| 成功 pose | **146**；失败 2（`DEC_YL5084_004/006`）；惩罚分无效 pose 1（`DEC_JNKIN8_017`） |
| 分析用 unique | **145**（4 active + 141 decoy） |
| Glide GScore AUC | **0.957** |
| Cohen's d | **0.09**（极弱，分数大量重叠） |
| decoy 压过最差 active | **14/141 ≈ 10%** |
| EF@1% / EF@5% | **18.1** / **9.1** |

### Active 的 GScore 全局排名

| 排名 | 化合物 | GScore | AF3 mPAE 排名 |
|------|--------|--------|---------------|
| #1 | JNK-IN-8 | −10.07 | #1 |
| #4 | YL2056 | −8.64 | #2（并列） |
| #11 | YL5084 | −7.55 | #2（并列） |
| #18 | 56d | −7.22 | #2（并列） |

Global Top 15 中 **decoy 占 13 个**；AF3 则 4 个 active 全部排在全部 decoy 之前。

### 分 anchor（Glide）

| Anchor | active 排名 | AUC | decoy 压过 active |
|--------|------------|-----|-------------------|
| JNK-IN-8 | #1 / 51 | **1.000** | 0 |
| YL5084 | #2 / 42 | 0.976 | 1 |
| **56d** | **#11 / 51** | **0.800** | **10** |

---

## 3. 两种方法排序能否“对应上”？

**有一定对应，但不强，也不能互换。**

### Spearman（144 化合物；方向对齐为越高越好）

| AF3 指标 | 与 Glide GScore 的 ρ | 对应程度 |
|----------|---------------------|----------|
| pTM | **0.56** | 中等（最强） |
| protein chain pTM | 0.54 | 中等 |
| ranking_score | 0.43 | 弱–中等 |
| ipTM | 0.40 | 弱–中等 |
| **mPAE** | **0.27** | **弱** |

要点：

- 和 Glide **最像**的是 **pTM**，不是 mPAE。
- **mPAE 与 Glide 最不像，但对 active/decoy 区分最好** → “像 Glide” ≠ “更适合筛选”。
- Top-15 与 AF3 Top-15 交集约 **5/15**（Jaccard ≈ 0.20）。

### 典型分歧

- **Glide 假阳性**：`DEC_56d_007/008` 等 Glide Top 内、mPAE ~2.3 Å。
- **分歧最大 active**：**56d**（AF3 很好，Glide #18）。
- **一致点**：**JNK-IN-8** 两种方法都最优。

---

## 4. 使用建议

| 用途 | 建议 |
|------|------|
| Phase 0 gate | **以 AF3 mPAE 为主** |
| Glide 角色 | 辅助看 pose / 相互作用；**不作唯一筛选标准** |
| 联合 | Glide Top ∩ mPAE≤阈值，去掉部分 Glide 假阳性 |
| 慎用 | **单独用 Glide 筛 56d 系列** |

---

## 5. 本目录文件

### `af3_analysis/`

| 文件 | 说明 |
|------|------|
| `af3_all_metrics.csv` | 每 job 主指标 |
| `af3_all_metrics_extended.csv` | 扩展置信度指标 |
| `af3_actives_ranking.csv` | 4 个 active 排名 |
| `af3_per_anchor_summary.csv` | 按锚点汇总 |
| `af3_summary.json` | 全局 AUC/EF 摘要 |
| `metric_*.csv` | 指标分离度 / effect size / 相关 |

### `schrodinger_analysis/`

| 文件 | 说明 |
|------|------|
| `glide_best_per_compound.csv` | 每化合物最优 GScore |
| `glide_all_jobs.csv` | 全部 job |
| `glide_af3_comparison.csv` | AF3 与 Glide 并排 |
| `glide_af3_correlation.csv` | Spearman/Kendall 相关表 |
| `glide_decoys_beating_actives.csv` | 14 个 Glide 假阳性 decoy |
| `glide_metric_ranking.csv` | Glide 各打分指标区分力 |

---

## 6. 与 PaperSpine / COValid 的关系

本对比支持立项叙事中的 **Cys116/8ELC 可复现筛选框架**：在 property-matched decoy 上，AF3+mPAE 达到 COValid 风格 gate（AUC=1、EF@1%≫2），而传统共价对接 enrichment 较弱——与 COValid 论文中物理对接 vs AF3 的相对表现一致。
