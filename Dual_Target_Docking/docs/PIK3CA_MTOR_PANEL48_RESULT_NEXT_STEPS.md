# PIK3CA/mTOR panel48 v0 — 结果解读与下一步

> 数据：`data/pik3ca_mtor_panel48_v0/`（E=16，seed=20260727，RTM best-of-9）  
> 对照：EGFR/HER2 panel40 v0  
> P1 已完成：协议回写 + `warning_flags.csv` + `analysis/decision_ablation_v0/`

---

## 1. 这一枪完成了什么

| 项 | 状态 |
|----|------|
| Cognate QC (PI-103) | **Go @ E=16**（4L23 mode1=0.62；4JT6 mode1≈7.1 但 best9≈0.3–1.4） |
| 全面板 | 96 作业完成（PM48_34@4JT6 仅 8 个有效 mode，已标注） |
| RTM + 消融表 | 已有；**并列报告** `vina_mean` 与 `rtm_min_z` |
| 失败分型 v0 | 已有（T1/T2/T5） |
| 化学型警告层 | `tables/warning_flags.csv`（诊断列；不进分） |
| 冻结决策消融 | **无法同时满足**（见 `decision_ablation_v0`） |
| Bootstrap CI | N=2000；见 `data/cross_pair_bootstrap_v0/BOOTSTRAP_CI_CONCLUSION_V0.md` |

**CI 要点：** 两对靶上 `rtm_min_z−vina_mean` 的 ΔAUROC **均未**在 95% 排除 0（EGFR 接近边界）；仅 EGFR 的 `consensus_rank_mean` ΔAUROC 显著；硬负 Top10 降幅均不显著。

**协议可迁移结论（弱正向）：** 朴素 Vina mean 不够；加 RTM/min_z 后 AUROC 略升。  
**但失败模式与 EGFR 不同，不能宣称“同一尺子已外推成功”（C4 未闭环）。**

---

## 2. 主表数字（Dual vs rest，N=48）

| 臂 | AUROC | Top10 dual | Top10 硬负 |
|----|-------|------------|------------|
| vina_mean | 0.633 | 5 | **5** |
| **rtm_min_z** | **0.685** | 6 | 4 |
| rtm_shortfall (λ=0.5) | 0.687 | 6 | 4 |
| consensus_rank_mean | 0.668 | 4 | **6** |
| consensus_and_top25 | 0.696 | 6 | 4 |

---

## 3. 关键发现（摘要）

- RTM：**B_only 可压**（WYE-132），但 **A_only T2 抬头**（PM48_26/20/21）+ **T5 误伤**（Torin1/Omipalisib）
- clash / shortfall / consensus（冻结阈值）：**不能**同时降低硬负 Top10 并保护 T5
- 主文卖点：诊断尺子 + 跨对靶失败分型 + 化学型警告；不卖“规则已闭环”

---

## 4. 下一步

### 已完成（P1）
1. ✅ 协议定稿回写（E 对靶特异；双报告臂；Limitations）
2. ✅ 化学型警告层
3. ✅ 冻结决策消融 + 诚实结论

### 可选 / P2
- Torin1/Omipalisib 全 9-mode hinge 一页纸（采样 vs RTM 盲区）
- 扩面板 120–200；第三对靶（**分型语言进主文后再开对接**）

### 明确不做
- 全面板重跑 / E=32 化妆；拧 clash 打 PM48_26；乘客主线；宣称 C4 成功外推

---

## 5. 给本地 / 下一任 agent

> P1 已冻结：并列读 `vina_mean`+`rtm_min_z`，用 warning flags 标 T2/T5；决策规则未额外增益。下一步若动手，优先可选 T5 9-mode 深挖或扩面板设计——**不要**重跑 Vina / 升 E / 开第三对对接除非用户明确要求。
