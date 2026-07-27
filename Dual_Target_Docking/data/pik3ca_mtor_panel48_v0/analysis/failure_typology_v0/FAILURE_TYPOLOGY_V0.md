# PIK3CA/mTOR failure typology v0

> 输入：`analysis/rtm_best_pose_export_v1/`（7 配体 × 两端 RTM-best）  
> 对照：EGFR/HER2 `EH40_23` / `EH40_18`  
> 表：`tables/interaction_summary.csv`, `score_asymmetry.csv`, `chem_similarity*.csv`, `scaffold_flags.csv`

---

## 0. 一句话结论

第二对靶上，RTM/min **不是** EGFR 故事的简单复制：

| 类型 | 代表 | 机制 | clash 门控 | 简单 shortfall |
|------|------|------|------------|----------------|
| **T2 化学型同源假双靶** | PM48_21；PM48_26/20 | 弱端也给出干净 ATP/hinge pose，两端 RTM 都高 | 无效（clash=0） | **无效**（Top1 仍是 PM48_26） |
| **T1 单端虚高被救** | PM48_34 WYE-132 | Vina 把 B_only 抬进 Top10；RTM 打穿 PIK3CA 弱端 | 不需要 | min 已够 |
| **T5 重打分误伤真 dual** | PM48_10 Torin1；PM48_02 Omipalisib | Vina 很强，但 PIK3CA 端 RTM-best 离 hinge/共晶位差 | 无效 | 会更伤 |
| **金标准 OK** | PM48_01 PI-103 | 两端 hinge+高占用；4JT6 需 mode3 | — | — |

**主文可卖：** 诊断尺子 + 跨对靶失败分型。  
**不能卖：** “同一 RTM/min 尺子已在第二对闭环。”

---

## 1. 跨对靶对照（诚实版）

| | EGFR/HER2 panel40 | PIK3CA/mTOR panel48 |
|--|-------------------|---------------------|
| Cognate / E | E=8 够 | **需 E=16**（4JT6 mode1 常 ~7 Å，best9 才过） |
| vina_mean AUROC | ~0.55 | ~0.63 |
| 最佳臂 AUROC | rtm_min_z ~0.71 | rtm_min_z ~0.69 |
| RTM 对硬负 | 整体下降（Top10 硬负 6→3） | **B↓ 但 A↑**（Top10 含 4 个 A_only） |
| 典型被救 | EH40_18 类（RTM 压下） | WYE-132（#10→#40） |
| 典型顽固 | EH40_23（anilinoquinazoline） | PM48_26/20/21（ATP 化学型交叉） |
| 典型误伤 | 较少强调 | Torin1 / Omipalisib 从 Vina Top 掉出 |
| clash 门控 | 打不掉 EH40_23 | 打不掉本对全部个案（均为 0） |

---

## 2. 个案摘要

### 2.1 金标准 — PM48_01 PI-103

| 端 | RTM mode | hinge (VAL) | vs X6K occupancy | MCS RMSD | RTM z |
|----|----------|-------------|------------------|----------|-------|
| 4L23 | 1 | yes 2.80 Å | 1.00 | 1.09 Å | +1.02 |
| 4JT6 | **3** | yes 2.80 Å | 1.00 | **0.45 Å** | +1.21 |

- 排名：vina_mean #9 → rtm_min_z **#4**  
- 含义：ensemble/RTM-best 对 cognate 必要且有效；协议保留 best-of-9。

### 2.2 顽固 A_only — PM48_21（最像 EH40_23）

- 标签：pChEMBL PIK3CA 8.70 / mTOR **5.92**（硬负）  
- 两端 **hinge=yes、clash=0**；对 X6K 占用 **1.00 / 0.97**，MCS RMSD **1.58 / 0.92 Å**  
- 含 morpholine；对 PI-103 Tanimoto 0.27、MCS 13 原子 → **ATP 位点化学型同源**  
- rtm_min_z **#5**；几何门控与 shortfall 都压不住  
- **分型：T2 chemotype_homolog**（对标 EH40_23）

### 2.3 顽固 A_only — PM48_26 / PM48_20

- 互相 Tanimoto **0.46**、MCS 18（同系列 amino-triazine 样 PI3K 化学型）  
- 两端 hinge=yes、clash=0；弱端（mTOR）RTM 仍高（z +1.53 / +1.13）  
- rtm_min_z **#1 / #2** — 比多数真 dual 还靠前  
- **分型：T2**（inactive 端假阳性分数，不是 pose 脏）

### 2.4 被救 B_only — PM48_34 WYE-132

- vina_mean **#10** → rtm_min_z **#40**  
- 4L23 RTM z **−1.38**（短板端被打穿）；4JT6 占用尚可  
- **分型：T1 score_artifact_rescued**（对标 EH40_18）  
- 注：4JT6 仅 8 个有效 mode（已记 Limitations）

### 2.5 误伤 dual — PM48_10 Torin1 / PM48_02 Omipalisib

| | Torin1 | Omipalisib |
|--|--------|------------|
| vina_mean | **#1** | **#3** |
| rtm_min_z | **#31** | **#30** |
| 弱端 | 4L23 mode7 | 4L23 mode7 |
| 4L23 hinge | **no** (5.59 Å) | **no** (3.88 Å) |
| 4L23 vs X6K | occ 0.44, cent 4.2 Å | occ 0.64, cent 2.8 Å |
| 4JT6 | hinge yes；RTM 中高 | hinge yes；RTM 很高 |

- Vina mode1 在 4L23 上分数更好，但 RTM 更偏好 mode7，且 mode7 几何偏离 hinge/共晶位  
- **分型：T5 rescoring_injury / pose-family mismatch**  
- 含义：只报 RTM/min 会误伤经典 dual；主文必须同时报 Vina 与 RTM，并讨论 chemotype 覆盖。

---

## 3. 对决策规则的直接含义（消融前预判）

试算：`score = min(zA,zB) − λ|zA−zB|`，λ∈{0,0.25,0.5,1}  

- AUROC 几乎不动（~0.685→0.691）  
- Top10 仍 **4 个 A_only**，**PM48_26 仍是 #1**  

因为顽固例的问题是 **弱端绝对分也高**，不是两端不平衡。

| 候选规则 | 预期 |
|----------|------|
| clash 门控 | 已证伪（本对个案全 0） |
| 简单 shortfall | 已弱证伪（压不住 T2） |
| hinge 门控 | 压不住 T2（两端都有 hinge）；可能误伤已有 hinge 的 dual |
| 化学型警告层（不改分） | **推荐**：T2 诊断旗标，对标 EH40_23 |
| Vina∩RTM 共识 / 双读数 | **推荐** 处理 T5，避免只信 RTM |
| 端特异校准 / 选择性先验 | P1 试验；需冻结阈值防刷榜 |

**不要**为打掉 PM48_26 去拧 clash 阈值。

---

## 4. 与 EGFR 分型对照表

| 分型 | EGFR/HER2 | PIK3CA/mTOR |
|------|-----------|-------------|
| T1 被救单端虚高 | EH40_18 | PM48_34 WYE-132 |
| T2 化学型同源顽固 | EH40_23 | PM48_21；PM48_26/20 |
| T5 重打分误伤 dual | （次要） | Torin1；Omipalisib |
| Cognate OK | TAK-285 | PI-103（4JT6 用 RTM mode3） |

---

## 5. 下一步操作顺序（执行清单）

### 现在就做（本包已完成 / 可审阅）
1. ✅ RTM-best 姿态导出与口袋级分型  
2. ✅ 跨对靶诚实对照表  
3. ✅ shortfall 预实验（阴性结果写入 Methods）

### 下一步（P1，按序）— **已完成 2026-07-27**
1. ✅ **协议定稿回写**：主读数 `rtm_min_z` + **并列报告 `vina_mean`**；写明 T2/T5 边界；E=16 仅本对；PM48_34 8-mode  
2. ✅ **化学型警告层（不改排序分）**：`tables/warning_flags.csv`（panel48 + panel40）  
3. ✅ **有限决策消融**：vina / rtm_min_z / shortfall / consensus → 见 `../decision_ablation_v0/DECISION_ABLATION_V0.md`（**无法同时满足**）  
4. **Torin1/Omipalisib 深挖（可选 1 页）**：4L23 全 9 mode 的 hinge/占用表  

### 明确后做（P2）
- 扩面板 120–200；第三对靶  
- **不要**在 T2 边界未写入主文前开第三对对接  

### 明确不做
- 全面板重跑 / E=32 化妆  
- 乘客/moiety 主线  
- 宣称 C4 已成功外推  
- 拧 clash 打掉 PM48_26；把 warning flags 写进 gated score

---

## 6. 文件索引

| 路径 | 内容 |
|------|------|
| `../rtm_best_pose_export_v1/` | Maestro 可读 pose 包 |
| `tables/interaction_summary.csv` | hinge / clash / vs X6K |
| `tables/score_asymmetry.csv` | 排名与 shortfall |
| `tables/chem_similarity*.csv` | 化学型相似 |
| `CASE_*.md` | 各案一页纸 |
| `CROSS_PAIR_COMPARISON.md` | 跨对靶表 |
| `NEXT_ACTIONS.md` | 给本地/下一任 agent 的动作 |
