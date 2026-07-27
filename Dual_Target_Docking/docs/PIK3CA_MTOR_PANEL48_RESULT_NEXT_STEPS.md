# PIK3CA/mTOR panel48 v0 — 结果解读与下一步

> 数据：`data/pik3ca_mtor_panel48_v0/`（E=16，seed=20260727，RTM best-of-9）  
> 对照：EGFR/HER2 panel40 v0

---

## 1. 这一枪完成了什么

| 项 | 状态 |
|----|------|
| Cognate QC (PI-103) | **Go @ E=16**（4L23 mode1=0.62；4JT6 mode1≈7.1 但 best9≈0.3–1.4） |
| 全面板 | 96 作业完成（PM48_34@4JT6 仅 8 个有效 mode，已标注） |
| RTM + 消融表 | 已有；主臂 `rtm_min_z` |

**协议可迁移结论（弱正向）：** 朴素 Vina mean 不够；加 RTM/min_z 后 AUROC 略升。  
**但失败模式与 EGFR 不同，不能直接宣称“同一尺子已外推成功”。**

---

## 2. 主表数字（Dual vs rest，N=48）

| 臂 | AUROC | Top10 dual | Top10 硬负 |
|----|-------|------------|------------|
| vina_mean | 0.633 | 5 | **5** |
| vina_min | 0.606 | 6 | 4 |
| rtm_mean | 0.669 | 6 | 4 |
| rtm_min | 0.681 | 6 | 4 |
| **rtm_min_z** | **0.685** | 6 | 4 |
| gated_rtm_min | 0.681 | 6 | 4 |

对照 EGFR/HER2（Dual vs A∪B）：vina_mean≈0.55 → rtm_min_z≈0.71；硬负 Top10 从 6→3。

本对：提升更小（0.63→0.69），且 **RTM Top10 里 A_only 反而更多**。

姿态换位：RTM-best ≠ Vina mode1 — 4L23 **32/48**，4JT6 **30/48**（与 EGFR 同量级，支持 ensemble 必要）。

---

## 3. 关键新发现（决定下一步）

### 3.1 Vina 的污染类型

Vina mean Top10：Torin1/Dactolisib/Omipalisib 等 dual 靠前，但混入 **OSI-027、WYE-132（B_only）、neither、A_only**。

### 3.2 RTM 的“双向效应”

| 现象 | 例子 | 含义 |
|------|------|------|
| 压下部分 B_only | WYE-132 #10→#40；若干 B_only 出 Top10 | 短板/重打分对 mTOR 单端虚高有效 |
| **抬高 A_only（新顽固型）** | **PM48_26 #1**、PM48_20 #2、PM48_21 #5、PM48_23 #9（rtm_min_z） | 类似 EH40_23：两端 RTM 都高，min 也压不住 |
| 误伤真 dual | Torin1 #1→#31；Omipalisib #3→#30；Dactolisib #2→#18 | 需查 pose/化学型，不能只报 AUROC |
| 金标准尚可 | PI-103 #9→**#4** | 重打分对 cognate 有利 |
| 命名硬负大多正确靠后 | Alpelisib/Taselisib/AZD-8055/Ku-0063794 均很后 | 标签锚点大体可信 |

**一句话：** 这不是“RTM 全面复制 EGFR 成功”，而是 **B 端假阳性可压、A 端假双靶抬头 + 部分经典 dual 被打落**。

---

## 4. 接下来该干什么（按优先级）

### P0 — 失败分型个案（本周，发文 sharpness）

写 4–6 个一页纸 case（对接 pose + 两端分数 + 是否 hinge 氢键/化学型同源）：

1. **顽固 A_only：** PM48_26、PM48_20、PM48_21（RTM Top）  
2. **被救 B_only：** WYE-132（PM48_34）  
3. **被误伤 dual：** Torin1（PM48_10）、Omipalisib（PM48_02）  
4. **金标准：** PI-103（PM48_01）mode1 vs RTM-best @4JT6  

产出：`analysis/failure_typology_v0/` + 与 EH40_18/23 对照表。

### P0 — 跨对靶对照表（C4 诚实版）

一张表写清：

| | EGFR/HER2 | PIK3CA/mTOR |
|--|-----------|-------------|
| Cognate | E=8 够 | **需 E=16** |
| Vina AUROC | ~0.55 | ~0.63 |
| 最佳 RTM AUROC | ~0.71 | ~0.69 |
| RTM 对硬负 | 整体下降 | **B↓ A↑** |
| 结论 | 重打分必要且偏成功 | 重打分必要但**规则未闭环** |

### P1 — 决策规则补强（C3，专治 A_only Top）

在已有姿态上试（**不要**为打掉假阳性改 clash 阈值乱调）：

- shortfall / 非对称门槛（A 端分数通胀校正）  
- 相互作用指纹门控（hinge 匹配但弱端无关键接触则降权）  
- 对照臂：只 RTM vs RTM+规则（证明增益来自决策不是换打分）

若仍压不住 PM48_26 类 → Methods 诚实写边界（ATP 化学型交叉），主文卖**诊断+分型**。

### P1 — 协议定稿回写

- `exhaustiveness_v0_1=16` 仅本对（已有 freeze 文件）  
- 主读数臂暂定 `rtm_min_z`，但注明 **跨对靶 Top10 组成不同**  
- PM48_34 8-mode 例外写入 Limitations  

### P2 — 扩规模 / 第三对（发文必需，但排在个案之后）

- 两对靶面板扩到 **120–200**  
- 第三对（Mcl-1/Bcl-xL 或 AChE/BChE）  
- **不要**在没做完 A_only 顽固分型前盲目开第三对对接

### 明确先不做

- 为抬 AUROC 全面板重跑 / 乱升 E=32  
- 乘客/moiety 主线  
- 宣称“尺子已在第二对完美外推”

---

## 5. 给本地 / 下一任 agent 的一句话

> 第二对已跑通 E=16+RTM；下一步优先做 **A_only 顽固假阳性 + Torin1/Omipalisib 误伤** 的 pose 级分型，并完成与 EGFR 的对照表；再决定是否加 shortfall/门控。扩面板与第三对排在分型之后。
