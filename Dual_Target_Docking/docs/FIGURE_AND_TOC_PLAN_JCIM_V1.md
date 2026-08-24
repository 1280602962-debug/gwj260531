# 主图、对接展示与 TOC graphic（JCIM Article）

> 取代过时的 `FIGURE_PLAN_V1.md`（未含 holdout / 换晶 / S12 / GNINA best9）。  
> 主张天花板仍服从 `CLAIM_CEILING.md`：图不能比正文更满。

**现状（2026-08-24）：** 主文 Fig 1–7、SI Fig S1–S2、TOC graphic 已由 `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`（扩展图：`plot_jcim_si_composites_v1.py`）从冻结 CSV 绘制，输出 `figures/jcim_article/`。数值锁定见 `plotted_values.json`。未做主面板姿态图（git 无主面板 pose）。

JCIM Articles 无硬图数上限。主文 **7 张图 + 1 张 TOC**：Fig 6 是错口袋方向在 holdout 上反转（改写故事）；Fig 7 是 Fig 4 未展开的指纹/全描述符/协变量/匹配子集。协议旋钮与 `as_is` 供给规则留在 SI，避免与 Fig 3、Fig 2 重复。

---

## 1. 哪些结果值得可视化（按是否值得占主文版面）

| 优先级 | 结果 | 做成图？ | 原因 |
|---|---|---|---|
| 必做 | 四类任务 + 口袋匹配定义 | **Fig 1 示意** | 指标是本文贡献的核心；表说不清 |
| 必做 | 49 对硬负供给 | **Fig 2** | “为什么只有 K=4”必须一眼看见；EGFR=7 要标出来 |
| 必做 | Table 2：口袋匹配 summary_min ± CI，并排最强描述符 | **Fig 3 森林图** | 全文主结果图 |
| 必做 | 弱臂不对称 + 对接 vs 描述符 | **Fig 4** | 否则读者只记住 0.692 |
| 必做 | 主面板 vs holdout，以及换晶 | **Fig 5** | 这是 3.9–3.10 真正改写故事的两块；比 E=8 对照重要得多 |
| 必做 | 错口袋：主面板 matched>wrong，holdout 反转，匹配后仍不恢复 | **Fig 6** | 开放结果；比协议旋钮更改故事 |
| 必做 | 指纹 GroupKFold、全部描述符、协变量、匹配子集 | **Fig 7** | Fig 4 只给最强描述符与 TPSA；这张是检验 |
| 可做（主文或 SI） | PI-103 共晶回收 + 一个 T2 硬负两端都“看起来对” | **SI 姿态图** | 对接论文需要一张姿态图，但不要假装做了 PLIF |
| SI | 阈值网格、GNINA mode01/best9、PM110、E8/E16、单靶 enrichment | **Fig S1** | 协议旋钮不改排序 |
| SI | S12 as_is vs equal_only、holdout 抽样偏移 | **Fig S2** | 供给规则 + 抽样位移；错口袋机制已在 Fig 6 |
| 不要画 | 1000-panel、median 全面板、PLIF 热图、Framework Step 1–5、LigPrep 对比 | — | 没做或禁止写入 |
| 不要画 | `pocket_matched_size_strata_v1.csv` | — | 多层 underpowered，容易误读成稳健分层 |
| 不要画 | `asymmetry_pooled_vs_directional_v1.csv` 的 vina_mean | — | 不是 Table 2 的 θ=6.0 口袋匹配指标（EGFR 0.2824 ≠ 0.4297） |
| 不要画 | `pocket_specificity_gap_v1.csv` 单独成图 | — | 与 Fig 6 的 matched−wrong 重复 |
| 不要画成主图 | 全链序列一致性 vs AUROC（n=4 四个点） | 若画只放 SI | 看起来像相关，实际禁止当相关 |

**旧计划里的 Fig 3 baseline-gate 单独成图：建议并进 Fig 3 或 Fig 4，不要占一张主图。** Δ CI 含 0 用误差条颜色（灰=跨 0）就能说清。

---

## 2. 建议的五张主图（画什么、不画什么）

### Fig 1 — 任务与指标（新画，示意）

两块并排，**不要**写成 Framework Step 1–5。

- **左：** 两个口袋卡通 + 四类配体（dual / A_only / B_only / neither）。A_only、B_only 明确标 *experimental hard-negative*。
- **右：** 口袋匹配读出：dual vs A_only 用**口袋 B** 的分数；dual vs B_only 用**口袋 A**；`summary_min` = 两臂较小值。可加一个浅色“pooled mean”虚影，表示弱臂被掩盖。
- 图注关键词：*evaluation protocol* / DualFourClass-Bench；禁止 D-DRAF、novel framework named。

### Fig 2 — 供给（数据已有）

横轴：49 对的 `min_strict_hardneg`（对数或截断在 200）；纵轴排序。

- 三条水平线：≥50 厚面板、≥20 薄面板。
- 高亮三对厚面板；EGFR/HER2 用箭头标 **7**。
- HDAC1/HDAC6 灰色（金属酶，排除）。
- 可选极小插图：同一四对上 ChEMBL vs BindingDB `equal_only` min HN（S12），标题写 *count-level, no docking*。

### Fig 3 — 主森林图（投稿级重画，≥300 dpi）

每对一组：

- Vina 口袋匹配 summary_min + 95% CI（粗）
- 该对最强平凡描述符（细，另一种标记）
- 可选：RTM、GNINA best-of-9 作更小的点，**不要**用“冠军色”

竖线 AUROC=0.5。按 summary_min 排序或固定 Table 2 顺序。标题/轴必须写 *pocket-matched*，不要用旧的 pooled 标注。

读图应得到：三对贴线或低于描述符；PM 点最高但误差条碰到 0.5。

### Fig 4 — 混淆与弱臂（2 或 3 面板）

- **(A)** 分组柱：每对 D/A（口袋 B）与 D/B（口袋 A）。EGFR 的 B 臂是视觉锚。
- **(B)** 对接 summary_min 对最强描述符，或 Δ ± CI（跨 0 用灰）。
- **(C，可选）** AChE：dual vs 硬负的 TPSA 或重原子分布（小提琴/箱线）。一句图注：TPSA AUROC > Vina。

不要把指纹 0.78–0.91 画成“对接失败证明”；若画，标 *chemotype–label association*。

### Fig 5 — 稳健性：配体层 vs 受体层（这张比旧 Fig 5 的 E8/E16 更值）

- **(A)** 三对：主面板 vs holdout 的 summary_min ± CI（EGFR 不在 holdout，不要空出误导）。
- **(B)** PM48：4L23/4JT6（主）→ 换 4JPS / 5DXT / 4JSX。点估计掉到 ~0.5 必须看得见；CI 含 0.5 用灰色。

E8/E16、PM110、单靶 enrichment 改 SI Fig S1：它们不改变故事。

### Fig 6 — 错口袋（主文；取代旧 SI Fig S1）

四面板，**7.00 in**：

- **(A)** 主面板 K=4：pocket-matched vs wrong-pocket。四对都是 matched > wrong。
- **(B)** unused-pool holdout：三对 wrong ≥ matched（EGFR 无 holdout）。
- **(C)** holdout 效价/尺寸匹配后，九个格子仍 wrong ≥ matched。
- **(D)** scoring-free `contact_count`（不是 PLIF）。B 臂高于随机；幅度不能复制 PIK3CA/mTOR 的 Vina 错口袋。

### Fig 7 — 混淆检验（主文；Fig 4 的展开，不是重复）

- **(A)** ECFP4 GroupKFold，标 *chemotype–label association*。
- **(B)** 四个平凡描述符全画，相对 Fig 4B 只画最强描述符。
- **(C)** 协变量 logistic；EGFR 0.5703 ≠ Table 2 的 0.4297。
- **(D)** 效价/尺寸匹配后的弱臂。

---

## 3. 对接部分可以展示什么

本文是**评测**，不是发现结合模式。姿态图最多一张主文或一张 SI，作用是让“硬负也会给出干净姿态”可看见。

**可以展示（本地有姿态，或 git 已有）：**

| 画面 | 材料 | 图注必须说的话 |
|---|---|---|
| PI-103 / X6K 共晶 vs 对接 | 4L23 + 4JT6 cognate QC；failure typology 金标准 | 协议能回收近晶姿态；4JT6 可能 best-of-9 才像晶体 |
| 一个 T2 硬负两端都 hinge 阳性 | **PM48_21**（A_only；弱端 mTOR 占用 0.97、hinge yes） | 几何干净 ≠ 实验 dual；不要写成新结合模式 |
| 4L23 vs 5DXT 口袋 Cα 叠合 | `jcim_structure_robust_v0` 蛋白 PDB；局域 RMSD 0.343 Å | 骨架几乎不动，AUROC 仍崩；**不是**机制已解决 |
| 可选 SI：holdout 一个 dual 与一个 B_only 的接触表面 | holdout `out.pdbqt` 在 git | 对应 contact_count 的 B 臂混淆，粗粒度埋藏，非 PLIF |

**不要展示：**

- 主面板 8 个“漂亮”姿态拼盘（像结构生物学发现）
- PLIF 指纹热图、残基守恒网络（没做）
- LigPrep 与 RDKit 姿态对比
- 换晶后的对接姿态若要上主文，必须标明 *illustrative of receptor dependence*，不要暗示已用 PLIF 解释崩盘
- GNINA 最小化姿态当“更好的真姿态”

软件：PyMOL/ChimeraX 即可。共晶配体一种颜色、对接另一种；氢债券只画 hinge 1–2 根，避免花哨。

---

## 4. Abstract / TOC graphic 怎么做、怎么写

ACS 要求（JCIM 跟同一套 TOC 规范）：

- **尺寸：** 实际印刷尺寸不超过 **3.25 × 1.75 in**（约 8.25 × 4.45 cm），按这个尺寸交稿，不要先画很大再缩小到看不清。
- **格式：** 彩色 TIF 300 dpi，或 EPS（RGB、字体转曲线）。无衬线字体 Helvetica/Arial，**8 pt，不小于 6 pt**。
- **位置：** 稿件最后一页，标题 **For Table of Contents Only**；也可另传 Graphics for manuscript。
- **内容规则：** 原创、未发表；**不要复用正文里的某一张图**；不要写长句；**不要放具体结果数字**（规范原文：*without providing specific results*）。禁止人像、商标、货币。

### 构图（推荐一屏三块，从左到右）

因为横条很扁，只留 **三个视觉词**：四类硬负、分臂打分、对靶且有限。

1. **左（任务）：** 两个简单蛋白表面 + 四种小分子图标（双色 dual，单色 A_only / B_only，灰色 neither）。标签最多 *dual / hard-negatives*。
2. **中（指标）：** 一根箭头表示“看对端口袋”：A_only 去口袋 B 评分。不要写 AUROC 公式长句；最多 `min(arm)` 一个符号。
3. **右（定性结局，不要数字）：** 四个小靶对点：三个靠近一条虚线 *chance*，一个略高于虚线但旁边一个小警告标记（size / crystal）。**禁止写 0.692、0.430。**

配色建议：dual 与硬负强对比；chance 线中性灰；不要用金/银暗示“PM 赢了”。

### 不要用的 TOC

- 森林图缩小版（有具体 AUROC，且复用 Fig 3）
- Step 1→2→3→4→5 流程图
- 两个结合口袋的美图 + “dual-target drug design”
- DualDiff 生成分子（未评测）
- 大标题 *Docking identifies dual-target ligands*

### 稿里怎么描述（图注 + cover letter）

**最后一页图注（短，几乎无句）：**

> **For Table of Contents Only.** DualFourClass-Bench evaluates whether docking ranks dual-target ligands above experimental single-end hard negatives in both pockets.

**Cover letter 一句（可选）：**

> The TOC graphic depicts the four-class hard-negative task and the pocket-matched readout; it does not reproduce a manuscript figure or report numerical AUROCs.

**Abstract 正文仍然 3–4 句、不描述 TOC。** TOC 是独立视觉摘要，不要在 Abstract 里写 “as shown in the TOC graphic”。

### 中文设计说明（给你自己画图用，不要印在 TOC 上）

> 左：两端口袋与四类配体（真 dual 对实验硬负）。中：对端口袋打分，弱臂决定 summary。右：多数靶对贴近随机，一对略高但受分子属性与晶体选择限制。不标注具体 AUROC。

---

## 5. 执行顺序

1. 重画 Fig 3 森林图（脚本已有，改标注 + ≥300 dpi PDF/PNG）。这张决定审稿人 30 秒印象。
2. 新画 Fig 1 示意 + Fig 2 供给（J0 CSV）。
3. Fig 4 弱臂/描述符；Fig 5 holdout+换晶（S8/S9 CSV）。
4. 主文 Fig 6–7（错口袋；指纹/描述符/协变量/匹配子集）。SI 只留协议旋钮（S1）与 equal-relation 供给 + 抽样偏移（S2）。姿态图若本地有 pose，另做 SI，不占 Fig 6。
5. TOC 单独画，**不要从 Fig 1 裁一块交差**（ACS 不鼓励复用正文图）。按 3.25×1.75 in 交 TIF。

核对：脚本对每个 plotted 值回读冻结 CSV，失败即退出。S1 图注不得写“PM 在全网格都是最高点”（θ=5.5 的 PM 是 underpowered 0.5017，低于 AChE 0.6058）。S1B 的 −0.04..+0.08 是 best9 相对 mode01，不是相对 Vina。Fig 7C EGFR 0.5703 是 logistic AUROC，不是 Table 2 的 0.4297。主文与 SI 不得复用同一张图。
