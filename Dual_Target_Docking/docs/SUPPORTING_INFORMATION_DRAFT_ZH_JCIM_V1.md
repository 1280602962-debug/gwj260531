# Supporting Information

本文件对应正文 `MANUSCRIPT_JCIM_ZH.md`。Supporting Tables 按 Table S1–S10 排列；Supporting Figures 按正文首次引用顺序编号为 Figure S1–S4。历史坐标匹配 RMSD、五种子长表、固定成员交集、分层 bootstrap 备选、样本量情景模拟、BindingDB/PubChem 原始供给清点、发表年份子集，以及逐配体或探索性切片保留在公开仓库，不在投稿 SI 重复排版。

正文主表不在 SI 重复：Table 1 评价集组成与对接条件；Table 2 方向性主结果；Table 3 dual–neither 对照比较。

---

## Table S1. 计算设置与统计定义

| 项目 | 取值 |
|------|------|
| 生产对接 Python / RDKit / Meeko | 本机冻结环境；RDKit 2026.3.1；Meeko 0.7.1 |
| 零对接复分析 | Python 3.12.13；RDKit 2026.3.5；NumPy 2.5.2；pandas 3.0.5；SciPy 1.18.1；scikit-learn 1.9.0（`requirements-analysis.txt`） |
| AutoDock Vina | 1.2.7；默认 `vina` 打分；n_modes = 9；energy_range = 3 kcal mol⁻¹ |
| GNINA | 1.3.2；CPU `--no_gpu` |
| RTMScore | `rtmscore_model1` |
| 配体准备 | 去盐（最大有机片段）→ AddHs → ETKDGv3（seed 20260727）→ MMFF 最多 200 步 → Meeko PDBQT |
| 面板抽样 / holdout 种子 | 20260729 / 20260731 |
| Vina 生产种子 | 20260727；五种子归档另加 20260811–20260814 |
| exhaustiveness | PIK3CA/mTOR = 16；其余主面板与 holdout = 8 |
| Bootstrap | B = 2000。Table 2 报告锁定的配体层非分层百分位 bootstrap 下两条方向性区间及逐次 `summary_min` 区间。dual–neither 为两类分层百分位区间；对应/非对应口袋差值为配对 bootstrap。类别分层 `summary_min` 敏感性见归档表 `summary_min_stratified_sensitivity_review_v1.csv`，不替换 Table 2 |
| GroupKFold / 逻辑回归 | 折数 = min(5, 骨架数, 两类样本量)；C = 1.0；最多 4000 次迭代；默认不打乱 GroupKFold |
| 对接盒子 | 共晶配体 AABB 外扩 5 Å，边长下限 20 Å |
| 共晶 QC 门槛 | 全部保存姿态中的最低重原子 RMSD < 2.0 Å（证明搜索覆盖，不证明 top-1 排序） |

两端均有有效分数的配体进入方向性 AUROC；n_scored 可低于 n_panel（Table 1）。源：`ENV_PIN.md`；`docking_failure_census_v1.csv`。

---

## Table S2. 受体、对接盒子与统一共晶重对接 QC

八个靶对使用 14 个主受体结构（JAK1 的 6N7A 与 PPARA 的 6LXA 被两对共用）。PIK3CA/mTOR 因 4JT6 在 E = 8 时未过门槛而采用 E = 16。EGFR 3POZ 与 HER2 3RCD 按规范重原子盒重对接。9V8H 为 PPARγ LBD–BRL–PG08-NL 三元复合物，对接保留肽链。Figure S4 与近天然判定对 14 个受体统一使用化学对应 CalcRMS 表（S2b）。AChE 4EY7 与 TYK2 3LXP 实际保存 8 个姿态，其最低 RMSD 不是统一的 best-of-9。早期坐标匈牙利匹配结果仅留在仓库，不再作为第二套正式 RMSD 表。

**S2a. 对接盒子（Å）与分辨率**

| 蛋白 | PDB | 共晶配体 | 分辨率 (Å) | center (x, y, z) | size (x, y, z) |
|------|-----|----------|-----------:|------------------|----------------|
| PIK3CA | 4L23 | X6K | 2.50 | 32.443, 45.431, 42.139 | 20.000, 20.000, 20.000 |
| mTOR | 4JT6 | X6K | 3.60 | 51.949, 0.065, −47.707 | 20.332, 20.000, 20.000 |
| AChE | 4EY7 | E20 | 2.35 | −13.988, −43.906, 27.108 | 23.341, 20.000, 20.355 |
| BChE | 4BDS | THA | 2.10 | 133.076, 116.113, 41.335 | 20.000, 20.000, 20.000 |
| EGFR | 3POZ | 03P | 1.50 | 18.680, 32.127, 11.865 | 22.189, 20.000, 22.836 |
| HER2 | 3RCD | 03P | 3.21 | 12.463, 3.371, 27.619 | 23.222, 23.155, 20.000 |
| F2 | 4UDW | N6L | 1.16 | 129.739, −14.402, 164.946 | 20.600, 20.000, 20.000 |
| F10 | 2JKH | BI7 | 1.25 | −3.643, 10.770, 22.771 | 20.000, 20.000, 21.118 |
| JAK1 | 6N7A | KEV | 1.33 | 5.911, −24.317, 21.495 | 20.000, 20.597, 20.000 |
| TYK2 | 3LXP | IZA | 1.65 | −4.498, 25.578, −31.064 | 20.000, 20.000, 20.000 |
| JAK2 | 8BXH | C87 | 1.30 | −12.589, 14.594, 21.263 | 24.069, 20.000, 20.000 |
| PPARG | 9V8H | BRL | 1.39 | −5.454, −5.373, −23.591 | 20.000, 20.140, 20.000 |
| PPARA | 6LXA | EPA | 1.23 | 11.998, 5.742, −7.443 | 20.000, 20.200, 23.432 |
| PPARD | 5U3Q | 7UJ | 1.50 | 40.599, 0.533, 135.353 | 20.223, 20.000, 21.346 |

**S2b. 14 个主受体的化学对应共晶 RMSD（RDKit CalcRMS）**

RMSD 使用 RDKit 对称感知 CalcRMS，不做叠合。原子对应使用同一套流程：准备配体可映射到晶体时用 Meeko 拓扑或 SDF/CCD 图 CalcRMS；PIK3CA/mTOR 的准备配体不在晶体坐标系中，因此使用图自同构 CalcRMS。`2JKH/BI7` 因 OpenBabel SDF 氮原子价态无效而改用 CCD SMILES。EGFR 3POZ 与 HER2 3RCD 为规范重原子盒重对接。

| 蛋白 | PDB | E | top-1 (Å) | top-3 (Å) | 全部保存姿态最低 (Å) | 覆盖 | top-1 < 2 Å |
|------|-----|--:|----------:|----------:|-----------------:|:--------:|:-----------:|
| EGFR | 3POZ | 8 | 1.019 | 1.019 | 1.019 | 通过 | 是 |
| HER2 | 3RCD | 8 | 1.947 | 1.947 | 1.947 | 通过 | 是 |
| JAK1 | 6N7A | 8 | 0.459 | 0.459 | 0.459 | 通过 | 是 |
| JAK2 | 8BXH | 8 | 10.596 | 0.807 | 0.807 | 通过 | 否 |
| TYK2 | 3LXP | 8 | 0.196 | 0.196 | 0.196 | 通过 | 是 |
| PIK3CA | 4L23 | 16 | 0.624 | 0.624 | 0.624 | 通过 | 是 |
| mTOR | 4JT6 | 16 | 7.118 | 0.445 | 0.445 | 通过 | 否 |
| AChE | 4EY7 | 8 | 0.339 | 0.339 | 0.339 | 通过 | 是 |
| BChE | 4BDS | 8 | 4.794 | 0.386 | 0.386 | 通过 | 否 |
| F2 | 4UDW | 8 | 0.382 | 0.382 | 0.382 | 通过 | 是 |
| F10 | 2JKH | 8 | 0.658 | 0.658 | 0.658 | 通过 | 是 |
| PPARG | 9V8H | 8 | 7.085 | 1.636 | 1.636 | 通过 | 否 |
| PPARA | 6LXA | 8 | 7.857 | 7.848 | 1.401 | 通过 | 否 |
| PPARD | 5U3Q | 8 | 1.510 | 1.510 | 1.510 | 通过 | 是 |

源：`all14_cognate_rmsd_calcrrms_v1.csv`；逐姿态 `cognate_rank_rmsd_reaudit_v1.csv`。14 个受体均通过搜索覆盖。JAK2、PPARG、PPARA、BChE 与 mTOR 未过 2 Å top-1 门槛。EGFR 3POZ 与 HER2 3RCD 均通过 top-1 < 2 Å。

---

## Table S3. 活性标签与聚合敏感性

在冻结 Vina 分数上重标。主分析为 θ = 6.0（Table 2）。本表列出 EGFR/HER2 与 PIK3CA/mTOR 上改变类别组成的全部阈值，以及其他六对的锁定 θ = 6.0 行。这六对在 θ = 6.5 和严格 6.5/5.5 规则下保持相同的主要类别组成；完整阈值 × 靶对展开留在仓库。

| 靶对 | 标签规则 | n (D / A / B) | summary_min | 95% CI |
|------|----------|--------------:|------------:|--------|
| EGFR/HER2 | θ = 5.5 | 69 / 22 / 10 | 0.425 | [0.242, 0.626] |
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 | 0.324 | [0.195, 0.471] |
| EGFR/HER2 | θ = 6.5 | 26 / 29 / 29 | 0.460 | [0.304, 0.609] |
| EGFR/HER2 | 严格 6.5/5.5 | 26 / 17 / 7 | 0.324 | [0.138, 0.525] |
| PIK3CA/mTOR | θ = 5.5 | 33 / 9 / 5 | 0.502 | [0.257, 0.625] |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 | 0.692 | [0.470, 0.813] |
| PIK3CA/mTOR | θ = 6.5 | 17 / 15 / 12 | 0.674 | [0.438, 0.797] |
| PIK3CA/mTOR | 严格 6.5/5.5 | 17 / 7 / 4 | 0.639 | [0.317, 0.792] |
| AChE/BChE | θ = 6.0 | 27 / 26 / 28 | 0.606 | [0.443, 0.735] |
| F2/F10 | θ = 6.0 | 31 / 32 / 32 | 0.345 | [0.211, 0.477] |
| JAK1/TYK2 | θ = 6.0 | 31 / 32 / 32 | 0.365 | [0.231, 0.503] |
| JAK1/JAK2 | θ = 6.0 | 32 / 32 / 32 | 0.588 | [0.444, 0.725] |
| PPARG/PPARA | θ = 6.0 | 32 / 31 / 32 | 0.649 | [0.504, 0.751] |
| PPARA/PPARD | θ = 6.0 | 32 / 32 / 32 | 0.446 | [0.296, 0.584] |

**pChEMBL 最大值对中位数（ChEMBL 37；冻结 Vina 分数；不替代 Table 2）：** Table 2 所用的同一套 ChEMBL 37 记录在 θ = 6.0 下按中位数重新汇总。最大 pChEMBL 与全部已评分配体一致（缺失端 0；不一致 0）。最大值到中位数的类别翻转：EGFR/HER2 6/110（主分析 `summary_min` 0.324）；AChE/BChE 1/96（CHEMBL659；0.606 → 0.629）；PPARA/PPARD 1/110（CHEMBL121；A-only 32→31；dual–A-only 0.646 → 0.636；`summary_min` 仍为 0.446）。PIK3CA/mTOR 与其余四对保持类别组成和 `summary_min` 点估计。源：`eight_pair_dump_gated_v1/max_vs_median_auroc_v1.csv`；`eight_pair_dump_gated_v1/parity_v1.csv`。

---

## Table S4. 固定评分通道对照类别比较及两个最大效应的簇重采样

口袋评分保持不变。dual 配体重采样一次；选择性和 neither 阴性独立重采样。Δ = AUROC(dual vs neither) − AUROC(dual vs 对应选择性类别)。正 Δ 表示 neither 对照看起来更容易。PIK3CA/mTOR 的 neither n = 4 效能不足。采用双口袋平均分的 dual 对全部非 dual 比较仅为混合库描述性参考，随分数表归档，不支撑正文结论。

| 靶对 | 评分通道 | dual vs 选择性 | dual vs neither | Δ | 95% CI | neither 效能不足 |
|------|----------|---------------:|----------------:|--:|--------|:----------------:|
| EGFR/HER2 | 口袋 A（对 B-only） | 0.324 | 0.786 | 0.462 | [0.262, 0.651] | 否 |
| EGFR/HER2 | 口袋 B（对 A-only） | 0.666 | 0.720 | 0.054 | [−0.157, 0.246] | 否 |
| AChE/BChE | 口袋 A | 0.606 | 0.590 | −0.016 | [−0.154, 0.120] | 否 |
| AChE/BChE | 口袋 B | 0.650 | 0.709 | 0.058 | [−0.085, 0.196] | 否 |
| PIK3CA/mTOR | 口袋 A | 0.692 | 0.472 | −0.220 | [−0.514, 0.019] | 是 |
| PIK3CA/mTOR | 口袋 B | 0.714 | 0.583 | −0.131 | [−0.433, 0.179] | 是 |
| F2/F10 | 口袋 A | 0.345 | 0.508 | 0.163 | [−0.005, 0.339] | 否 |
| F2/F10 | 口袋 B | 0.413 | 0.528 | 0.115 | [−0.018, 0.254] | 否 |
| JAK1/TYK2 | 口袋 A | 0.365 | 0.809 | 0.444 | [0.263, 0.620] | 否 |
| JAK1/TYK2 | 口袋 B | 0.575 | 0.705 | 0.130 | [−0.060, 0.295] | 否 |
| JAK1/JAK2 | 口袋 A | 0.728 | 0.723 | −0.004 | [−0.177, 0.170] | 否 |
| JAK1/JAK2 | 口袋 B | 0.588 | 0.730 | 0.142 | [−0.058, 0.334] | 否 |
| PPARG/PPARA | 口袋 A | 0.706 | 0.759 | 0.053 | [−0.139, 0.228] | 否 |
| PPARG/PPARA | 口袋 B | 0.649 | 0.644 | −0.005 | [−0.208, 0.205] | 否 |
| PPARA/PPARD | 口袋 A | 0.446 | 0.484 | 0.038 | [−0.139, 0.216] | 否 |
| PPARA/PPARD | 口袋 B | 0.646 | 0.665 | 0.019 | [−0.172, 0.198] | 否 |

**两个口袋 A 差值的簇重采样（dual–neither 减 dual–B-only）。** 上表配体层区间不把同一文献或同一 Bemis–Murcko 骨架的配体视为独立。簇区间不替代 Table 2。EGFR/HER2 的簇 bootstrap 未在规范重原子盒上重算；正式配体层差值为 0.462 [0.262, 0.651]。JAK1/TYK2 簇行未改，因为该对未重做。

| 靶对 | 重采样单位 | Δ 点估计 | 95% CI | CI 排除 0 |
|------|-----------------|--------:|--------|:-------------:|
| EGFR/HER2 | 配体层（正式） | 0.462 | [0.262, 0.651] | 是 |
| JAK1/TYK2 | 骨架簇 | 0.444 | [0.234, 0.633] | 是 |
| JAK1/TYK2 | 文献簇 | 0.444 | [−0.034, 0.682] | 否 |

源：`formulation_equal_score_negative_v1.csv`；`equal_score_negative_s34_v1.csv`；`equal_score_cluster_bootstrap_v1.csv`（仅 JAK1/TYK2）。双口袋平均分的 dual–neither 比较见正文 Table 3。

---

## Table S5. 配体化学基线与对接增量信息

ECFP4 与 ECFP4+对接评分 AUROC 为同一骨架分组交叉验证下的折外预测。最后一列为 Table 2 的原始 Vina 排序 AUROC，仅作描述性参照。Δ = (ECFP4+对接) − ECFP4。16 个方向上最大 |Δ| 为 0.023。主结果使用未缩放逻辑回归。

| 靶对 | 方向 | ECFP4 | ECFP4+对接 | Δ | Vina 排序 AUROC（Table 2） |
|------|------|------:|--------------:|--:|--------------------------:|
| EGFR/HER2 | D vs A | 0.745 | 0.751 | +0.006 | 0.666 |
| EGFR/HER2 | D vs B | 0.890 | 0.887 | −0.002 | 0.324 |
| AChE/BChE | D vs A | 0.895 | 0.893 | −0.002 | 0.652 |
| AChE/BChE | D vs B | 0.821 | 0.808 | −0.013 | 0.606 |
| PIK3CA/mTOR | D vs A | 0.762 | 0.742 | −0.020 | 0.714 |
| PIK3CA/mTOR | D vs B | 0.889 | 0.898 | +0.009 | 0.692 |
| F2/F10 | D vs A | 0.937 | 0.929 | −0.007 | 0.413 |
| F2/F10 | D vs B | 0.665 | 0.687 | +0.021 | 0.345 |
| JAK1/TYK2 | D vs A | 0.846 | 0.840 | −0.006 | 0.575 |
| JAK1/TYK2 | D vs B | 0.902 | 0.903 | +0.001 | 0.365 |
| JAK1/JAK2 | D vs A | 0.915 | 0.916 | +0.001 | 0.588 |
| JAK1/JAK2 | D vs B | 0.968 | 0.969 | +0.001 | 0.728 |
| PPARG/PPARA | D vs A | 0.833 | 0.820 | −0.013 | 0.649 |
| PPARG/PPARA | D vs B | 0.668 | 0.676 | +0.008 | 0.706 |
| PPARA/PPARD | D vs A | 0.932 | 0.928 | −0.004 | 0.646 |
| PPARA/PPARD | D vs B | 0.858 | 0.835 | −0.023 | 0.446 |

**Vina `summary_min` 相对最佳单一描述符的配对 Δ。** 四个单描述符完整矩阵留在仓库，正式 SI 只排每对最佳描述符。AChE/BChE 的 TPSA 方向性 AUROC 为 0.733 / 0.801。八对中六对 95% CI 包含 0；F2/F10 与 JAK1/TYK2 不包含 0。Figure S2 绘制 Vina 区间和描述符点估计，不绘制差值区间。

| 靶对 | 最佳描述符 | 描述符 summary_min | Δ | 95% CI | CI 排除 0 |
|------|-----------------|-----------------------:|--:|--------|:-------------:|
| EGFR/HER2 | cLogP | 0.482 | −0.052 | [−0.200, 0.116] | 否 |
| AChE/BChE | TPSA | 0.733 | −0.128 | [−0.304, 0.049] | 否 |
| PIK3CA/mTOR | 重原子数 | 0.463 | 0.229 | [−0.011, 0.435] | 否 |
| F2/F10 | cLogP | 0.515 | −0.170 | [−0.324, −0.012] | 是 |
| JAK1/TYK2 | cLogP | 0.580 | −0.215 | [−0.374, −0.001] | 是 |
| JAK1/JAK2 | 重原子数 | 0.578 | 0.010 | [−0.084, 0.171] | 否 |
| PPARG/PPARA | TPSA | 0.627 | 0.022 | [−0.166, 0.193] | 否 |
| PPARA/PPARD | cLogP | 0.564 | −0.117 | [−0.331, 0.107] | 否 |

源：`pocket_matched_vs_best_descriptor_delta_v1.csv`；`descriptor_paired_delta_s19_v1.csv`；`incremental_information_v1.csv`；`ecfp4_incremental_s20s24_v1.csv`；`ligand_ml_baseline_scaffold_cv_v1.csv`。

**特征缩放敏感性。** 同一 GroupKFold 划分在各训练折上拟合 `StandardScaler`。16 个方向上最大 |Δ| 为 0.008（PIK3CA/mTOR D vs A）。缩放不替代未缩放的 0.023 主结果。源：`ecfp4_docking_scaler_sensitivity_v1.csv`。

---

## Table S6. 对应与非对应口袋及未使用池留出集

Δ = 对应 `summary_min` − 非对应 `summary_min`。正 Δ 表示对应口袋的较弱方向更高。主评价集中仅 AChE/BChE 的 95% CI 排除 0；EGFR/HER2 包含 0，不解释为具有 matched-pocket 优势。七个已评分留出集均包含 0。区间包含 0 并不证明不存在优势。EGFR/HER2 无留出集。`较弱方向切换 = 是` 表示对应与非对应评分的较弱方向不同，因此 Δ`summary_min` 不能同时代表两个方向。

| 靶对 | 集合 | Δ | 95% CI | CI 排除 0 | 较弱方向切换 |
|------|------|--:|--------|:---------:|:-------------------:|
| EGFR/HER2 | 主集 | 0.056 | [−0.044, 0.160] | 否 | 否 |
| AChE/BChE | 主集 | 0.177 | [0.050, 0.297] | 是 | 是 |
| PIK3CA/mTOR | 主集 | 0.090 | [−0.122, 0.263] | 否 | 否 |
| F2/F10 | 主集 | −0.031 | [−0.117, 0.040] | 否 | 否 |
| JAK1/TYK2 | 主集 | −0.065 | [−0.152, 0.038] | 否 | 否 |
| JAK1/JAK2 | 主集 | −0.019 | [−0.097, 0.054] | 否 | 否 |
| PPARG/PPARA | 主集 | 0.030 | [−0.081, 0.163] | 否 | 是 |
| PPARA/PPARD | 主集 | 0.012 | [−0.092, 0.145] | 否 | 是 |
| AChE/BChE | 留出 | 0.025 | [−0.090, 0.108] | 否 | 是 |
| PIK3CA/mTOR | 留出 | −0.023 | [−0.117, 0.079] | 否 | 是 |
| F2/F10 | 留出 | −0.079 | [−0.251, 0.075] | 否 | 是 |
| JAK1/TYK2 | 留出 | 0.025 | [−0.087, 0.133] | 否 | 否 |
| JAK1/JAK2 | 留出 | 0.008 | [−0.073, 0.111] | 否 | 否 |
| PPARG/PPARA | 留出 | 0.006 | [−0.168, 0.183] | 否 | 是 |
| PPARA/PPARD | 留出 | 0.150 | [−0.056, 0.294] | 否 | 否 |

留出配体来自同一套 ChEMBL 37 记录、排除主评价集成员后按冻结配额抽取。JAK1/JAK2 抽得 20 / 20 / 18。留出集为方向性面板（dual / A-only / B-only），不含 neither。该分析是内部成员敏感性，不作为外部验证。

| 靶对 | 主集 summary_min [95% CI] | 留出 n (D / A / B) | 留出 summary_min [95% CI] |
|------|----------------------------:|----------------------:|------------------------------|
| AChE/BChE | 0.606 [0.437, 0.730] | 20 / 20 / 20 | 0.615 [0.407, 0.764] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 20 / 20 / 20 | 0.765 [0.603, 0.891] |
| F2/F10 | 0.345 [0.211, 0.477] | 19 / 20 / 20 | 0.392 [0.214, 0.573] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 20 / 20 / 20 | 0.475 [0.282, 0.660] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 20 / 20 / 18 | 0.619 [0.420, 0.749] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 20 / 19 / 20 | 0.535 [0.350, 0.717] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 20 / 20 / 20 | 0.445 [0.241, 0.559] |

源：`wrong_pocket_paired_delta_bootstrap_v1.csv`（`set=main_panel` / `unused_pool_holdout`）；`wrong_pocket_by_channel_v1.csv`；`pocket_unidirectional_delta_v1.csv`；`holdout_pocket_matched_v1.csv`；`table2_comparable_by_channel_v1.csv`（`holdout_vina_20260727`）。

---

## Table S7. 计算实现敏感性：受体替换、独立 GNINA 与 PPARG 重评分

独立 GNINA 生成新姿态，不是对 Vina 姿态重评分。范围是 EGFR/HER2、PIK3CA/mTOR 与 JAK1/TYK2。独立 GNINA 并未对每个配体都返回双端评分（EGFR/HER2 缺 EH120_109；PIK3CA/mTOR 缺 PM48_19；JAK1/TYK2 缺 1 个 dual 和 3 个 B-only）。因此 EGFR/HER2 的 dual–neither 使用 n_neither = 11，而主要 Vina Table 3 为 12。EGFR/HER2 与 PIK3CA/mTOR 独立 GNINA 的 `summary_min` 行在源文件中 CI 列为空；下表区间标在对应单臂上，不是对两臂逐次取最小值的区间。五种子 Vina 范围、固定成员交集和完整病例种子表留在仓库。这些种子上数值有波动，但主要任务和靶对模式未变。EGFR/HER2 的设定差距在五个 Vina 种子上均为正。

**S7a. 独立 GNINA 姿态生成**

| 靶对 | 引擎 | n_dual / n_A / n_B / n_neither | summary_min | 较弱臂 AUROC [95% CI] | Dual vs neither |
|------|------|------|------------:|---------------------------|----------------:|
| EGFR/HER2 | Vina 主分析 | 28 / 38 / 32 / 12 | 0.324 [0.195, 0.471] | dual–B-only（口袋 A） 0.324 [0.188, 0.464] | 0.759 [0.557, 0.923] |
| EGFR/HER2 | GNINA 独立 | 28 / 38 / 32 / 11 | 0.265 | dual–B-only（口袋 A） 0.265 [0.148, 0.394] | 0.737 [0.536, 0.903] |
| PIK3CA/mTOR | Vina 主分析 | 18 / 14 / 12 / 4 | 0.692 [0.470, 0.813] | dual–B-only（口袋 A） 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] |
| PIK3CA/mTOR | GNINA 独立 | 18 / 13 / 12 / 4 | 0.633 | dual–A-only（口袋 B） 0.633 [0.427, 0.825] | 0.569 [0.222, 0.889] |
| JAK1/TYK2 | Vina 主分析 | 31 / 32 / 32 / 14 | 0.365 [0.231, 0.503] | dual–B-only（口袋 A） 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] |
| JAK1/TYK2 | GNINA 独立 | 30 / 32 / 29 / 14 | 0.317 [0.183, 0.463] | dual–B-only（口袋 A） 0.317 [0.183, 0.463] | 0.705 [0.517, 0.876] |

**S7b. PPARG/PPARA 同姿态重评分（主优势不稳定）**

| 通道 | summary_min [95% CI] | Dual vs neither |
|------|----------------------|----------------:|
| Vina 主分析 | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] |
| RTMScore（全部保存姿态） | 0.369 [0.233, 0.475] | 0.817 |
| GNINA CNN affinity | 0.500 [0.356, 0.623] | 0.884 |
| 未使用池留出集 | 0.535 [0.350, 0.717] | — |

**S7c. PIK3CA/mTOR 晶体替换**

一次只替换一个口袋：另一口袋保持冻结主集分数。仅此身份已核实的靶对接受了预先规定的替代晶体对接。

| 替换 | 被替换口袋 | D vs A | D vs B | summary_min [95% CI] |
|------|------------|-------:|-------:|----------------------|
| 主集 4L23 / 4JT6 | — | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| PIK3CA → 4JPS | A | 0.714 | 0.486 | 0.486 [0.259, 0.692] |
| PIK3CA → 5DXT | A | 0.714 | 0.505 | 0.505 [0.292, 0.696] |
| mTOR → 4JSX | B | 0.639 | 0.692 | 0.639 [0.418, 0.776] |

源：`independent_dock_formulation_v1.csv`；`table2_comparable_by_channel_v1.csv`；`pocket_matched_PM48_alt4JPS_v1.csv`，`..._alt5DXT_v1.csv`，`..._alt4JSX_v1.csv`。刚性 Cα 叠合为探索性分析，已归档于仓库。

---

## Table S8. 独立性过滤后的外部数据准入

对全部八个靶对检索了 BindingDB 与 PubChem。正式表和 Figure 6C,D 是 BindingDB 按 \(\theta=6.0\) 标记的独立性过滤剩余，不是原始供给清点。PubChem 作为额外的成对数据可用性检查，未并入这些计数。外部对接还要求剔除共享文献来源、重复结构以及 ECFP4 Tanimoto ≥ 0.70 的分子，并要求 dual / A-only / B-only 各类 n ≥ 20 且每类至少 3 个来源。开发分子集合包括主评价集、扩大的 PIK3CA/mTOR PM110 面板和内部留出集。没有靶对满足独立外部评价准入标准，因此未进行外部对接。BindingDB/PubChem 原始供给清点和已建面板上的发表年份子集留在仓库；二者均不作为外部验证。

| 靶对 | 过滤后 dual / A-only / B-only | n_sources (D / A / B) | 门槛 |
|------|------------------------------:|-------------------:|------|
| EGFR/HER2 | 180 / 10 / 20 | 16 / 5 / 4 | 未过（A-only n = 10 < 20） |
| AChE/BChE | 4 / 8 / 14 | 2 / 6 / 3 | 未过 |
| PIK3CA/mTOR | 91 / 4 / 1 | 9 / 2 / 1 | 未过 |
| F2/F10 | 46 / 15 / 16 | 4 / 1 / 3 | 未过（A/B n = 15/16；A-only 仅 1 个来源） |
| JAK1/TYK2 | 323 / 7 / 103 | 20 / 3 / 6 | 未过（A-only n = 7） |
| JAK1/JAK2 | 928 / 40 / 14 | 28 / 7 / 4 | 未过（B-only n = 14） |
| PPARG/PPARA | 0 / 0 / 1 | 0 / 0 / 1 | 未过 |
| PPARA/PPARD | 0 / 1 / 0 | 0 / 1 / 0 | 未过 |

源：`external_slice_summary_v1.csv`。本表是准入筛查，不作为外部验证。

---

## Table S9. 最终靶对准入审计

七对达到 G5 流程兼容性门槛，并因具有可在统一刚性受体 Vina 流程下表示的常规非共价口袋而纳入：PIK3CA/mTOR、AChE/BChE、F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD。本表列出供给受限的 EGFR/HER2 例外，以及未过最终结构或流程兼容性门槛的靶对。这是靶对选择审计，不是对接性能表。源：`pair_eligibility_audit_s14_v1.csv`；`FEASIBLE_PAIR_LADDER_V1.md`；`TIER1_DOCKING_ROSTER_V1.md`；`pair_ligand_identity_qc_v1.csv`。

| 靶对 | 最后到达门槛 | 纳入/排除 | 理由 | 所用证据 |
|------|-------------------|-------------------|--------|---------------|
| EGFR/HER2 | 供给受限例外 | 纳入 | 未满足严格 6.5/5.5 选择性供给标准（最小选择性类别计数为 7），但具有合适的人源全配体结构、可由共晶定义的对接位点，并在 \(\theta=6.0\) 下有足够 dual、A-only 和 B-only 配体做方向性评价。 | `TIER1_DOCKING_ROSTER_V1.md`；Table 1 |
| CTSK/CTSS | G4 配体身份 | 排除 | 两端共晶均为可逆共价半胱氨酸蛋白酶复合物，需超出统一非共价刚性受体 Vina 流程处理。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| CREBBP/BRD4 | G4 配体身份 | 排除 | CREBBP 同时具有 HAT 催化位点和溴结构域，拟对接结构域在统一流程下不能唯一确定。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| F2/PRSS1 | G4 配体身份 | 排除 | 胰蛋白酶（PRSS1）是药理学抗靶，而不是统一流程下的设计双靶伙伴。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| CNR1/CNR2 | G4 配体身份 | 排除 | 膜蛋白 GPCR 对，需超出统一可溶刚性受体 Vina 流程处理构建体与构象状态。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| HCRTR1/HCRTR2 | G4 配体身份 | 排除 | 膜蛋白 GPCR 对，需超出统一可溶刚性受体 Vina 流程处理构建体与构象状态。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| OPRM1/OPRD1 | G4 配体身份 | 排除 | 膜蛋白 GPCR 对，需超出统一可溶刚性受体 Vina 流程处理构建体与构象状态。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| OPRD1/OPRK1 | G4 配体身份 | 排除 | 膜蛋白 GPCR 对，需超出统一可溶刚性受体 Vina 流程处理构建体与构象状态。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| S1PR3/S1PR1 | G4 配体身份 | 排除 | 膜蛋白 GPCR 对，需超出统一可溶刚性受体 Vina 流程处理构建体与构象状态。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| SLC6A4/SLC6A3 | G4 配体身份 | 排除 | 膜蛋白 SLC6 转运体对，需超出统一可溶刚性受体 Vina 流程处理。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| SLC6A2/SLC6A4 | G4 配体身份 | 排除 | 膜蛋白 SLC6 转运体对，需超出统一可溶刚性受体 Vina 流程处理。 | `TIER1_DOCKING_ROSTER_V1.md`；`FEASIBLE_PAIR_LADDER_V1.md` |
| OPRM1/OPRK1 | G3 人源全配体供给 | 排除 | 药物类小分子过滤后最小严格选择性类别计数由 56 降至 46，因此未过 G4 配体身份门槛。 | `pair_ligand_identity_qc_v1.csv`；`FEASIBLE_PAIR_LADDER_V1.md` |
| JAK3/TYK2 | G3 人源全配体供给 | 排除 | 药物类小分子过滤后最小严格选择性类别计数由 51 降至 48，因此未过 G4 配体身份门槛。 | `pair_ligand_identity_qc_v1.csv`；`FEASIBLE_PAIR_LADDER_V1.md` |

---

## Table S10. 八对候选排序操作点

八个主评价集均按双口袋平均 Vina 评分 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\) 排序，并列时按配体 ID 升序。排序读数为完整四状态面板（含 neither）的 Top 10%，\(k=\lceil 0.10\,n\rceil\)。用固定筛选比例而不是固定条数，才能比较不同规模的靶对。Top 10% dual 比例为 \(\mathrm{dual}/k\)。相对全面板 dual 基线，\(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\)；小于 1 表示 Top 10% 中 dual 比例低于全面板。AND 过滤作用于 Dual+A-only+B-only、排除 neither，并保留 \(S_{\mathrm{worst}}\geq\) dual 中位 \(S_{\mathrm{worst}}\) 的配体。这些行是描述性操作点，不是对接质量的八对排行。JAK1/TYK2 Top 10% 见 Figure 2D，AND 过滤见 Figure S1。源：`eight_pair_ranking_operating_point_v1.csv`；`review_scored_membership_v1.csv`。

| 靶对 | n_ranked (D / A / B / N) | k | Top 10% D / A / B / N | dual / k | n_dual / n | EF_dual,10% | AND 输入 → 通过 (D / A / B) | AND dual precision |
|------|-------------------------:|--:|----------------------:|---------:|-----------:|------------:|------------------------------:|-------------------:|
| EGFR/HER2 | 110 (28 / 38 / 32 / 12) | 11 | 1 / 5 / 5 / 0 | 0.091 | 0.255 | 0.357 | 98 → 47 (14 / 9 / 24) | 0.298 |
| JAK1/JAK2 | 110 (32 / 32 / 32 / 14) | 11 | 6 / 5 / 0 / 0 | 0.545 | 0.291 | 1.875 | 96 → 35 (16 / 13 / 6) | 0.457 |
| JAK1/TYK2 | 109 (31 / 32 / 32 / 14) | 11 | 1 / 3 / 7 / 0 | 0.091 | 0.284 | 0.320 | 95 → 50 (16 / 12 / 22) | 0.320 |
| PIK3CA/mTOR | 48 (18 / 14 / 12 / 4) | 5 | 4 / 1 / 0 / 0 | 0.800 | 0.375 | 2.133 | 44 → 17 (9 / 4 / 4) | 0.529 |
| AChE/BChE | 96 (27 / 26 / 28 / 15) | 10 | 5 / 3 / 1 / 1 | 0.500 | 0.281 | 1.778 | 81 → 32 (14 / 7 / 11) | 0.438 |
| F2/F10 | 107 (31 / 32 / 32 / 12) | 11 | 4 / 1 / 6 / 0 | 0.364 | 0.290 | 1.255 | 95 → 59 (16 / 20 / 23) | 0.271 |
| PPARG/PPARA | 109 (32 / 31 / 32 / 14) | 11 | 7 / 3 / 0 / 1 | 0.636 | 0.294 | 2.168 | 95 → 31 (16 / 9 / 6) | 0.516 |
| PPARA/PPARD | 110 (32 / 32 / 32 / 14) | 11 | 5 / 2 / 3 / 1 | 0.455 | 0.291 | 1.562 | 96 → 58 (16 / 22 / 20) | 0.276 |

---

## 与上一版 S1–S14 的对应

上一版正式 SI 的 S2b、S6b、S9c–S9e、S10 分层 bootstrap 与样本量情景、S11a、S12 表格以及 Figure S5 均改为仓库归档。簇重采样的两个最大固定通道差值并入 Table S4；留出集并入 Table S6；受体替换、独立 GNINA 与 PPARG 重评分并入 Table S7；独立性过滤剩余为 Table S8；靶对审计压缩为 Table S9；八对候选排序操作点为 Table S10。

---

## Supporting Figures

补充图按正文首次引用顺序编号为 Figure S1–S4。原 holdout、BindingDB 矩阵、簇重采样图以及可检测效应情景模拟保留在仓库 `figures/jcim_article/`，不在投稿 SI 排版。

### Figure S1. JAK1/TYK2 的 AND 型双口袋过滤

![Figure S1](../figures/jcim_article/FigS1_posthoc_diagnostics.png)

**Figure S1.** JAK1/TYK2 操作点：通过 dual 中位最差靶点评分阈值的化合物类别组成，分母为 Dual+A-only+B-only（n = 95），排除 neither。对应 Top 10% 排序见 Figure 2D；八对同一规则见 Table S10。属探索性分析；数字已在正文 Results 3.2 报告。

### Figure S2. 口袋匹配 summary_min 森林图

![Figure S2](../figures/jcim_article/FigS2_pocket_matched_forest.png)

**Figure S2.** 八个主评价靶对的 Vina \(\mathrm{summary}_{\min}\) 配体水平 bootstrap 95% 置信区间，以及最佳单一描述符的点估计。该图不展示 Vina−描述符差值的置信区间；差值及区间见 Table S5。

### Figure S3. PIK3CA/mTOR 的协议敏感性

![Figure S3](../figures/jcim_article/FigS3_protocol_sensitivity.png)

**Figure S3.** (A) PM48 与 PM110 面板的 Vina \(\mathrm{summary}_{\min}\)。PM48 为主评价集（配额 18/14/12/4，n_scored dual/A/B = 18/14/12，exhaustiveness = 16）；PM110 为同一靶对的更大协议敏感性面板。(B) 在 PM48 上 exhaustiveness 16 与 8。均为描述性点估计，不是 Table 2 的区间。

### Figure S4. 14 个主受体的共晶重对接 RMSD

![Figure S4](../figures/jcim_article/FigS4_cognate_rmsd.png)

**Figure S4.** 各主受体共晶配体重对接的重原子 RMSD，为搜索覆盖检查。圆点为 top-1，菱形为全部保存姿态中的最低重原子 RMSD。AChE 4EY7 与 TYK2 3LXP 保存 8 个姿态，其余主受体保存 9 个姿态。虚线为 2 Å。EGFR 3POZ 的 top-1 为 1.019 Å，HER2 3RCD 的 top-1 为 1.947 Å，均低于 2 Å。14 个受体统一使用化学对应 CalcRMS 表（Table S2）。

---

## Data provenance

表中数值来自冻结 CSV；完整清单、脚本与 SHA-256 记录见仓库 `REVISION_CHECKSUM_MANIFEST_v1.csv`。逐配体长表、探索性切片与已退出体系不在本 SI 排版。
