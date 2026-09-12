# 主张—证据对应清单（唯一工作清单）

日期：2026-09-12  
本轮基准：PR #32 head `cursor/jcim-submission-pack-0b1a` @ `381bfe1692e0a1b93d8105194bcdbb7401ea75a8`  
修订分支：`cursor/jcim-manuscript-consolidation-0b1a`  
正式文稿：`docs/MANUSCRIPT_JCIM_ZH.md`、`docs/MANUSCRIPT_JCIM_EN.md`（由章节工作稿组装）  
正式 SI：`docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`、`docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`  
投稿包副本：`submission_pack/manuscript/`、`submission_pack/figures/`

定位：基于双端实验活性数据的双靶对接评价设计与判别来源分析。`summary_min` 是两个方向 AUROC 较低值的描述性汇总，不是新的配体评分函数。

历史审计文件、旧 lock 备忘和对话数字不能自动作为当前主分析依据。本页只指向现行正式稿使用的证据路径。

| 核心结论 | 稿件位置 | 表/图 | 机器可读结果 | 生成或锁定脚本 | 数据版本 | 当前证据允许的解释 |
|---|---|---|---|---|---|---|
| 固定同一评分通道后，更换对照类别及对应化合物集合可改变性能判断，且具有靶对和方向依赖性 | 摘要；Results 3.2；Discussion 4.1；Conclusions 1 | Figure 2A；Table S4 | `formulation_equal_score_negative_v1.csv`；`equal_score_negative_s34_v1.csv` | `build` 路径见 `STATISTICAL_LOCK_V1.md` | 主评价 Vina mode-1；θ=6.0 | 描述性设定对比。更换对照同时更换分子。区间含 0 不能写成无差异或等效。 |
| EGFR/HER2 口袋 A：dual–B-only 0.430 vs dual–neither 0.808，Δ 0.378 [0.205, 0.547] | 摘要；Results 3.2 | Figure 2A；Table S4 | `formulation_equal_score_negative_v1.csv` | 同上 | 冻结生产面板 / 20260727 | 该方向的固定通道差值。不是双口袋平均分效应，也不是因果“活性类别本身”的分离。 |
| JAK1/TYK2 口袋 A：Δ 0.444 [0.263, 0.620] | 摘要；Results 3.2 | Figure 2A；Table S4 | `equal_score_negative_s34_v1.csv` | 五对 stack | 五对生产 Vina / θ=6.0 | 第二对固定通道差值。文献簇区间含 0（Table S10）。 |
| Table 2 八对方向性 AUROC 与 `summary_min` | Results 3.2 Table 2 | Figure 2B；Figure S2 | 原三对：`unified_threshold_sensitivity_v2.csv`（`label_rule=theta_6.0`）；五对：`table2_comparable_theta6_v1.csv` | `scripts/primary/bootstrap_primary.py` | θ=6.0；max pChEMBL；生产 Vina rank-1 affinity | 非分层配体 bootstrap，B=2000；每次重算两方向再取 min。区间含 0.5 ≠ 证明等同随机。八对不平均。 |
| Table 3 dual–neither（`vina_mean`） | Results 3.2 Table 3 | Figure 2C | 原三对：`formulation_conventional_vs_directional_v1.csv`；五对：`table2_comparable_theta6_v1.csv` | Table 3 锁定 | 同上 | 同时改变聚合形式和对照。对照组成的单独影响以 Table S4 为准。 |
| 配体化学可不使用受体结构获得判别；加入对接后 16 个方向 |Δ| ≤ 0.023 | 摘要；Results 3.3；Conclusions 2 | Figure 3；Table S5 | `incremental_information_v1.csv`；`ecfp4_incremental_s20s24_v1.csv` | ECFP4 GroupKFold 脚本 | 主评价集；骨架分组折外预测 | 折外预测 AUROC ≠ Table 2 原始分数 AUROC。不能写成对接不含三维信息，也不能写成提高了前瞻命中率。 |
| 对应口袋优势仅主集两对区间排除 0；七个留出集均含 0 | 摘要；Results 3.4；Discussion 4.3 | Figure 4；Tables S6, S7 | `wrong_pocket_paired_delta_bootstrap_v1.csv`；`wrong_pocket_by_channel_v1.csv`；`holdout_pocket_matched_v1.csv` | 口袋对应脚本 | 主评价与未使用池留出 | Δ`summary_min` ≠ 单方向差值。较弱方向切换时不得按原方向解释。区间含 0 ≠ 证明无结构信息。 |
| 独立 GNINA 姿态生成下仍见类似任务差异 | 摘要；Results 3.4；Methods 2.4.5 | Figure 5A；Table S9a | `independent_dock_formulation_v1.csv`；`table2_comparable_by_channel_v1.csv`（JAK1/TYK2） | 独立 GNINA 分析 | EGFR/HER2、PIK3CA/mTOR、JAK1/TYK2；rank-1 `minimizedAffinity` | 另一姿态生成流程，不是固定通道效应的严格重复，也不是独立外部验证。EGFR `summary_min` 点估计 0.220，无 min-of-two CI；[0.109, 0.343] 属于 dual–B-only。 |
| 同姿态 RTM / GNINA CNN 重评分 | Results 3.4；Methods 2.4.5 | Table S9b | 五对 channel 表 | 同姿态重评分 | PPARG/PPARA 等 | 与独立姿态生成不同。CNN affinity 为主要 CNN 读数。不比较引擎总体性能。 |
| 五种子与固定成员 | Results 3.4；Table S9 | Figure 5C；Table S9c/e | `multiseed_auroc_aggregate_v2.csv`；`fiveseed_summary_min_aggregate_v1.csv`；`multiseed_fixed_membership_v1.csv` | `analyze_multiseed_vina_v2.py`；`multiseed_fixed_membership_v1.py` | 20260727 + 20260811–20260814 | 完整病例 ≠ 固定交集。固定成员只覆盖新增五对。PPARG/PPARA 三种子较弱方向切换。脚本排除空字符串，未另检 NaN/inf 字面量。 |
| API max/median（原三对） | Results 3.5；Discussion 4.5；Table S3 | Table S3 | `assay_max_vs_median_agreement_v1.csv`；`assay_max_vs_median_auroc_v1.csv` | `assay_aggregation_max_vs_median_v1.py` | 2026-08-26 API 快照，不是冻结缓存 | EGFR：一致率 93.6%；冻结 0.430 / API-max 0.417 / median 0.424，含 EH120_060 快照差异。AChE：median 0.629 ≠ 主结果 0.606。不能写成“标签和结果均未改变”。 |
| 高置信字段筛查（原三对） | Results 3.5；Table S3 | Table S3 | `high_confidence_summary_v1.csv` | `high_confidence_label_rebuild_v1.py` | 同日 API；自动字段规则 | 253/253 类别一致，方向性结果不变。不是逐篇读原文，也不覆盖八对。 |
| 五对转储 max/median | Results 3.5；Table S3 | Table S3 | `five_pair_dump_gated_v1/max_vs_median_auroc_v1.csv` | 五对 dump-gated 脚本 | ChEMBL 37 dump | 仓库已完成范围。仅 PPARA/PPARD 1 个翻转（CHEMBL121）；五对 `summary_min` 点估计三位小数不变。区间不进入正文。 |
| 2018 年时间切分 | Results 3.6；Table S12 | Table S12 | `time_split_class_counts_v1.csv`；`five_pair_dump_gated_v1/time_split_v1.csv` | 年份切分脚本 | 已构建评价面板 | JAK1/TYK2、JAK1/JAK2 达报告门槛。内部敏感性，不是外部验证。不得沿用“全部时间切分均失败”。 |
| BindingDB / PubChem | Results 3.6；Table S11 | Figure 6C,D | `external_slice_summary_v1.csv`；`crossdb_strict_supply_v1.csv` | 独立性过滤脚本 | 所用过滤规则与准入门槛 | 未形成外部评价集。不能写成客观上不存在外部数据。S11a 供给清点 ≠ S11b 独立性过滤。 |
| 可检测效应模拟 | Discussion 4.5；Figure S5；Table S10 | Figure S5（磁盘文件 `FigS6_detectable_effect`） | `detectable_effect_simulation_v1.csv` | 独立情景模拟脚本 | 原三对类别样本量 | 与 Table 2 非分层 bootstrap 不同。真实 AUROC=0.50 时事件概率 0.109/0.119/0.128，存在覆盖问题。不能证明主结果只是 n 不足，也不能校准 Table 2。 |
| 共晶重对接 QC | Methods 2.4.4；Results 3.4；Table S2 | Figure S4 | `layer3_cognate_rmsd_calcrrms_v1.csv`；原有受体既有 QC 表 | 不重对接 | 保存姿态 | 最低 RMSD < 2 Å 只说明搜索覆盖。8 姿态体系不是 best-of-9。EGFR 3POZ 为重建 QC。S2b 匈牙利匹配为历史对照。 |
| 数据可用性 | 数据与软件可用性 | — | 分数长表、回放表、checksum | `build_checksum_manifest_v1.py` | git 已提交文件 | 共晶姿态在仓库；五对生产姿态不在 git。checksum ≠ 原始文件。未签发 DOI。 |

本轮未开展高可比性确认、新对接、新种子、新评分模型、新外部集或新模拟。
