# Supporting Information（中文工作稿 · JCIM Articles）

> 与 [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md)（2.1–2.13 协议稿）、[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。  
> **原则：** 仅收录仓库中已有实验/冻结产物；数字均可追溯到下列源文件。未做的分析不填空表。

**溯源（机器可读原料）：**

| SI 表 | 源文件 |
|-------|--------|
| Table S1 | `data/jcim_strengthen_t0t1_v0/ENV_PIN.md`；各 panel `protocol.yaml` |
| Table S2 | `data/*/boxes/*.json`（冻结面板所用条目） |
| Table S3 | PM：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT*.md`；AChE/BChE 与 PIK3CB：各 `cognate_qc/`；ranked re-audit：`data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv`；EGFR/HER2：历史 `protocol/protocol.yaml` + `analysis/exhaustiveness_sensitivity_v1/SENSITIVITY_VERDICT.md` |
| Table S4 | `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` |
| Table S5 | `data/jcim_strengthen_t0t1_v0/tables/matched_subset_directional_v1.csv` |
| Table S6 | `data/jcim_bench_v0/tables/pocket_matched_directional_v1.csv`；worst-pocket 亦见 `data/jcim_strengthen_t0t1_v0/tables/aggregation_sensitivity_v1.csv` |
| Table S7 | `data/jcim_bench_v0/analysis/structural_context_v1/full_chain_identity_v1.py`（脚本）与 `full_chain_identity_v1_output.tsv`（输出） |
| Table S8 | `data/jcim_holdout_v0/tables/holdout_panel_*.csv` + `holdout_ligand_scores_v1.csv` + `holdout_pocket_matched_v1.csv`；结论见 `analysis/HOLDOUT_VERDICT.md` |
| Table S9 | `data/jcim_structure_robust_v0/analysis/STRUCTURE_ROBUSTNESS_QC_V1.md` / `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`；`tables/pocket_matched_PM48_alt*_v1.csv`；`tables/pocket_matched_PAB_alt*_v1.csv`；`tables/receptor_realization_two_pair_v1.csv` |
| Table S10 | `data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/POCKET_MECHANISM_VERDICT_V1.md` + `pocket_superposition_v1.py`（脚本，零新对接，仅用已冻结晶体坐标） |
| Table S11 | `data/jcim_holdout_v0/analysis/WRONG_POCKET_MECHANISM_VERDICT_V1.md` + `scripts/wrong_pocket_contact_v1.py`（脚本，零新对接，仅用已冻结姿态坐标） |
| Table S12 | `data/jcim_supply_crossdb_v0/tables/crossdb_strict_supply_v1.csv`；结论见 `analysis/SUPPLY_CROSSDB_VERDICT_V1.md`（BindingDB REST + PubChem PUG REST 计数，零对接） |
| Table S13 | `data/jcim_holdout_v0/tables/holdout_matched_wrong_pocket_summary_v1.csv` + `holdout_vs_main_potency_size_v1.csv`；结论见 `analysis/HOLDOUT_WRONG_POCKET_POTENCY_VERDICT_V1.md` |
| Table S14 | `data/jcim_bench_v0/tables/gnina_mode01_vs_best9_auroc.csv`；结论见 `analysis/GNINA_BEST9_STATUS.md`（worst-pocket 敏感性，零新对接，仅重打分） |
| Table S15 | `data/jcim_bench_v0/tables/gnina_pocket_matched_mode01_vs_best9_k4_v1.csv` + `..._stability_v1.csv`；结论见 `analysis/GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`（真口袋匹配 GNINA，零新对接，仅重打分） |
| Table S16 | `data/jcim_strengthen_t0t1_v0/tables/endpoint_hierarchy_v1.csv` + `frozen_vs_holdout_v1.csv` |
| Table S17 | `data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv` |
| Table S18 | `data/jcim_strengthen_t0t1_v0/tables/pose_fairness_channels_v1.csv` |
| Table S19 | `data/jcim_strengthen_t0t1_v0/tables/pocket_matched_vs_best_descriptor_delta_v1.csv` |
| Table S20 | `data/jcim_strengthen_t0t1_v0/tables/ligand_ml_scaffold_vs_random_v1.csv` |
| Table S21 | `data/jcim_strengthen_t0t1_v0/tables/ranking_top10_vina_mean_exploratory_v1.csv` |
| Table S22 | `data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv` |
| Table S23 | `data/jcim_novelty_v0/tables/chemotype_matched_hardneg_v1.csv` |
| Table S24 | `data/jcim_novelty_v0/tables/incremental_information_v1.csv` |
| Table S25 | `data/jcim_novelty_v0/tables/mixed_library_enrichment_v1.csv` |
| Table S26 | `data/jcim_novelty_v0/tables/aggregation_min_mean_geometric_harmonic_v1.csv` |
| Table S27 | `docking_failure_census_v1.csv` + `docking_failed_ligand_properties_v1.csv` + `docking_failure_rank_extreme_v1.csv` |
| Table S28 | `data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv` |
| Table S29 | `data/jcim_novelty_v0/tables/assay_max_vs_median_agreement_v1.csv` + `assay_max_vs_median_{summary,auroc,flips}_v1.csv` |
| Table S30 | `data/jcim_structure_robust_v0/tables/receptor_realization_two_pair_v1.csv` |
| Table S31 | `data/jcim_novelty_v0/tables/detectable_effect_simulation_v1.csv` |
| Table S32 | `data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv` + `independent_dock_summary_v1.csv` + `independent_dock_enrichment_v1.csv`；结论见 `analysis/INDEPENDENT_DOCK_VERDICT_V1.md` |
| Table S33 | `data/jcim_structure_robust_v0/analysis/plif_v1/plif_residue_shift_top10_v1.csv`；结论见 `PLIF_VERDICT_V1.md`（几何占有率，非 ProLIF 因果） |
| Table S34 | `data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv`（固定口袋分数，仅替换负类） |
| Table S35 | `measurement_frequency_by_class_v1.csv` + `measurement_frequency_max_median_v1.csv`（测量频次诊断） |
| Table S36 | `high_confidence_summary_v1.csv` + `high_confidence_activity_audit_v1.csv`（当前 ChEMBL 高置信标签视图） |
| Table S37 | `complete_case_usable_pchembl_overlap_v1.csv` + `source_document_concentration_v1.csv`（完整病例覆盖与来源集中度） |
| Table S38 | `class_chemistry_summary_v1.csv`（类别化学空间与骨架诊断） |
| Table S39 | `document_blocked_cv_summary_v1.csv` + `document_cluster_bootstrap_v1.csv`（文献阻断 CV） |
| Table S40 | `document_blocked_cv_methods_v1.csv` + `document_blocked_cv_folds_v1.csv`（同一折上的 ECFP4/物化/docking） |
| Table S41 | `time_split_class_counts_v1.csv` + `time_split_auroc_v1.csv`（冻结文献年份分割） |
| Table S42 | `assay_context_priority_ligands_v1.csv` + `assay_context_audit.csv`（assay-context 机器审计 + 元数据纳入/排除） |
| Table S43 | `bindingdb_independence_summary_v1.csv`（BindingDB REST 结构/文献独立性计数；历史供给核对；未对接，不包装为外部验证） |
| Table S44 | `theta6_pair_census_v1.csv`（θ = 6.0 四状态标签普查；不是对接扩面） |
| Table S45 | `property_caliper_match_v1.csv`（物化 caliper 1:1 匹配） |
| Table S46 | `and_filter_operating_point_v1.csv`（AND 式双口袋工作点） |
| Table S47 | `ligand_only_fullmap_auroc_v1.csv`（四对完整 ChEMBL 图的配体层 AUROC；不是对接） |
| Table S48 | `external_candidate_flow.csv`（BindingDB 202608 原生切片分层计数；未对接） |
| Table S49 | `external_slice_summary_v1.csv`（原生切片主门槛摘要；packaged=0） |
| Table S50 | `mcl1_bclxl_panel_freeze_v1.csv` + `mcl1_bclxl_chembl_panel96_v1.csv`（MCL1/Bcl-xL 面板；LC6 gate 失败后仅作 applicability stress-test） |
| Table S51 | `mcl1_bclxl_receptor_freeze_v1.csv` + `data/mcl1_bclxl_panel_v0/tables/cognate_qc_lc6_v1.csv`（3WIY/3WIZ；LC6 gate FAIL） |
| Table S52 | `benchmark_literature_comparator_v1.csv`（文献对照；不是 bake-off） |
| Table S53 | `data/mcl1_bclxl_panel_v0/tables/formulation_auroc_MBX_v1.csv`（MCL1 Vina AUROC；不进 Table 2） |
| Table S54 | `data/jcim_multiseed_v0/tables/multiseed_auroc_aggregate_v2.csv` + `multiseed_auroc_by_seed_v2.csv`（五种子 Vina；Dual-versus-neither = `AUC(vina_mean)`，与 Table 3 对齐） |
| Figure S7 | `figures/jcim_article/FigS7_posthoc_diagnostics.png` |
| Figure S8 | `figures/jcim_article/FigS_bindingdb_native_slice_v1.png` |
| Master index | `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` |
| Figure S4 | `figures/jcim_article/FigS_pocket_matched_forest.png`（原主文森林图；GNINA 为 Vina 姿态的 CNN 再打分） |
| Figure S5 | `figures/jcim_article/FigS_unused_pool_holdout.png` |
| Figure S6 | `figures/jcim_article/FigS_detectable_effect.png` |
| Figure S9 | `figures/jcim_article/FigS9_ligand_controls.png` |
| Supporting Note S1 | `data/pik3ca_mtor_panel48_v0/analysis/failure_typology_v0/CASE_PM48_21_Aonly.md` + `CASE_PM48_10_02_injured_duals.md` |
| ChEMBL 聚合 | `data/jcim_novelty_v0/tables/assay_max_vs_median_agreement_v1.csv`；审计 `analysis/A4_B5_STATISTICAL_AUDIT_V1.md` |

---

## Table S1. 计算参数与软件版本

| 项目 | 取值 |
|------|------|
| Python | 原始对接：3.x（本机 conda）；当前零对接复分析：3.12.13 |
| RDKit | 原始对接/配体准备：2026.3.1；当前复分析：2026.3.5 |
| meeko | 0.7.1 |
| AutoDock Vina | 1.2.7；`scoring_function = vina` |
| GNINA | 1.3.2（CPU，`--no_gpu`） |
| RTMScore 权重 | `rtmscore_model1` |
| Open Babel | Vina PDBQT → SDF（GNINA 前） |
| 分析库 | 当前复分析：NumPy 2.5.2 / SciPy 1.18.1 / scikit-learn 1.9.0 / pandas 3.0.5（见根目录 `requirements-analysis.txt`） |
| Biopython | `PDBParser` / `Superimposer` / `PairwiseAligner`（序列一致性与 Cα 叠合） |
| ChEMBL | Web API；靶对审计锁定 2026-07-23（未记录 release 编号） |
| 配体准备 | 去盐（最大有机片段）→ RDKit AddHs → ETKDGv3（seed 20260727）→ MMFFOptimizeMolecule（maxIters=200）→ meeko 默认 PDBQT |
| 描述符 | RDKit GetNumHeavyAtoms / MolWt / MolLogP / TPSA |
| 逻辑回归 | scikit-learn `LogisticRegression`（C = 1.0，max_iter = 2000）；支架 `GroupKFold`，折数 min(5, n_pos, n_neg, n_groups) |
| 面板抽样种子 | 20260729 |
| Holdout 抽样种子 | 20260731；仅三对（不含 EGFR/HER2）；PM 排除 PM110 超集 |
| Bootstrap | B = 2000；seed 20260729；SHA-256 稳定子种子；配体层 2.5%–97.5% 百分位区间；不做多重比较校正 |
| Vina `n_modes` | 9 |
| Vina `energy_range` | 3 |
| Vina 随机种子 | 20260727 |
| Vina exhaustiveness | PIK3CA/mTOR 主面板/扩面/单靶对照/换晶 = 16；其余靶对与对应 holdout = 8；PM 另报 E = 8 对照 |
| GNINA | `--cnn_scoring rescore --minimize --seed 20260727` |
| 盒子定义 | 共晶配体 AABB + 5 Å；每边下限 20 Å |
| Cognate 通过门槛 | 重原子 RMSD，`best_of_9` &lt; 2.0 Å（同坐标系，不叠合） |
| 受体 PDBQT | PIK3CA/mTOR/EGFR/HER2：含氢蛋白坐标 + `mk_prepare_receptor.py --read_pdb`；AChE/BChE/PIK3CB：沉积 ATOM + `mk_prepare_receptor`（default altloc A） |
| 换晶 | 一次只换一个口袋；未换端保留冻结主面板分数；4JPS/5DXT 用于 PM48 与 PIK3CA/PIK3CB；4JSX 仅 PM48；3T8M 排除 |
| contact_count | mode-1；配体–受体 重原子距离 ≤4.0 Å；非 PLIF |

采集快照见 `ENV_PIN.md`（2026-07-29）。

---

## Table S2. 对接盒子坐标（冻结受体）

坐标单位：Å。中心与边长来自各面板冻结的 `boxes/*.json`。分辨率：4L23/4JT6/3POZ/3RCD 取自冻结 prepared/protein PDB 头；4EY7/4BDS/2WXF 取自 RCSB entry metadata（与所用 PDB ID 对应）。

| 靶标 | PDB | 分辨率 (Å) | 共晶配体（HET） | center_x | center_y | center_z | size_x | size_y | size_z | 源面板 |
|------|-----|----------:|-----------------|---------:|---------:|---------:|-------:|-------:|-------:|--------|
| PIK3CA | 4L23 | 2.50 | X6K (PI-103) | 32.443 | 45.431 | 42.139 | 20.000 | 20.000 | 20.000 | PM48 / PIK3CA–PIK3CB |
| mTOR | 4JT6 | 3.60 | X6K (PI-103) | 51.949 | 0.065 | −47.707 | 20.332 | 20.000 | 20.000 | PM48 |
| AChE | 4EY7 | 2.35 | E20 (donepezil) | −13.988 | −43.906 | 27.108 | 23.341 | 20.000 | 20.355 | AChE–BChE |
| BChE | 4BDS | 2.10 | THA (tacrine) | 133.076 | 116.113 | 41.335 | 20.000 | 20.000 | 20.000 | AChE–BChE |
| PIK3CB | 2WXF | 1.90 | 039 | −5.454 | −0.547 | 22.243 | 20.000 | 20.000 | 20.000 | PIK3CA–PIK3CB |
| EGFR | 3POZ | 1.50 | 03P (TAK-285) | 18.680 | 32.127 | 11.865 | 22.189 | 20.000 | 22.836 | EGFR–HER2 |
| HER2 | 3RCD | 3.21 | 03P (TAK-285) | 12.463 | 3.371 | 27.619 | 23.222 | 23.155 | 20.000 | EGFR–HER2 |

说明：PIK3CA 在 PIK3CA/mTOR 与 PIK3CA/PIK3CB 两套面板中共用同一 4L23 盒子与受体冻结。

---

## Table S3. 共晶配体重对接（八个靶对口袋槽位，七个唯一受体）

**各靶对口袋均做过 cognate redock。** PIK3CA 的 4L23 被两个靶对共用，因此是八个靶对槽位、七个唯一受体。协议门槛为 `best_of_9` 重原子 RMSD &lt; 2.0 Å。主 seed = 20260727（EGFR 历史 as-run 与敏感性表同 seed）。

**RMSD 定义：** 对接坐标系、不做蛋白叠合。PIK3CA/mTOR 与 EGFR/HER2 历史值：meeko `REMARK SMILES IDX` 映射 + 模板自同构上最小 CalcRMS（见各 panel `rmsd_definition.md`）。AChE/BChE 与 PIK3CB 的旧冻结 QC 使用无元素约束的重原子坐标匈牙利匹配；本次 re-audit 不再直接沿用该方法。对含 Meeko 拓扑备注的 4EY7，直接重建 RDKit 分子；对旧式 4BDS/2WXF PDBQT，先用元素约束 + 晶体坐标（最大映射误差 ≤0.001 Å）映射到参考 SDF 拓扑，再用 symmetry-aware `CalcRMS`。三者均与旧 best 值一致。脚本 `cognate_rank_qc_v1.py`。

### S3a. 冻结受体在协议 exhaustiveness 下的结果

| 靶标 | PDB | 配体 | E | top-1 RMSD (Å) | top-3 best (Å) | all deposited best (Å) | best mode | 面板 E |
|------|-----|------|--:|-----------------:|----------------:|-------------------------:|----------:|-------:|
| PIK3CA | 4L23 | X6K | 16 | 0.624 | 0.624 | 0.624 | 1 | 16 |
| mTOR | 4JT6 | X6K | 8 | 7.118 | 5.003* | 5.003 | 2 | — |
| mTOR | 4JT6 | X6K | 16 | **7.118** | **0.445** | **0.445** | 3 | 16 |
| AChE | 4EY7 | E20 | 8 | 0.339 | 0.339 | 0.339 (8 poses) | 1 | 8 |
| BChE | 4BDS | THA | 8 | **4.794** | **0.386** | **0.386** | 3 | 8 |
| PIK3CB | 2WXF | 039 | 8 | 0.405 | 0.405 | 0.405 | 1 | 8 |
| EGFR | 3POZ | 03P | 8 | **9.505** | **6.227** | **0.760** | 7 | 8 |
| HER2 | 3RCD | 03P | 8 | 1.855 | 1.394 | 1.394 | 3 | 8 |

\* E = 8 历史表仅给 mode1、best mode = 2 与 best-all；因此 top-3 等于已知 best-all。

EGFR/HER2 行为 **reconstructed QC**（原始九姿态生产 PDBQT 未找回；Vina seed 20260727、E = 8、9 modes）。稿件数字取 `cognate_rank_rmsd_reaudit_v1.py` 的 RDKit symmetry-aware `CalcRMS`，不是重建脚本中的匈牙利匹配 RMSD。历史 as-run summary（3POZ top-1 9.483 / best 0.955；3RCD top-1 1.941）仅作量级对照，不以推断值填 top-3。

补充（EGFR 面板历史 as-run，与上表敏感性诊断一致量级）：

| 靶标 | PDB | Vina top1 RMSD (Å) | RTM 选中 mode RMSD (Å) |
|------|-----|-------------------:|----------------------:|
| EGFR | 3POZ | 9.514 | 1.015（mode 2） |
| HER2 | 3RCD | 1.869 | 1.974（mode 3） |

### S3b. 判读（与正文 Methods 2.5 对齐）

1. **八个靶对口袋槽位（七个唯一受体）都做过共晶重对接**；另有未入选候选结构的失败记录（见 S3c），不纳入主协议。
2. **在 E = 8、门槛 = best_of_9 &lt; 2 Å 时：** 4L23、4EY7、4BDS、2WXF、3POZ、3RCD 通过；**仅 4JT6 未通过**（5.003 Å）。升至 E = 16 后 4JT6 的 best_of_9 = 0.445 Å，故 PIK3CA/mTOR 全面板采用 E = 16。
3. **不能把“best_of_9 &lt; 2 Å”等同于“Vina mode1 &lt; 2 Å”。** 4JT6 与 4BDS 都是 top-1 失败、top-3 成功；重建的 3POZ QC 在 E = 8 时 top-1 = 9.505 Å、top-3 = 6.227 Å，近晶构象出现在 mode 7（0.760 Å）。因此本 QC 证明搜索覆盖，不证明 Vina 正确排序 pose。
4. EGFR/HER2 面板仍用 E = 8：敏感性显示升 E 不能修复 3POZ 的 mode1 排序失败，且 E = 8 时 best_of_9 已 &lt; 2 Å。重建的 3RCD QC 通过 top-1（1.855 Å）。

### S3c. 受体筛选中试过但未冻结的候选（实验记录，非主结果）

| 拟用端 | PDB | 配体 | 结果摘要 |
|--------|-----|------|----------|
| BChE | 6ZWI | QRH | best_of_9 ≈ 2.3–2.5 Å @ E8/16，未过门槛 |
| BChE | 6QAA / 5DYW | HUN / 5HF | PDBQT 解析失败 |
| PIK3CB | 2Y3A | GD9 | best_of_9 ≈ 3.85 Å @ E8/16，未过 |
| PIK3CB | 4BFR | J82 | Vina PDBQT 解析失败 |

---

## Table S4. 统一 θ = 6.0 主表述的阈值敏感性网格（支持性）

来源：`unified_threshold_sensitivity_v2.csv`。正文 Table 2 已采用本表 θ = 6.0 行作为四对统一主结果（Results 3.2）；本表 θ = 5.5/6.5 与严格 6.5/5.5 行为支持性阈值敏感性分析，用于说明排序不随阈值网格翻转，不是与 Table 2 竞争的第二套主标准。阈值网格图见 Figure S1A。

| 靶对 | 标签规则 | n (D / A / B) | AUROC D vs A | AUROC D vs B | summary_min | 95% CI | underpowered |
|------|----------|--------------:|-------------:|-------------:|------------:|--------|:------------:|
| EGFR/HER2 | θ = 5.5 | 69 / 22 / 10 | 0.773 | 0.425 | 0.425 | [0.242, 0.626] | 否 |
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 | [0.282, 0.578] | 否 |
| EGFR/HER2 | θ = 6.5 | 26 / 29 / 29 | 0.735 | 0.460 | 0.460 | [0.304, 0.609] | 否 |
| EGFR/HER2 | 严格 6.5/5.5 | 26 / 17 / 7 | 0.799 | 0.324 | 0.324 | [0.138, 0.525] | **是** |
| AChE/BChE | θ = 5.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.442, 0.735] | 否 |
| AChE/BChE | θ = 6.0 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.437, 0.730] | 否 |
| AChE/BChE | θ = 6.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.438, 0.742] | 否 |
| AChE/BChE | 严格 6.5/5.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.442, 0.734] | 否 |
| PIK3CA/PIK3CB | θ = 5.5 | 30 / 25 / 28 | 0.729 | 0.522 | 0.522 | [0.363, 0.667] | 否 |
| PIK3CA/PIK3CB | θ = 6.0 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.350, 0.650] | 否 |
| PIK3CA/PIK3CB | θ = 6.5 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.350, 0.646] | 否 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.342, 0.652] | 否 |
| PIK3CA/mTOR | θ = 5.5 | 33 / 9 / 5 | 0.502 | 0.506 | 0.502 | [0.257, 0.625] | **是** |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 | [0.470, 0.813] | 否 |
| PIK3CA/mTOR | θ = 6.5 | 17 / 15 / 12 | 0.710 | 0.674 | 0.674 | [0.438, 0.797] | 否 |
| PIK3CA/mTOR | 严格 6.5/5.5 | 17 / 7 / 4 | 0.639 | 0.669 | 0.639 | [0.317, 0.792] | **是** |

说明：AChE/BChE 与 PIK3CA/PIK3CB 建造时已按严格配额冻结，重标后计数基本不变。EGFR/HER2 与 PIK3CA/mTOR 在严格规则下 B_only 过少，仅作稳健性描述。早期 `threshold_sensitivity_v1.csv`（vina_mean 通道）保留在仓库作内部对照，不进正文主敏感性表。

---

## Table S5. 效价 / 尺寸匹配子集（口袋匹配 Vina）

来源：`matched_subset_directional_v1.csv`。效价匹配：\|ΔpChEMBL\| ≤ 0.5；尺寸匹配：\|Δheavy atoms\| ≤ 2。单对比 AUROC 与配体层 bootstrap 95% CI。

| 靶对 | 子集 | n_dual / n_other | AUROC | 95% CI |
|------|------|-----------------:|------:|--------|
| EGFR/HER2 | potency D vs A | 17 / 17 | 0.747 | [0.561, 0.914] |
| EGFR/HER2 | potency D vs B | 14 / 14 | 0.469 | [0.260, 0.704] |
| EGFR/HER2 | size D vs A | 24 / 24 | 0.547 | [0.396, 0.707] |
| EGFR/HER2 | size D vs B | 22 / 22 | 0.519 | [0.339, 0.692] |
| AChE/BChE | potency D vs A | 20 / 20 | 0.593 | [0.395, 0.765] |
| AChE/BChE | potency D vs B | 24 / 24 | 0.601 | [0.443, 0.754] |
| AChE/BChE | size D vs A | 16 / 16 | 0.484 | [0.285, 0.680] |
| AChE/BChE | size D vs B | 15 / 15 | 0.596 | [0.373, 0.800] |
| PIK3CA/PIK3CB | potency D vs A | 20 / 20 | 0.680 | [0.483, 0.853] |
| PIK3CA/PIK3CB | potency D vs B | 20 / 20 | 0.458 | [0.280, 0.645] |
| PIK3CA/PIK3CB | size D vs A | 20 / 20 | 0.630 | [0.435, 0.805] |
| PIK3CA/PIK3CB | size D vs B | 17 / 17 | 0.450 | [0.246, 0.654] |
| PIK3CA/mTOR | potency D vs A | 13 / 13 | 0.710 | [0.497, 0.900] |
| PIK3CA/mTOR | potency D vs B | 9 / 9 | 0.728 | [0.407, 0.975] |
| PIK3CA/mTOR | size D vs A | 12 / 12 | 0.778 | [0.569, 0.951] |
| PIK3CA/mTOR | size D vs B | 12 / 12 | 0.722 | [0.493, 0.917] |

PIK3CA/mTOR 若干臂 n &lt; 15，区间宽；正文仅作方向是否同向的描述，不以子集点估计作主主张。

---

## Table S6. 分数聚合对照（Vina）

来源：池化 / 口袋匹配 / 错口袋取自 `pocket_matched_directional_v1.csv`；worst-pocket 取自 `aggregation_sensitivity_v1.csv`。同一面板、同一标签下的四种聚合。

| 靶对 | 池化 summary_min | 口袋匹配 summary_min | 错口袋 min | worst-pocket min |
|------|-----------------:|---------------------:|-----------:|-----------------:|
| EGFR/HER2 | 0.311 | 0.430 | 0.260 | 0.271 |
| AChE/BChE | 0.530 | 0.606 | 0.444 | 0.579 |
| PIK3CA/PIK3CB | 0.412 | 0.500 | 0.349 | 0.439 |
| PIK3CA/mTOR | 0.671 | 0.692 | 0.602 | 0.627 |

正文主指标为口袋匹配；池化、错口袋与 worst-pocket 仅作对照。

---

## Table S7. 靶对内全链序列一致性（结构决定因素，探索性）

来源：`data/jcim_bench_v0/analysis/structural_context_v1/full_chain_identity_v1.py`（脚本）与同目录 `full_chain_identity_v1_output.tsv`（输出）。方法：Biopython `PDBParser` 从各冻结受体 `*_protein.pdb` 提取最长蛋白链的一级序列（仅计入标准氨基酸 ATOM 残基），用 `Bio.Align.PairwiseAligner`（BLOSUM62，全局比对，gap open = −11，gap extend = −1）两两比对靶对内两条链。一致性两种归一化：以比对总长为分母，以较短链长度为分母。

| 靶对 | 受体 PDB (A / B) | 链长 (A / B，残基数) | 匹配残基数 | 比对长度 | 一致性（/比对长度，%） | 一致性（/较短链，%） |
|------|------------------|----------------------:|-----------:|---------:|----------------------:|----------------------:|
| PIK3CA/mTOR | 4L23 / 4JT6 | 1014 / 1054 | 213 | 1175 | 18.1 | 21.0 |
| PIK3CA/PIK3CB | 4L23 / 2Y3A | 1014 / 976 | 423 | 1045 | 40.5 | 43.3 |
| AChE/BChE | 4EY7 / 4BDS | 535 / 524 | 278 | 536 | 51.9 | 53.1 |
| EGFR/HER2 | 3POZ / 3RCD | 275 / 274 | 210 | 294 | 71.4 | 76.6 |

说明：这是全链一级序列一致性，不是经结构叠合的口袋残基级 RMSD 或 PLIF 相似度；后者需 TM-align/PyMOL 等经验证的结构叠合工具与手工核对口袋残基对应关系，本轮未做，不虚构该数值。n = 4，仅作描述性对照，不做正式相关性检验。详见 Results 3.6。

---

## Table S8. 未参与建面的 unused-pool holdout 配体与口袋匹配分数

来源：`data/jcim_holdout_v0/`。抽样种子 20260731；标签规则 strict 6.5/5.5；对接协议与主面板相同。完整配体级表见 `holdout_ligand_scores_v1.csv`（含 SMILES、类别、两端 Vina mode-1、physchem）。下表为口袋匹配汇总。

| 靶对 | n (D/A/B) | pocket_matched summary_min | 95% CI | D vs A | D vs B | 主面板 | 最强平凡基线 | Δ(dock−基线) |
|------|-----------|---------------------------:|-------:|-------:|-------:|-------:|-------------:|-------------:|
| PIK3CA/mTOR | 60 (20/20/20) | 0.765 | 0.603–0.891 | 0.860 | 0.765 | 0.692 | heavy 0.555 | +0.210 |
| AChE/BChE | 60 (20/20/20) | 0.618 | 0.422–0.759 | 0.635 | 0.618 | 0.606 | cLogP 0.575 | +0.043 |
| PIK3CA/PIK3CB | 59 (20/19/20) | 0.425 | 0.241–0.618 | 0.766 | 0.425 | 0.500 | heavy 0.691 | −0.266 |

失败：`HOAP_028`（含硼，AutoDock 原子类型 `B` 不支持）两端未得分，已剔除。错口袋对照见 `holdout_pocket_matched_v1.csv` 中 `wrong_pocket_control_vina` 行。

---

## Table S9. Receptor-realization sensitivity: cognate QC and PM48 one-pocket swap

来源：`data/jcim_structure_robust_v0/analysis/STRUCTURE_ROBUSTNESS_QC_V1.md`。协议：Vina E=16（PM48），seed=20260727，best_of_9 &lt; 2 Å。这是 receptor-realization sensitivity，不是稳健性证明。

| 靶标 | PDB | 共晶配体 | mode1 RMSD (Å) | best_of_9 (Å) | 结论 |
|------|-----|----------|---------------:|--------------:|------|
| PIK3CA | 4JPS | 1LT | 0.607 | 0.607 | PASS |
| PIK3CA | 5DXT | 5H5 | 0.624 | 0.624 | PASS |
| mTOR | 4JSX | 17G (Torin2) | 0.515 | 0.515 | PASS |

PASS 受体已写入 `receptors/`。PIK3CA/mTOR PM48 单端替换后的口袋匹配 summary_min：

| 替代结构 | 替换口袋 | 保留口袋 | summary_min [95% CI] | 主面板 | Δ |
|----------|----------|----------|---------------------:|-------:|--:|
| 4JPS | A | 4JT6 | 0.486 [0.259, 0.692] | 0.692 | −0.206 |
| 5DXT | A | 4JT6 | 0.505 [0.292, 0.696] | 0.692 | −0.187 |
| 4JSX | B | 4L23 | 0.639 [0.418, 0.776] | 0.692 | −0.053 |

PIK3CA/PIK3CB 的同一 PIK3CA 晶体替换见 Table S30（方向相反）。详见 `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`。

---

## Table S10. 受体依赖的探索性结构对照：晶体间 Cα 叠合

来源：`data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/POCKET_MECHANISM_VERDICT_V1.md`。方法：Biopython `PDBParser` + `Superimposer`；口袋残基由参考结构自身共晶配体重原子 ≤5 Å 界定，按残基编号+残基名精确匹配（匹配位点零错配）；全域与口袋局域 RMSD 共用同一次刚体拟合。PIK3CA 替代结构 n = 2，mTOR n = 1。

| 参考（主面板） | 替代 | 匹配 Cα 数 | 全域 Cα RMSD (Å) | 口袋残基数 | 口袋局域 Cα RMSD (Å) | 共晶配体质心距离 (Å) |
|------|------|-----------:|------------------:|-----------:|----------------------:|----------------------:|
| 4L23（PIK3CA） | 4JPS | 982 | 1.486 | 20 | 0.867 | 2.566 |
| 4L23（PIK3CA） | 5DXT | 862 | 1.441 | 20 | 0.343 | 2.072 |
| 4JT6（mTOR） | 4JSX | 1054 | 0.454 | 18 | 0.467 | 2.196 |

PIK3CA 口袋残基（4L23，共晶 X6K）：Met772、Trp780、Ile800、Lys802、Leu807、Asp810、Leu814、Tyr836、Cys838、Ile848、Glu849、Val850、Val851、Ser854、Thr856、Gln859、Met922、Phe930、Ile932、Asp933。mTOR 口袋残基（4JT6，共晶 X6K/PI-103）：Ile2163、Pro2169、Leu2185、Lys2187、Glu2190、Leu2192、Asp2195、Tyr2225、Val2227、Ile2237、Gly2238、Trp2239、Val2240、Met2345、Leu2354、Ile2356、Asp2357、Phe2358。5DXT 匹配 862 个 Cα，少于 4JPS 的 982 个，全域 RMSD 不是等覆盖比较。数字与换晶后不对称**方向一致，不作定量因果解释**。详见 Results 3.4。

---

## Table S11. Holdout 错口袋对照的几何对照（scoring-free contact_count）

来源：`data/jcim_holdout_v0/analysis/WRONG_POCKET_MECHANISM_VERDICT_V1.md`；脚本 `data/jcim_holdout_v0/scripts/wrong_pocket_contact_v1.py`。定义：配体重原子中与受体重原子距离 ≤4.0 Å 的原子数（`contact_count`），直接取自已冻结的 mode-1 姿态坐标，不涉及 Vina 能量函数。Vina 错口袋分臂取自 `holdout_pocket_matched_v1.csv`；重原子均值取自 `holdout_ligand_scores_v1.csv`。

| 靶对 | Vina 错口袋 summary_min（D/A，D/B） | contact_count AUROC（口袋 A / 口袋 B） | contact_count 的 min | dual / A_only / B_only 重原子均值 |
|------|------------------------------------:|----------------------------------------:|---------------------:|----------------------------------:|
| AChE/BChE | 0.643（0.643 / 0.653） | 0.581 / 0.706 | 0.581 | 35.1 / 34.0 / 29.5 |
| PIK3CA/mTOR | 0.788（0.788 / 0.858） | 0.552 / 0.698 | 0.552 | 33.5 / 32.3 / 31.0 |
| PIK3CA/PIK3CB | 0.520（0.640 / 0.520） | 0.622 / 0.714 | 0.622 | 34.5 / 31.6 / 28.3 |

B 臂 contact_count 高于随机（0.698–0.714），与 dual 对 B_only 尺寸差更大一致；A 臂接近随机（0.552–0.622），与 dual 对 A_only 尺寸差很小一致。contact_count **不能按幅度复现** Vina 错口袋（尤其 PM：0.788 对 0.552）。详见 Results 3.5。

---

## Table S12. 冻结 K=4 的 BindingDB / PubChem 严格硬负计数核对（零对接）

来源：`data/jcim_supply_crossdb_v0/tables/crossdb_strict_supply_v1.csv`；脚本 `scripts/bindingdb_pubchem_strict_count_v1.py`；结论 `analysis/SUPPLY_CROSSDB_VERDICT_V1.md`。规则与 J0 相同（dual 两端 ≥ 6.5；A_only A ≥ 6.5 且 B ≤ 5.5；B_only 对称）。BindingDB：REST `getLigandsByUniprots`，cutoff = 1 mM，按 monomerid 配对。PubChem：PUG REST `protein/accession/…/concise`，按 CID 配对。`equal_only` 只保留等式（或无修饰）测定，作为与 ChEMBL pChEMBL 的主比较；`as_is` 把 `>`/`<` 的数值当作点估计（敏感性）。**不做**跨库 InChIKey 合并，不重建面板，不对接。

**主比较（ChEMBL pChEMBL vs BindingDB/PubChem `equal_only`）**

| 靶对 | ChEMBL both / dual / A/B / min HN | BindingDB equal_only both / dual / A/B / min HN | PubChem equal_only both / dual / A/B / min HN | ≥50 厚面板门槛是否翻转 |
|------|----------------------------------:|-----------------------------------------------:|---------------------------------------------:|:----------------------:|
| PIK3CA/mTOR | 2713 / 1552 / 80/81 / **80** | 2739 / 1579 / 76/96 / **76** | 2955 / 1602 / 86/93 / **86** | 否（仍过） |
| AChE/BChE | 2537 / 687 / 189/78 / **78** | 2711 / 698 / 181/92 / **92** | 2916 / 742 / 214/97 / **97** | 否（仍过） |
| PIK3CA/PIK3CB | 1990 / 602 / 56/67 / **56** | 2545 / 855 / 58/75 / **58** | 2860 / 908 / 61/74 / **61** | 否（仍过） |
| EGFR/HER2 | 1751 / 951 / 39/7 / **7** | 2269 / 1336 / 34/31 / **31** | 2068 / 1121 / 43/30 / **30** | 否（仍不过 ≥50；升至薄面板 ≥20） |

**敏感性（`as_is`，含 `>` 截尾）**

| 靶对 | BindingDB as_is A/B / min HN | PubChem as_is A/B / min HN |
|------|-----------------------------:|---------------------------:|
| PIK3CA/mTOR | 389/151 / 151 | 405/153 / 153 |
| AChE/BChE | 228/141 / 141 | 275/153 / 153 |
| PIK3CA/PIK3CB | 208/129 / 129 | 212/144 / 144 |
| EGFR/HER2 | 85/92 / **85**（过 ≥50） | 88/92 / **88**（过 ≥50） |

EGFR/HER2 的 as_is 抬升不可直接当成“ChEMBL 漏检”：BindingDB 92 个 as-is B_only 中，**49** 个在 EGFR 端只有 `>` 记录（典型 IC50 > 10 µM 选择性面板），**43** 个至少有一条 EGFR 等式记录。PubChem 与 BindingDB 数量接近，符合沉积重叠，不是两次独立普查。详见 Results 3.1。

---

## Table S13. Holdout 错口袋：效价/尺寸匹配诊断（零新对接）

来源：`data/jcim_holdout_v0/tables/holdout_matched_wrong_pocket_summary_v1.csv`、`holdout_vs_main_potency_size_v1.csv`；脚本 `scripts/wrong_pocket_potency_match_v1.py`。匹配规则与 Table S5 相同（效价：共享活性端 \|ΔpChEMBL\| ≤ 0.5；尺寸：\|Δheavy\| ≤ 2）。口袋匹配：D/A 用 vina_B，D/B 用 vina_A；错口袋为对调。主 holdout 数字仍以 Table S8 未匹配全样本为准。

**相对主面板的均值偏移（holdout − 主面板）**

| 靶对 | dual pA / pB | A_only pA / pB | B_only pA / pB |
|------|-------------:|---------------:|---------------:|
| AChE/BChE | −0.16 / −0.46 | +0.25 / +0.12 | 0.00 / −0.15 |
| PIK3CA/PIK3CB | +0.61 / −0.20 | +0.18 / −0.03 | −0.03 / −0.15 |
| PIK3CA/mTOR | **−1.07 / −0.34** | **−1.26 / −0.30** | −0.43 / **−1.76** |

**匹配后口袋匹配 vs 错口袋 summary_min**

| 靶对 | 家族 | n_min | 口袋匹配 | 错口袋 | 错口袋 ≥ 匹配？ |
|------|------|------:|---------:|-------:|:---------------:|
| AChE/BChE | unmatched | 20 | 0.618 | 0.642 | 是 |
| AChE/BChE | potency_matched | 18 | 0.593 | 0.642 | 是 |
| AChE/BChE | size_matched | 9 | 0.407 | 0.432 | 是 |
| PIK3CA/PIK3CB | unmatched | 19 | 0.425 | 0.520 | 是 |
| PIK3CA/PIK3CB | potency_matched | 11 | 0.363 | 0.562 | 是 |
| PIK3CA/PIK3CB | size_matched | 13 | 0.302 | 0.426 | 是 |
| PIK3CA/mTOR | unmatched | 20 | 0.765 | 0.788 | 是 |
| PIK3CA/mTOR | potency_matched | 12 | 0.715 | 0.734 | 是 |
| PIK3CA/mTOR | size_matched | 12 | 0.715 | 0.818 | 是 |

效价与尺寸匹配均不翻转“错口袋 ≥ 口袋匹配”。PIK3CA/mTOR holdout 比主面板更弱（不是更强），抽样偏移存在但不足以解释悖论。详见 Results 3.5。

---

## Table S14. GNINA mode_01 与全 9 姿态公平重打的敏感性（worst-pocket，零新对接）

来源：`data/jcim_bench_v0/tables/gnina_mode01_vs_best9_auroc.csv`；脚本 `scripts/compare_gnina_mode01_vs_best9.py`；结论 `analysis/GNINA_BEST9_STATUS.md`。2026-08-24：用户本地对已冻结 K=4 面板的全部 9 个 Vina 姿态分别做 GNINA CNN 重打分（`--cnn_scoring rescore --minimize`），取每端最高 CNNscore，与 RTM 的 best-of-9 覆盖对齐；mode_01 结果保留为历史备份。本表用 `min(score_A, score_B)` 同时代入 dual 对 A_only 与 dual 对 B_only 两个对比（**worst-pocket**，与 `gnina_cnn_min`/`vina_worst`/`rtm_worst` 同一约定），**不是** Methods 2.8 的方向性口袋匹配定义；后者见 Table S15。

| 靶对 | n | worst-pocket mode01 | worst-pocket best9 | Δ | mode_01 是最佳姿态的比例 |
|------|---:|---------------------:|---------------------:|-----:|-----------------------:|
| AChE/BChE | 84 | 0.372 | 0.359 | −0.013 | 0.194 |
| PIK3CA/PIK3CB | 84 | 0.506 | 0.434 | −0.073 | 0.271 |
| PIK3CA/mTOR | 44 | 0.564 | 0.595 | +0.032 | 0.292 |
| EGFR/HER2 | 98 | 0.263 | 0.265 | +0.001 | 0.286 |

mode_01 是 9 个姿态中 CNNscore 最高的比例仅 19–29%，即多数配体的最佳姿态并非 Vina 排名第一姿态；但汇总的 worst-pocket AUROC 变化很小（−0.07 至 +0.03），说明姿态覆盖不对称本身不是此前 GNINA 结论偏弱的主要原因。

---

## Table S15. GNINA 真口袋匹配（Methods 2.8 定义），mode_01 与全 9 姿态对照（零新对接）

来源：`data/jcim_bench_v0/tables/gnina_pocket_matched_mode01_vs_best9_k4_v1.csv` + `..._stability_v1.csv`；脚本 `scripts/gnina_pocket_matched_best9_v1.py`；结论 `analysis/GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`。定义与 Vina/RTM 的 `pocket_matched_vina`/`pocket_matched_rtm` 完全一致：dual 对 A_only 用口袋 B 的 GNINA CNNscore；dual 对 B_only 用口袋 A 的 CNNscore；取两者较小值为 summary_min，配体层 bootstrap 95% CI（B = 2000，种子 20260729）。mode_01 数值与建面时冻结的历史值精确一致（PM48 0.5794、PM110 0.5222），确认备份数据与脚本口径无误。

**K = 4（正文冻结集）**

| 靶对 | 通道 | n (D/A/B) | D vs A（口袋B） | D vs B（口袋A） | summary_min [95% CI] | 同面板 Vina 参考 |
|------|------|-----------|----------------:|----------------:|----------------------:|------------------:|
| EGFR/HER2 | mode01 | 28/38/32 | 0.490 | 0.327 | 0.327 [0.19, 0.46] | 0.430 |
| EGFR/HER2 | best9 | 28/38/32 | 0.471 | 0.290 | **0.290** [0.16, 0.42] | 0.430 |
| AChE/BChE | mode01 | 27/25/28 | 0.486 | 0.442 | 0.442 [0.29, 0.54] | 0.606 |
| AChE/BChE | best9 | 27/25/28 | 0.553 | 0.413 | **0.413** [0.25, 0.55] | 0.606 |
| PIK3CA/PIK3CB | mode01 | 28/27/28 | 0.607 | 0.554 | 0.554 [0.39, 0.68] | 0.500 |
| PIK3CA/PIK3CB | best9 | 28/27/28 | 0.570 | 0.533 | **0.533** [0.37, 0.64] | 0.500 |
| PIK3CA/mTOR | mode01 | 18/14/12 | 0.579 | 0.671 | 0.579 [0.36, 0.75] | 0.692 |
| PIK3CA/mTOR | best9 | 18/14/12 | 0.655 | 0.685 | **0.655** [0.44, 0.81] | 0.692 |

**稳定性核对面板（PM48 / PM110，非独立验证）**

| 面板 | 通道 | n (D/A/B) | D vs A（口袋B） | D vs B（口袋A） | summary_min [95% CI] |
|------|------|-----------|----------------:|----------------:|----------------------:|
| PM48 | mode01 | 18/14/12 | 0.579 | 0.671 | 0.579 [0.36, 0.75] |
| PM48 | best9 | 18/14/12 | 0.655 | 0.685 | **0.655** [0.43, 0.81] |
| PM110 | mode01 | 30/30/30 | 0.522 | 0.713 | 0.522 [0.38, 0.67] |
| PM110 | best9 | 30/30/30 | 0.613 | 0.682 | **0.613** [0.46, 0.74] |

全 9 姿态公平重打后，GNINA 真口袋匹配 summary_min 相对 mode_01 的变化为 −0.04（EGFR/HER2）、−0.03（AChE/BChE）、−0.02（PIK3CA/PIK3CB）、+0.08（PIK3CA/mTOR，含 PM48/PM110）。EGFR/HER2、AChE/BChE、PIK3CA/mTOR 三对上 **GNINA best-of-9 均未超过同面板 Vina 口袋匹配**；EGFR/HER2 与 AChE/BChE 上 GNINA best-of-9 仍低于随机（<0.5）。**PIK3CA/PIK3CB 是例外**：GNINA 口袋匹配（mode01 0.554、best9 0.533）略高于同面板 Vina（0.500），但该关系在 mode-1 时已存在，不是 best-of-9 新产生的现象，且二者均接近随机、bootstrap 区间明显重叠，不构成“GNINA 优于 Vina”的主张，仅表示二者在该对上统计不可分。PIK3CA/mTOR（含 PM48、PM110）在 GNINA 通道下点估计上升但仍不超过 Vina。姿态覆盖对齐后，`RTMScore 与 GNINA 未改变这一格局`（Results 3.2）这一表述继续成立（三对不超过 Vina、一对与 Vina 统计不可分），且现在有方向性（非仅池化/worst-pocket）GNINA 数字支持。PM48/PM110 稳定性核对文本（Results 3.4）与 `PM110_VS_PM48.md`/`B_GROUP_VERDICT.md` 的 GNINA 引用值已改为本表 best9 数值；mode_01 保留仅作追溯校验。

---

## Table S16. 终点层级与主面板 / holdout 对照（写作冻结）

来源：`endpoint_hierarchy_v1.csv`、`frozen_vs_holdout_v1.csv`、`unified_threshold_sensitivity_v2.csv`（θ = 6.0）、`holdout_pocket_matched_v1.csv`。主终点只有一个：统一 θ = 6.0 口袋匹配 Vina `summary_min`。PM48 是 PIK3CA/mTOR 的主面板，不是第二套主指标。EGFR/HER2 holdout = not eligible。

**A. 终点层级（摘录）**

| 角色 | 终点 | 报告位置 |
|------|------|----------|
| primary | θ = 6.0 口袋匹配 Vina summary_min | Table 2；Figure 4A；Figure S4 |
| pre-specified secondary | 方向臂 D/A、D/B；RTM；GNINA best-of-9；最强描述符 | Table 2；Figure 4；Figure S4 |
| robustness | θ 网格、PM110、E=8、holdout、换晶、错口袋配对 Δ | Figure 5–6；S1；S3；S5；Table S4/S8/S9/S17 |
| exploratory | ECFP4、contact_count（非 PLIF）、vina_mean Top-10 | Figure 7；S3D；Table S11/S20/S21 |
| 不作主指标 | pooled vina_mean（EGFR 0.2824 ≠ 0.4297） | Table S6 |

**B. 主面板 vs unused-pool holdout（口袋匹配 Vina）**

| 靶对 | 集合 | n (D/A/B) | D vs A | D vs B | summary_min [95% CI] |
|------|------|-----------|-------:|-------:|----------------------|
| EGFR/HER2 | 主面板 | 28/38/32 | 0.6664 | 0.4297 | 0.4297 [0.2818, 0.5775] |
| EGFR/HER2 | holdout | — | — | — | not eligible |
| AChE/BChE | 主面板 | 27/25/28 | 0.6504 | 0.6058 | 0.6058 [0.4370, 0.7303] |
| AChE/BChE | holdout | 20/20/20 | 0.635 | 0.6175 | 0.6175 [0.4216, 0.7593] |
| PIK3CA/PIK3CB | 主面板 | 28/27/28 | 0.6905 | 0.5 | 0.5 [0.3502, 0.6495] |
| PIK3CA/PIK3CB | holdout | 20/19/20 | 0.7658 | 0.425 | 0.425 [0.2406, 0.6184] |
| PIK3CA/mTOR | 主面板 | 18/14/12 | 0.7143 | 0.6921 | 0.6921 [0.4702, 0.8133] |
| PIK3CA/mTOR | holdout | 20/20/20 | 0.86 | 0.765 | 0.765 [0.6025, 0.8911] |

主面板 CI 来自 `unified_threshold_sensitivity_v2.csv`（Table 2）；holdout CI 来自 `holdout_pocket_matched_v1.csv`。不得把 holdout 写成跨库外部验证。

---

## Table S17. 错口袋配对 Δ bootstrap（同一配体样本）

来源：`wrong_pocket_paired_delta_bootstrap_v1.csv`。B = 2000，seed 20260729。点估计 Δ = 四位小数的 matched − wrong，与 Table 2 / Figure 6 算术一致。Figure S3A–B 绘本表。

| 集合 | 靶对 | matched | wrong | Δ | 95% CI | 不含 0？ |
|------|------|--------:|------:|--:|--------|---------|
| 主面板 | EGFR/HER2 | 0.4297 | 0.26 | 0.1697 | [0.06, 0.2803] | yes |
| 主面板 | AChE/BChE | 0.6058 | 0.4444 | 0.1614 | [0.037, 0.269] | yes |
| 主面板 | PIK3CA/PIK3CB | 0.5 | 0.3489 | 0.1511 | [−0.0215, 0.3105] | no |
| 主面板 | PIK3CA/mTOR | 0.6921 | 0.6019 | 0.0902 | [−0.1222, 0.2626] | no |
| holdout | AChE/BChE | 0.6175 | 0.6425 | −0.025 | [−0.1119, 0.0714] | no |
| holdout | PIK3CA/PIK3CB | 0.425 | 0.52 | −0.095 | [−0.2814, 0.1143] | no |
| holdout | PIK3CA/mTOR | 0.765 | 0.7875 | −0.0225 | [−0.1165, 0.079] | no |

不得把 holdout 点估计反转写成“CI 已排除 0”。EGFR/HER2 无 holdout。

---

## Table S18. 打分通道姿态覆盖（pose fairness）

来源：`pose_fairness_channels_v1.csv`。同一受体、同一盒子、同一组 9 个 Vina 姿态。

| 通道 | 生成姿态 | 实际打分 | 每口袋聚合 | Table 2？ |
|------|----------|----------|------------|----------|
| Vina 1.2.7 | 9 | mode 1（最负 E） | \(S=-E\) | 是（主终点） |
| RTMScore | 9 | 全部 9 | max RTM | 否（次级） |
| GNINA CNN mode-1 | 9 | 仅 mode 1 | CNNscore（minimize） | 否（历史） |
| GNINA CNN best-of-9 | 9 | 全部 9 | max CNNscore | 否（次级；与 RTM 对齐） |

best9 − mode01 为 −0.04 至 +0.08，**不是**相对 Vina。

---

## Table S19. 口袋匹配 Vina 对最强描述符的配对 Δ

来源：`pocket_matched_vs_best_descriptor_delta_v1.csv`。**不是** `baseline_gate_bootstrap_v1.csv`（后者用 pooled `vina_mean`）。Figure S3C。

| 靶对 | 描述符 | Vina | 描述符 | Δ | 95% CI | 不含 0？ |
|------|--------|-----:|-------:|--:|--------|---------|
| EGFR/HER2 | cLogP | 0.4297 | 0.4821 | −0.0524 | [−0.2, 0.1155] | no |
| AChE/BChE | TPSA | 0.6058 | 0.7333 | −0.1275 | [−0.3039, 0.0493] | no |
| PIK3CA/PIK3CB | heavy | 0.5 | 0.6217 | −0.1217 | [−0.3197, 0.0891] | no |
| PIK3CA/mTOR | heavy | 0.6921 | 0.463 | 0.2291 | [−0.0105, 0.4352] | no |

---

## Table S20. ECFP4 支架 GroupKFold 对随机 StratifiedKFold（泄漏核对）

来源：`ligand_ml_scaffold_vs_random_v1.csv`。支架折是主 ML 读出；随机折不是为了找更大 gap。八个方向对比的 mean(random − scaffold) = 0.0258。Figure S3D。

| 靶对 | 对比 | scaffold | random | Δ(random−scaffold) | 对接口袋匹配 |
|------|------|--------:|-------:|-------------------:|-------------:|
| EGFR/HER2 | D vs A | 0.7453 | 0.7961 | +0.0508 | 0.6664 |
| EGFR/HER2 | D vs B | 0.8895 | 0.8884 | −0.0011 | 0.4297 |
| AChE/BChE | D vs A | 0.8948 | 0.9096 | +0.0148 | 0.6504 |
| AChE/BChE | D vs B | 0.8214 | 0.8241 | +0.0027 | 0.6058 |
| PIK3CA/PIK3CB | D vs A | 0.7817 | 0.8042 | +0.0225 | 0.6905 |
| PIK3CA/PIK3CB | D vs B | 0.7691 | 0.8890 | +0.1199 | 0.5000 |
| PIK3CA/mTOR | D vs A | 0.7619 | 0.7262 | −0.0357 | 0.7143 |
| PIK3CA/mTOR | D vs B | 0.8889 | 0.9213 | 0.0324 | 0.6921 |

---

## Table S21. Top-10 硬负计数（探索性；pooled vina_mean，非 Table 2）

来源：`ranking_top10_vina_mean_exploratory_v1.csv` ← `top10_hardneg_bootstrap_v1.csv` 的 `vina_mean` 行。这是排序读出，**不是**口袋匹配主指标。EGFR/HER2 的 Top-10 中有 9 个硬负。

| 靶对 | n_top10 dual | A_only | B_only | hardneg | hardneg bootstrap mean [95% CI] |
|------|-------------:|-------:|-------:|--------:|--------------------------------|
| EGFR/HER2 | 1 | 5 | 4 | 9 | 8.9215 [7, 10] |
| AChE/BChE | 6 | 3 | 1 | 4 | 4.374 [1, 8] |
| PIK3CA/PIK3CB | 3 | 4 | 3 | 7 | 6.519 [3, 9] |
| PIK3CA/mTOR | 6 | 2 | 2 | 4 | 4.0775 [1, 8] |

---

## Table S22. Benchmark-formulation comparison on the same frozen Vina scores

来源：`formulation_conventional_vs_directional_v1.csv`。方向性两臂使用口袋匹配分数；Dual-versus-neither 与 Dual versus all non-duals 使用 pooled `vina_mean`。这些任务的负样本集合不同，均为描述性对照，不是配对显著性检验。`vina_worst` 是较差口袋分数的 AND-like 辅助读出。单靶式类比行保留在源 CSV 中。

| 靶对 | D/A (pocket B) | D/B (pocket A) | directional min | D vs neither mean | D vs neither worst | D vs all non-duals | n neither |
|------|---------------:|---------------:|----------------:|------------------:|-------------------:|-------------------:|----------:|
| EGFR/HER2 | 0.6664 | 0.4297 | 0.4297 | 0.7560 | 0.7440 | 0.5514 | 12 |
| AChE/BChE | 0.6504 | 0.6058 | 0.6058 | 0.6494 | 0.6765 | 0.5792 | 15 |
| PIK3CA/PIK3CB | 0.6905 | 0.5000 | 0.5000 | 0.5592 | 0.6384 | 0.5558 | 16 |
| PIK3CA/mTOR | 0.7143 | 0.6921 | 0.6921 | 0.5139* | 0.4028* | 0.6741 | 4 |

\* neither n = 4，underpowered；不解释为反向效应。

---

## Table S23. Chemotype-constrained selectivity hard negatives

来源：`chemotype_matched_hardneg_v1.csv`。表中 constrained 使用每个硬负相对任一 dual 的最大 ECFP4 Tanimoto ≥ 0.3；distant 为 < 0.3。T ≥ 0.3 只是 similarity-constrained subset，不是 chemically matched analogue set。T ≥ 0.4/0.5 的完整结果见源 CSV，许多格子 n_neg ≤ 7；T ≥ 0.7 的匹配集合为空。

| 靶对 | 对比 | all AUROC (n_neg) | T ≥ 0.3 AUROC (n_neg) | T < 0.3 AUROC (n_neg) |
|------|------|------------------:|-----------------------:|-----------------------:|
| EGFR/HER2 | D vs A | 0.6664 (38) | 0.6548 (27) | 0.6948 (11) |
| EGFR/HER2 | D vs B | 0.4297 (32) | 0.4257 (25) | 0.4439 (7) |
| AChE/BChE | D vs A | 0.6504 (25) | 0.5714 (7) | 0.6811 (18) |
| AChE/BChE | D vs B | 0.6058 (28) | 0.5320 (11) | 0.6536 (17) |
| PIK3CA/PIK3CB | D vs A | 0.6905 (27) | 0.5032 (11) | 0.8192 (16) |
| PIK3CA/PIK3CB | D vs B | 0.5000 (28) | 0.5114 (22) | 0.4583 (6) |
| PIK3CA/mTOR | D vs A | 0.7143 (14) | 0.4815 (3) | 0.7778 (11) |
| PIK3CA/mTOR | D vs B | 0.6921 (12) | 0.6667 (6) | 0.7176 (6) |

---

## Table S24. Incremental information from docking beyond ECFP4

来源：`incremental_information_v1.csv`。模型均使用相同 Bemis–Murcko scaffold GroupKFold。Δ = AUROC(ECFP4+docking) − AUROC(ECFP4)。该 logistic docking AUROC 与 Table 2 的 rank AUROC 不是同一估计量。最大绝对变化为 0.0198，若干方向为负。

| 靶对 | 对比 | ECFP4 | ECFP4+docking | Δ | Table 2 rank docking |
|------|------|------:|--------------:|---:|---------------------:|
| EGFR/HER2 | D vs A | 0.7453 | 0.7509 | +0.0056 | 0.6664 |
| EGFR/HER2 | D vs B | 0.8895 | 0.8873 | −0.0022 | 0.4297 |
| AChE/BChE | D vs A | 0.8948 | 0.8933 | −0.0015 | 0.6504 |
| AChE/BChE | D vs B | 0.8214 | 0.8082 | −0.0132 | 0.6058 |
| PIK3CA/PIK3CB | D vs A | 0.7817 | 0.7857 | +0.0040 | 0.6905 |
| PIK3CA/PIK3CB | D vs B | 0.7691 | 0.7717 | +0.0026 | 0.5000 |
| PIK3CA/mTOR | D vs A | 0.7619 | 0.7421 | −0.0198 | 0.7143 |
| PIK3CA/mTOR | D vs B | 0.8889 | 0.8981 | +0.0092 | 0.6921 |

---

## Table S25. Mixed-library Top-10 composition under pooled `vina_mean`

来源：`mixed_library_enrichment_v1.csv`。这是探索性排序读出，不是 Table 2 的方向主终点。EF5、EF10、`vina_worst` 与单口袋排序见源 CSV。

| 靶对 | library n (dual n) | Top-10 dual | A_only | B_only | neither | EF Top-10 | hard-negative fraction |
|------|-------------------:|------------:|-------:|-------:|--------:|----------:|-----------------------:|
| EGFR/HER2 | 110 (28) | 1 | 5 | 4 | 0 | 0.393 | 0.900 |
| AChE/BChE | 95 (27) | 5 | 3 | 1 | 1 | 1.759 | 0.400 |
| PIK3CA/PIK3CB | 99 (28) | 3 | 4 | 2 | 1 | 1.061 | 0.600 |
| PIK3CA/mTOR | 48 (18) | 6 | 2 | 1 | 1 | 1.600 | 0.300 |

---

## Table S26. summary_min vs arithmetic, geometric, and harmonic means

来源：`aggregation_min_mean_geometric_harmonic_v1.csv`。主终点仍为 min。四对排序在四种聚合下完全相同（PM > AChE > PIK3CB > EGFR）。EGFR Dual-versus-neither（0.756）相对 min / arithmetic / geometric / harmonic 的差分别为 +0.326 / +0.208 / +0.221 / +0.234，方向不变。该对照不是配对显著性检验。

| 靶对 | D/A | D/B | min | arithmetic | geometric | harmonic | Dual vs neither | rank (all four) |
|------|----:|----:|----:|-----------:|----------:|---------:|----------------:|----------------:|
| EGFR/HER2 | 0.6664 | 0.4297 | 0.4297 | 0.5481 | 0.5351 | 0.5225 | 0.756 | 4 |
| AChE/BChE | 0.6504 | 0.6058 | 0.6058 | 0.6281 | 0.6277 | 0.6273 | 0.6494 | 2 |
| PIK3CA/PIK3CB | 0.6905 | 0.5000 | 0.5000 | 0.5953 | 0.5876 | 0.5800 | 0.5592 | 3 |
| PIK3CA/mTOR | 0.7143 | 0.6921 | 0.6921 | 0.7032 | 0.7031 | 0.7030 | 0.5139* | 1 |

\* neither n = 4，underpowered。

---

## Table S27. Docking attempted / successful / failed

来源：`docking_failure_census_v1.csv`、`docking_failed_ligand_properties_v1.csv`、`docking_failure_rank_extreme_v1.csv`；脚本 `docking_failure_sensitivity_v1.py`。主 AUROC 以两端均得分的配体为条件。HOAP_028 为 AutoDock 原子类型 `B`（硼）覆盖失败，不是 silent missingness。

| 集合 | 靶对 | attempted | both-end success | fail either | fail A | fail B |
|------|------|----------:|-----------------:|------------:|-------:|-------:|
| 主面板 | EGFR/HER2 | 110 | 110 | 0 | 0 | 0 |
| 主面板 | AChE/BChE | 100 | 95 | 5 | 4 | 5 |
| 主面板 | PIK3CA/PIK3CB | 100 | 99 | 1 | 1 | 0 |
| 主面板 | PIK3CA/mTOR | 48 | 48 | 0 | 0 | 0 |
| holdout | AChE/BChE | 60 | 60 | 0 | 0 | 0 |
| holdout | PIK3CA/PIK3CB | 60 | 59 | 1 | 1 | 1 |
| holdout | PIK3CA/mTOR | 60 | 60 | 0 | 0 | 0 |

AChE 主面板失败：AB_001（dual，两端）、AB_053/054/056（A_only，两端）、AB_097（neither，B 端）。PIK3CB 主面板：PAB_034（A_only，A 端，`timeout_900s_torsdof=23`）。holdout：HOAP_028 两端硼原子类型失败。PAB_034 在 4JPS/5DXT 替换中同样超时（Table S30），不是标签过滤。

### S27b. 主面板失败配体的化学覆盖

描述符按去盐规则保留最大重原子片段。四个 AChE/BChE 失败是协议在 torsdof ≥25 时主动跳过；另两个是超时。因此缺失不是随机计算噪声。

| 靶对 | ligand | 类别 | A/B score | heavy | MW | cLogP | TPSA | charge | rotatable | 原因 |
|------|--------|------|-----------|------:|---:|------:|-----:|-------:|----------:|------|
| AChE/BChE | AB_001 | dual | 0/0 | 85 | 1151.4 | 8.88 | 234.1 | 0 | 29 | torsdof=31，协议跳过 |
| AChE/BChE | AB_053 | A-only | 0/0 | 56 | 801.1 | 10.74 | 40.6 | +2 | 27 | torsdof=31，协议跳过 |
| AChE/BChE | AB_054 | A-only | 0/0 | 48 | 667.0 | 7.65 | 65.6 | 0 | 27 | torsdof=29，协议跳过 |
| AChE/BChE | AB_056 | A-only | 0/0 | 50 | 693.1 | 9.32 | 40.6 | +2 | 27 | torsdof=29，协议跳过 |
| AChE/BChE | AB_097 | neither | 1/0 | 54 | 729.0 | 9.39 | 75.4 | 0 | 17 | B 端 600 s timeout |
| PIK3CA/PIK3CB | PAB_034 | A-only | 0/1 | 60 | 957.7 | 4.96 | 170.5 | 0 | 22 | A 端 900 s timeout |

### S27c. 方向终点的 arm-available 与 rank-extreme 缺失敏感性

arm-available 使用该方向所需口袋的全部已有分数，不要求另一端也成功。rank-extreme lower/upper 把每个涉及缺失所需口袋分数的比较全部判为逆向/顺向；这是确定性边界，不是插补模型，也没有 bootstrap CI。

| 靶对 | 方向 | complete n+/n− | complete AUROC | arm-available n+/n− | arm AUROC | rank-extreme [lower, upper] |
|------|------|----------------:|---------------:|--------------------:|----------:|----------------------------:|
| AChE/BChE | D vs A, pocket B | 27/25 | 0.6504 | 27/25 | 0.6504 | [0.5599, 0.6990] |
| AChE/BChE | D vs B, pocket A | 27/28 | 0.6058 | 27/28 | 0.6058 | [0.5842, 0.6199] |
| PIK3CA/PIK3CB | D vs A, pocket B | 28/27 | 0.6905 | 28/28 | 0.6952 | [0.6952, 0.6952] |
| PIK3CA/PIK3CB | D vs B, pocket A | 28/28 | 0.5000 | 28/28 | 0.5000 | [0.5000, 0.5000] |

这些敏感性不改变当前 pair-level `summary_min` 判读，但只支持“可被该协议处理的化学空间”。

---

## Table S28. 四个预先指定描述符的方向 AUROC

来源：`descriptor_all_four_directional_v1.csv`。全部报告；最高者为 best single-descriptor reference，不是 confirmatory competitor。

| 靶对 | heavy min | MW min | cLogP min | TPSA min | best reference |
|------|----------:|-------:|----------:|---------:|----------------|
| EGFR/HER2 | 0.3694 | 0.4163 | **0.4821** | 0.4275 | cLogP |
| AChE/BChE | 0.5820 | 0.5785 | 0.4669 | **0.7333** | TPSA |
| PIK3CA/PIK3CB | **0.6217** | 0.6204 | 0.5952 | 0.4180 | heavy |
| PIK3CA/mTOR | **0.4630** | 0.4484 | 0.3102 | 0.2599 | heavy |

---

## Table S29. Max vs median pChEMBL (full scored panels; θ = 6.0)

来源：`assay_max_vs_median_agreement_v1.csv`。冻结 Vina 分数不重算。分母是 scored n。正文 A4 只报告标签一致率与 pair-level Δsummary_min。EGFR 冻结 Table 2（0.430）与 API-max（0.417）差 1 个配体（EH120_060 / CHEMBL24828）；配体层缓存/API 不一致只记录于本表。数值 max≠median 计数只放本表。

| 靶对 | n scored | 类别翻转 | 标签一致率 | 数值 max≠median | 冻结 summary_min | API-max min | API-median min |
|------|--------:|--------:|----------:|----------------:|-----------------:|------------:|---------------:|
| EGFR/HER2 | 110 | 7 | 103/110 = 93.6% | 40/110 | 0.430 (28/38/32) | 0.417 (29/37/32) | 0.424 (26/35/33) |
| AChE/BChE | 95 | 1 | 94/95 = 98.9% | 13/95 | 0.606 | 0.606 | 0.629 |
| PIK3CA/PIK3CB | 99 | 1 | 98/99 = 99.0% | 25/99 | 0.500 | 0.500 | 0.500 |
| PIK3CA/mTOR | 48 | 0 | 48/48 = 100% | 27/48 | 0.692 | 0.692 | 0.692 |

翻转清单见 `assay_max_vs_median_flips_v1.csv`。EGFR：EH40_08 dual→A_only；EH40_27、EH120_070/076/077 A_only→neither；EH120_041 dual→B_only；EH120_060 API-max dual vs median A_only。AChE：AB_018 dual→neither。PIK3CB：PAB_053 A_only→neither。PM：无。

---

## Table S30. Two-pair PIK3CA receptor-realization effect (B pocket frozen)

来源：`receptor_realization_two_pair_v1.csv`。CI 用沉积 CSV，不用临时重算。PAB_034：100 尝试 / 99 成功 / 1 超时（原始 4L23 与 4JPS、5DXT 均失败）。

| 靶对 | PIK3CA | 保留 B | attempted / success / fail | D/A | D/B | summary_min [95% CI] | Δ |
|------|--------|--------|---------------------------:|----:|----:|----------------------|--:|
| PIK3CA/mTOR | 4L23 | 4JT6 | 48/48/0 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | — |
| PIK3CA/mTOR | 4JPS | 4JT6 | 48/48/0 | 0.714 | 0.486 | 0.486 [0.259, 0.692] | −0.206 |
| PIK3CA/mTOR | 5DXT | 4JT6 | 48/48/0 | 0.714 | 0.505 | 0.505 [0.292, 0.696] | −0.187 |
| PIK3CA/mTOR | 4L23 | 4JSX | 48/48/0 | 0.639 | 0.692 | 0.639 [0.418, 0.776] | −0.053 |
| PIK3CA/PIK3CB | 4L23 | 2WXF | 100/99/1 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | — |
| PIK3CA/PIK3CB | 4JPS | 2WXF | 100/99/1 | 0.691 | 0.707 | 0.691 [0.516, 0.779] | +0.191 |
| PIK3CA/PIK3CB | 5DXT | 2WXF | 100/99/1 | 0.691 | 0.685 | 0.685 [0.506, 0.768] | +0.185 |

PIK3CA/PIK3CB 弱臂：原始 D/B = 0.500；4JPS 后弱臂切到冻结 D/A = 0.691；5DXT 两臂接近平衡。同一 PIK3CA 扰动、方向相反。不是普遍定律。

---

## Table S31. Detectable-effect simulation for `summary_min`

来源：`detectable_effect_simulation_v1.csv`；脚本 `scripts/detectable_effect_simulation_v1.py`。双正态分数模型；配体层 bootstrap 与 Methods 2.4 相同（B = 2000，seed 20260729）；N_MC = 1000。单元格为 **P(95% CI 排除 0.5)**。这不是观察后功效。`summary_min` 在真实 AUROC = 0.50 时点估计偏低（min 的偏倚），因此空假设下排除 0.5 的概率可略高于两臂单独的 ~0.05；完整四对照、含 Dual versus neither，见源 CSV。Figure S6。

| Pair | n_scored (dual / A / B) | true 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |
|------|------------------------:|----------:|-----:|-----:|-----:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.025 | 0.065 | 0.268 | 0.621 | 0.907 |
| AChE/BChE | 27 / 25 / 28 | 0.020 | 0.049 | 0.225 | 0.504 | 0.828 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.032 | 0.041 | 0.226 | 0.564 | 0.849 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.037 | 0.025 | 0.072 | 0.219 | 0.452 |

当前样本量更容易分辨较大的方向性效应。CI 未能排除 0.5 不能写成与随机等价。

---

## Table S32. Independent GNINA pose generation versus frozen Vina (θ = 6.0)

来源：`data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv`、`independent_dock_summary_v1.csv`、`independent_dock_enrichment_v1.csv`；脚本 `scripts/run_gnina_independent_dock_v1.py`、`analyze_independent_dock_v1.py`；结论 `analysis/INDEPENDENT_DOCK_VERDICT_V1.md`。GNINA 1.3.2 对接搜索（mode-1 `minimizedAffinity`，种子 20260727），**不是**对冻结 Vina 姿态的 CNN 重打分（后者仍见 Tables S14–S15）。复用冻结 Meeko 配体、受体坐标、对接盒与 exhaustiveness。失败配体：EGFR/HER2 neither `EH120_109`（两端）；PIK3CA/mTOR A-only `PM48_19`（两端）。分析只用两端均有分数的配体。主文四舍五入至三位；未四舍五入值见源 CSV。

不允许写成“GNINA 优于/劣于 Vina”。允许句：EGFR/HER2 上设定效应在该独立姿态生成协议下不是 Vina 特有的。

| Pair | Engine | n_scored (dual / A-only / B-only) | Dual vs neither (`mean`) [95% CI] | n_neither | D/A (pocket B) | D/B (pocket A) | summary_min [95% CI] |
|------|--------|----------------------------------:|----------------------------------:|----------:|---------------:|---------------:|----------------------|
| EGFR/HER2 | Vina (frozen) | 28 / 38 / 32 | 0.756 [0.562, 0.920] | 12 | 0.666 | 0.430 | 0.430 [0.282, 0.578] |
| EGFR/HER2 | GNINA dock | 28 / 38 / 32 | 0.783 [0.610, 0.922] | 11 | 0.660 | 0.220 | 0.220 [0.109, 0.343] |
| PIK3CA/mTOR | Vina (frozen) | 18 / 14 / 12 | 0.514 [0.222, 0.806] | 4 | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| PIK3CA/mTOR | GNINA dock | 18 / 13 / 12 | 0.569 [0.222, 0.889] | 4 | 0.633 | 0.704 | 0.633 [0.427, 0.825] |

EGFR/HER2 mixed-library Top-10 by GNINA mean pocket score: 1 dual, 4 A-only, 5 B-only, 0 neither（EF10 = 0.389）。Vina `vina_mean` Top-10 为 1 / 5 / 4 / 0（Table S25；EF10 = 0.393）。PIK3CA/mTOR Dual versus neither 仍因 n_neither = 4 而效能不足，不作为设定对照主读出。

---

## Table S33. PIK3CA pocket-contact occupancy shifts (4JPS / 5DXT vs 4L23)

来源：`data/jcim_structure_robust_v0/analysis/plif_v1/plif_residue_shift_top10_v1.csv`；脚本 `scripts/run_receptor_plif_v1.py`；结论 `analysis/plif_v1/PLIF_VERDICT_V1.md`。同一套 PM48 配体（n = 48）在 4L23 / 4JPS / 5DXT 的 mode-1 姿态上，对 20 个冻结口袋残基做重原子距离 ≤ 4.5 Å 的几何占有率。这是 SOP 允许的 ProLIF 等价快照，不是完整 ProLIF 指纹，也不是残基层因果。占有率为有接触配体比例。下表为 |Δ| 最大的 10 个残基。

不允许写成“残基 X 导致 AUROC 变化”或“PLIF 解释了 PIK3CA/PIK3CB 的相反位移”。

| Residue | 4L23 | 4JPS | 5DXT | Δ 4JPS−4L23 | Δ 5DXT−4L23 |
|---------|-----:|-----:|-----:|------------:|------------:|
| Met772 | 1.000 | 0.771 | 0.958 | −0.229 | −0.042 |
| Leu807 | 0.583 | 0.396 | 0.396 | −0.188 | −0.188 |
| Gln859 | 0.625 | 0.646 | 0.812 | +0.021 | +0.188 |
| Thr856 | 0.688 | 0.500 | 0.542 | −0.188 | −0.146 |
| Cys838 | 0.062 | 0.229 | 0.000 | +0.167 | −0.062 |
| Glu849 | 0.917 | 0.854 | 0.750 | −0.062 | −0.167 |
| Phe930 | 0.438 | 0.604 | 0.479 | +0.167 | +0.042 |
| Asp933 | 1.000 | 0.854 | 0.938 | −0.146 | −0.062 |
| Asp810 | 0.729 | 0.583 | 0.667 | −0.146 | −0.062 |
| Trp780 | 0.958 | 0.833 | 0.917 | −0.125 | −0.042 |

---

## Table S34. 固定口袋分数后的负类定义对照

来源：`data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv`；脚本 `scripts/benchmark_formulation_v1.py`。每个对照保持口袋分数不变，仅将方向性 selective 负类替换为 neither。每次 bootstrap 共享 dual 抽样，两个负类分层独立抽样（B = 2000）。这是两个 AUROC 估计量差值的描述性区间，不是负配体层面的配对检验；8 个方向未做多重性校正。

| 靶对 | 固定口袋 | selective 负类 | n dual/selective/neither | AUROC D vs selective | AUROC D vs neither | Δ neither−selective [95% CI] |
|------|----------|----------------|-------------------------:|---------------------:|-------------------:|-------------------------------:|
| EGFR/HER2 | B | A-only | 28/38/12 | 0.666 | 0.720 | 0.054 [−0.157, 0.246] |
| EGFR/HER2 | A | B-only | 28/32/12 | 0.430 | 0.808 | 0.378 [0.205, 0.547] |
| AChE/BChE | B | A-only | 27/25/15 | 0.650 | 0.709 | 0.058 [−0.085, 0.196] |
| AChE/BChE | A | B-only | 27/28/15 | 0.606 | 0.590 | −0.016 [−0.154, 0.120] |
| PIK3CA/PIK3CB | B | A-only | 28/27/16 | 0.691 | 0.501 | −0.189 [−0.421, 0.052] |
| PIK3CA/PIK3CB | A | B-only | 28/28/16 | 0.500 | 0.659 | 0.159 [0.001, 0.321] |
| PIK3CA/mTOR | B | A-only | 18/14/4 | 0.714 | 0.583 | −0.131 [−0.433, 0.179]* |
| PIK3CA/mTOR | A | B-only | 18/12/4 | 0.692 | 0.472 | −0.220 [−0.514, 0.019]* |

\* n_neither = 4，underpowered。

---

## Table S35. 按实验状态分组的测量频次诊断

来源：`data/jcim_novelty_v0/tables/measurement_frequency_by_class_v1.csv` 和 `measurement_frequency_max_median_v1.csv`；脚本 `measurement_frequency_audit_v1.py`。activity records 是 ChEMBL API 返回的记录数，不是独立实验次数，也不等于独立文献数。本表不包含 assay confidence、species 或 document-level 去重，因此是 profiling-intensity 诊断，不是高置信标签重建。

| 靶对 | 类别 | n | 双端 activity records 中位数 [IQR] | 最大值 | A 端重复 pChEMBL 比例 | B 端重复 pChEMBL 比例 |
|------|------|--:|--------------------------------------:|-------:|--------------------------:|--------------------------:|
| EGFR/HER2 | dual | 28 | 4 [2–6.5] | 318 | 0.536 | 0.429 |
| EGFR/HER2 | A-only | 38 | 4 [2–8] | 62 | 0.474 | 0.158 |
| EGFR/HER2 | B-only | 32 | 4 [2–5.25] | 11 | 0.125 | 0.188 |
| EGFR/HER2 | neither | 12 | 3 [2–3] | 4 | 0.333 | 0.000 |
| AChE/BChE | dual | 27 | 2 [2–3] | 207 | 0.185 | 0.074 |
| AChE/BChE | A-only | 25 | 2 [2–4] | 11 | 0.280 | 0.000 |
| AChE/BChE | B-only | 28 | 2 [2–4] | 6 | 0.107 | 0.143 |
| AChE/BChE | neither | 15 | 2 [2–2] | 3 | 0.000 | 0.067 |
| PIK3CA/PIK3CB | dual | 28 | 2.5 [2–5] | 46 | 0.286 | 0.286 |
| PIK3CA/PIK3CB | A-only | 27 | 2 [2–4] | 21 | 0.222 | 0.037 |
| PIK3CA/PIK3CB | B-only | 28 | 2 [2–3] | 15 | 0.036 | 0.429 |
| PIK3CA/PIK3CB | neither | 16 | 2 [2–2] | 6 | 0.000 | 0.000 |
| PIK3CA/mTOR | dual | 18 | 22 [10.75–34] | 90 | 0.944 | 0.944 |
| PIK3CA/mTOR | A-only | 14 | 3 [2–4.75] | 28 | 0.500 | 0.143 |
| PIK3CA/mTOR | B-only | 12 | 3 [2–6.75] | 40 | 0.167 | 0.500 |
| PIK3CA/mTOR | neither | 4 | 2 [2–2.5] | 4 | 0.250 | 0.250 |

API-max 与 median 差值和重复记录数有预期的正相关，但该诊断包含差值必然为 0 的单次记录，因此不将相关系数作为因果证据。主要结论是类别间 profiling intensity 不均衡，尤其是 PIK3CA/mTOR dual 类。

---

## Table S36. 当前 ChEMBL 快照的高置信标签稳健性

抓取时间：2026-08-26 UTC。来源：`high_confidence_activity_audit_v1.csv`、`high_confidence_labels_v1.csv`、`high_confidence_summary_v1.csv`；脚本 `high_confidence_label_rebuild_v1.py`；运行元数据和 API 缓存见 `cache/high_confidence_v1/`。

保留规则：Homo sapiens `SINGLE PROTEIN` target；assay confidence score ≥ 8；`standard_relation = "="`；endpoint 属于 IC50/Ki/Kd/EC50/Potency；无 `data_validity_comment`；`potential_duplicate = 0`。共审计 2748 条 activity records，保留 1546 条；513 条因 potential duplicate 标记被排除。完整的逐记录排除理由见源 CSV。

| 靶对 | 冻结 scored n | 高置信双端完整 n | 与冻结类别一致 | dual/A-only/B-only/neither | 高置信 summary_min |
|------|--------------:|----------------------:|------------------:|---------------------------:|-----------------------:|
| EGFR/HER2 | 110 | 110 | 110/110 | 28/38/32/12 | 0.430 |
| AChE/BChE | 95 | 95 | 95/95 | 27/25/28/15 | 0.606 |
| PIK3CA/PIK3CB | 99 | 99 | 99/99 | 28/27/28/16 | 0.500 |
| PIK3CA/mTOR | 48 | 48 | 48/48 | 18/14/12/4 | 0.692 |

这一 post-hoc current-database 稳健性视图说明上述显式记录过滤不改变当前面板分类。它不是 2026-07-23 数据库状态的重建，也未统一 assay condition、protein construct、mutation context 或 source-document sampling；不得称为 assay-harmonized ground truth。

---

## Table S37. 完整病例覆盖与来源文献集中度

来源：`complete_case_usable_pchembl_overlap_v1.csv`、`source_document_concentration_v1.csv`；脚本 `complete_case_document_audit_v1.py`。覆盖率分母为至少一个靶点存在可用 pChEMBL 的结构并集；“A-only measured / B-only measured”表示另一端在冻结可用值映射中缺失，**不表示无活性**。来源集中度使用 2026-08-26 高置信视图中保留的 activity records；record 不是独立实验重复。

| 靶对 | A 有值 | B 有值 | 两端有值 | 仅 A measured | 仅 B measured | 两端/并集 |
|------|-------:|-------:|---------:|--------------:|--------------:|----------:|
| EGFR/HER2 | 11198 | 2619 | 1751 | 9447 | 868 | 0.145 |
| AChE/BChE | 6197 | 3798 | 2537 | 3660 | 1261 | 0.340 |
| PIK3CA/PIK3CB | 7732 | 2786 | 1990 | 5742 | 796 | 0.233 |
| PIK3CA/mTOR | 7732 | 5209 | 2713 | 5019 | 2496 | 0.265 |

| 靶对 | 类别 | n ligands | retained records | unique documents | top-document record fraction | top-document ligand fraction |
|------|------|----------:|-----------------:|-----------------:|-----------------------------:|-----------------------------:|
| EGFR/HER2 | dual | 28 | 335 | 102 | 0.134 | 0.071 |
| EGFR/HER2 | A-only | 38 | 185 | 65 | 0.141 | 0.105 |
| EGFR/HER2 | B-only | 32 | 74 | 16 | 0.162 | 0.188 |
| EGFR/HER2 | neither | 12 | 28 | 8 | 0.429 | 0.333 |
| AChE/BChE | dual | 27 | 129 | 77 | 0.047 | 0.074 |
| AChE/BChE | A-only | 25 | 58 | 20 | 0.121 | 0.120 |
| AChE/BChE | B-only | 28 | 64 | 20 | 0.125 | 0.107 |
| AChE/BChE | neither | 15 | 31 | 12 | 0.194 | 0.200 |
| PIK3CA/PIK3CB | dual | 28 | 84 | 25 | 0.119 | 0.071 |
| PIK3CA/PIK3CB | A-only | 27 | 66 | 20 | 0.136 | 0.185 |
| PIK3CA/PIK3CB | B-only | 28 | 75 | 15 | 0.160 | 0.143 |
| PIK3CA/PIK3CB | neither | 16 | 32 | 11 | 0.313 | 0.313 |
| PIK3CA/mTOR | dual | 18 | 264 | 118 | 0.083 | 0.278 |
| PIK3CA/mTOR | A-only | 14 | 59 | 30 | 0.119 | 0.286 |
| PIK3CA/mTOR | B-only | 12 | 58 | 19 | 0.293 | 0.083 |
| PIK3CA/mTOR | neither | 4 | 8 | 1 | 1.000 | 1.000 |

该审计量化了完整病例选择与来源集中风险，不能推断未测结构的真实活性。文献阻断与时间分割见 Tables S39–S41。

---

## Table S38. 类别化学空间与骨架诊断

来源：`class_chemistry_summary_v1.csv`；脚本 `class_chemistry_audit_v1.py`。数值为中位数 [IQR]；NN-dual 为 radius-2、2048-bit ECFP4 到同靶对 dual 类的最大 Tanimoto（dual 行为最近的另一个 dual）。scaffolds 为类别内 Bemis–Murcko 数；singleton fraction 为落在仅出现一次骨架中的配体比例。完整 CSV 另含 heavy atoms、formal charge、rotatable bonds 的中位数与 IQR；四对所有类别的 formal-charge 中位数与 IQR 均为 0。该表是 post-hoc 描述性混杂审计，不做多重未校正检验。

| 靶对 | 类别 | n | MW [IQR] | cLogP [IQR] | TPSA [IQR] | scaffolds | singleton fraction | median NN-dual |
|------|------|--:|----------|--------------|------------|----------:|-------------------:|---------------:|
| EGFR/HER2 | dual | 28 | 481.6 [439.5–539.3] | 4.88 [3.73–5.75] | 95.2 [86.2–106.3] | 26 | 0.857 | 0.421 |
| EGFR/HER2 | A-only | 38 | 433.8 [352.8–508.8] | 4.65 [3.42–5.54] | 80.2 [65.9–96.4] | 30 | 0.658 | 0.392 |
| EGFR/HER2 | B-only | 32 | 512.8 [475.4–548.7] | 4.90 [4.19–5.46] | 97.1 [91.0–109.7] | 26 | 0.688 | 0.367 |
| EGFR/HER2 | neither | 12 | 373.0 [323.5–416.7] | 3.76 [2.82–5.04] | 67.6 [52.1–99.5] | 10 | 0.667 | 0.154 |
| AChE/BChE | dual | 27 | 471.6 [366.4–559.9] | 5.64 [3.26–6.98] | 76.0 [51.4–95.7] | 26 | 0.926 | 0.474 |
| AChE/BChE | A-only | 25 | 409.5 [380.5–490.4] | 4.50 [3.83–5.54] | 38.8 [32.8–71.9] | 19 | 0.560 | 0.236 |
| AChE/BChE | B-only | 28 | 425.6 [389.5–459.4] | 5.18 [4.36–6.69] | 48.1 [39.2–59.1] | 25 | 0.821 | 0.253 |
| AChE/BChE | neither | 15 | 422.5 [358.6–443.6] | 4.64 [4.04–5.47] | 50.8 [30.5–59.1] | 12 | 0.667 | 0.236 |
| PIK3CA/PIK3CB | dual | 28 | 477.7 [408.2–548.9] | 3.69 [2.64–4.32] | 97.5 [78.8–123.4] | 24 | 0.750 | 0.607 |
| PIK3CA/PIK3CB | A-only | 27 | 417.6 [395.0–482.1] | 3.14 [2.52–3.87] | 105.1 [100.2–121.6] | 22 | 0.704 | 0.250 |
| PIK3CA/PIK3CB | B-only | 28 | 396.2 [371.4–442.7] | 2.13 [1.65–2.55] | 88.0 [76.0–98.1] | 20 | 0.464 | 0.419 |
| PIK3CA/PIK3CB | neither | 16 | 377.1 [320.1–497.0] | 3.43 [2.62–4.28] | 92.1 [81.6–100.1] | 14 | 0.750 | 0.177 |
| PIK3CA/mTOR | dual | 18 | 421.5 [388.4–491.3] | 2.31 [1.65–3.16] | 102.5 [82.8–114.8] | 18 | 1.000 | 0.367 |
| PIK3CA/mTOR | A-only | 14 | 427.5 [400.4–513.0] | 2.52 [1.31–3.50] | 124.2 [105.2–134.5] | 14 | 1.000 | 0.226 |
| PIK3CA/mTOR | B-only | 12 | 443.0 [403.5–467.5] | 3.27 [2.64–3.97] | 97.2 [93.1–114.7] | 11 | 0.833 | 0.306 |
| PIK3CA/mTOR | neither | 4 | 328.6 [292.6–369.7] | 2.60 [2.39–2.98] | 103.4 [93.8–107.6] | 4 | 1.000 | 0.395 |

---

## Table S39. 文献阻断交叉验证与 document-cluster bootstrap

来源：`document_blocked_cv_summary_v1.csv`、`document_cluster_bootstrap_v1.csv`；脚本 `document_blocked_cv_v1.py`。共享任一保留高置信 `document_id` 的配体为同一组。分组规则在看到 AUROC 后未更改。rank AUROC 为冻结 Vina 口袋匹配分数；OOF 模型使用相同 `GroupKFold`。PIK3CA/mTOR Dual versus B-only 只有 1 个同时含两类的折，报告为无法稳定估计。

| 靶对 | 对比 | n_pos/n_neg | groups | documents | valid folds | rank full | rank mean-fold | ECFP4 OOF | physchem OOF | dock logistic OOF | doc-cluster 95% CI | status |
|------|------|-------------|-------:|----------:|------------:|----------:|---------------:|----------:|-------------:|------------------:|--------------------|--------|
| EGFR/HER2 | D_vs_A | 28/38 | 28 | 156 | 5 | 0.666 | 0.734 | 0.627 | 0.726 | 0.582 | [0.464, 0.795] | ok |
| EGFR/HER2 | D_vs_B | 28/32 | 23 | 113 | 5 | 0.430 | 0.497 | 0.623 | 0.336 | 0.469 | [0.321, 0.617] | ok |
| AChE/BChE | D_vs_A | 27/25 | 39 | 97 | 5 | 0.650 | 0.662 | 0.847 | 0.596 | 0.615 | [0.497, 0.809] | ok |
| AChE/BChE | D_vs_B | 27/28 | 41 | 95 | 5 | 0.606 | 0.654 | 0.767 | 0.726 | 0.540 | [0.429, 0.753] | ok |
| PIK3CA/PIK3CB | D_vs_A | 28/27 | 33 | 43 | 5 | 0.691 | 0.734 | 0.652 | 0.497 | 0.536 | [0.460, 0.868] | ok |
| PIK3CA/PIK3CB | D_vs_B | 28/28 | 30 | 38 | 5 | 0.500 | 0.501 | 0.834 | 0.755 | 0.337 | [0.310, 0.698] | ok |
| PIK3CA/mTOR | D_vs_A | 18/14 | 8 | 138 | 3 | 0.714 | 0.604 | 0.514 | 0.521 | 0.590 | [0.400, 0.886] | ok |
| PIK3CA/mTOR | D_vs_B | 18/12 | 9 | 129 | 1 | 0.692 | — | — | — | — | [0.000, 0.818] | cannot_stably_estimate |

---

## Table S40. 文献阻断折明细

来源：`document_blocked_cv_folds_v1.csv`。每个测试折必须同时含正负两类才会进入 OOF AUROC。训练与测试的 group overlap 必须为 0。

完整折级 n_pos/n_neg/n_documents 见源 CSV。

---

## Table S41. 冻结文献年份分割

来源：`time_split_class_counts_v1.csv`、`time_split_auroc_v1.csv`、`document_year_lookup_v1.csv`；协议见 `docs/TIME_SPLIT_PROTOCOL_FREEZE.md`。配体年份 = 保留高置信记录中最早的 `document.year`。主截止年 2018 在计算 AUROC 之前冻结。n < 10 的格子不报告 AUROC，也不为凑类别改截止年。主截止年可评估靶对为 0，**不包装为外部验证**。

| cutoff | 靶对 | 测试 dual/A-only/B-only/neither | gate | D_vs_A | D_vs_B |
|-------:|------|--------------------------------:|------|--------|--------|
| 2015 | EGFR/HER2 | 9/8/15/6 | descriptive_only | — | — |
| 2015 | AChE/BChE | 11/11/24/10 | underpowered_report | 0.603 | 0.636 |
| 2015 | PIK3CA/PIK3CB | 17/20/4/8 | descriptive_only | — | — |
| 2015 | PIK3CA/mTOR | 6/2/1/0 | descriptive_only | — | — |
| 2018 | EGFR/HER2 | 6/3/14/2 | descriptive_only | — | — |
| 2018 | AChE/BChE | 8/5/15/6 | descriptive_only | — | — |
| 2018 | PIK3CA/PIK3CB | 12/11/0/3 | unevaluable | — | — |
| 2018 | PIK3CA/mTOR | 2/0/1/0 | unevaluable | — | — |
| 2020 | EGFR/HER2 | 6/0/14/2 | unevaluable | — | — |
| 2020 | AChE/BChE | 5/3/14/6 | descriptive_only | — | — |
| 2020 | PIK3CA/PIK3CB | 6/4/0/1 | unevaluable | — | — |
| 2020 | PIK3CA/mTOR | 1/0/0/0 | unevaluable | — | — |

---

## Table S42. Assay-context 优先分子审计与元数据纳入/排除

来源：`assay_context_priority_ligands_v1.csv`、`assay_context_audit.csv`；脚本 `assay_context_audit_v1.py` 与 `assay_context_human_review_v1.py`。机器提取做风险分层；随后一次元数据审核在 ChEMBL assay 自由文本不可用时填写纳入/排除。`protein_construct` 与 `wildtype_or_mutant` 全部为 `unknown`。冻结标签未改。优先集合：EGFR/HER2 全部 dual/A-only/B-only（98）、PIK3CA/mTOR neither（4）、混合端点（40）、生化+功能（22）、对主 AUROC 影响最大的分子（20）及高集中文献系列（57）。人工 SOP：`docs/ASSAY_CONTEXT_HUMAN_REVIEW_SOP.md`。审核摘要：`analysis/ASSAY_CONTEXT_HUMAN_REVIEW_SUMMARY_V1.md`。

| 标记 | n ligands |
|------|----------:|
| egfr_directional_priority | 98 |
| top_document_series | 57 |
| mixed_endpoint | 40 |
| borderline_pA | 27 |
| biochem_and_functional | 22 |
| high_auroc_influence | 20 |
| borderline_pB | 20 |
| pm_neither_single_document | 4 |
| non_human_assay_organism | 3 |
| 优先配体合计（去重） | 186 |
| include / uncertain / exclude | 179 / 7 / 0 |
| 相对冻结类别翻转 | 0 |

7 个 uncertain 配体（EGFR/HER2 EH120_045；PIK3CA/mTOR PM48_04、PM48_05、PM48_22；AChE/BChE AB_089、AB_091、AB_094）保留冻结类别。该步骤不是文献级构建体/突变核查。

---

## Table S43. BindingDB REST 结构与文献独立性计数（历史供给核对）

来源：`bindingdb_independence_summary_v1.csv`；脚本 `bindingdb_independence_audit_v1.py`；SOP `docs/BINDINGDB_EXTERNAL_SOP.md`；结论 `analysis/BINDINGDB_INDEPENDENCE_VERDICT.md`。规则在查看独立类别计数前冻结：UniProt 锁定、等式关系 IC50/Ki/Kd/EC50、θ = 6.0。相对已打分面板的 InChIKey 去重在本地完成。本表保留为 REST pmax JSON 的历史供给计数，不是文献独立切片。原生归档重建见 Tables S48–S49。**未对接，不计算 AUROC，不包装为外部验证。**

| 靶对 | BindingDB 两端 dual/A/B/neither | 与已打分面板 InChIKey 重叠 n | 去掉面板结构后 dual/A/B/neither | 面板结构门槛 | 文献+ChEMBL 图门槛 | 包装为外部验证 |
|------|--------------------------------:|-----------------------------:|--------------------------------:|--------------|--------------------|:--------------:|
| EGFR/HER2 | 1621/186/92/370 | 96 | 1589/161/62/361 | supply_enough_to_dock | unevaluable | 0 |
| AChE/BChE | 988/467/253/1003 | 76 | 966/450/230/989 | supply_enough_to_dock | unevaluable | 0 |
| PIK3CA/PIK3CB | 1371/286/283/605 | 90 | 1341/261/257/596 | supply_enough_to_dock | unevaluable | 0 |
| PIK3CA/mTOR | 2027/251/279/182 | 44 | 2008/238/266/182 | supply_enough_to_dock | unevaluable | 0 |

原生归档重建（Tables S48–S49）已替代“完成本表的 UniChem + PMID 去重后再对接”这一后续步骤。不得为追逐 AUROC 而改阈值或只对接有利子集。

---

## Table S44. θ = 6.0 四状态标签普查

来源：`theta6_pair_census_v1.csv`。复用 J0 49 对（去掉顺序别名）与缓存 `mols_*.json`。方向门槛：dual/A-only/B-only 均 n ≥ 10。设定门槛：另加 neither n ≥ 10。**17 对同时通过两个门槛；其中 16 对非金属酶。已对接 4 对。不得写成 17 对对接基准。**

## Table S45. 物化 caliper 匹配方向 AUROC

来源：`property_caliper_match_v1.csv`。特征为 z 标准化 MW/cLogP/TPSA/重原子；1:1 贪心；caliper 0.5/1.0 SD。n_matched < 8 为效能不足。1.0 SD Dual-versus-B-only：EGFR/HER2 0.566（n = 16）、AChE/BChE 0.462（n = 13）、PIK3CA/PIK3CB 0.556（n = 9）；PIK3CA/mTOR n = 5 效能不足。不是因果调整。

## Table S46. AND 式双口袋工作点

来源：`and_filter_operating_point_v1.csv`。库 = Dual+A-only+B-only。EGFR/HER2 `vina_worst` Dual 中位截断：precision 0.298，硬负比例 0.702；第 90 百分位 precision 0.130。不是对 DualDiff/FuseDiff 的重打分。

## Table S47. 四对完整 ChEMBL 图的配体层 AUROC

来源：`ligand_only_fullmap_auroc_v1.csv`。每类最多 120 个分子；支架 GroupKFold。EGFR/HER2 Dual versus neither 0.921，Dual versus B-only 0.864，ECFP4 `summary_min` 0.801。**不是对接结果，不替换 Table 2。**

---

## Table S48. BindingDB-native 候选分层计数

来源：`external_candidate_flow.csv`；脚本 `bindingdb_native_slice_v1.py`；合约 `protocol/external_slice_contract.yaml`；归档锁 `bindingdb_archive_lock_v1.csv`（release 202608）。规则在对接前冻结。ChEMBL 文献解析 519/680，故后续剩余计数是完全文献独立集的上界。**未对接，不计算 AUROC。**

| 靶对 | 层 | n | dual/A/B/neither |
|------|----|--:|-----------------:|
| EGFR/HER2 | native_paired_theta6 | 468 | 371/30/54/13 |
| EGFR/HER2 | drop_shared_literature | 380 | 303/16/51/10 |
| EGFR/HER2 | drop_shared_structure | 346 | 284/11/44/7 |
| EGFR/HER2 | drop_neighbors_lt_0.70 | 216 | 180/10/20/6 |
| EGFR/HER2 | lt_0.70_plus_0.50_0.70_sensitivity | 327 | 269/11/40/7 |
| AChE/BChE | native_paired_theta6 | 318 | 159/46/37/76 |
| AChE/BChE | drop_shared_literature | 225 | 93/32/27/73 |
| AChE/BChE | drop_shared_structure | 170 | 62/19/23/66 |
| AChE/BChE | drop_neighbors_lt_0.70 | 85 | 4/8/14/59 |
| AChE/BChE | lt_0.70_plus_0.50_0.70_sensitivity | 136 | 39/13/20/64 |
| PIK3CA/PIK3CB | native_paired_theta6 | 296 | 114/29/4/149 |
| PIK3CA/PIK3CB | drop_shared_literature | 296 | 114/29/4/149 |
| PIK3CA/PIK3CB | drop_shared_structure | 199 | 81/11/4/103 |
| PIK3CA/PIK3CB | drop_neighbors_lt_0.70 | 112 | 9/0/3/100 |
| PIK3CA/PIK3CB | lt_0.70_plus_0.50_0.70_sensitivity | 147 | 40/4/3/100 |
| PIK3CA/mTOR | native_paired_theta6 | 1193 | 1000/115/30/48 |
| PIK3CA/mTOR | drop_shared_literature | 1192 | 999/115/30/48 |
| PIK3CA/mTOR | drop_shared_structure | 958 | 888/41/16/13 |
| PIK3CA/mTOR | drop_neighbors_lt_0.70 | 98 | 91/4/1/2 |
| PIK3CA/mTOR | lt_0.70_plus_0.50_0.70_sensitivity | 659 | 621/25/6/7 |
| MCL1/Bcl-xL | native_paired_theta6 | 90 | 32/22/6/30 |
| MCL1/Bcl-xL | drop_shared_literature | 90 | 32/22/6/30 |
| MCL1/Bcl-xL | drop_shared_structure | 16 | 11/1/3/1 |
| MCL1/Bcl-xL | drop_neighbors_lt_0.70 | 3 | 1/0/2/0 |
| MCL1/Bcl-xL | lt_0.70_plus_0.50_0.70_sensitivity | 9 | 5/1/3/0 |

---

## Table S49. BindingDB-native 切片主门槛摘要

来源：`external_slice_summary_v1.csv`；结论 `analysis/EXTERNAL_SLICE_FREEZE.md`。主门槛：dual/A-only/B-only 各 n ≥ 20、每类至少 3 个来源、最大单文献配体份额 ≤ 50%。薄复制门槛不把任何一对标为可对接。`packaged_as_external_evaluation = 0` 对全部五对。剩余 n 是上界。

| 靶对 | 过滤后 dual/A/B/neither | 来源 dual/A/B | gate | packaged |
|------|------------------------:|--------------:|------|:--------:|
| EGFR/HER2 | 180/10/20/6 | 16/5/4 | insufficient | 0 |
| AChE/BChE | 4/8/14/59 | 2/6/3 | insufficient | 0 |
| PIK3CA/PIK3CB | 9/0/3/100 | 4/0/1 | insufficient | 0 |
| PIK3CA/mTOR | 91/4/1/2 | 9/2/1 | insufficient | 0 |
| MCL1/Bcl-xL | 1/0/2/0 | 1/0/1 | insufficient | 0 |

主门槛对数 **0**；薄复制 **0**；未对接。停止规则：少于两对主外部靶对，保持四靶对评价设定审计。不得把剩余配体称为数据库外部集。

---

## Table S50. MCL1/Bcl-xL 面板冻结

来源：`mcl1_bclxl_panel_freeze_v1.csv`、`mcl1_bclxl_chembl_panel96_v1.csv`。ChEMBL θ = 6.0 图 305 个配对结构：dual/A/B/neither 82/77/24/122。冻结面板 24/24/24/24（种子 20260729）；B-only 穷尽。正式 LC6 pose-gold 未有效完成后，按预先声明对接为 **applicability stress-test**（`docked=1`）。同源 BCL-2 折叠上的 PPI/BH3 槽域位移，不是异质折叠对，也不是“首次非激酶对”，**不是第五个 Table 2 靶对**。

---

## Table S51. MCL1/Bcl-xL 受体与 LC6 初步坐标诊断（非正式 pose-gold 验证）

来源：`mcl1_bclxl_receptor_freeze_v1.csv`、`data/mcl1_bclxl_panel_v0/tables/cognate_qc_lc6_v1.csv`。链来自 RCSB polymer entity 1。选择依据为分辨率与野生型占位，不是 AUROC。Gate：两端 best-of-top3 RMSD < 2.0 Å（Vina seed 20260727，exhaustiveness 8，九个姿态）。

| 靶标 | PDB | 角色 | 分辨率 (Å) | 链 | top1 Å | best-top3 Å | 2 Å 数值筛选 |
|------|-----|------|----------:|----|-------:|------------:|------|
| MCL1 | 3WIY | primary | 2.15 | A | 1.689 | 1.689 | < 2 Å |
| BCL2L1 | 3WIZ | primary | 2.45 | A | 4.17 | 2.011 | ≥ 2 Å |
| MCL1 | 6UDV | alternate | 1.35 | A | — | — | 未跑 |
| BCL2L1 | 3SP7 | alternate | 1.4 | A | — | — | 未跑 |

上述数值来自按元素全局 Hungarian 坐标匹配；该映射没有保持分子图同构，不能称为 topology-aware symmetry-corrected RMSD。预声明的 physical-validity、关键相互作用恢复和第二随机种子也未完成。因此正式 pose-gold gate 视为 **未满足/未有效完成**，表中的 2 Å 筛选不构成姿态验证证据。不得包装为标准筛选性能或 domain extension。

---

## Table S52. 文献对照（不是 bake-off）

来源：`benchmark_literature_comparator_v1.csv`。把 DualFourClass 的负类定义、配对测定、化学对照、受体敏感性与 Zhou 2013、DUD-E、LIT-PCBA、CASF-2016、DOCKSTRING 并置。不得写成引擎竞赛或“首次非激酶基准”。

---

## Table S53. MCL1/Bcl-xL Vina applicability stress-test（不进 Table 2）

来源：`data/mcl1_bclxl_panel_v0/tables/formulation_auroc_MBX_v1.csv`。LC6 gate FAIL 后按预先声明对接。93/96 配体两端有分；`panel_role=applicability_stress_test`。

| contrast | n_pos / n_neg | AUROC | 95% CI |
|----------|-------------:|------:|--------|
| Dual vs neither (mean) | 23 / 24 | 0.628 | [0.462, 0.786] |
| Dual vs A-only @ Bcl-xL | 23 / 24 | 0.793 | [0.655, 0.915] |
| Dual vs B-only @ MCL1 | 23 / 22 | 0.609 | [0.439, 0.776] |
| `summary_min` | 23 / — | 0.609 | — |

Dual-versus-B-only 区间包含 0.5。不得写成 PPI 筛选成功、第五对主结果或外部验证。

---

## Table S54. 五种子 Vina 敏感性（与 Table 3 同一 Dual-versus-neither 估计量）

来源：`data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv`、`multiseed_auroc_aggregate_v2.csv`。生产种子 20260727 复用冻结分数，另四个冻结种子为 20260811–20260814。配体、受体、对接盒、exhaustiveness、姿态数与分析规则不变。Dual versus neither 为逐配体 `vina_mean` 再算一个 AUROC，与 Table 3 相同；主种子恢复 0.756 / 0.649 / 0.559 / 0.514。v1 的 `mean(AUC_A, AUC_B)` 列只作敏感性，不引用。

方向性 `summary_min` 中位数（范围）：EGFR/HER2 0.373（0.321–0.430），AChE/BChE 0.599（0.553–0.606），PIK3CA/PIK3CB 0.478（0.468–0.502），PIK3CA/mTOR 0.704（0.676–0.726）。EGFR/HER2 设定差距在 5/5 种子上为正。PIK3CA/mTOR Dual versus neither 仍因 n_neither = 4 效能不足，不解释为反转。不得用多种子均值替换 Table 2，也不得按观察到的性能挑选种子。

---

## Supporting Note S1. Exploratory PIK3CA/mTOR pose-level diagnostics

来源：`data/pik3ca_mtor_panel48_v0/analysis/failure_typology_v0/`，仅为代表性案例，不是全面板 PLIF 或机制分析。

- **T2 / PM48_21（A_only）**：pChEMBL PIK3CA/mTOR = 8.70/5.92。RTM-best pose 在 4L23 与 4JT6 均 hinge-positive、无 clash，共晶位占用分别为 1.00/0.97；其弱端仍形成几何上干净的 ATP-site pose。因此“看似合理的双端 pose”不等于实验 dual。
- **T5 / Torin1 与 omipalisib（dual）**：Vina pooled 排名分别为 1/3，但 RTM-best PIK3CA pose 偏离 hinge/cognate family，RTM `min_z` 排名降至 31/30。该案例只说明 pose selection / rescoring 可改变排序，不支持 residue-level mechanism。

完整字段与限制见 `CASE_PM48_21_Aonly.md` 和 `CASE_PM48_10_02_injured_duals.md`。

---


- 本文件是**已有数据的汇编**，不是新实验。若某分析尚无机器可读表，宁缺毋填。
- 投稿英文 SI 时：Table 编号可按期刊习惯重排；数字不得改动。
- Cognate 表必须同时报告 mode1 与 best_of_9，避免审稿人误读“全部 &lt; 2 Å”。
- 早期借用 Schrodinger 处理过的姿态对照**不写入投稿稿**（无正式使用权限；主协议已统一为 RDKit/meeko）。仓库内 `pm48_directional_by_prep_v1.csv` 仅作内部记录。
- ChEMBL median：全面板 A4 已完成（Table S29）。不得把 27 配体诊断样当成 SI 表。当前数据库的 confidence≥8 / 物种 / relation / validity / duplicate 过滤见 Table S36，但不是冻结日期重建或 assay harmonization。
- Table S12 是计数核对（BindingDB REST + PubChem PUG REST），不是对接结果；不得把 `as_is` 的 EGFR ≥50 写成已建成 BindingDB 厚面板。
- Table S13 是 holdout 效价/尺寸匹配诊断，不替换 Table S8；不得写成错口袋悖论已解决。
- Table S16–S21 是冻结分数上的补表（零新对接）。S17 的 holdout Δ CI 均含 0；S19 四对描述符 Δ CI 均含 0；S21 是 vina_mean Top-10，不是 Table 2。
- Table S22–S54 来自 `data/jcim_novelty_v0/`、`data/jcim_structure_robust_v0/`、`data/jcim_independent_dock_v0/` 与 `data/jcim_multiseed_v0/`：S22 formulation comparison（主文 Figure 3）；S23 chemotype-constrained hard-negatives（T ≥ 0.7 为空；T ≥ 0.3 不是 analogue matching）；S24 incremental ECFP/docking；S25 mixed-library EF；S26 min/arithmetic/geometric/harmonic 聚合敏感性（四对排序不变）；S27 docking N_attempted/success/fail；S28 四个描述符全报；S29 max vs median；S30 两对 PIK3CA receptor-realization；S31 detectable-effect simulation；S32 独立 GNINA 姿态生成；S33 PIK3CA 几何占有率位移；S34 固定口袋负类对照；S35 测量频次；S36 当前 ChEMBL 高置信视图；S37 完整病例覆盖与来源集中度；S38 类别化学空间；S39–S40 文献阻断 CV；S41 冻结时间分割（主截止年不可包装为外部验证）；S42 assay-context 元数据审核（无冻结标签翻转）；S43 BindingDB REST 独立性计数（历史供给；未对接）；S44 θ = 6.0 标签普查（不是对接扩面）；S45 物化 caliper 匹配；S46 AND 过滤工作点；S47 配体层全图 AUROC（不是对接）；S48 BindingDB 原生候选分层；S49 原生切片门槛摘要（0 对通过主门槛；未对接）；S50–S51 MCL1/Bcl-xL 面板与 LC6 gate（FAIL；仅 stress-test）；S52 文献对照（不是 bake-off）；S53 MCL1 Vina AUROC（不进 Table 2）；S54 五种子 Vina 敏感性（Dual-versus-neither = `AUC(vina_mean)`，不能替换主种子结果或用于挑选有利种子；v1 错误口径不引用）。Figure S4 = 口袋匹配森林图；Figure S5 = unused-pool holdout；Figure S6 = detectable-effect heatmap；Figure S7 = 普查/AND/配体层全图；Figure S8 = BindingDB 原生切片过滤漏斗；Figure 8 = diagnostic workflow。不得把 Dual-vs-neither 写成 “conventional benchmark”；不得把 EGFR 0.756 vs 0.430 写成配对显著性；不得把 PIK3CA/mTOR Dual-vs-neither（n = 4）写成反转；不得把受体替换写成单向 collapse 或 robustness。不得把 2015 AChE/BChE 时间分割写成全文外部验证。不得把 BindingDB REST 计数或原生切片写成已对接外部验证。不得把 EGFR/HER2 reconstructed QC 写成历史生产 pose。不得把 17 对标签普查写成 17 对对接基准。不得把 Table S47 写成 Table 2 的替代。不得把 MCL1/Bcl-xL 写成第五个 Table 2 靶对、异质折叠对、首次非激酶对或已通过 pose-gold 的 domain extension。
- Figure S3 不得复用 Figure 6 的 AUROC 柱；它只画配对 Δ ± CI。
