# Supporting Information（中文工作稿 · JCIM Articles）

> 与 [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md)、[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。  
> **原则：** 仅收录仓库中已有实验/冻结产物；数字均可追溯到下列源文件。未做的分析不填空表。

**溯源（机器可读原料）：**

| SI 表 | 源文件 |
|-------|--------|
| Table S1 | `data/jcim_strengthen_t0t1_v0/ENV_PIN.md`；各 panel `protocol.yaml` |
| Table S2 | `data/*/boxes/*.json`（冻结面板所用条目） |
| Table S3 | PM：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT*.md`；AChE/BChE：`cognate_qc/COGNATE_QC.md`；PIK3CB：`cognate_qc/COGNATE_QC.md`；EGFR/HER2：`protocol/protocol.yaml` + `analysis/exhaustiveness_sensitivity_v1/SENSITIVITY_VERDICT.md` |
| Table S4 | `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` |
| Table S5 | `data/jcim_strengthen_t0t1_v0/tables/matched_subset_directional_v1.csv` |
| Table S6 | `data/jcim_bench_v0/tables/pocket_matched_directional_v1.csv`；worst-pocket 亦见 `data/jcim_strengthen_t0t1_v0/tables/aggregation_sensitivity_v1.csv` |
| Table S7 | `data/jcim_bench_v0/analysis/structural_context_v1/full_chain_identity_v1.py`（脚本）与 `full_chain_identity_v1_output.tsv`（输出） |
| Table S8 | `data/jcim_holdout_v0/tables/holdout_panel_*.csv` + `holdout_ligand_scores_v1.csv` + `holdout_pocket_matched_v1.csv`；结论见 `analysis/HOLDOUT_VERDICT.md` |
| Table S9 | `data/jcim_structure_robust_v0/analysis/STRUCTURE_ROBUSTNESS_QC_V1.md` / `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`；`tables/pocket_matched_PM48_alt*_v1.csv` |
| Table S10 | `data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/POCKET_MECHANISM_VERDICT_V1.md` + `pocket_superposition_v1.py`（脚本，零新对接，仅用已冻结晶体坐标） |
| Table S11 | `data/jcim_holdout_v0/analysis/WRONG_POCKET_MECHANISM_VERDICT_V1.md` + `scripts/wrong_pocket_contact_v1.py`（脚本，零新对接，仅用已冻结姿态坐标） |
| Table S12 | `data/jcim_supply_crossdb_v0/tables/crossdb_strict_supply_v1.csv`；结论见 `analysis/SUPPLY_CROSSDB_VERDICT_V1.md`（BindingDB REST + PubChem PUG REST 计数，零对接） |
| ChEMBL 聚合局限 | `data/jcim_strengthen_t0t1_v0/analysis/T0_SKIPS.md` |

---

## Table S1. 计算参数与软件版本

| 项目 | 取值 |
|------|------|
| Python | 3.x（本机 conda） |
| RDKit | 2026.3.1 |
| meeko | 0.7.1 |
| AutoDock Vina | 1.2.7；`scoring_function = vina` |
| GNINA | 1.3.2（CPU，`--no_gpu`） |
| RTMScore 权重 | `rtmscore_model1` |
| Open Babel | Vina PDBQT → SDF（GNINA 前） |
| 分析库 | NumPy / SciPy / scikit-learn / pandas（版本随公开复现环境；仓库 ENV_PIN 未钉死具体小版本） |
| Biopython | `PDBParser` / `Superimposer` / `PairwiseAligner`（序列一致性与 Cα 叠合） |
| ChEMBL | Web API；靶对审计锁定 2026-07-23（未记录 release 编号） |
| 配体准备 | 去盐（最大有机片段）→ RDKit AddHs → ETKDGv3（seed 20260727）→ MMFFOptimizeMolecule（maxIters=200）→ meeko 默认 PDBQT |
| 描述符 | RDKit GetNumHeavyAtoms / MolWt / MolLogP / TPSA |
| 逻辑回归 | scikit-learn `LogisticRegression`（C = 1.0，max_iter = 2000）；支架 `GroupKFold`，折数 min(5, n_pos, n_neg, n_groups) |
| 面板抽样种子 | 20260729 |
| Holdout 抽样种子 | 20260731；仅三对（不含 EGFR/HER2）；PM 排除 PM110 超集 |
| Bootstrap | B = 2000；seed 20260729；配体层 2.5%–97.5% 百分位区间；不做多重比较校正 |
| Vina `n_modes` | 9 |
| Vina `energy_range` | 3 |
| Vina 随机种子 | 20260727 |
| Vina exhaustiveness | PIK3CA/mTOR 主面板/扩面/单靶对照/换晶 = 16；其余靶对与对应 holdout = 8；PM 另报 E = 8 对照 |
| GNINA | `--cnn_scoring rescore --minimize --seed 20260727` |
| 盒子定义 | 共晶配体 AABB + 5 Å；每边下限 20 Å |
| Cognate 通过门槛 | 重原子 RMSD，`best_of_9` &lt; 2.0 Å（同坐标系，不叠合） |
| 受体 PDBQT | PIK3CA/mTOR/EGFR/HER2：含氢蛋白坐标 + `mk_prepare_receptor.py --read_pdb`；AChE/BChE/PIK3CB：沉积 ATOM + `mk_prepare_receptor`（default altloc A） |
| 换晶 | 一次只换一个口袋；未换端保留冻结主面板分数；4JPS/5DXT/4JSX 过 QC，3T8M 排除 |
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

## Table S3. 共晶配体重对接（八个冻结受体）

**均已做 cognate redock。** 协议门槛为 `best_of_9` 重原子 RMSD &lt; 2.0 Å。主 seed = 20260727（EGFR 历史 as-run 与敏感性表同 seed）。

**RMSD 定义：** 对接坐标系、不做蛋白叠合。PIK3CA/mTOR 与 EGFR/HER2：meeko `REMARK SMILES IDX` 映射 + 模板自同构上最小 CalcRMS（见各 panel `rmsd_definition.md`）。AChE/BChE 与 PIK3CB 冻结 QC：重原子坐标匈牙利匹配（`linear_sum_assignment`）。

### S3a. 冻结受体在协议 exhaustiveness 下的结果

| 靶标 | PDB | 配体 | E | RMSD mode1 (Å) | RMSD best_of_9 (Å) | best mode | best_of_9 &lt; 2 Å？ | 面板 E |
|------|-----|------|--:|---------------:|------------------:|----------:|:-------------------:|-------:|
| PIK3CA | 4L23 | X6K | 8 | 0.624 | 0.624 | 1 | 是 | 16（随 mTOR 端升 E） |
| mTOR | 4JT6 | X6K | 8 | 7.118 | **5.003** | 2 | **否** | — |
| mTOR | 4JT6 | X6K | 16 | 7.118 | **0.445** | 3 | 是 | 16 |
| AChE | 4EY7 | E20 | 8 | （mode1 = best） | **0.339** | 1 | 是 | 8 |
| BChE | 4BDS | THA | 8 | （记录为 best_of_9） | **0.386** | — | 是 | 8 |
| PIK3CB | 2WXF | 039 | 8/16 | （记录为 best_of_9） | **0.405** | — | 是 | 8 |
| EGFR | 3POZ | 03P | 8 | **9.483** | **0.955** | — | 是（best_of_9） | 8 |
| HER2 | 3RCD | 03P | 8 | **1.941** | **1.941** | 1 | 是 | 8 |

补充（EGFR 面板历史 as-run，与上表敏感性诊断一致量级）：

| 靶标 | PDB | Vina top1 RMSD (Å) | RTM 选中 mode RMSD (Å) |
|------|-----|-------------------:|----------------------:|
| EGFR | 3POZ | 9.514 | 1.015（mode 2） |
| HER2 | 3RCD | 1.869 | 1.974（mode 3） |

### S3b. 判读（与正文 Methods 2.4 对齐）

1. **八个冻结受体都做过共晶重对接**；另有未入选候选结构的失败记录（见 S3c），不纳入主协议。
2. **在 E = 8、门槛 = best_of_9 &lt; 2 Å 时：** 4L23、4EY7、4BDS、2WXF、3POZ、3RCD 通过；**仅 4JT6 未通过**（5.003 Å）。升至 E = 16 后 4JT6 的 best_of_9 = 0.445 Å，故 PIK3CA/mTOR 全面板采用 E = 16。
3. **不能把“best_of_9 &lt; 2 Å”等同于“Vina mode1 &lt; 2 Å”。** 4JT6 在 E = 16 时 mode1 仍约 7.1 Å；3POZ 在 E = 8/16/32 时 mode1 均约 9.5 Å，近晶构象出现在非 top1 mode。因此协议保留输出 9 个 mode，并在打分对照中使用 best-of-9 / 重打分。
4. EGFR/HER2 面板仍用 E = 8：敏感性显示升 E 不能修复 3POZ 的 mode1 排序失败，且 E = 8 时 best_of_9 已 &lt; 2 Å。

### S3c. 受体筛选中试过但未冻结的候选（实验记录，非主结果）

| 拟用端 | PDB | 配体 | 结果摘要 |
|--------|-----|------|----------|
| BChE | 6ZWI | QRH | best_of_9 ≈ 2.3–2.5 Å @ E8/16，未过门槛 |
| BChE | 6QAA / 5DYW | HUN / 5HF | PDBQT 解析失败 |
| PIK3CB | 2Y3A | GD9 | best_of_9 ≈ 3.85 Å @ E8/16，未过 |
| PIK3CB | 4BFR | J82 | Vina PDBQT 解析失败 |

---

## Table S4. 统一 θ = 6.0 主表述的阈值敏感性网格（支持性）

来源：`unified_threshold_sensitivity_v2.csv`。正文 Table 2 已采用本表 θ = 6.0 行作为四对统一主结果（Results 3.2）；本表 θ = 5.5/6.5 与严格 6.5/5.5 行为支持性阈值敏感性分析，用于说明排序不随阈值网格翻转，不是与 Table 2 竞争的第二套主标准。

| 靶对 | 标签规则 | n (D / A / B) | AUROC D vs A | AUROC D vs B | summary_min | 95% CI | underpowered |
|------|----------|--------------:|-------------:|-------------:|------------:|--------|:------------:|
| EGFR/HER2 | θ = 5.5 | 69 / 22 / 10 | 0.773 | 0.425 | 0.425 | [0.238, 0.622] | 否 |
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 | [0.284, 0.576] | 否 |
| EGFR/HER2 | θ = 6.5 | 26 / 29 / 29 | 0.735 | 0.460 | 0.460 | [0.305, 0.623] | 否 |
| EGFR/HER2 | 严格 6.5/5.5 | 26 / 17 / 7 | 0.799 | 0.324 | 0.324 | [0.130, 0.519] | **是** |
| AChE/BChE | θ = 5.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.446, 0.730] | 否 |
| AChE/BChE | θ = 6.0 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.440, 0.740] | 否 |
| AChE/BChE | θ = 6.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.450, 0.735] | 否 |
| AChE/BChE | 严格 6.5/5.5 | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | [0.449, 0.738] | 否 |
| PIK3CA/PIK3CB | θ = 5.5 | 30 / 25 / 28 | 0.729 | 0.522 | 0.522 | [0.365, 0.676] | 否 |
| PIK3CA/PIK3CB | θ = 6.0 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.347, 0.648] | 否 |
| PIK3CA/PIK3CB | θ = 6.5 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.331, 0.653] | 否 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 | [0.346, 0.647] | 否 |
| PIK3CA/mTOR | θ = 5.5 | 33 / 9 / 5 | 0.502 | 0.506 | 0.502 | [0.248, 0.635] | **是** |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 | [0.464, 0.802] | 否 |
| PIK3CA/mTOR | θ = 6.5 | 17 / 15 / 12 | 0.710 | 0.674 | 0.674 | [0.444, 0.791] | 否 |
| PIK3CA/mTOR | 严格 6.5/5.5 | 17 / 7 / 4 | 0.639 | 0.669 | 0.639 | [0.321, 0.798] | **是** |

说明：AChE/BChE 与 PIK3CA/PIK3CB 建造时已按严格配额冻结，重标后计数基本不变。EGFR/HER2 与 PIK3CA/mTOR 在严格规则下 B_only 过少，仅作稳健性描述。早期 `threshold_sensitivity_v1.csv`（vina_mean 通道）保留在仓库作内部对照，不进正文主敏感性表。

---

## Table S5. 效价 / 尺寸匹配子集（口袋匹配 Vina）

来源：`matched_subset_directional_v1.csv`。效价匹配：\|ΔpChEMBL\| ≤ 0.5；尺寸匹配：\|Δheavy atoms\| ≤ 2。单对比 AUROC 与配体层 bootstrap 95% CI。

| 靶对 | 子集 | n_dual / n_other | AUROC | 95% CI |
|------|------|-----------------:|------:|--------|
| EGFR/HER2 | potency D vs A | 17 / 17 | 0.747 | [0.564, 0.914] |
| EGFR/HER2 | potency D vs B | 14 / 14 | 0.469 | [0.235, 0.699] |
| EGFR/HER2 | size D vs A | 24 / 24 | 0.547 | [0.380, 0.717] |
| EGFR/HER2 | size D vs B | 22 / 22 | 0.519 | [0.353, 0.696] |
| AChE/BChE | potency D vs A | 20 / 20 | 0.593 | [0.403, 0.770] |
| AChE/BChE | potency D vs B | 24 / 24 | 0.601 | [0.443, 0.757] |
| AChE/BChE | size D vs A | 16 / 16 | 0.484 | [0.273, 0.691] |
| AChE/BChE | size D vs B | 15 / 15 | 0.596 | [0.378, 0.796] |
| PIK3CA/PIK3CB | potency D vs A | 20 / 20 | 0.680 | [0.497, 0.855] |
| PIK3CA/PIK3CB | potency D vs B | 20 / 20 | 0.458 | [0.283, 0.653] |
| PIK3CA/PIK3CB | size D vs A | 20 / 20 | 0.630 | [0.443, 0.803] |
| PIK3CA/PIK3CB | size D vs B | 17 / 17 | 0.450 | [0.259, 0.661] |
| PIK3CA/mTOR | potency D vs A | 13 / 13 | 0.710 | [0.485, 0.905] |
| PIK3CA/mTOR | potency D vs B | 9 / 9 | 0.728 | [0.444, 0.975] |
| PIK3CA/mTOR | size D vs A | 12 / 12 | 0.778 | [0.576, 0.938] |
| PIK3CA/mTOR | size D vs B | 12 / 12 | 0.722 | [0.479, 0.924] |

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

## Table S8. 面板外冻结验证集（holdout）配体与口袋匹配分数

来源：`data/jcim_holdout_v0/`。抽样种子 20260731；标签规则 strict 6.5/5.5；对接协议与主面板相同。完整配体级表见 `holdout_ligand_scores_v1.csv`（含 SMILES、类别、两端 Vina mode-1、physchem）。下表为口袋匹配汇总。

| 靶对 | n (D/A/B) | pocket_matched summary_min | 95% CI | D vs A | D vs B | 主面板 | 最强平凡基线 | Δ(dock−基线) |
|------|-----------|---------------------------:|-------:|-------:|-------:|-------:|-------------:|-------------:|
| PIK3CA/mTOR | 60 (20/20/20) | 0.765 | 0.603–0.891 | 0.860 | 0.765 | 0.692 | heavy 0.555 | +0.210 |
| AChE/BChE | 60 (20/20/20) | 0.618 | 0.422–0.759 | 0.635 | 0.618 | 0.606 | cLogP 0.575 | +0.043 |
| PIK3CA/PIK3CB | 59 (20/19/20) | 0.425 | 0.241–0.618 | 0.766 | 0.425 | 0.500 | heavy 0.691 | −0.266 |

失败：`HOAP_028`（含硼，AutoDock 原子类型 `B` 不支持）两端未得分，已剔除。错口袋对照见 `holdout_pocket_matched_v1.csv` 中 `wrong_pocket_control_vina` 行。

---

## Table S9. 结构稳健性：替代受体 cognate QC

来源：`data/jcim_structure_robust_v0/analysis/STRUCTURE_ROBUSTNESS_QC_V1.md`。协议：Vina E=16，seed=20260727，best_of_9 &lt; 2 Å。

| 靶标 | PDB | 共晶配体 | mode1 RMSD (Å) | best_of_9 (Å) | 结论 |
|------|-----|----------|---------------:|--------------:|------|
| PIK3CA | 4JPS | 1LT | 0.607 | 0.607 | PASS |
| PIK3CA | 5DXT | 5H5 | 0.624 | 0.624 | PASS |
| mTOR | 4JSX | 17G (Torin2) | 0.515 | 0.515 | PASS |

PASS 受体已写入 `receptors/`。PM48 单端替换后的口袋匹配 summary_min：

| 替代结构 | 替换口袋 | summary_min [95% CI] | 主面板 | Δ |
|----------|----------|---------------------:|-------:|--:|
| 4JPS | A | 0.486 [0.259, 0.692] | 0.692 | −0.206 |
| 5DXT | A | 0.505 [0.292, 0.696] | 0.692 | −0.187 |
| 4JSX | B | 0.639 [0.418, 0.776] | 0.692 | −0.053 |

详见 `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`。

---

## Table S10. 受体依赖的探索性结构对照：晶体间 Cα 叠合

来源：`data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/POCKET_MECHANISM_VERDICT_V1.md`。方法：Biopython `PDBParser` + `Superimposer`；口袋残基由参考结构自身共晶配体重原子 ≤5 Å 界定，按残基编号+残基名精确匹配（匹配位点零错配）；全域与口袋局域 RMSD 共用同一次刚体拟合。PIK3CA 替代结构 n = 2，mTOR n = 1。

| 参考（主面板） | 替代 | 匹配 Cα 数 | 全域 Cα RMSD (Å) | 口袋残基数 | 口袋局域 Cα RMSD (Å) | 共晶配体质心距离 (Å) |
|------|------|-----------:|------------------:|-----------:|----------------------:|----------------------:|
| 4L23（PIK3CA） | 4JPS | 982 | 1.486 | 20 | 0.867 | 2.566 |
| 4L23（PIK3CA） | 5DXT | 862 | 1.441 | 20 | 0.343 | 2.072 |
| 4JT6（mTOR） | 4JSX | 1054 | 0.454 | 18 | 0.467 | 2.196 |

PIK3CA 口袋残基（4L23，共晶 X6K）：Met772、Trp780、Ile800、Lys802、Leu807、Asp810、Leu814、Tyr836、Cys838、Ile848、Glu849、Val850、Val851、Ser854、Thr856、Gln859、Met922、Phe930、Ile932、Asp933。mTOR 口袋残基（4JT6，共晶 X6K/PI-103）：Ile2163、Pro2169、Leu2185、Lys2187、Glu2190、Leu2192、Asp2195、Tyr2225、Val2227、Ile2237、Gly2238、Trp2239、Val2240、Met2345、Leu2354、Ile2356、Asp2357、Phe2358。5DXT 匹配 862 个 Cα，少于 4JPS 的 982 个，全域 RMSD 不是等覆盖比较。数字与换晶后不对称**方向一致，不作定量因果解释**。详见 Results 3.11。

---

## Table S11. Holdout 错口袋对照的几何对照（scoring-free contact_count）

来源：`data/jcim_holdout_v0/analysis/WRONG_POCKET_MECHANISM_VERDICT_V1.md`；脚本 `data/jcim_holdout_v0/scripts/wrong_pocket_contact_v1.py`。定义：配体重原子中与受体重原子距离 ≤4.0 Å 的原子数（`contact_count`），直接取自已冻结的 mode-1 姿态坐标，不涉及 Vina 能量函数。Vina 错口袋分臂取自 `holdout_pocket_matched_v1.csv`；重原子均值取自 `holdout_ligand_scores_v1.csv`。

| 靶对 | Vina 错口袋 summary_min（D/A，D/B） | contact_count AUROC（口袋 A / 口袋 B） | contact_count 的 min | dual / A_only / B_only 重原子均值 |
|------|------------------------------------:|----------------------------------------:|---------------------:|----------------------------------:|
| AChE/BChE | 0.643（0.643 / 0.653） | 0.581 / 0.706 | 0.581 | 35.1 / 34.0 / 29.5 |
| PIK3CA/mTOR | 0.788（0.788 / 0.858） | 0.552 / 0.698 | 0.552 | 33.5 / 32.3 / 31.0 |
| PIK3CA/PIK3CB | 0.520（0.640 / 0.520） | 0.622 / 0.714 | 0.622 | 34.5 / 31.6 / 28.3 |

B 臂 contact_count 高于随机（0.698–0.714），与 dual 对 B_only 尺寸差更大一致；A 臂接近随机（0.552–0.622），与 dual 对 A_only 尺寸差很小一致。contact_count **不能按幅度复现** Vina 错口袋（尤其 PM：0.788 对 0.552）。详见 Results 3.9。

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

## 写法说明（不进投稿 SI 正文）

- 本文件是**已有数据的汇编**，不是新实验。若某分析尚无机器可读表，宁缺毋填。
- 投稿英文 SI 时：Table 编号可按期刊习惯重排；数字不得改动。
- Cognate 表必须同时报告 mode1 与 best_of_9，避免审稿人误读“全部 &lt; 2 Å”。
- 早期借用 Schrodinger 处理过的姿态对照**不写入投稿稿**（无正式使用权限；主协议已统一为 RDKit/meeko）。仓库内 `pm48_directional_by_prep_v1.csv` 仅作内部记录。
- ChEMBL median / confidence≥8 / 物种过滤：本地缓存无字段（见 `T0_SKIPS.md`），不得编造；写入 Limitations。
- Table S12 是计数核对（BindingDB REST + PubChem PUG REST），不是对接结果；不得把 `as_is` 的 EGFR ≥50 写成已建成 BindingDB 厚面板。
