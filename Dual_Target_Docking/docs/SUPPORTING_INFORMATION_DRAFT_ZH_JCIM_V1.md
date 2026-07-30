# Supporting Information（中文工作稿 · JCIM Articles）

> 与 [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md)、[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。  
> **原则：** 仅收录仓库中已有实验/冻结产物；数字均可追溯到下列源文件。未做的分析不填空表。

**溯源（机器可读原料）：**

| SI 表 | 源文件 |
|-------|--------|
| Table S1 | `data/jcim_strengthen_t0t1_v0/ENV_PIN.md`；各 panel `protocol.yaml` |
| Table S2 | `data/*/boxes/*.json`（冻结面板所用条目） |
| Table S3 | PM：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT*.md`；AChE/BChE：`cognate_qc/COGNATE_QC.md`；PIK3CB：`cognate_qc/COGNATE_QC.md`；EGFR/HER2：`protocol/protocol.yaml` + `analysis/exhaustiveness_sensitivity_v1/SENSITIVITY_VERDICT.md` |
| Table S4 | `data/jcim_bench_v0/tables/threshold_sensitivity_v1.csv` |
| Table S5 | `data/jcim_strengthen_t0t1_v0/tables/matched_subset_directional_v1.csv` |

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
| ChEMBL | Web API；靶对审计锁定 2026-07-23（未记录 release 编号） |
| 配体准备 | 去盐（最大有机片段）→ RDKit AddHs → ETKDGv3（seed 20260727）→ MMFFOptimizeMolecule（maxIters=200）→ meeko 默认 PDBQT |
| 面板抽样种子 | 20260729 |
| Bootstrap | B = 2000；seed 20260729；配体层 2.5%–97.5% 百分位区间；不做多重比较校正 |
| Vina `n_modes` | 9 |
| Vina `energy_range` | 3 |
| Vina exhaustiveness | PIK3CA/mTOR 主面板/扩面/单靶对照 = 16；其余靶对 = 8；PM 另报 E = 8 对照 |
| GNINA | `--cnn_scoring rescore --minimize --seed 20260727` |
| 盒子定义 | 共晶配体 AABB + 5 Å；每边下限 20 Å |
| Cognate 通过门槛 | 重原子 RMSD，`best_of_9` &lt; 2.0 Å（同坐标系，不叠合） |
| 受体 PDBQT | PIK3CA/mTOR/EGFR/HER2：含氢蛋白坐标 + `mk_prepare_receptor.py --read_pdb`；AChE/BChE/PIK3CB：沉积 ATOM + `mk_prepare_receptor`（default altloc A） |

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

## Table S4. 标签阈值敏感性（vina_mean）

在既有面板配体上，将两端 pChEMBL 按 θ ∈ {5.5, 6.0, 6.5} 重标四类后，用**两端 Vina 分数均值**（`vina_mean`）计算方向 AUROC 与 summary_min。  
**注意：** 此表不是正文主指标“口袋匹配”分数；口袋匹配主结果见正文 Table 2。完整多通道表见 `threshold_sensitivity_v1.csv`。

| 靶对 | θ | n (D / A / B / N) | AUROC D vs A | AUROC D vs B | summary_min |
|------|--:|------------------:|-------------:|-------------:|------------:|
| EGFR/HER2 | 5.5 | 69 / 22 / 10 / 9 | 0.790 | 0.371 | 0.371 |
| EGFR/HER2 | 6.0 | 28 / 38 / 32 / 12 | 0.681 | 0.282 | 0.282 |
| EGFR/HER2 | 6.5 | 26 / 29 / 29 / 26 | 0.740 | 0.304 | 0.304 |
| AChE/BChE | 5.5 | 27 / 25 / 28 / 15 | 0.530 | 0.585 | 0.530 |
| AChE/BChE | 6.0 | 27 / 25 / 28 / 15 | 0.530 | 0.585 | 0.530 |
| AChE/BChE | 6.5 | 27 / 25 / 28 / 15 | 0.530 | 0.585 | 0.530 |
| PIK3CA/PIK3CB | 5.5 | 30 / 25 / 28 / 16 | 0.746 | 0.430 | 0.430 |
| PIK3CA/PIK3CB | 6.0 | 28 / 27 / 28 / 16 | 0.703 | 0.412 | 0.412 |
| PIK3CA/PIK3CB | 6.5 | 28 / 27 / 28 / 16 | 0.703 | 0.412 | 0.412 |
| PIK3CA/mTOR | 5.5 | 33 / 9 / 5 / 1 | 0.488 | 0.630 | 0.488† |
| PIK3CA/mTOR | 6.0 | 18 / 14 / 12 / 4 | 0.722 | 0.671 | 0.671 |
| PIK3CA/mTOR | 6.5 | 17 / 15 / 12 / 4 | 0.702 | 0.667 | 0.667 |

† θ = 5.5 时 B_only 仅 5，标记为 underpowered（源表 `underpowered=1`）。  
AChE/BChE 与 PIK3CA/PIK3CB 面板在建造时已按严格配额冻结，简单 θ 重标后类别计数变化很小或不变。

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

## 写法说明（不进投稿 SI 正文）

- 本文件是**已有数据的汇编**，不是新实验。若某分析尚无机器可读表，宁缺毋填。
- 投稿英文 SI 时：Table 编号可按期刊习惯重排；数字不得改动。
- Cognate 表必须同时报告 mode1 与 best_of_9，避免审稿人误读“全部 &lt; 2 Å”。
- 早期借用 Schrodinger 处理过的姿态对照**不写入投稿稿**（无正式使用权限；主协议已统一为 RDKit/meeko）。仓库内 `pm48_directional_by_prep_v1.csv` 仅作内部记录。
