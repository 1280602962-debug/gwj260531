# SI 证据汇总清单（SI Evidence Register）

- 汇总日期：2026-10-17（Asia/Shanghai）
- 基底分支：`cursor/v4-2-final-audited-results-20260930`（HEAD = `7777e42e8f7f204d4ba275c7f1cbefc5f8f6bd5b`）
- 本分支：`cursor/si-evidence-consolidation-20261017`
- 性质：**只做加法**。本次只在新目录 `Dual_Target_Docking/si_evidence_consolidation_20261017/` 下新增文件，**未修改、未删除 v4-2 的任何已有文件**（尤其未动 `ANALYSIS_FREEZE_JCIM.md`、`Results_CN.md`、`Supporting_Information_CN.md`、`PRIMARY_DIRECTIONAL_METRICS.csv`、`results/canonical/` 等）。
- 复制原则：所有文件均用 `git archive` / `git show` 从各源分支按原始 blob 字节原样搬运，**未改动任何内容或数值**。
- 读不到的文件标注「未获取」。

---

## 一、证据文件总表

| 证据名称 | 原始分支 | 原始路径 | 本次复制后路径（相对仓库根） | 关键结论数值 | 是否可直接引用 | 备注 |
|---|---|---|---|---|---|---|
| 阶段 A 全套（A01–A10 + 快照 + supporting + 生成脚本） | `cursor/jcim-hardnegative-feasibility-20261012` | `Dual_Target_Docking/extension_jcim_validation/stage_A/` | `si_evidence_consolidation_20261017/stage_A/` | A09 判定 `STOP_NO_QUALIFYING_EXTERNAL_DATA`；A06：因无 E2/E3，描述符重叠与困难负类匹配可行性未运行 | 否（不在正式 v4-2 分支，属外部验证可行性的否定性结论） | 含 `supporting/build_stage_a.py`、`supporting/audit_summary.json`、`snapshots/pkis2_S4_pone.0181585.s004.xlsx` |
| FIVE_PAIR_CROSSDB_V1.md | `cursor/five-pair-local-channel-stats-0b1a` | `Dual_Target_Docking/data/jcim_chembl_universe_v0/local_track_b_v0/analysis/FIVE_PAIR_CROSSDB_V1.md` | `si_evidence_consolidation_20261017/bindingdb_feasibility/FIVE_PAIR_CROSSDB_V1.md` | 五对 equal_only 下 min HN 129/95/54/84/71，全部 >50；原文「Zero docking. Not Table 2. Do not hard-dock BindingDB as external validation.」 | 否（计数级一致性核对，非对接、非 Table 2、非外部验证） | 附支撑表见下 |
| BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md | `cursor/five-pair-local-channel-stats-0b1a` | `Dual_Target_Docking/data/jcim_novelty_v0/analysis/BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md` | `si_evidence_consolidation_20261017/bindingdb_feasibility/BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md` | 预过滤后 0 对合格，`packaged_as_external_evaluation = 0`，未对接 | 否（否定性可行性结论，不能作为外部验证） | 附 flow/funnel 等表 |
| BINDINGDB_EXTERNAL_SOP.md | `cursor/five-pair-local-channel-stats-0b1a` | `Dual_Target_Docking/docs/BINDINGDB_EXTERNAL_SOP.md` | `si_evidence_consolidation_20261017/bindingdb_feasibility/BINDINGDB_EXTERNAL_SOP.md` | BindingDB 外部切片的操作规程（SOP） | 作为方法 SOP 可引为附录 | 不在 v4-2 |
| BindingDB 支撑表（flow / funnel / external_candidate_flow / external_slice_summary / crossdb_strict_supply / fetch_log） | `cursor/five-pair-local-channel-stats-0b1a` | `Dual_Target_Docking/data/jcim_novelty_v0/tables/*` 与 `.../local_track_b_v0/tables/five_pair_crossdb_v1/*` | `si_evidence_consolidation_20261017/bindingdb_feasibility/tables/` | 漏斗各层计数（见上两份 md） | 随上述 md 一并作为 SI 支撑 | 原样搬运 |
| class_chemistry_summary.csv（8 对完整版） | `cursor/scientific-freeze-c7cc` | `Dual_Target_Docking/results/canonical/class_chemistry_summary.csv` | `si_evidence_consolidation_20261017/chemistry_neighbor/class_chemistry_summary.csv` | 32 行（8对×4类）；`nearest_dual_ecfp4` 中位数区间 0.1543–0.7386 | 作为描述量可引用（定义见 v4-2 的 CHEMISTRY_BASELINE_FREEZE.md） | **未放入 v4-2 的 canonical 目录**，仅放入本汇总目录，避免覆盖冻结结果 |
| source_document_concentration_v1.csv | `cursor/scientific-freeze-c7cc` | `Dual_Target_Docking/data/jcim_novelty_v0/tables/source_document_concentration_v1.csv` | `si_evidence_consolidation_20261017/source_sharing/source_document_concentration_v1.csv` | 文献来源集中度；比例对不上「25.5%–73.0%」 | 描述性审计可引用 | 原样搬运 |
| scaffold_cluster_bootstrap_v1.csv | `cursor/scientific-freeze-c7cc` | `Dual_Target_Docking/data/jcim_novelty_v0/tables/scaffold_cluster_bootstrap_v1.csv` | `si_evidence_consolidation_20261017/source_sharing/scaffold_cluster_bootstrap_v1.csv` | 骨架簇 bootstrap；比例对不上「25.5%–73.0%」 | 描述性审计可引用 | 原样搬运 |
| document_blocked_cv（folds / methods / summary） | `cursor/scientific-freeze-c7cc` | `Dual_Target_Docking/data/jcim_novelty_v0/tables/document_blocked_cv_*_v1.csv` | `si_evidence_consolidation_20261017/source_sharing/document_blocked_cv_*_v1.csv` | 文献连通簇阻断交叉验证 | 描述性审计可引用 | 原样搬运 |
| **H05_CHEMICAL_SPACE_OVERLAP.csv** | 未知 | （用户所述）`H05_CHEMICAL_SPACE_OVERLAP.csv` | **未获取** | — | 否 | 在 GitHub 上无法解析（见下「未获取项」） |
| **holdout_history_H0/H03_SOURCE_PROVENANCE.csv** | 未知（提交 `f4352f414`） | （用户所述）`extension_jcim_validation/holdout_history_H0/H03_SOURCE_PROVENANCE.csv` | **未获取** | — | 否 | 提交 `f4352f414` 无法解析（见下「未获取项」） |

---

## 二、已核实的关键事实（原样记录，未改数）

### 1. 权威指标表（在 v4-2 上，可直接引用）
- 路径：`Dual_Target_Docking/results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv`，共 120 行（8 对 × 5 方法 × 3 对比）。
- 列定义：`score_used = min(AUC_B, AUC_A)`、`auroc_tie = 0.5`、`universe = shared_dual`、`boot_B = 10000`、`boot_seed = 271828`、`ci_method = linear`。
- M0 八对 `summary_min` 与 95% CI（已逐条核实，与权威表一致）：

| 靶对 | summary_min | 95% CI | n_dual |
|---|---|---|---|
| EGFR/HER2 | 0.2467 | 0.1328–0.3783 | 28 |
| JAK1/JAK2 | 0.5908 | 0.4512–0.7246 | 32 |
| JAK1/TYK2 | 0.4219 | 0.2812–0.5635 | 32 |
| PIK3CA/mTOR | 0.5000 | 0.2778–0.6944 | 18 |
| AChE/BChE | 0.6341 | 0.4607–0.7540 | 27 |
| F2/F10 | 0.3623 | 0.2305–0.4805 | 32 |
| PPARG/PPARA | 0.6021 | 0.4561–0.7051 | 32 |
| PPARA/PPARD | 0.5107 | 0.3652–0.6475 | 32 |

- 聚合结论：精确中位数 **0.5054**；CI 下界最高者为 0.4607（AChE/BChE），**无一对的 95% CI 下界超过 0.5**。

### 2. 近邻描述量
- `class_chemistry_summary.csv` 实读 32 行（8对×4类），`nearest_dual_ecfp4` 中位数区间 **0.1543–0.7386**（含 q1/q3 时 0.1483–0.8634）；dual–dual 近邻 0.37–0.74；绝大多数中位数 <0.5，仅 PPARA/PPARD (0.7386) 与 PPARG/PPARA (0.566) 的 dual–dual 偏高。
- 最近邻定义出处：`analysis_plan_jcim/CHEMISTRY_BASELINE_FREEZE.md` D1 节，**在 v4-2 上可直接引用**；原文限定其为**描述量、非第三条化学基线 AUROC**；指纹参数 ECFP4, radius=2, nBits=2048, useChirality=false, useFeatures=false, binary。

### 3. 不可引用结论（重点记录）
- 曾流传的「25.5%–73.0%」跨折近邻 / 来源共享比例，**在现有任何表中均无可复现来源**：
  - 表中存在 ≈0.255 的点（JAK1/JAK2 A_only 中位数 0.2549），但**没有任何 0.730 的点**（最接近为 PPARA/PPARD dual 中位数 0.7386、PPARG/PPARA dual q3 0.7391，均为 73.9%，非 73.0%）；
  - `source_document_concentration_v1.csv`、`scaffold_cluster_bootstrap_v1.csv` 的比例也对不上。
- **结论：「25.5%–73.0%」应从所有稿件中删除，或改引上述实读区间。**

### 4. 骨架分组折（在 v4-2 上，可直接引用）
- `results/canonical/model_fold_assignments.csv`：929 条折记录、含 `scaffold` 列，按（靶对、方向、骨架）聚合为 767 组、**同组跨折 0**。
- 正确表述：这是 **Murcko 骨架分组交叉验证**，非普通冻结五折；其局限是「骨架分组 ≠ 来源独立的外部泛化验证」。

### 5. 阶段 A（外部验证可行性，否定性结论）
- A09 判定 `STOP_NO_QUALIFYING_EXTERNAL_DATA`：无来源达 E2，无一套数据同时满足三方向类别各≥20，未进入外部对接。
- A06 原文：因无 E2/E3，描述符重叠与困难负类匹配可行性未运行。

### 6. BindingDB 两份可行性文档
- `FIVE_PAIR_CROSSDB_V1.md`：五对（F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA、PPARA/PPARD）equal_only 下 min HN 分别 129/95/54/84/71，全部 >50；原文「Zero docking. Not Table 2. Do not hard-dock BindingDB as external validation.」
- `BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md`：五对（EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB、PIK3CA/mTOR、MCL1/Bcl-xL）预过滤后 0 对合格，`packaged_as_external_evaluation = 0`，未对接；各对 ECFP 后 n 与 dual/A/B/neither：216(180/10/20/6)、85(4/8/14/59)、112(9/0/3/100)、98(91/4/1/2)、3(1/0/2/0)。

### 7. 活性来源恢复（在 v4-2 上，可引用）
- `activity_label_source_audit_20261008/ChEMBL37来源恢复说明.md`。
- ChEMBL_37，SHA-256 `33c2037…`（2026-10-08 下载）；805 条合格观察、783 不重复分子、3555 原始活性行（A端 2076 / B端 1479，3506 个 activity_id）；未取整 MAX 与历史面板 805/805 一致，θ=6.0 回放 805/805；两端并列最大值均有文献标识 794/805；11 条为 DATASET 文档（PubChem BioAssay / DrugMatrix）缺 DOI/PMID/专利号；3 条历史审核差异（PM48_04/05/22）理由码仅 `source_adjudication`、实质理由待核实、冻结类别不变；`ACTIVITY_PROVENANCE_VERIFIED=PARTIAL`。

---

## 三、未获取项（待用户从本地推送后才能引用）

- `H05_CHEMICAL_SPACE_OVERLAP.csv` 与 `extension_jcim_validation/holdout_history_H0/H03_SOURCE_PROVENANCE.csv`（提交 `f4352f414`）：
  - 在 GitHub 上无法解析（contents API 与 get_commit 均返回 404/422「No commit found」）；各已取分支（v4-2、hardnegative、five-pair、scientific-freeze-c7cc 等）的树中也搜不到这两个文件名。
  - 推断：疑**仅存于本地未推送工作区**。**需用户将相应提交/文件推送到远端后方可获取并引用。**
