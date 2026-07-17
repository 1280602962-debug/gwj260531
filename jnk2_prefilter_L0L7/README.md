# JNK2 对接前预筛选 L0–L7 结果（ARS×PaperSpine）

按本分支 `paper_rewriting_output/ars_deep_research/phase4_plan/prefilter_replan_v1.md`
对合并后的 **527,779** 个唯一伯/仲胺（ChEMBL ∪ 陶术）跑通对接前 chemotype triage。

## 漏斗计数

| 阶段 | 结果 |
|------|------|
| L0 胺库体检 | 527,779 |
| L1 胺预筛（MW 120–450；sites=1 优先；sites=2 且 MW≤350；≥3 丢） | 345,626 |
| L2–L4 标准丙烯酰胺枚举 + QC + 特征 | 239,990 产物（ok 80,174 / watch 159,816；bad 硬删） |
| L5 校准 | core↔full Tc ρ=0.66；Sim 窗 Jaccard=0.27；ErG p75=0.734 |
| L6 分仓 | Sim 8,063 / Novel 103,187 / Pan 9,204 |
| **L7 dock_ready** | **10,000**：sim_yl 4200 + sim_56d 1800 + novel 3500 + pan 500（仓间 ID 交集=0） |

## 文件

- `L7/L7_dock_ready_{sim_yl,sim_56d,novel,pan}.csv` — 下游 Glide 松阈值粗筛输入（分仓分文件，禁止合并成单一总分榜）
- `config/thresholds.json` — L5 校准后的工作阈值
- `config/anchors.json` — 四锚点 full/core SMILES（去弹头）
- `L5_calibration.md` · `L5/tc_grid.json` · `L5/loo.json` — 校准、Tc 网格敏感性、留一锚点
- `L0_report.md` · `L1_stats.json` · `L2L4_stats.json` · `FUNNEL_SUMMARY.md` · `ACCEPTANCE.md`
- `scripts/run_prefilter_L0L7.py` — 可重启流水线（L0/L1/L2L4/L5/L6L7）

## 阈值如何确定（L5）

1. 全库算 `max_tc_core`（去弹头）与 `max_tc_full`，Spearman ρ≈0.66、Sim 窗 Jaccard≈0.27 → 证明 Track-Sim 应用 core-Tc。
2. Tc 网格（lo∈{0.15…0.25}, hi∈{0.50…0.70}）看 Sim 规模稳定性 → 采用 **0.22–0.55**。
3. ErG 分布分位数 → Novel 门槛取 **p75=0.734**（或 hinge_hits≥1）。
4. `max_tc_full>0.70` 判近重复剔除；watch 在 Novel ≤20%；L7 配额为 ARS 默认。



## Step A–C（L7b 收紧 + Glide/AF3 交接包）

在 Novel 50 面板启发式 keep 审阅后执行（**未跑对接**）：

| 步骤 | 内容 |
|------|------|
| Step A | `stepA_novel_review/`：50 分子分层抽查面板 + checklist + 自动 keep |
| Step B | 拆分 hard_bad / soft_watch；重建 `L7b/`（ok > soft_watch） |
| Step C | `handoff_glide_af3/`：`.smi`、`af3_manifest_*.csv`、`GLIDE_AF3_HANDOFF.md` |

**L7b 交付（合计 8543，chemotype triage，非活性命中）**

| 仓 | n | ok | soft_watch |
|----|---|----|------------|
| sim_yl | 4198 | 456 | 3742 |
| sim_56d | 1800 | 610 | 1190 |
| novel | 2245 | 2245 | 0（仅 Step A keep=yes） |
| pan | 300 | 11 | 289 |

- 仓间 ID 交集 = 0
- hard_bad 主因：epoxide / β-lactam / aldehyde（共 5）
- soft_watch 主因：`aromatic_acrylamide`
- novel `unsure` 另存 `novel_unsure_hold.csv`，未进对接包
- 脚本：`scripts/run_stepBC.py`；审计：`watch_reason_counts.json`

## 下游

L7b → Schrödinger 共价对接（松阈值粗筛，仅几何/相互作用 QC）→ AF3 mPAE 最终门控与采购/合成名单。Sim/Novel 分榜；56d 仓提高 AF3 配额。

> 注：大体积中间产物（L1/L3/L4 特征表约 167MB、L6 分仓表约 88MB）保留在本地
> `D:\CADD paper exercise\JNK2\chembl_amine_pipeline\prefilter_L0L7\`，未纳入仓库。
