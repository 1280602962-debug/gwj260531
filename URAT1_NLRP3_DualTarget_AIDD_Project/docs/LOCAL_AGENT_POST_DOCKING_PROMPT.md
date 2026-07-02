# 本地 Agent：XP 对接打分后的数据处理

**适用**：你已在 Maestro/Glide 完成 **XP 对接**，需要合并分数、跑 Pareto、更新图表与 Git 可追踪数据。

**先分清两条线（不要混用）**：

| 路线 | 化合物数 | 靶点 | 用途 | 主脚本 |
|------|----------|------|------|--------|
| **A. 重定位主流程** | ~1588（P≥0.5） | 9DKB + 8ETR 双靶 | 主文 Fig4 Pareto、短名单 | `merge_docking_pareto.py` |
| **B. 8973 回顾** | 8973 | 仅 9DKB | 主文 Fig3 URAT1 回顾 | `merge_8973_docking_results.py` |

---

## 一、人工准备（对接刚算完时做）

### 1. 从 Maestro 导出 CSV

每个靶点导出 **XP 分数表**（可含多 pose；脚本会按 SMILES 保留最佳 pose）。

**必须能映射到列**（列名可不同，脚本有别名表）：

- SMILES：`canonical_smiles` / `ligprep_smiles` / `r_m_chemaxon_smiles` …
- XP 分：`r_glide_xp_gscore` / `glide_score_xp` / `r_i_glide xp` …
- 可选状态：`docking_status` / `r_i_glide_pose` …

### 2. 放到建议目录

**路线 A（1588 双靶）**：

```text
URAT1_NLRP3_DualTarget_AIDD_Project/results/repurposing/docking_raw/
  urat1_9dkb_p05.csv      # 1588 @ 9DKB XP
  nlrp3_8etr_p05.csv      # 1588 @ 8ETR XP
```

**路线 B（8973 回顾，若尚未合并）**：

```text
URAT1_NLRP3_DualTarget_AIDD_Project/results/docking/raw/
  9DKB_glide-dock_XP_*.csv
```

### 3. 确认对接池 manifest 存在

路线 A 需要（仓库已有 Git 副本）：

- `data/repurposing/screening/docking_pool_p05.csv`（1588 条）
- `data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv`（8319 条 ML 分）

---

## 二、直接复制给本地 Agent 的整段命令（路线 A：主流程 Pareto）

```
你在本地离线/半离线环境工作。项目根目录：URAT1_NLRP3_DualTarget_AIDD_Project/

背景（必读 docs/WORKFLOW_CURRENT.md）：
- NLRP3：8319 ML 预筛 → P≥0.5 约 1588 → 双靶 XP 对接 → Pareto 短名单
- 8973 仅用于 URAT1 回顾，不参与 Pareto
- 不要复用旧 TAPE-GATE / 8973 双靶 Pareto 逻辑

【我已完成】
Maestro Glide XP 对接，原始导出在：
  results/repurposing/docking_raw/urat1_9dkb_p05.csv
  results/repurposing/docking_raw/nlrp3_8etr_p05.csv
（若路径/文件名不同，请先 ls 并改参数，不要猜列名）

【任务 1】检查 CSV 列名是否与脚本别名兼容
- 读 scripts/merge_docking_pareto.py 中 SMILES_ALIASES / SCORE_ALIASES
- 打印两个 CSV 的 columns 前 20 项；若无法识别 SMILES 或 XP 分列，写一个小脚本 normalize_glide_export.py 统一为：
  canonical_smiles, glide_score_xp, docking_status, pdb_id
- 规范化后再跑 merge

【任务 2】合并双靶 + Pareto
python3 scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_8etr_p05.csv \
  --pool data/repurposing/screening/docking_pool_p05.csv \
  --sn-mode both \
  --output-dir results/repurposing

【任务 3】验收 merge 输出
- results/repurposing/pareto_merged_scores.csv
- results/repurposing/pareto_shortlist.csv
- results/repurposing/pareto_summary.json
检查 pareto_summary.json：
  - n_merged_dual_dock 应接近 1588（允许略少：对接失败/缺 SMILES）
  - n_pareto_front > 0
  - known_controls 含 lesinurad / colchicine 等（若在池内）

【任务 4】同步可提交 Git 的 data/ 副本（若 .gitignore 忽略 results/）
mkdir -p data/repurposing/pareto
cp results/repurposing/pareto_*.csv data/repurposing/pareto/
cp results/repurposing/pareto_summary.json data/repurposing/pareto/

【任务 5】更新 NLRP3 漏斗统计（Fig02d「pending」→ 实数）
- 编辑或生成 nlrp3_screening_summary 中 dual docking 完成数
- 重跑：python3 scripts/plot_available_figures.py
- 确认 figures/generated/main/fig02_* 漏斗第三级不再是 pending

【任务 6】（可选）Fig4 Pareto 散点图
- 若 plot_available_figures.py 尚无 Pareto 面板，在 scripts/ 中新增 plot_pareto_shortlist.py：
  - 读 pareto_merged_scores.csv
  - x = s_u_percentile（URAT1 对接），y = s_n_percentile（NLRP3 ML∪dock）
  - 标注 Pareto 前沿、对照药（灰/粉）
  - 输出 figures/generated/main/fig04_pareto_dual_evidence.pdf|.png
- 风格遵循 figures/jmm_style.py（Arial 8 pt，无网格）

【禁止】
- 用 8973 对接分跑 merge_docking_pareto
- 在 Pareto 纵轴仅用 NLRP3 ML 却不说明（默认 --sn-mode both 取 ML 与 NLRP3 对接百分位较大值）

【若失败】
- SMILES 对不上：对 pool 与 dock CSV 各抽样 5 行打印 canonical_smiles（RDKit 规范化后）比对
- 合并条数过少：检查 Maestro 导出是否缺盐剥离/是否与 docking_pool_p05 同一 SMILES 规范
- 把报错、列名、n_merged 写入 results/repurposing/docking_merge_qc.txt
```

---

## 三、路线 B：8973 回顾（仅当你又跑了一批 8973 或需重合并）

```
项目根：URAT1_NLRP3_DualTarget_AIDD_Project/

【输入】results/docking/raw/<你的 9DKB XP 导出>.csv

【执行】
python3 scripts/merge_8973_docking_results.py \
  --glide-csv results/docking/raw/<文件名>.csv \
  --pdb 9DKB

python3 scripts/analyze_urat1_docking_vs_ml.py \
  --merged data/docking/8973_9DKB_with_manifest.csv

【输出】
data/docking/8973_9DKB_with_manifest.csv
data/docking/urat1_docking_vs_ml_summary.json
data/docking/urat1_benchmark_rankings_docking.csv

【重画图】
python3 scripts/plot_available_figures.py
```

---

## 四、Agent 可能需要「新写」脚本的情形

| 情况 | 建议 |
|------|------|
| Maestro 列名特殊 | 写 `scripts/normalize_glide_export.py` 转成标准四列 |
| 一个 CSV 混了 SP+XP | 过滤 `docking_stage==XP` 或只保留 `r_glide_xp_gscore` 非空行 |
| 双文件分批对接 | 先 `pd.concat` + `groupby(smiles).first()`（XP 分越低越好） |
| 需要 7ALV SI | 另存 `nlrp3_7alv_p05.csv`，复制 merge 逻辑或加 `--nlrp3-dock-7alv` 参数 |
| QC 报告 | 写 `scripts/qc_docking_coverage.py`：池内 1588 中有多少双靶都有分 |

**优先复用**：`merge_docking_pareto.py`、`merge_8973_docking_results.py` 已含列别名与按 SMILES 取最佳 pose，不要重写 Pareto 数学。

---

## 五、本地一键命令（列名已标准时）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project

# 主流程 Pareto
python3 scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_8etr_p05.csv \
  --pool data/repurposing/screening/docking_pool_p05.csv

# 看图
python3 scripts/plot_available_figures.py
```

---

## 六、提交 Git 时建议包含

```text
data/repurposing/pareto/pareto_merged_scores.csv
data/repurposing/pareto/pareto_shortlist.csv
data/repurposing/pareto/pareto_summary.json
figures/generated/main/fig04_*   # 若已生成
```

`results/` 通常 gitignore，以 `data/` 副本为准（与当前仓库惯例一致）。

---

*Workflow: `docs/WORKFLOW_CURRENT.md` · 脚本：`scripts/merge_docking_pareto.py`*
