# 本地 Agent：ChEMBL 三库清洗 / 合并 / 打标签

**适用**：你已从 ChEMBL 下载三个 Excel：

| 文件 | 含义 |
|------|------|
| `Phase1_2_3.xls` | 临床阶段库 |
| `Level1 ATC.xls` | Level 1 ATC 大类库 |
| `Level 2 ATC.xls` | Level 2 ATC 精细库 |

---

## 一、人工准备（1 分钟）

1. 在项目里建目录并把三个文件放进去：

```text
URAT1_NLRP3_DualTarget_AIDD_Project/data/repurposing/raw/
  Phase1_2_3.xls
  Level1 ATC.xls
  Level 2 ATC.xls
```

2. 安装 Excel 依赖（若未装）：

```bash
pip install openpyxl xlrd
```

---

## 二、直接复制给本地 Agent 的整段命令

```
你在离线环境工作，不要 git clone / push。

项目根目录：URAT1_NLRP3_DualTarget_AIDD_Project/

任务：清洗、合并、分类标签化我从 ChEMBL 下载的三个药物库，生成双靶重定位筛选用的 manifest。

【输入】
data/repurposing/raw/Phase1_2_3.xls
data/repurposing/raw/Level1 ATC.xls
data/repurposing/raw/Level 2 ATC.xls

【执行脚本】（项目已提供）
python3 scripts/build_repurposing_library.py \
  --input-dir data/repurposing/raw \
  --phase-file "Phase1_2_3.xls" \
  --atc-l1-file "Level1 ATC.xls" \
  --atc-l2-file "Level 2 ATC.xls" \
  --primary-mode atc_phase \
  --min-phase-primary 3

【脚本应完成】
1. 读三个 ChEMBL Excel，自动识别 SMILES / ChEMBL ID / Max Phase / ATC 列
2. 去盐（最大有机片段）+ RDKit 规范化 SMILES
3. MW 150–800；剔除肽/蛋白类
4. InChIKey 第一块去重合并
5. 打来源标签：source_phase / source_atc_l1 / source_atc_l2
6. 打库面板标签 library_panel：
   - primary_atc_phase = ATC(L1或L2) 且 max_phase>=3（主筛选用）
   - atc_only / phase_only / benchmark_forced
7. 从 data/benchmarks/literature_benchmarks.csv 强制并入 benchmark（防 ATC 缺码漏药）
8. 输出：
   - data/repurposing/repurposing_manifest.csv（全量合并表）
   - data/repurposing/repurposing_primary.csv（主筛选面板，默认 atc_phase）
   - data/repurposing/repurposing_build_summary.json（规模统计）

【验收标准】
- repurposing_build_summary.json 中 n_manifest_unique_inchikey 约在 800–3000
- primary_atc_phase 子集约 500–2500（若过大可把 primary-mode 改为 atc_only 再导出一份对比）
- manifest 必含列：canonical_smiles, library_panel, source_phase, source_atc_l1, source_atc_l2, include_screen
- lesinurad / allopurinol 等 benchmark 在 manifest 中 benchmark_ref 非空或 library_panel=benchmark_forced

【若脚本报错】
- 列名不匹配：打印三个文件的前 5 列名，补 COL_ALIASES 后重跑
- .xls 读失败：pip install xlrd；或另存为 .xlsx
- 若文件名不同：改 --phase-file / --atc-l1-file / --atc-l2-file 参数

【不要做】
- 不要把 manifest 与 data/distill/distill_manifest.csv（8973）合并
- 不要用 8973 做 NLRP3 筛选

【完成后汇报】
贴出 repurposing_build_summary.json 全文 + primary 面板各 library_panel 计数。
```

---

## 三、输出字段说明（给后续筛选用）

| 列 | 用途 |
|----|------|
| `canonical_smiles` | Glide 对接 + NLRP3 ML 输入 |
| `library_panel` | 论文 Methods 分层（ATC vs 全临床对照） |
| `source_atc_l1` / `source_atc_l2` | ATC 来源追溯 |
| `source_phase` | 是否来自 Phase 导出 |
| `max_phase` | 临床阶段（3=III 期，4=上市） |
| `benchmark_ref` | 已知对照药 ID（如 URAT1_POS_01） |
| `include_screen` | 是否纳入广义筛选 union |

**主筛选文件**：`repurposing_primary.csv`（默认 = ATC 富集 ∩ max_phase≥3）

---

## 四、后续筛选命令（库建好后）

完整流程见 [`WORKFLOW_CURRENT.md`](WORKFLOW_CURRENT.md)。

```bash
# Step 1: NLRP3 ML 全临床库 + 导出 P(active)≥0.5 对接池
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all \
  --export-p05-pool \
  --skip-tanimoto
# → results/repurposing/docking_pool_p05.csv（n≈1588）

# Step 2: Maestro 双靶对接（本地 Schrödinger）
# 输入：docking_pool_p05.csv 的 canonical_smiles
# URAT1 @ 9DKB XP；NLRP3 @ 8ETR XP

# Step 3: Pareto 整合
python3 scripts/merge_docking_pareto.py \
  --ml-scores results/repurposing/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_8etr_p05.csv \
  --pool results/repurposing/docking_pool_p05.csv
```

**主筛选库**：全临床 manifest（`clinical_all`，n=8319），不是仅 `repurposing_primary.csv`。  
**SI 敏感性**：`--panel phase_ge3`（III 期+上市子集）。

**不要**：把 manifest 与 8973 distill 合并；不要在 8973 上跑 NLRP3 ML。

---

## 五、可选：导出三个面板作审稿对照

```bash
# 全临床 Phase 库（对照）
python3 scripts/build_repurposing_library.py --primary-mode phase_only \
  --output-dir data/repurposing/panel_phase_only

# 仅 ATC 富集（不设 phase 门槛）
python3 scripts/build_repurposing_library.py --primary-mode atc_only \
  --output-dir data/repurposing/panel_atc_only
```

论文 Supplementary 可报告三面板 hit rate 差异，降低 ATC 偏倚质疑。
