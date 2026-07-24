# URAT1 / NLRP3 痛风双节点 — 临床药物重定位计算项目

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

面向 **高尿酸血症/痛风** 的 **URAT1（代谢）+ NLRP3（炎症）** 双节点：先在 **TrueDecoy/RandomDecoy** 上为 URAT1（9DKB）选定开源对接协议 Π\*，再嵌入 **NLRP3 ML 缩库 → 双靶对接 → Pareto → 成药性审计** 的不对称临床库漏斗。

> **当前论文路线（V2 · 2026-07-21）**：[`docs/MANUSCRIPT_OUTLINE_V2.md`](docs/MANUSCRIPT_OUTLINE_V2.md)  
> **中文写作**：引言 [`docs/INTRO_DRAFT_CN.md`](docs/INTRO_DRAFT_CN.md) · Methods [`docs/METHODS_DRAFT_CN.md`](docs/METHODS_DRAFT_CN.md) · 大纲 [`docs/MANUSCRIPT_OUTLINE_V2_CN_DRAFT.md`](docs/MANUSCRIPT_OUTLINE_V2_CN_DRAFT.md)  
> **本机重对接烟雾表（P0–P5）**：[`docs/REDOCK_SMOKE_TEST_SHEET.md`](docs/REDOCK_SMOKE_TEST_SHEET.md) · `bash scripts/run_redock_smoke_local.sh`  
> **已过时（含 Glide XP 主叙事）**：[`docs/MANUSCRIPT_DRAFT_CN.md`](docs/MANUSCRIPT_DRAFT_CN.md) — 勿再当正文  
> **更旧归档**：[`docs/LEGACY_ARCHIVE.md`](docs/LEGACY_ARCHIVE.md)

**目标期刊（首投）**：*Journal of Computer-Aided Molecular Design*（Hybrid / 可选非 OA）；备选 *Molecular Diversity* 等。  
**不以** Schrödinger Glide XP 为默认对接答案（有许可仅可作 SI 对照）。

---

## 科学定位（一句话）

> 先按双诱饵框架锁定 URAT1 开源对接排序协议，再以 **NLRP3 ML** 缩临床库并对 **P(active)≥0.5** 池做 **9DKB + 7ALV** 双靶对接与 Pareto，经模块 A–F 审计后给出可检验假说（如 canagliflozin 类），**不声称**已验证双口袋抑制剂。

**不是**：双靶新药发现、默认 Glide XP 漏斗、Teacher M-CPDL、百万库虚筛、OAT 迁移创新。

---

## 三套数据（禁止混用）

| 数据集 | 规模 | 用途 |
|--------|------|------|
| **`data/repurposing/repurposing_manifest.csv`** | 8319 | **主筛选**：NLRP3 ML → 对接 → Pareto |
| **`data/distill/distill_manifest.csv`** | 8973 | **仅 URAT1 回顾**：A vs D 富集（Vina 重对接） |
| **Benchmark 六药** | 4 URAT1 + 2 NLRP3 | 对照定位 + MD |

---

## 当前计算流程

```
[协议筛选 · URAT1@9DKB]
  TrueDecoy + RandomDecoy 基准
    → P0–P5（Vina / gnina / RTMScore）富集选优 → 锁定 Π*

[主漏斗 · 临床库]
  ChEMBL 临床药物库 (8319)
    → NLRP3 ML 全库打分                    [screen_repurposing_library.py]
    → P(active) ≥ 0.5  (n≈1588)           [对接池]
    → URAT1 @ 9DKB + NLRP3 @ 7ALV（Π*）  [run_vina_batch / run_gnina_batch]
    → Pareto + 模块 A–F 审计提名           [merge_docking_pareto.py]
    → 代表药 MD（可选主文）

并行（可选 SI，不替代 TrueDecoy）：
    8973 A vs D 敏感性 / 历史 Glide 表对照
```

详见 [**当前工作流**](docs/WORKFLOW_CURRENT.md)（部分命令仍在更新中；论文叙事以 V2 为准）。

---

## 实现状态

| 模块 | 状态 | 脚本 / 输出 |
|------|------|-------------|
| 数据清洗 URAT1/NLRP3 | ✅ | `00_prepare_data.py` |
| NLRP3 + URAT1 模型训练 | ✅ | `02_train_asymmetric_models.py` |
| Benchmark 回测 | ✅ | `07_benchmark_backtest.py` |
| ChEMBL 重定位库 manifest | ✅ | `data/repurposing/repurposing_manifest.csv` |
| **NLRP3 ML 全库筛选** | ✅ | `screen_repurposing_library.py` |
| 8973 对接合并 + 回顾分析 | ✅ | `merge_8973_docking_results.py`, `analyze_urat1_docking_vs_ml.py` |
| 重对接烟雾（lesinurad@9DKB） | ✅ 表+脚本 | `docs/REDOCK_SMOKE_TEST_SHEET.md`；`run_redock_smoke_local.sh` |
| 重定位库双靶对接 | ⏳ Vina 批量 | `scripts/run_vina_batch.py`；输入：`docking_pool_p05.csv` |
| Pareto 整合 | ✅ 脚本就绪 | `merge_docking_pareto.py`（待本地对接 CSV） |
| 代表药 MD | ⏳ | 2+2 benchmark |

事实数字见 [`docs/DATA_FACT_CHECK.md`](docs/DATA_FACT_CHECK.md)。

---

## 文档导航（按当前路线）

| 文档 | 内容 |
|------|------|
| [**论文大纲 V2**](docs/MANUSCRIPT_OUTLINE_V2.md) | **主规划**：协议筛选 + 不对称漏斗 |
| [**中文 Methods**](docs/METHODS_DRAFT_CN.md) | 大纲、公式、§2 正文草稿 |
| [**中文引言**](docs/INTRO_DRAFT_CN.md) | 引言草稿（原创表述） |
| [**中文先行大纲**](docs/MANUSCRIPT_OUTLINE_V2_CN_DRAFT.md) | ✅/⏳ 可写章节 |
| [**当前工作流**](docs/WORKFLOW_CURRENT.md) | 命令、路径、阶段 |
| [**重定位库指南**](docs/REPURPOSING_DRUG_LIBRARY_GUIDE.md) | ChEMBL manifest、对接池 |
| [**模型质量报告**](docs/MODEL_QUALITY_REPORT.md) | URAT1_NO_GO / NLRP3 可用 |
| [**数据事实核验**](docs/DATA_FACT_CHECK.md) | 投稿前必读 |
| [**旧路线归档**](docs/LEGACY_ARCHIVE.md) | TAPE-GATE、Glide XP 旧稿等（勿再执行） |

---

## 快速开始

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
pip install -r requirements.txt

# 1) 训练模型（若尚无 results/training/*.joblib）
python3 scripts/00_prepare_data.py
python3 scripts/02_train_asymmetric_models.py --no-oat-transfer

# 2) NLRP3 ML 筛临床库 + 导出 P≥0.5 对接池
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto

# 3) URAT1 TrueDecoy / RandomDecoy 基准（协议筛选用；现行集已在 data/benchmarks/urat1_true_decoy/）
# 重建需 taosu 预过滤大池，勿再用默认 distill_subset_d 当作正式 VS 基准：
# python3 scripts/build_urat1_true_decoy.py --pool /path/to/taosu_pool_prefiltered.csv \
#   --ratio 10 --inactive-pactivity-max 5 --seed 42
# 对接请用 unique_docking_pool.csv（9849），不要对 true/random 两份 CSV 各跑一遍
# 随后在服务器按 docs/MANUSCRIPT_OUTLINE_V2.md 跑 P0–P5 + RTMScore，锁定 Π*
```

---

## 目录结构（当前相关）

```
URAT1_NLRP3_DualTarget_AIDD_Project/
├── README.md
├── config/
│   ├── targets.yaml
│   └── docking_ensemble.yaml      # 9DKB, 7ALV, 8ETR
├── data/
│   ├── repurposing/               # ChEMBL 临床药物 manifest + screening 输出
│   │   └── screening/             # docking_pool_p05.csv（1588，已提交 Git）
│   ├── models/                    # 训练模型副本（已提交 Git）
│   ├── benchmarks/backtest/       # 六药回测（已提交 Git）
│   ├── processed/                 # 训练集（已提交 Git）
│   ├── distill/                   # 8973（仅 URAT1 回顾）
│   ├── docking/                   # 8973 合并分 + 回顾分析
│   └── benchmarks/
├── docs/
│   ├── WORKFLOW_CURRENT.md        # ★ 主流程
│   ├── MANUSCRIPT_OUTLINE_CURRENT.md
│   └── LEGACY_ARCHIVE.md
└── scripts/
    ├── screen_repurposing_library.py
    ├── merge_docking_pareto.py
    ├── merge_8973_docking_results.py
    ├── analyze_urat1_docking_vs_ml.py
    ├── build_repurposing_library.py
    └── 02_train_asymmetric_models.py
```

---

## 推荐论文题目

*Clinical drug repurposing for gout-related URAT1 and NLRP3 targets: NLRP3 machine-learning prescreening, dual-target docking, and molecular dynamics of benchmark inhibitors*

---

## 许可

MIT License — 见 [LICENSE](LICENSE)。
