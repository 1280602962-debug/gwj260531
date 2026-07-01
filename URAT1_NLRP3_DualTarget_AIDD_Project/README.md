# URAT1 / NLRP3 痛风双节点 — 临床药物重定位计算项目

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

面向 **高尿酸血症/痛风** 的 **URAT1（代谢）+ NLRP3（炎症）** 双节点，在 ChEMBL **临床药物库** 上做 **NLRP3 ML 预筛 → 双靶对接 → Pareto 短名单**，并用 **8973 蒸馏集** 单独完成 URAT1 对接回顾验证。

> **当前论文路线（2026-07）**：见 [`docs/MANUSCRIPT_OUTLINE_CURRENT.md`](docs/MANUSCRIPT_OUTLINE_CURRENT.md)  
> **旧版 TAPE-GATE / MASFL / 8973 双靶 Pareto / OAT 迁移** 已归档，**不再按该路线执行**：[`docs/LEGACY_ARCHIVE.md`](docs/LEGACY_ARCHIVE.md)

**目标期刊**：*Journal of Molecular Modeling*（计算药理学 + 重定位，非湿实验 hit 发现）

---

## 科学定位（一句话）

> 在 ChEMBL 训练 **0 SMILES 重叠** 条件下，用 **NLRP3 ML** 筛临床药物库，对 **P(active)≥0.5** 命中分子做 **URAT1@9DKB + NLRP3@8ETR** 双靶对接并 Pareto 整合；用 **8973** 仅证明 URAT1 应对接而非 ML；用 **代表药 MD（2+2）** 解释机制。

**不是**：双靶新药发现、Teacher M-CPDL、百万库虚筛、OAT 迁移创新。

---

## 三套数据（禁止混用）

| 数据集 | 规模 | 用途 |
|--------|------|------|
| **`data/repurposing/repurposing_manifest.csv`** | 8319 | **主筛选**：NLRP3 ML → 对接 → Pareto |
| **`data/distill/distill_manifest.csv`** | 8973 | **仅 URAT1 回顾**：A vs D 富集（已 9DKB XP） |
| **Benchmark 六药** | 4 URAT1 + 2 NLRP3 | 对照定位 + MD |

---

## 当前计算流程

```
ChEMBL 临床药物库 (8319)
    → NLRP3 ML 全库打分                    [screen_repurposing_library.py]
    → P(active) ≥ 0.5  (n≈1588)           [对接池]
    → URAT1 @ 9DKB XP + NLRP3 @ 8ETR XP   [Maestro / 本地]
    → Pareto 双证据短名单                  [merge_docking_pareto.py]
    → 代表药 MD 2+2                       [benzbromarone, dotinurad, MCC950, GDC]

并行（独立一节，不用于 NLRP3 筛选）：
    8973 @ 9DKB XP → URAT1 ML vs 对接回顾  [merge_8973_docking_results.py]
```

详见 [**当前工作流**](docs/WORKFLOW_CURRENT.md)。

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
| 重定位库双靶对接 | ⏳ 本地 Maestro | 输入：`docking_pool_p05.csv` |
| Pareto 整合 | ✅ 脚本就绪 | `merge_docking_pareto.py`（待本地对接 CSV） |
| 代表药 MD | ⏳ | 2+2 benchmark |

事实数字见 [`docs/DATA_FACT_CHECK.md`](docs/DATA_FACT_CHECK.md)。

---

## 文档导航（按当前路线）

| 文档 | 内容 |
|------|------|
| [**当前工作流**](docs/WORKFLOW_CURRENT.md) | **主索引**：命令、路径、阶段 |
| [**论文定稿思路**](docs/MANUSCRIPT_OUTLINE_CURRENT.md) | Results 结构、主图、不写清单 |
| [**重定位库指南**](docs/REPURPOSING_DRUG_LIBRARY_GUIDE.md) | ChEMBL manifest、对接池 |
| [**8973 对接整理**](docs/LOCAL_AGENT_8973_DOCKING_PROMPT.md) | URAT1 回顾验证 only |
| [**模型质量报告**](docs/MODEL_QUALITY_REPORT.md) | URAT1_NO_GO / NLRP3 可用 |
| [**数据事实核验**](docs/DATA_FACT_CHECK.md) | 投稿前必读 |
| [**旧路线归档**](docs/LEGACY_ARCHIVE.md) | TAPE-GATE、MASFL 等（勿再执行） |

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

# 3) 8973 URAT1 回顾（可选，已完成则跳过）
python3 scripts/merge_8973_docking_results.py \
  --glide-csv results/docking/raw/9DKB_glide-dock_XP_8000_343e.csv
python3 scripts/analyze_urat1_docking_vs_ml.py
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
│   ├── repurposing/               # ChEMBL 临床药物 manifest
│   ├── distill/                   # 8973（仅 URAT1 回顾）
│   ├── docking/                   # 8973 合并分（已提交）
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
