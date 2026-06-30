# URAT1/NLRP3 双靶点 AI 辅助药物发现项目（TAPE-GATE）

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**TAPE-GATE**（**T**ransporter-**A**ware **P**aired-path **E**vidence fusion with **G**enerative **A**nd library screening for dual-**T**arget **E**valuation）是一套面向 **高尿酸血症/痛风** 的 **URAT1（代谢层）+ NLRP3（炎症层）** 双靶点小分子发现的 **纯计算** 方法学框架。

> **v2.0 升级**：相对 STAD-AIDD v1.0 与 PLK1/NLRP3 类不对称框架，引入 **库筛 + 生成式双路径**、**assay-conditioned NLRP3 建模**、**可靠性加权 Pareto 融合**，并刻意避开锚点相似性 + 固定 0.5/0.5 融合等雷同设计。

---

## 核心问题

| 挑战 | TAPE-GATE 应对 |
|------|---------------|
| URAT1/NLRP3 **0 重叠** SMILES | 独立双模型 + 证据融合（MTL 仅作消融） |
| URAT1 是转运体而非激酶 | **$S_{\text{trap}}$ 构象系综**（非 PLK1 式单结构对接） |
| NLRP3 **7.2% 跨 assay 活性离散（curated 39 assays）** | **Assay-conditioned 分类**（非锚点 ECFP 相似性） |
| 候选化学空间局限 | **Path A 库筛** + **Path B CLM+RL 生成** |
| 与 PLK1/NLRP3 方法撞车 | 差异化模块 + **PLK1-style 阴性对照**消融 |
| 无湿实验 | 分路径 benchmark 回收 + 7 组消融（消融脚本为骨架） |

---

## 实现状态（GitHub 可核验）

| 模块 | 状态 | 证据 |
|------|------|------|
| 数据清洗 `00_prepare_data.py` | ✅ 已实现 | 822/513/0 重叠 |
| 双模型训练 + conformal | ✅ 已实现 | `02_train_asymmetric_models.py` |
| Benchmark 回测 | ✅ 已实现 | `07_benchmark_backtest.py` → URAT1_NO_GO |
| $S_{\text{trap}}$ / 对接系综 | ☐ 设计+配置 | `03_structure_screening.py` 骨架 |
| Path A/B 库筛/生成 | ☐ 骨架 | `03_library_screening.py` 等 |
| OAT 迁移训练 | ☐ 配置已写 | 辅助 CSV 待 ChEMBL 导出 |
| MASFL v3.1 | ☐ 设计稿 | 见 `MASFL_V3_WORKFLOW.md` |

事实与 ID 黑名单见 [**数据事实核验**](docs/DATA_FACT_CHECK.md)。

## 文档导航

| 文档 | 内容 |
|------|------|
| [**完整流程与文件清单**](docs/COMPLETE_WORKFLOW_AND_FILES.md) | **端到端流程、数据库、文件树（主索引）** |
| [**TAPE-GATE 框架总览**](docs/TAPE_GATE_FRAMEWORK.md) | v2.0 架构、双路径、融合策略 |
| [**MASFL v3.1 完整流程**](docs/MASFL_V3_WORKFLOW.md) | v3.1 扩展路线（**设计稿**，多数脚本未实现） |
| [**算法框架详解**](docs/ALGORITHM_FRAMEWORK.md) | 公式、伪代码、各 Stage 技术细节 |
| [**与 PLK1/NLRP3 差异化**](docs/DIFFERENTIATION_VS_PLK1_NLRP3.md) | 模块对照、避雷同清单（**重要**） |
| [**项目总体设计**](docs/PROJECT_DESIGN.md) | 科学逻辑、实施计划、期刊策略 |
| [**创新点与差异化**](docs/INNOVATION_POINTS.md) | 论文 Contribution、对比表 |
| [**URAT1 转运体验证**](docs/URAT1_TRANSPORTER_VALIDATION.md) | 转运体 vs 酶验证要求 |
| [**论文大纲**](docs/MANUSCRIPT_OUTLINE.md) | SCI 稿件结构、图表清单 |
| [**准备清单**](docs/PREPARATION_CHECKLIST.md) | 数据、软件、结构 |
| [**模型质量报告**](docs/MODEL_QUALITY_REPORT.md) | CV 指标 + benchmark 回测结论（**已运行**） |
| [**Benchmark 选择标准**](docs/BENCHMARK_SELECTION_CRITERIA.md) | 化合物合理性、分层考试、文献来源 |
| [**数据事实核验**](docs/DATA_FACT_CHECK.md) | ChEMBL/PDB/PMID 与规模数字（**投稿前必读**） |
| [**SLC22 辅助库逻辑**](docs/SLC22_AUXILIARY_RATIONALE.md) | OAT 主迁移 vs OCT 脱靶；错误来源与可信边界 |
| [**参考文献**](docs/REFERENCES.md) | 可核验文献列表 |

---

## 双路径流水线概览

```
Path A (库筛)                    Path B (生成式)
Enamine ~10⁶                     CLM cross-fine-tune + RL
    │                                │
    ├─ URAT1 conformal 过滤           ├─ 双靶 ML 奖励
    ├─ NLRP3 assay-conditioned       ├─ S_trap + NLRP3 对接奖励
    └─ 构象系综对接                   └─ QED/SA/新颖性
              │                                │
              └──────────┬─────────────────────┘
                         ▼
              可靠性加权 + Pareto 融合
                         ▼
              回顾性 benchmark 验证
```

---

## 目录结构

```
URAT1_NLRP3_DualTarget_AIDD_Project/
├── README.md
├── config/
│   ├── targets.yaml
│   ├── docking_ensemble.yaml
│   ├── model_hierarchy.yaml      # 不对称双证据模型配置
│   └── dual_path.yaml            # 库筛 + 生成式双路径配置
├── data/
│   ├── structures/
│   └── benchmarks/
├── docs/
└── scripts/
    ├── 00_prepare_data.py
    ├── 01_dataset_analysis.py
    ├── 02_train_asymmetric_models.py   # URAT1 + NLRP3 独立模型
    ├── 03_library_screening.py         # Path A
    ├── 04_generative_optimization.py   # Path B
    ├── 05_fusion_and_ranking.py        # 可靠性 Pareto 融合
    ├── 06_retrospective_validation.py  # 含 PLK1-style baseline
    ├── 07_benchmark_backtest.py          # benchmark 回测
    ├── run_model_build_and_validate.py   # 数据+训练+回测一键脚本
    ├── run_tape_gate_pipeline.py
    └── run_stad_pipeline.py            # v1.0 兼容入口
```

---

## 快速开始

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
pip install -r requirements.txt

# 1) 建模 + 质量评估 + benchmark 回测（推荐先跑）
python3 scripts/run_model_build_and_validate.py

# 2) 端到端 TAPE-GATE 流水线
python3 scripts/run_tape_gate_pipeline.py

# 仅库筛路径（算力有限时）
python3 scripts/run_tape_gate_pipeline.py --skip-generative
```

---

## 推荐论文题目

1. *TAPE-GATE: Transporter-aware paired-path evidence fusion for URAT1/NLRP3 dual-target discovery under assay-heterogeneous conditions*
2. *Assay-conditioned and conformation-ensemble dual evidence for hyperuricemia dual-target screening with generative augmentation*

---

## 许可

MIT License — 见 [LICENSE](LICENSE)。
