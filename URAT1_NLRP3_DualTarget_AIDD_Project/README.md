# URAT1/NLRP3 双靶点 AI 辅助药物发现项目（STAD-AIDD）

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**STAD-AIDD**（**S**tructure-constrained **T**ransporter-aware **A**daptive **D**ual-target AIDD）是一套面向 **高尿酸血症/痛风** 的 **URAT1（代谢层）+ NLRP3（炎症层）** 双靶点小分子发现的 **纯计算** 方法学框架。

> **定位**：稳健型方法论文（偏 AI/生成模型），**不依赖湿实验**，通过严格回顾性验证、结构约束虚拟筛选与可复现基准支撑结论。

---

## 核心问题

| 挑战 | 本项目应对 |
|------|-----------|
| URAT1/NLRP3 公开活性数据仅千级 | 分子基础模型预训练 + 多任务迁移 + 辅助靶点数据 |
| URAT1 是转运体而非酶 | **构象系综对接** + 转运循环阻断评分（非静态口袋对接） |
| NLRP3 构象动态、结合位点不确定 | NACHT 域晶体/冷冻电镜结构 + 变构锁定机制约束 |
| 缺乏 URAT1/NLRP3 双靶系统研究 | 双任务学习 + 结构约束生成式优化 + 协同评分函数 |
| 无生物测试条件 | 文献 benchmark 回顾、骨架分组 CV、外部专利集验证 |

---

## 文档导航

| 文档 | 内容 |
|------|------|
| [**项目总体设计**](docs/PROJECT_DESIGN.md) | 科学逻辑、论文故事线、分阶段实施、目标期刊 |
| [**算法框架详解**](docs/ALGORITHM_FRAMEWORK.md) | 预训练、MTL、对接系综、生成式 RL、评分函数（**重点**） |
| [**URAT1 转运体验证**](docs/URAT1_TRANSPORTER_VALIDATION.md) | 转运体 vs 酶：必须验证什么、如何验证 |
| [**创新点与差异化**](docs/INNOVATION_POINTS.md) | 方法学创新、与现有工作对比 |
| [**准备清单**](docs/PREPARATION_CHECKLIST.md) | 数据、软件、结构、文献 benchmark |
| [**算力需求**](docs/COMPUTE_REQUIREMENTS.md) | GPU/CPU、存储、时间估算 |
| [**论文大纲**](docs/MANUSCRIPT_OUTLINE.md) | SCI 稿件结构、图表清单、写作要点 |
| [**参考文献**](docs/REFERENCES.md) | 可核验文献列表（DOI/PMID） |

---

## 目录结构

```
URAT1_NLRP3_DualTarget_AIDD_Project/
├── README.md
├── requirements.txt
├── LICENSE
├── config/
│   ├── targets.yaml              # ChEMBL 靶点 ID、数据清洗阈值
│   ├── docking_ensemble.yaml     # URAT1/NLRP3 对接系综与评分
│   └── model_hierarchy.yaml      # 模型层级与超参
├── data/
│   ├── structures/
│   │   ├── docking_ensemble_pdb.csv
│   │   └── README.md
│   └── benchmarks/
│       └── literature_benchmarks.csv
├── docs/                         # 全部设计文档
└── scripts/
    ├── 00_prepare_data.py        # ChEMBL + 专利数据清洗
    ├── 01_dataset_analysis.py    # 化学空间 / 活性分布 / 重叠分析
    ├── 02_train_mtl_models.py    # Chemprop MTL + XGBoost baseline
    ├── 03_structure_screening.py # 构象系综对接 + 转运体评分
    ├── 04_generative_optimization.py  # 双靶生成式优化（RL）
    ├── 05_retrospective_validation.py # 文献 benchmark 回顾验证
    └── run_stad_pipeline.py      # 端到端流程
```

---

## 快速开始（规划阶段）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
pip install -r requirements.txt

# 1. 数据准备（ChEMBL API + 专利活性）
python3 scripts/00_prepare_data.py

# 2. 数据集分析
python3 scripts/01_dataset_analysis.py

# 3. 多任务活性模型训练
python3 scripts/02_train_mtl_models.py

# 4. 结构约束虚拟筛选
python3 scripts/03_structure_screening.py

# 5. 生成式双靶优化（可选，算力需求较高）
python3 scripts/04_generative_optimization.py

# 6. 回顾性验证
python3 scripts/05_retrospective_validation.py
```

> 当前仓库以 **方法学设计文档 + 配置骨架** 为主；脚本为可执行框架，需按 `docs/PREPARATION_CHECKLIST.md` 准备数据后运行。

---

## 推荐论文题目（备选）

1. *STAD-AIDD: A structure-constrained transporter-aware framework for dual-target URAT1/NLRP3 inhibitor discovery under small-data regimes*
2. *Integrating molecular foundation models and conformational ensemble docking for synergistic URAT1 and NLRP3 inhibition in hyperuricemia*
3. *From metabolic urate clearance to inflammasome blockade: an AI-assisted dual-target design pipeline for gout*

---

## 许可

MIT License — 见 [LICENSE](LICENSE)。
