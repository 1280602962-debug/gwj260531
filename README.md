# JNK1 Selectivity Inhibitor Screening (CADD/AIDD Pipeline)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

基于 **ChEMBL JNK1/2/3 数据** 的 **JNK1 亚型选择性抑制剂** 计算机辅助药物设计（CADD）与 AI 药物发现（AIDD）完整流程。

## 项目目标

1. 分别下载并清洗 JNK1、JNK2、JNK3 生物活性数据  
2. **比较三个数据集的结构相似性与活性分布差异**  
3. 逐步整合 **5 类结构–活性关系（SAR）学习方法**，构建选择性预测模型  
4. 对 **百万级化合物库** 进行虚拟筛选  
5. 使用 **SHAP** 及其他可解释性方法评估模型  

## 文档

| 文档 | 说明 |
|------|------|
| [**完整流程思路**](docs/JNK1_selectivity_screening_workflow.md) | 科学设计、分步操作、评分函数、参考文献 |
| [**模型对比报告**](results/model_comparison/MODEL_COMPARISON_REPORT.md) | Chemprop vs XGBoost 初步结果 |
| [**参考文献**](docs/REFERENCES.md) | 格式化文献列表 |
| [**配置说明**](config/targets.yaml) | ChEMBL 靶点 ID 与阈值 |

## 目录结构

```
.
├── README.md
├── requirements.txt
├── config/
│   └── targets.yaml
├── docs/
│   ├── JNK1_selectivity_screening_workflow.md
│   └── REFERENCES.md
└── scripts/
    ├── 00_prepare_user_data.py       # 解析 docs/JNK*.csv → processed 数据
    ├── 01_download_chembl_data.py    # ChEMBL API 下载（可选）
    ├── 02_dataset_similarity.py      # JNK1/2/3 数据集相似性比较
    ├── 03_sar_analysis.py            # SAR / MMP / 药效团分析
    ├── 04_train_selectivity_model.py # XGBoost 多任务 + 选择性模型
    ├── 04b_train_chemprop_mtl.py     # Chemprop 2.0 原生 MTL
    ├── 05_model_interpretation.py    # SHAP / 子结构归因
    ├── 06_virtual_screening.py       # Virtual screening funnel
    ├── 07_compare_models.py          # Chemprop vs XGBoost comparison
    ├── build_demo_library.py         # Build demo SMILES library
    ├── plot_style.py                 # Journal figure style (Arial, 300 dpi)
    └── run_selectivity_pipeline.py   # End-to-end 04→05→06 pipeline
```

## Selectivity pipeline (04 → 05 → 06)

```bash
# Full pipeline: data prep + train + SHAP + screening
python3 scripts/run_selectivity_pipeline.py

# Or skip data prep if processed CSVs already exist
python3 scripts/run_selectivity_pipeline.py --skip-data-prep
```

All pipeline figures use **Arial**, **English labels**, and **300 dpi** (see `scripts/plot_style.py`).
For Arial on Linux: `sudo apt-get install ttf-mscorefonts-installer`

## 快速开始

```bash
# 0. 数据已在 docs/JNK1.csv, JNK2.csv, JNK3.csv
pip install -r requirements.txt

# 1. 一键：数据预处理 + 相似性分析 + 模型对比
python3 scripts/07_compare_models.py

# 或分步运行：
python3 scripts/00_prepare_user_data.py
python3 scripts/02_dataset_similarity.py --input data/processed --output results/similarity
python3 scripts/07_compare_models.py --skip-prepare
```

**初步结论（v2 改进后）**：**XGBoost** 平均 holdout R² = **0.699**（JNK1=0.703, JNK3=0.775 均超过 0.7），详见 [`results/model_comparison/MODEL_COMPARISON_REPORT.md`](results/model_comparison/MODEL_COMPARISON_REPORT.md)。

## ChEMBL 靶点 ID

| 亚型 | ChEMBL ID | 基因 | UniProt |
|------|-----------|------|---------|
| JNK1 | CHEMBL2276 | MAPK8 | P45983 |
| JNK2 | CHEMBL4179 | MAPK9 | P45984 |
| JNK3 | CHEMBL2637 | MAPK10 | P53779 |

## 流程概览

```
ChEMBL JNK1/2/3 下载
        ↓
数据集相似性比较 (化学空间 / 骨架 / 配对化合物)
        ↓
5 类 SAR 方法逐步分析
  ① 经典 SAR/MMP  ② 单靶点 ML  ③ 多任务 MTL
  ④ 选择性专用模型  ⑤ SBDD+ML 联合
        ↓
选择性预测模型 + SHAP 解释
        ↓
百万分子库漏斗筛选 → Top 候选
        ↓
酶活验证 (JNK1/2/3 IC50 + SI)
```

## 引用

如使用本流程，请引用 ChEMBL、RDKit、Chemprop、SHAP 等核心工具文献，详见 [docs/REFERENCES.md](docs/REFERENCES.md)。

## License

MIT
