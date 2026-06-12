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
    ├── 01_download_chembl_data.py      # ChEMBL 数据下载与清洗
    ├── 02_dataset_similarity.py        # JNK1/2/3 数据集相似性比较
    ├── 03_sar_analysis.py              # SAR / MMP / 药效团分析
    ├── 04_train_selectivity_model.py   # 多任务 + 选择性模型训练
    ├── 05_model_interpretation.py      # SHAP / 子结构归因
    └── 06_virtual_screening.py         # 百万分子库筛选
```

## 快速开始

```bash
# 1. 安装依赖（建议 conda 环境）
pip install -r requirements.txt

# 2. 下载 ChEMBL 数据
python scripts/01_download_chembl_data.py --output data/raw

# 3. 数据集相似性分析
python scripts/02_dataset_similarity.py --input data/processed --output results/similarity

# 4. SAR 分析
python scripts/03_sar_analysis.py --input data/processed --output results/sar

# 5. 训练选择性模型
python scripts/04_train_selectivity_model.py --input data/processed --output models/

# 6. SHAP 可解释性分析
python scripts/05_model_interpretation.py --model models/best_model.joblib --output results/shap

# 7. 百万分子库虚拟筛选
python scripts/06_virtual_screening.py \
    --model models/best_model.joblib \
    --library data/libraries/enamine_real.smi \
    --output results/screening
```

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
