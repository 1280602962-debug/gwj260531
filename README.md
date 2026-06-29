# gwj260531

## 项目目录

| 项目 | 说明 |
|------|------|
| [**JNK1_Selectivity_Project/**](JNK1_Selectivity_Project/README.md) | JNK1/2/3 亚型选择性抑制剂 CADD/AIDD 完整流程 |
| [**URAT1_NLRP3_DualTarget_AIDD_Project/**](URAT1_NLRP3_DualTarget_AIDD_Project/README.md) | URAT1/NLRP3 双靶 AI 药物发现（**TAPE-GATE v2.0**：库筛+生成式双路径） |

### JNK1 快速开始

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py --library data/libraries/your_library.csv --output results/screening_v2
```

### URAT1/NLRP3 双靶项目（TAPE-GATE）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
# 框架总览与 PLK1/NLRP3 差异化
cat docs/TAPE_GATE_FRAMEWORK.md
cat docs/DIFFERENTIATION_VS_PLK1_NLRP3.md
# 端到端流水线
python3 scripts/run_tape_gate_pipeline.py --skip-generative  # Path A only
```
