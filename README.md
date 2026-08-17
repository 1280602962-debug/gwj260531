# gwj260531

## 项目目录

| 目录 | 说明 |
|------|------|
| [`JNK1_Selectivity_Project/`](JNK1_Selectivity_Project/) | JNK1 选择性虚拟筛选与机器学习流程 |
| [`NLRP3_JNK1_Dual_Target_Design/`](NLRP3_JNK1_Dual_Target_Design/) | NLRP3/JNK1 双靶点药物设计调研、考量与实验验证清单 |

### JNK1 选择性筛选（快速开始）

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py --library data/libraries/your_library.csv --output results/screening_v2
```

详见 [JNK1_Selectivity_Project/README.md](JNK1_Selectivity_Project/README.md)。

### NLRP3/JNK1 双靶点药物设计调研

详见 [NLRP3_JNK1_Dual_Target_Design/README.md](NLRP3_JNK1_Dual_Target_Design/README.md)。
