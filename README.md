# gwj260531

本仓库包含两个**相互独立**的课题目录：

| 目录 | 说明 |
|------|------|
| [`JNK1_Selectivity_Project/`](JNK1_Selectivity_Project/README.md) | JNK1 选择性筛选与建模 |
| [`Dual_Target_Docking/`](Dual_Target_Docking/README.md) | 双靶对接方法调研、打分设计与共晶编目 |

---

## JNK1 选择性

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py --library data/libraries/your_library.csv --output results/screening_v2
```

详见 [JNK1_Selectivity_Project/README.md](JNK1_Selectivity_Project/README.md)。

## 双靶对接

文档与数据均在 `Dual_Target_Docking/`，从 [Dual_Target_Docking/README.md](Dual_Target_Docking/README.md) 进入。
