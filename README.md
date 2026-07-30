# gwj260531

本仓库包含两个**相互独立**的课题目录：

| 目录 | 说明 |
|------|------|
| [`Dual_Target_Docking/`](Dual_Target_Docking/README.md) | **DualFourClass-Bench**：双靶对接四类硬负评测（冲 JCIM） |
| [`JNK1_Selectivity_Project/`](JNK1_Selectivity_Project/README.md) | JNK1 选择性筛选与建模（独立课题） |

---

## DualFourClass（当前主线）

入口：[`Dual_Target_Docking/README.md`](Dual_Target_Docking/README.md)  

下一步（P0）：[`Dual_Target_Docking/docs/JCIM_P0_COMPLETION_GUIDE.md`](Dual_Target_Docking/docs/JCIM_P0_COMPLETION_GUIDE.md)

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
python3 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py
```

## JNK1 选择性

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
# 见该目录 README
```
