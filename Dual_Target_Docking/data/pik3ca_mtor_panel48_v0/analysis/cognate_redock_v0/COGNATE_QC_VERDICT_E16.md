# COGNATE QC VERDICT E16 — PM48_01 (PI-103 / X6K)

**Verdict: Go**

## Go 条件核对
1. 4L23 @ 20260727 best9=0.624 (<2? True)
2. 4JT6 @ 20260727 best9=0.445 (<2? True)
3. 4JT6 噪声 seeds{20260727,7,42} pass count = 3/3（需≥2）

## 全表

| target | seed | rmsd_mode1 | rmsd_best_of_9 | best_mode | mode1<2 | best9<2 |
|--------|------|------------|----------------|-----------|---------|---------|
| 4L23 | 20260727 | 0.624 | 0.624 | 1 | True | True |
| 4JT6 | 20260727 | 7.118 | 0.445 | 3 | False | True |
| 4JT6 | 7 | 7.123 | 0.314 | 3 | False | True |
| 4JT6 | 42 | 7.123 | 1.374 | 5 | False | True |

## 备注

- mode1 失败但 best_of_9 成功：全面板必须输出 9 modes，并计划 RTM best-of-9。

## 参数
- exhaustiveness: 16
- n_modes: 9
- RMSD: heavy-atom, meeko SMILES IDX map, template automorphism min CalcRMS, no superposition
- poses: `poses/cognate_E16/`（未覆盖 E=8）
- table: `analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16.csv`

**允许启动全面板 48×2 @ E=16 @ seed=20260727。**
