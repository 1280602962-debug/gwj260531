# Track B 共晶 RMSD 化学对应复核

对象：已保存的 8 个 `*_out_E8.pdbqt`，不重对接。

主数字：RDKit `CalcRMS`（考虑对称、不做蛋白叠合）。匈牙利坐标匹配只作对照。`2JKH/BI7` 的 OpenBabel SDF 键级无效，改用 CCD SMILES。

| 蛋白 | PDB | Hungarian top-1 / best | CalcRMS top-1 / top-3 / best | best&lt;2 | top-1&lt;2 |
|------|-----|------------------------:|-----------------------------:|:---:|:---:|
| F2 | 4UDW | 0.382 / 0.382 | 0.382 / 0.382 / 0.382 | 是 | 是 |
| F10 | 2JKH | 0.658 / 0.658 | 0.658 / 0.658 / 0.658 | 是 | 是 |
| JAK1 | 6N7A | 0.459 / 0.459 | 0.459 / 0.459 / 0.459 | 是 | 是 |
| TYK2 | 3LXP | 0.196 / 0.196 | 0.196 / 0.196 / 0.196 | 是 | 是 |
| JAK2 | 8BXH | 4.064 / 0.807 | 10.596 / 0.807 / 0.807 | 是 | 否 |
| PPARG | 9V8H | 6.493 / 1.459 | 7.085 / 1.636 / 1.636 | 是 | 否 |
| PPARA | 6LXA | 7.508 / 1.098 | 7.857 / 7.848 / 1.401 | 是 | 否 |
| PPARD | 5U3Q | 1.452 / 1.452 | 1.510 / 1.510 / 1.510 | 是 | 是 |

搜索覆盖门槛（全部保存姿态最低 RMSD < 2.0 Å）在化学对应复核后仍全部通过。JAK2、PPARG、PPARA 的第一姿态未回到 2 Å 以内；PPARA 前三姿态也未回到 2 Å 以内。

表：`local_track_b_v0/tables/layer3_cognate_rmsd_calcrrms_v1.csv`；逐姿态：`local_track_b_v0/tables/layer3_cognate_rmsd_calcrrms_modes_v1.csv`。
