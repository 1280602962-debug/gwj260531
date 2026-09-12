# 14 个主受体共晶 RMSD 统一化学对应复核

不重对接。原子对应优先级与 `reaudit_layer3_cognate_rmsd_v1.py` 相同。匈牙利坐标匹配只作对照。EGFR 3POZ 仍为重建 QC。

| 蛋白 | PDB | E | Hungarian top-1 / best | CalcRMS top-1 / top-3 / best | best&lt;2 | top-1&lt;2 | 方法 |
|------|-----|--:|------------------------:|-----------------------------:|:---:|:---:|------|
| AChE | 4EY7 | 8 | 0.339 / 0.339 | 0.339 / 0.339 / 0.339 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |
| BChE | 4BDS | 8 | 0.474 / 0.386 | 4.794 / 0.386 / 0.386 | 是 | 否 | SDF graph + element-constrained map + CalcRMS |
| EGFR | 3POZ | 8 | 7.780 / 0.680 | 9.505 / 6.227 / 0.760 | 是 | 否 | Meeko topology + RDKit CalcRMS; no superposition |
| HER2 | 3RCD | 8 | 1.729 / 1.113 | 1.855 / 1.394 / 1.394 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |
| PIK3CA | 4L23 | 16 | 0.624 / 0.624 | 0.624 / 0.624 / 0.624 | 是 | 是 | Meeko pose vs CCD/SMILES + bond-aware crystal coords + graph automorphism CalcRMS |
| mTOR | 4JT6 | 16 | 1.648 / 0.445 | 7.118 / 0.445 / 0.445 | 是 | 否 | Meeko pose vs CCD/SMILES + bond-aware crystal coords + graph automorphism CalcRMS |
| F2 | 4UDW | 8 | 0.382 / 0.382 | 0.382 / 0.382 / 0.382 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |
| F10 | 2JKH | 8 | 0.658 / 0.658 | 0.658 / 0.658 / 0.658 | 是 | 是 | CCD SMILES + bond-aware crystal coords + CalcRMS |
| JAK1 | 6N7A | 8 | 0.459 / 0.459 | 0.459 / 0.459 / 0.459 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |
| TYK2 | 3LXP | 8 | 0.196 / 0.196 | 0.196 / 0.196 / 0.196 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |
| JAK2 | 8BXH | 8 | 4.064 / 0.807 | 10.596 / 0.807 / 0.807 | 是 | 否 | Meeko topology + RDKit CalcRMS; no superposition |
| PPARG | 9V8H | 8 | 6.493 / 1.459 | 7.085 / 1.636 / 1.636 | 是 | 否 | Meeko topology + RDKit CalcRMS; no superposition |
| PPARA | 6LXA | 8 | 7.508 / 1.098 | 7.857 / 7.848 / 1.401 | 是 | 否 | Meeko topology + RDKit CalcRMS; no superposition |
| PPARD | 5U3Q | 8 | 1.452 / 1.452 | 1.510 / 1.510 / 1.510 | 是 | 是 | Meeko topology + RDKit CalcRMS; no superposition |

搜索覆盖（全部保存姿态最低 RMSD < 2.0 Å）：14/14。

表：`data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv`；逐姿态：`data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_modes_v1.csv`。
