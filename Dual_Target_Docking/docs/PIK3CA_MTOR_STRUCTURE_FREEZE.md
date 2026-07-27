# PIK3CA / mTOR：对接结构冻结（2026-07-27 RCSB 复核）

> 回答：第二对靶到底用哪两个蛋白、哪两块晶体。

## 蛋白对（靶点）

| 端 | 蛋白 | UniProt | ChEMBL |
|----|------|---------|--------|
| A | **PIK3CA**（PI3Kα / p110α） | P42336 | CHEMBL4005 |
| B | **mTOR** | P42345 | CHEMBL2842 |

这是公开主对 DTPAIR-01，不是 PIK3Cγ、不是同工酶对 PIK3CA/PIK3CB。

## 对接主结构（冻结）

| 端 | PDB | 分辨率 | 配体 | chem_comp | 身份核对 |
|----|-----|--------|------|-----------|----------|
| PIK3Cα | **4L23** | 2.5 Å | **PI-103** | **X6K** | 真 α（P42336）+ ATP 位点抑制剂 |
| mTOR | **4JT6** | 3.6 Å | **PI-103** | **X6K** | 真 mTOR（P42345）+ 同一配体 |

**姿态金标准：** 与 EGFR/HER2 上 TAK-285（3POZ/3RCD）同逻辑 — 同一双靶配体 **PI-103** 双端共晶；cognate/self-dock 要求双端 RMSD < 2 Å。

## 明确不要用的 PDB

| PDB | 为什么否 |
|-----|----------|
| **7L1C** | HLA-A\*03:01 + 突变 **PIK3CA 肽**，不是 p110α ATP 口袋 |
| **9CMK** | p110α **RBD** + molecular glue，不是标准 ATP 双抑制剂口袋 |
| **4DRI** | FKBP51 + rapamycin + mTOR **FRB**，不是 ATP-competitive 激酶位点 |
| **3ML9** 等 | **PI3Kγ**（P48736）代用，不是 CHEMBL4005 的 α |

## 可选对照（不替代主结构）

- **8EXL**：真 α + taselisib（更新、分辨率更好，但非 PI-103 双端金标准）
- **4JT5**：mTOR + pp242（ATP 位点对照，非同配体双端）

## ChEMBL 规模（已审计）

paired 2713；dual 2002 / A_only 266 / B_only 236 / dual_weak 209。
