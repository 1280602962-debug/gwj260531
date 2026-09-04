# 生产对接结构

权威口袋与搜索盒：`config/docking_production_p2.yaml`、`docs/OPEN_SOURCE_DOCKING.md`。
C5 W1 交叉对接盒子：`docs/C5_DOCKING_WORKLIST.md`、`config/docking_c5_w1.yaml`。

## 生产受体（仅此两处）

| 靶点 | PDB | 说明 |
|------|-----|------|
| URAT1 | [9DKB](https://www.rcsb.org/structure/9DKB) | lesinurad，inward-open（Suo/Fedor 2025） |
| NLRP3 | [7ALV](https://www.rcsb.org/structure/7ALV) | NACHT + NP3-146（MCC950 类类似物，非 MCC950 共晶） |

准备：去配体、pH 7.4 加氢、转 PDBQT；盒中心见开源对接文档。生产引擎是 **gnina P2**，不是 Glide。

## C5 W1 交叉对接受体（不进临床库主表）

同构建体四联体。坐标文件已放入 `pdb/`。

| PDB | 文件 | 说明 |
|------|------|------|
| 9DK9 | `9DK9.cif` | apo，2.68 Å；盒中心由 9DKB 锁定盒经 CA 叠合转入 |
| 9DKA | `9DKA.cif` | 苯溴马隆（R75），3.00 Å |
| 9DKC | `9DKC.cif` | TD-3（A1A45），2.55 Å；必须 mmCIF，`.pdb` 下载 404 |

## 不要当作生产口袋

| PDB | 说明 |
|------|------|
| 9JDZ | 也是 lesinurad inward-open，**不是** occluded/outward |
| 9B1K / 9B1L | Dai 2024 的 occluded / outward；**本文不做三态对接** |
| 9DKA / 9JDY / 9JE1 / 8ETR | 对照或 redock 用，不进临床库主表 |

详见 `docs/URAT1_THREE_STATE_DOCKING.md`。
