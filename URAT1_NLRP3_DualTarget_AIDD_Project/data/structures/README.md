# 生产对接结构

权威口袋与搜索盒：`config/docking_production_p2.yaml`、`docs/OPEN_SOURCE_DOCKING.md`。

## 生产受体（仅此两处）

| 靶点 | PDB | 说明 |
|------|-----|------|
| URAT1 | [9DKB](https://www.rcsb.org/structure/9DKB) | lesinurad，inward-open（Fedor/Suo 2025） |
| NLRP3 | [7ALV](https://www.rcsb.org/structure/7ALV) | NACHT + NP3-146（MCC950 类类似物，非 MCC950 共晶） |

准备：去配体、pH 7.4 加氢、转 PDBQT；盒中心见开源对接文档。生产引擎是 **gnina P2**，不是 Glide。

## 不要当作生产口袋

| PDB | 说明 |
|-----|------|
| 9JDZ | 也是 lesinurad inward-open，**不是** occluded/outward |
| 9B1K / 9B1L | Dai 2024 的 occluded / outward；**本文不做三态对接** |
| 9DKA / 9JDY / 9JE1 / 8ETR | 对照或红ock 用，不进临床库主表 |

详见 `docs/URAT1_THREE_STATE_DOCKING.md`。
