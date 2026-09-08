# URAT1 / NLRP3 结构选择

本文生产对接只用：

| 靶点 | PDB | 构象 / 口袋 | 共晶参考 |
|------|-----|-------------|----------|
| URAT1 | **9DKB** | inward-open | lesinurad |
| NLRP3 | **7ALV** | NACHT | MCC950 类类似物 NP3-146 / RM5 |

**不要用 9JDZ 当 occluded/outward。** RCSB 标注 9JDZ 也是 lesinurad 结合的 inward-open。occluded / outward 对应 Dai 等的 **9B1K / 9B1L**，本文不作为生产口袋。

搜索盒见 [`OPEN_SOURCE_DOCKING.md`](OPEN_SOURCE_DOCKING.md) 与 `config/docking_production_p2.yaml`。
