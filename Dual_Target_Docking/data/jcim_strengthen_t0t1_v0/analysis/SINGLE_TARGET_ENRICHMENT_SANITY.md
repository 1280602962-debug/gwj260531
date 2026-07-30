# SINGLE_TARGET_ENRICHMENT SANITY — PIK3CA (4L23) / mTOR (4JT6)

> 协议：actives pChEMBL≥6.5；property-matched decoys pChEMBL≤5.5；Vina E=16，seed 20260727。  
> 面板规模：每端 ≈50 actives + 150 decoys（并行 ChEMBL SMILES）。  
> 表：`tables/single_target_enrichment_v1.csv`

| Receptor | Target | n_active | n_decoy | n_docked | AUROC | EF1% | EF5% |
|----------|--------|----------|---------|----------|-------|------|------|
| 4L23 | PIK3CA | 49 | 149 | 198 | **0.603** | 2.04 | 1.22 |
| 4JT6 | MTOR | 50 | 147 | 197 | **0.629** | 2.00 | 3.20 |

## 结论

两端 AUROC 均略高于 0.5 但仍 **≈0.60–0.63**（未 ≪0.6，亦未达强 enrichment）。  
→ 单靶 docking 有弱辨别力，**不足以单独支撑强对接声明**；PM 双靶口袋匹配信号应与此对照并降调「靠对接赢」的措辞。
