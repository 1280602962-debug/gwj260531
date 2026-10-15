# 临床化合物库审计

已核实：`data/repurposing/repurposing_manifest.csv` 8319 行，与 `results/nlrp3_model_v2/stage2/clinical_library_scores.csv` 的 ID 集合相同。`docking_pool_p05.csv` 1588 行，是 8319 的子集，选择原因字段为 `nlrp3_p_active_ge_threshold`，池内最低历史概率大于 0.5。

来源脚本是 `scripts/build_repurposing_library.py`，读入 ChEMBL Explore Drugs 的分期表和 ATC 表。不是 Broad Drug Repurposing Hub。`data/repurposing/raw/` 只有说明，原始 Excel 不在仓库，导出日期无法恢复。

去重按规范 SMILES。每个化合物的 v1 概率、v2 的 `p_n392` 与 `p_n513`、历史排名和 v2 排名都在 `03_MACHINE_LEARNING/ALL_MODEL_SCORES_AND_RANKS.csv`，8319 行，没有删掉未过阈值的分子。

1588 之后：酸性等价 303，化学软过滤 156。v2 对历史池重数酸性等价是 304，差在名称排除的 Vecabrutinib 琥珀酸盐。303 是历史对接池的数字，304 是后来的重数，不是更正。

结构分布的已存汇总：酸性池 303 个有 273 个 Murcko 骨架；对接的 156 个有 148 个骨架；主候选 12 个有 12 个骨架。见 `data/manuscript/tables/chemical_space_summary.csv`。
