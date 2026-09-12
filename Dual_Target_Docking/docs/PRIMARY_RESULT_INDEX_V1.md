# 主结果来源索引

定位键：靶对—任务—评分通道—数据集—种子—标签规则。  
基线提交：`1ea95aa5`。本页只指向现行主分析，不把撤回体系或旧汇总当作当前数字。

| 稿件位置 | 任务 | 通道 / 数据集 / 种子 / 标签 | 分数或标签来源 | 统计脚本或锁定表 |
|---|---|---|---|---|
| Table 2 方向性 AUROC | dual–A-only、dual–B-only | Vina mode-1 / 主评价 / 20260727 / θ=6.0 | 原三对：`jcim_bench_v0` 分数长表；新增五对：`local_track_b_v0/tables/scores_vina_mode1_v1.csv` | `track_b_directional_auroc_v1.csv`；原三对锁定表见 `ARTICLE_ASSET_INDEX_v1.csv` |
| Table 2 `summary_min` | 两方向取较小值 | 同上 | 同上 | `track_b_summary_min_v1.csv` 及原三对锁定表 |
| Table 3 dual–neither | 同一评分通道换 neither | 同上 | 同上 | Table 3 锁定 CSV |
| Table S2a 对接盒 | — | 主受体 | `boxes/*.json` | SI Table S2a |
| Table S2b 共晶 RMSD | 搜索覆盖（坐标匹配） | 保存共晶姿态 | `cognate_qc/*_out_E8.pdbqt` 等 | `layer3_cognate_rmsd_v1.csv`；Figure S4 仍用此表 |
| Table S2c 共晶 RMSD | 化学对应复核 | 同一套保存姿态，不重对接 | 同上 | `layer3_cognate_rmsd_calcrrms_v1.csv` |
| 新增五对分数回放 | 生产 Vina | 本地 `poses/`（未入 git） | `scores_vina_mode1_pose_replay_v1.csv` | `replay_track_b_vina_mode1_v1.py` |
| Table S6 / S6b | 对应 vs 非对应口袋 | Vina 主评价与 holdout | 已锁定口袋分数 | `pocket_unidirectional_delta_v1.csv` |
| Table S9c / S9e | 多种子 | Vina / 主评价 / 五种子 / θ=6.0 | `multiseed/scores_vina_mode1_seed*.csv` | 完整病例：`fiveseed_summary_min_aggregate_v1.csv`；固定成员：`multiseed_fixed_membership_v1.csv` |

缺失项不得标 PASS。新增五对生产姿态目前只在本机磁盘核验过，公开仓库能重建的是分数长表上的统计，不是从 `REMARK VINA RESULT` 重新读出。
