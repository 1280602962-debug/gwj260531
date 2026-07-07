# 数据表格文件夹

本文件夹收录 JNK1/2/3 计算筛选项目报告中引用的全部可复现数据表。
所有数值均可回溯至仓库内原始文件；合成汇总表（00、16–24 号）的字段与
`docs/JNK1_PROJECT_REPORT.md`（v2.8）正文一致。

## 文件索引

| 文件名 | 说明 | 原始路径 |
|--------|------|----------|
| 01_文献benchmark化合物.csv | 见文件名 | `data/benchmarks/literature_benchmarks.csv` |
| 02_ML模型性能comparison.json | 见文件名 | `results/model_comparison/comparison.json` |
| 03_ML外部decoy验证指标.json | 见文件名 | `results/ml_external_validation/ml_external_validation_metrics_9bd8.json` |
| 04_再对接验证redocking.csv | 见文件名 | `results/docking_validation/redocking_summary_7725.csv` |
| 05_benchmark对接差值.csv | 见文件名 | `results/docking_validation/benchmark_deltas_51c1.csv` |
| 06_benchmark方向混淆矩阵.csv | 见文件名 | `results/docking_validation/direction_confusion_27c3.csv` |
| 07_benchmark_MMGBSA标定.csv | 见文件名 | `results/docking_validation/benchmark_mmgbsa_calibration.csv` |
| 08_Gly87自检.csv | 见文件名 | `results/docking_validation/gly87_selfcheck_16be.csv` |
| 09_各亚型排序相关.csv | 见文件名 | `results/docking_validation/isoform_rank_correlations_299a.csv` |
| 10_选择性标签统计.csv | 见文件名 | `results/similarity/sel_class_counts.csv` |
| 11_采购清单purchase_after_md.csv | 见文件名 | `data/purchase/purchase_after_md.csv` |
| 12_2231延伸MD_RMSD分位数.csv | 见文件名 | `results/md_2231_200ns/tables/09_production_rmsd_percentiles.csv` |
| 13_2231延伸MD_MMGBSA分量.csv | 见文件名 | `results/md_2231_200ns/tables/14_mm_gbsa_components.csv` |
| 14_ML阈值校准threshold.json | 见文件名 | `results/calibration/threshold_recommendation.json` |
| 15_ML虚拟筛选demo_screening.json | 见文件名 | `results/screening_v2/screening_report.json` |
| 25_shortlist_25to16.csv | 25→16 组内排序与 MD/落选明细 | `data/shortlist/md_shortlist_final.csv`；G1/G2 落选 ID 部分待归档 |
| 26_对接后筛选漏斗.csv | 对接后主线筛选各阶段数量 | 报告 §3.5 |
| 27_MD16_选择性排序与报价.csv | 16 人 **MD 后 JNK1 偏好**排序 + HIT 报价对照 | `purchase_after_md.csv` + 报告 §6.4 hinge/RMSD；`delta_sel_dock` 为参考列 |

| 00_端到端漏斗汇总.csv | 摘要漏斗各阶段数量 | 报告§摘要 |
| 16–24 号 CSV | 报告正文汇总表 | 报告对应章节 |

生成脚本：`scripts/build_popular_science_doc.py`
