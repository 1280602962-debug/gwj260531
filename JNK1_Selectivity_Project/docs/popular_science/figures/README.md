# ML / 筛选流程配图索引

样式遵循 `scripts/plot_style.py`（Arial、300 dpi、英文轴标签）；通俗解读 fig01–06 为中文标签（150 dpi）。

## 一键生成（推荐）

```bash
python3 scripts/plot_ml_pipeline_figures.py
python3 scripts/build_popular_science_doc.py --figures-only
```

## 文件对照

| 图片 | 对应段落主题 | 数据来源 |
|------|-------------|----------|
| `model_comparison_r2.png` | XGBoost vs Chemprop | `results/model_comparison/comparison.json` |
| `screening_funnel.png` | ML 虚拟筛选漏斗 | `results/screening_v2/screening_report.json` |
| `score_distribution.png` | `final_score` 加权排序分布 | `results/screening_v2/all_hits.csv` |
| `decoy_validation_metrics.png` | 外部 decoy 回测 | `results/ml_external_validation/ml_external_validation_metrics_9bd8.json` |
| `ml_benchmark_isoform_prediction.png` | ML 无法判断亚型方向 | `results/screening_v2/benchmark_validation.csv` |
| `selectivity_label_scarcity.png` | 选择性标签稀缺 | `results/similarity/sel_class_counts.csv` |
| `fig01_funnel.png` | 端到端漏斗 | `data_tables/00_端到端漏斗汇总.csv` |
| `fig03_benchmark_direction.png` | 对接选择性方向 | `data_tables/06_benchmark方向混淆矩阵.csv` |
| `fig06_direction_accuracy.png` | 方向准确率不达标 | `data_tables/19_benchmark定量结果.csv` |

原始高分辨率副本亦在 `results/model_comparison/` 与 `results/screening_v2/`（含 `.pdf`）。
