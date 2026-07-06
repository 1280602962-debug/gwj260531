# JNK1 筛选项目通俗解读材料

面向**完全不懂计算机辅助药物设计（CADD）**读者的交付物。

## 文件说明

| 文件/文件夹 | 说明 |
|-------------|------|
| `JNK1筛选项目通俗解读报告.docx` | **主文档**：全流程通俗解读、术语表、方法原理、图表与核心数据表、参考文献 |
| `data_tables/` | **全部数据表**（25 个 CSV/JSON + README），与报告数值一一对应 |
| `figures/` | Word 文档与 ML 流程配图（见 `figures/README.md`） |

## 如何重新生成

```bash
cd JNK1_Selectivity_Project
# ML 建模 / 筛选 / decoy 配图（plot_style.py，300 dpi）
python3 scripts/plot_ml_pipeline_figures.py
# 通俗解读 fig01–06 + 刷新 data_tables
python3 scripts/build_popular_science_doc.py --figures-only
# 完整 Word 文档
python3 scripts/build_popular_science_doc.py
```

## 数据原则

- 所有数字来自仓库内可复现文件或技术报告 `docs/JNK1_PROJECT_REPORT.md`（v2.7）
- 参考文献仅使用技术报告 §11 已列条目及 benchmark CSV 中 `source` 字段记载的来源
- 未在文献中登记的软件（Schrödinger Desmond/QikProp/Prime 等）在文中标注为商业软件，不编造 DOI

## 与技术报告的关系

本通俗版是 `JNK1_PROJECT_REPORT.md` 的**非专业读者改写**，不引入新结论；若数值冲突，以技术报告与 `data_tables/` 原始文件为准。
