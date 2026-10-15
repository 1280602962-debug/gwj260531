# URAT1/NLRP3 双靶点研究完整归档

本目录是 2026-10-15 从提交 `f3fb3d5a95079335f06c117272be6060f07f0dbc` 建立的只读整理。它不包含新的模型训练、新的对接或新的分子动力学。原始结果仍在原路径。这里给出索引、派生汇总和中文说明。

从项目根目录重建派生表：

```bash
python3 scripts/complete_research_archive_20261015/build_file_inventory.py
python3 scripts/complete_research_archive_20261015/build_master_tables.py
python3 scripts/complete_research_archive_20261015/build_docking_tables.py
python3 scripts/complete_research_archive_20261015/build_gate_tables.py
python3 scripts/complete_research_archive_20261015/build_remaining_tables.py
```

先读 `11_FINAL_REPORT/COMPLETE_PROJECT_RESULTS_zh.md`。审稿质疑先读 `10_SCIENTIFIC_AUDIT/CRITICAL_MAJOR_ISSUES.md`。

## 必须分开的三件事

1. 历史当时的决定。冻结的 12 个主候选和 21 个备用候选以 2026-09-05 的 C5 规则为准，没有被后来的统一规则改写。
2. 2026-10-15 的回顾性评价。同一批已存姿势用 A2 和完整 NLRP3 结构门重读。468 个历史种子与当时的判定逐项一致。
3. 第三阶段 B 的新增计算。32 个羧酸分子的 192 次对接，加上 REP_07837 在 NLRP3 上种子 43 和 44 的补算。这 13 个结构与化学都通过的分子不是冻结名单的成员。

## 不能从 GitHub 独立获得的材料

临床库原始 Excel、ChEMBL 导出的确切版本、Maestro 与 LigPrep 的完整命令、历史配体准备时的软件版本，以及 PF-03882845 的分子动力学轨迹，都不在本仓库。分子动力学分析包在本机 `D:` 盘的 `viz_pack_20260907`，完整溶剂化轨迹的符号链接已经失效。详见 `08_ADMET_AND_MD/MD_SOURCE_INVENTORY.csv`。
