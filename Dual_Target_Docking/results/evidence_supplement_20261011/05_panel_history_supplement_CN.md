# 尚未复现评价集的历史文件

本轮只核对已有 Git 对象和缓存，没有改种子，也没有重新设计抽样。

## AChE/BChE

上一轮用 `d824eb6b` 的 `build_strict_panels.py` 重放，集合与沉积表不同。脚本里的多样性上限是 `molecule_chembl_id` 前 8 位，种子是 `random.Random(20260729)`，类别顺序是 dual、A_only、B_only、neither，配额是 28/28/28/16。

`edf0bb2f8` 的说明是采用 no-ID-prefix 面板。同一次提交把旧表面归档到 `data/_legacy_archive/ache_bche_panel_pre_fix/`，README 写明该表面用 ChEMBL 编号前缀上限建成，不用于修复后的主分析。脚本文件的 blob 在 `c213c485f`、`edf0bb2f8` 和 `d824eb6b` 上相同，仓库里的脚本没有改成 no-prefix。沉积表面在 `a0ef2726` 与 `edf0bb2f8` 上是同一个 blob。

用 `a0ef2726` 的 `mols_ACHE.json` 和 `mols_BCHE.json`，保持种子、排序后的交集、配额和类别顺序：

- 保留前缀上限：与归档的 pre-fix 表面集合完全相同。
- 去掉前缀上限：与沉积 `panel_v0_strict.csv` 的成员和顺序完全相同。

两张表只差一个 A_only 位置，panel_id 都是 AB_056。前缀上限版是 CHEMBL173607（7.78/5.19）。沉积版是 CHEMBL3960861（6.94/4.30）。沉积表里这个编号前缀出现 6 次，所以它不可能由“每前缀最多 5 个”的脚本直接产生。

因此 AChE/BChE 的沉积成员和抽样顺序可以由这条历史规则恢复。恢复的是提交说明里的 no-prefix 抽样，不是把脚本原文又跑一遍。正式 98 条仍是这 100 条去掉后来不合格的两条，这一关系和上一轮相同。

## Track-B

`31aba6564` 和 `a0ef2726` 的 `jcim_chembl_universe_v0/tables/` 里只有五张已提交面板和 `track_b_panel_summary_v1.csv`。没有抽样前的有序候选名单。`extract_track_b_panels_v1.py` 写明：五张面板提交时交集还没有排序，后来补上的 `sorted` 不能重现已提交结果。按这个说明，重新从 ChEMBL 抽取也不能恢复当时的排列。此处停止，没有重跑。

## EGFR/HER2 与 PIK3CA/mTOR

`a0ef2726` 有 `egfr_her2_panel120_v0/tables/smiles_cache.json`，共 138 条，不是完整的有序候选池。PIK3CA/mTOR 仍只有 `mols_PIK3CA.json`、`mols_MTOR.json` 和沉积的 48 行。没有找到可以代替当时 API 遍历顺序的候选列表。没有为了凑齐成员而重放。
