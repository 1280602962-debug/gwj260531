# 对接验证

协议登记在 `DOCKING_PROTOCOL_REGISTRY.csv`。姿势级分数在 `ALL_POSE_LEVEL_RESULTS.csv`，覆盖酸性战役、第三阶段 B 和 7 次参照复现，共 10215 个姿势。phase-1 的 9839 个分子和 Rank 路线的 1583 个分子已有每分子一行的存档表，没有再拆姿势。

## 搜索和打分要分开

生产协议上的参照对接：

| 配体 | 受体 | 种子 | top-1 RMSD | 9 个姿势中的最好 RMSD | 搜索找到 ≤2 Å | 打分把它排在第一 |
| --- | --- | --- | --- | --- | --- | --- |
| lesinurad | 9DKB | 42 | 4.30 | 1.18 | 是 | 否 |
| lesinurad | 9DKB | 43 | 4.31 | 4.31 | 否 | 否 |
| lesinurad | 9DKB | 44 | 4.88 | 3.87 | 否 | 否 |
| benzbromarone | 9DKA | 42–44 | 约 3.59 | 约 1.10 | 是 | 否 |
| NP3-146 | 7ALV | 42 | 0.82 | 0.68 | 是 | 是 |
| NP3-146 | 7ALV | 43 | 0.67 | 0.67 | 是 | 是 |
| NP3-146 | 7ALV | 44 | 0.68 | 0.68 | 是 | 是 |

RMSD 是受体坐标系里的对称感知重原子 RMSD，没有再做配体叠合。lesinurad 种子 42 说明搜索可以产生接近晶体的姿势，CNNscore 没有选中它。种子 43 和 44 连搜索都没有产生 ≤2 Å 的姿势。NP3-146 三个种子都是打分和搜索一起通过。

更早的约 20 Å 盒子冒烟测试不是同一次实验。那次 Vina 的 top-1 RMSD 是 0.86 Å，GNINA CNNaffinity 的 top-1 是 4.16 Å、最好姿势 0.99 Å。不要把它和生产盒子的 4.30 Å 写成同一个数字。

## 7 次严格复现

`reference7/task_comparison.csv` 的 7 行都是 `strict_numerical_and_conclusion`。匹配姿势的 CNNscore 绝对差是 0，重原子 RMSD 是 0。这 7 次是 lesinurad 与 NP3-146 的种子 42，以及 REP_00015、REP_00300 的指定种子。种子 44 没有重跑。原始 SDF、命令和日志都在 `results/nlrp3_model_v2/stage3b_local/reference7/`。

在 A2 规则下，lesinurad 种子 42 选中的是模式 6，不是 CNNscore 最高的模式 1。这是规则差异，不是复现失败。

## 为什么有的评分没有用于最终排序

phase-1 真诱饵集 EF@1%：P5 2.80、P2 2.54、P0 1.90、P4 0.65、P1 0.43、P3 0.42。置信区间在 `RANKING_METHOD_ENRICHMENT_CI.csv`。P2 的区间约 1.31–4.07。P5 富集最高，但随机诱饵集 EF@1% 为 0，且覆盖不是完整的，没有进入生产。P1 和 P3 接近随机，被拒绝，表仍保留。对接评价没有计算 PR-AUC，只有表里的 AUC。

生产后来关闭的是“用这些分数找回 URAT1 活性”这件事，不是删掉了差的方法。最终候选不按 CNNaffinity 排序，而按几何门和化学标记。

## 缺失

Glide 分数已从仓库删除，没有幸存表。尿酸没有对接实验。phase-1 的 Vina 版本只有配置里的 1.2.5 建议，没有运行时版本字符串。酸性战役的临床 SDF 没有逐个保存 GNINA 版本日志；参照计算和第三阶段 B 的日志是 v1.3.2 `master:f23dd2b`。旧文稿写过 1.3.1，已存日志里没有看到这个版本。
