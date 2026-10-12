# 可行性判断

全研究状态：`STOP_INSUFFICIENT_MATCHES`

第三层化学共同人群上，没有任何靶对的单方向最大匹配达到 15，也没有任何共同三元组。因此既不是 `READY_FOR_HARD_NEGATIVE_DESIGN`，也不是 `LIMITED_DIRECTIONAL_FEASIBILITY`。人数门与 `CHEMISTRY_POPULATION_LAYERS.csv` 一致：938、934、928。第三层 928 条记录的 `canonical_smiles` 全部能解析。两个匹配算法的边数一致。没有读取评分列。

共同三元组不是单方向评价的前提。本轮没有任何方向达到 15 对，所以不存在“单方向已经够、只是三元组不够”的靶对。PIK3CA/mTOR 的三元组理论上限是 12，因为它只有 12 个 B-only；这一点如实保留，规则没有放宽。即使不考虑化学条件，它也不能进入 15 个三元组的设计门槛。它的单方向最大匹配也只有 1 和 1。

三层的最大匹配人数相同。M0 无效和缺少折号去掉的是别的分子，没有去掉仅有的候选边。下表的三个数字因此同时代表三层。`cap` 是不考虑化学条件时的三元组人数上限。

| 靶对 | dual–A-only | dual–B-only | 共同三元组 | 第三层三元组上限 | 第三层状态 |
| --- | ---: | ---: | ---: | ---: | --- |
| EGFR/HER2 | 2 | 0 | 0 | 28 | `STOP_INSUFFICIENT_MATCHES` |
| JAK1/JAK2 | 2 | 0 | 0 | 32 | `STOP_INSUFFICIENT_MATCHES` |
| JAK1/TYK2 | 1 | 1 | 0 | 31 | `STOP_INSUFFICIENT_MATCHES` |
| PIK3CA/mTOR | 1 | 1 | 0 | 12 | `STOP_INSUFFICIENT_MATCHES` |
| AChE/BChE | 0 | 0 | 0 | 25 | `STOP_INSUFFICIENT_MATCHES` |
| F2/F10 | 0 | 1 | 0 | 31 | `STOP_INSUFFICIENT_MATCHES` |
| PPARG/PPARA | 1 | 1 | 0 | 31 | `STOP_INSUFFICIENT_MATCHES` |
| PPARA/PPARD | 1 | 1 | 0 | 32 | `STOP_INSUFFICIENT_MATCHES` |

化学共同层相对标签层减少的是折号，不是匹配：JAK1/TYK2、F2/F10、PPARG/PPARA 的三元组上限各少 1。AChE/BChE 的 M0 层上限从 27 降到 25，因为正式 M0 有无效记录。这些减少没有改变最大匹配边数。

候选边所在的骨架很少。第三层每个方向的候选 dual 骨架数是 0、1 或 2。EGFR/HER2 对 A-only 的 2 条边来自 2 个 dual 骨架；其余有边的方向往往是 1 个骨架占全部候选边。即便人数门槛被满足，这样的集合也只能代表孤立结构对，不能代表该靶对的选择性化学系列。当前人数本身没有通过，所以这个集中性是附加限制，不是一条新的筛选规则。

## 各靶对的限制

EGFR/HER2 和 JAK1/JAK2 的 A-only 方向各能配 2 对，B-only 方向是 0。限制是 B-only 一侧没有同时满足指纹和描述符窗口的边，不是 dual 人数不足。

JAK1/TYK2、PPARG/PPARA、PPARA/PPARD 两个方向各有 1 对最大匹配，但没有任何一个 dual 能同时配上两侧，共同三元组是 0。

PIK3CA/mTOR 两个方向也各有 1 对。共同三元组为 0。就算化学条件全部放开，B-only 只有 12 人，三元组仍不能到 15。

AChE/BChE 两个方向的最大匹配都是 0。对 A-only 没有 Tanimoto ≥ 0.50 的对。对 B-only 有 6 对指纹过线，但没有一对同时落入三个描述符窗口。类别人数在化学条件之前仍有 25 个以上的三元组空位。

F2/F10 对 A-only 的最大匹配是 0，对 B-only 是 1。

## 是否值得做正式困难负类分析

以这一套预先固定的条件，不值得。最大的单方向供给是 2 对，共同三元组是 0。用 0 到 2 对去做方向性比较，结果会由个别分子对决定，不能实质增强八对基准的证据，只会把可分析样本缩到无法支撑结论的规模。

这不表示“困难负类”这个科学问题在任何别的、事先冻结的设计下都没有意义。本审计不允许根据本次人数去放宽 0.50 或取消描述符窗口。若要换条件，那是一次新的设计，不能从这次容量里挑选。

READY 没有成立。不进入对接，不计算 AUROC，不生成正式匹配成员。
