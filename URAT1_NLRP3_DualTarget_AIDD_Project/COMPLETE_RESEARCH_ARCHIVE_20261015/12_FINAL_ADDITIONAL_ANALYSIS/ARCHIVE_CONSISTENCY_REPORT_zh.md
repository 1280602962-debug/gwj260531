# 归档基准核对

本报告只记录核对结果。历史归档文件没有被修改。

## 与已存结果一致的计数

| 项目 | 核对值 | 来源 |
| --- | ---: | --- |
| 临床库 | 8319 | `ALL_COMPOUNDS_MASTER_INDEX.csv` |
| v1 冻结池 | 1588 | `in_nlrp3_pool_1588` |
| 酸性等价 | 303 | `in_acid_pool_303` |
| 化学软过滤 / 已对接 | 156 | `in_chemistry_soft_156`，并与 A2 对接成员一致 |
| 历史双靶结构合格，≥2/3 种子 | 54 | `acid_dual_keep_structural_seed{42,43,44}.csv` |
| 历史化学资格表 | 40 | `acid_a2_eligible_audited.csv` |
| 上述二者交集 | 37 | 与 `tier_pool_w2_annotated.csv` 的成员集合相同 |
| 冻结主候选 / 备用 / 可报告 | 12 / 21 / 33 | 冻结 CSV 与主索引标志一致 |
| 第三阶段 B 对接分子 | 32 | 主索引 |
| 第三阶段 B 对接运行 | 194 | `ALL_DOCKING_RUNS.csv`，`protocol_id = P_STAGE3B_NEW32` |
| 参照复现运行 | 7 | 同一文件，`P_REFERENCE7_REPRODUCTION` |
| 新增结构且化学通过 | 13 | `stage3b_new_structural_and_chemistry_pass` |
| 结构通过、化学未通过 | 8 | `stage3b_new_structural_pass_chemistry_fail` |
| 结构未通过 | 11 | `stage3b_new_not_structural_pass` |
| 规则不适用 | 7 | `stage3a_excluded_rule_not_applicable` |

`final_status` 十类人数之和为 8319。

## 37 为什么变成 33

40 个分子的化学资格来自 `scripts/build_c1_acid_shortlist_a2.py`，不是后来的 β-内酰胺规则。进入这张表需要同时满足：

- seed 42 的 `keep_dual_acid_geometry` 为真；
- 在可用种子上，同一宽松双靶几何标志至少 2/3 为真；
- RDKit PAINS 与 Brenk 均为阴性（`alert_clean`）；
- 名称不含当时列出的软排除片段。

β-内酰胺 SMARTS `[#6]1[#6][#7][#6]1=O` 不在这 40 人的入组条件里。对这 40 人用同一套 FilterCatalog 重算，PAINS 与 Brenk 和存档标志 40/40 一致，规则可以比较。

C5 再把这 40 人与“A2 结构门 ≥2/3 种子”相交，得到 37 人。这 37 人就是 tier1 的 13 人加 tier2 的 24 人。冻结时只做了一次排除：

- REP_07907 GSK-3008348 FREE BASE，tier1，理由是结构对照；
- REP_02721 CEFETRIZOLE、REP_04581 CEFAZEDONE、REP_04623 CEFOXAZOLE，tier2，理由是 β-内酰胺。

13 − 1 = 12，24 − 3 = 21。没有为了凑 33 而改成员。这四个人都在 37 人里面，而且 β-内酰胺没有在更早的资格表里先删过他们一次。

REP_04687 CEFAZAFLUR 同时命中名称软排除和 β-内酰胺 SMARTS，但它没有进入 40 人，因此也没有进入 37 人。它不是“先按 β-内酰胺删掉、冻结时再删一次”。

54 人里有 17 人不在 40 人表中。其中 REP_02366 AZD1981 和 REP_03238 SODELGLITAZAR 重算后 PAINS、Brenk、β-内酰胺和名称软排除都是阴性，结构门在 seed 43 和 44 通过，seed 42 的 NLRP3 没有通过。40 人表以 seed 42 几何通过为行入口，所以这两人当时进不去。本次不把他们补进冻结名单。

## 文字与文件不一致的一处

`FINAL_ARCHIVE_ACCEPTANCE.md` 写“四个生成脚本”。`scripts/complete_research_archive_20261015/` 里实际有 5 个 Python 文件。这是验收文字和目录不一致，原始验收文件按本次限制保持不动。
