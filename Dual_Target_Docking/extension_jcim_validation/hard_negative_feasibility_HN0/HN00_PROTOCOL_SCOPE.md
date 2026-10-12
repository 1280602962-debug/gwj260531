# HN-0 困难负类可行性审计

任务：`HARD_NEGATIVE_FEASIBILITY_AUDIT_HN0`

本阶段只估计冻结八对主面板里，dual 与单靶选择性配体在固定化学条件下最多能形成多少不重复匹配。它不是正式困难负类性能评价，也不称为外部验证。

## 允许与禁止

允许计算描述符、ECFP4 Tanimoto、候选边和最大匹配容量。

不修改正式标签、PRIMARY、评分主表、别名或交叉验证折。不读取任何 `*_score_*` 列。不运行对接，不计算 AUROC、AUPRC、富集或 Top-K，不训练模型，不做 bootstrap 或显著性检验，不按已知表现调整阈值，不生成正式困难负类成员名单，不重新查询 ChEMBL。

## 输入

正式路径以 `FORMAL_AUTHORITY_PATHS.yaml` 为准。人群来自 `PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv` 的成员、类别和 `M0_A_valid`、`M0_B_valid`。结构只来自 `pair_ligand_mapping.csv` 的 `canonical_smiles`。折号只来自 `model_fold_assignments.csv`。别名登记排除 `AB_040`。

三层人群用 `jcim_stage2_lib.py` 的 `build_layer_members()` 和 `directional_m0_ok()` 构建，并与 `CHEMISTRY_POPULATION_LAYERS.csv` 对账。该函数本身不解析 SMILES。第三层全部配体另外用 RDKit 解析 `canonical_smiles`。解析失败则阻断对应靶对和方向，不替换结构，也不悄悄减人。

## 化学定义

ECFP4：半径 2，2048 位，不用手性，不用特征，二值 Morgan，Tanimoto。描述符八项从同一分子重算。FormalCharge 是原子形式电荷之和。另计 Bemis–Murcko 骨架。不做去盐、中和、互变异构或三维构象。

候选边同时满足 Tanimoto ≥ 0.50，以及 MolWt、MolLogP、TPSA 的标准化绝对差都 ≤ 0.50。均值和总体标准差（`ddof = 0`）只在第一层该方向的 dual 与相应单靶类合并集合上计算，然后原样用于第二层和第三层。标准差为 0 时记 `ZERO_VARIANCE`，标准化差定义为 0。不尝试其他阈值。

## 容量

单方向容量是最大二分图匹配边数，Hopcroft–Karp 与最大流互相核对，不一致则该单元为 `MATCHING_CAPACITY_UNVERIFIED`。共同三元组用 `SOURCE → A-only → dual 入点 → dual 出点 → B-only → SINK` 的最大流，每个分子容量为 1。它必须不大于两个单方向容量的较小值，也不能用那个较小值代替。候选边标记为 `CANDIDATE_EDGES_ONLY`，不冻结为正式样本。

## 判定

第三层上，至少两个靶对各有不少于 15 个共同三元组，才是 `READY_FOR_HARD_NEGATIVE_DESIGN`。否则，若至少一个靶对有不少于 15 个共同三元组，或至少两个靶对在任一方向各有不少于 15 个最大匹配，为 `LIMITED_DIRECTIONAL_FEASIBILITY`。都不满足则为 `STOP_INSUFFICIENT_MATCHES`。核心人数或身份不符为 `BLOCKED_DATA_INTEGRITY`。

15 是预先写明的供给门槛，不是统计功效证明。共同三元组不足，不单独否定已经达到 15 对的单方向供给。本数据里没有任何方向达到 15 对。PIK3CA/mTOR 的共同三元组理论上限是 B-only 的 12 人，规则不因此放宽。
