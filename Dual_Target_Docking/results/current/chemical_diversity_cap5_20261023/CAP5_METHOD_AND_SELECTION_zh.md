# CAP5 方法与选择

这是正式八对评价集冻结之后的事后化学多样性敏感性分析。它不是预先登记的抽样规则，也不是八对分子当初共同使用的选样算法。本轮没有统一 GroupKFold，没有修改对接参数，也没有新增对接。

骨架调用现行 `jcim_stage2_compute._murcko`：`MurckoScaffold.GetScaffoldForMol` 之后 `MolToSmiles`。没有改用其他 Murcko 函数。结构来自现行 mapping 的 `canonical_smiles`，不做脱盐、互变异构枚举、三维构象、generic scaffold 或 Tanimoto 过滤。空骨架且环数为 0 记为无环；无环分子不参与“非空骨架最多 5 个”的限制。上限只在同一 pair、同一实验类别内计算。同一骨架出现在 dual 和 A-only 时不合并。不回补名额，不因药名或分数保留分子。

超额组内的排序键是 UTF-8 编码的 SHA-256，字段顺序为 `BM5_SENSITIVITY_20261023_V1|pair|class|global_ligand_entity_id|panel_id`，分隔符固定为 `|`。哈希升序保留前 5 个，哈希相同再按 `panel_id`。哈希只用于超额骨架内部，统计输入保持正式人群原行顺序。

与历史 fold 骨架字符串的核对：1082 行字符串完全一致；4 行是表示差异，其中 JAK1/JAK2 `J1J2_034`、JAK1/TYK2 `J1TYK2_065` 的两个方向是芳香写法不同但连接性 InChIKey 相同，PPARG/PPARA `PGPA_029` 的历史 fold 用单词 `acyclic` 表示空骨架，现行 `_murcko` 返回空字符串且环数为 0。AChE/BChE `AB_040` 是别名 fold 行，不是合格观察。没有无法解释的连接性或分组冲突，因此允许冻结。全部例外留在 `SCAFFOLD_CONCORDANCE_AUDIT.csv`。

合格观察 867 条全部计算。保留 861，排除 6，无环 4。排除只发生在：

- JAK1/JAK2 A-only：同一非空骨架 8 个，排除 `J1J2_031`、`J1J2_028`、`J1J2_013`
- PPARA/PPARD A-only：同一非空骨架 6 个，排除 `PAPD_027`
- PPARA/PPARD B-only：两个非空骨架各 6 个，排除 `PAPD_038`、`PAPD_055`

111 个 neither 的最大同骨架重复是 3，没有超过 5，因此 neither 没有排除。其余 pair 和类别也没有超额骨架。6 这个数字是全量计算后的结果，不是事先设定的排除人数。

描述性阈值在任何敏感性统计之前写入 `CAP5_SELECTION_FREEZE.json`。AUROC 的零效应是 0.5，差值的零效应是 0，点变化带宽是 0.02。`point_change_class` 取 STABLE、ATTENUATED、AMPLIFIED、REVERSED、NEAR_NULL、NOT_ESTIMABLE。`evidence_status` 单独记录正式区间和 CAP5 区间是否包含零效应。带宽只用于描述，不是显著性或等效性界值。完整点估计、区间和样本量都保留。

PRIMARY 使用 `q110_stats.single_method` 和 `paired_delta`。Stage 2 只调用 `module_b`、`module_c`、`module_d`、`module_i`。CAP5 的 Module I 先对八对 CAP5 行重新计算四组 `paired_delta`，再把带有 `pair`、`delta`、`delta_summary_min` 的结果传入 `module_i`。没有传入未过滤的正式 PRIMARY 差值，也没有调用 `execute()`、`run_stage2_modules()` 或 `run_primary_metrics_v42.py`。`q110_stats.py` 自身的 import 列表还包含其他 Stage 2 名称，本次运行没有调用它们。

Module D 沿用冻结 fold。原本缺 fold 的合格成员退出化学 OOF 并记录编号：AChE/BChE 的 `AB_001`、`AB_053`、`AB_054`，EGFR/HER2 的 `EH40_31`，F2/F10 的 `F2F10_080`，JAK1/TYK2 的 `J1TYK2_092`，PPARG/PPARA 的 `PGPA_030`。训练折或最终 OOF 单类别才标 NOT_ESTIMABLE；测试折单类别只作警告。本次没有新的 fold 连接失败，也没有单类别训练折或单类别 OOF。
