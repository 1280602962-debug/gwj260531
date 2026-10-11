# PPAR 药理作用可比较性

解决程度：PARTIAL。220 条观察的两端作用类型已经按既有语义表归组，并只在同一种作用定义下计算探索性 AUROC。原始文献中的作用方式、以及“缺少某种作用”是否代表无活性，都没有被重新写成标签。

## 证据

输入是 `activity_semantics_audit_20261009/PPAR_activity_semantics.csv` 的全部并列最大值，共 561 行、220 条观察、440 个端点。没有只用 `decisive_for_threshold=YES` 的 315 行。输出：

- `01_PPAR_pharmacology_review.csv`：561 行
- `01b_PPAR_observation_summary.csv`：220 条
- `01c_PPAR_same_type_counts.json`：同型人数与 AUROC 状态

两端作用组相同的观察有 161 条，两端不同的有 59 条。相同的定义来自语义表已有类别：激动或转录激活、拮抗或抑制功能、直接或竞争结合、部分激动或其他调节、作用方向不能确认。论文题名没有当作主要证据。FRET 仍记为作用方向不能确认。缺少某一作用类型的实验没有被记成无活性，拮抗也没有被自动记成低活性。

PAPD_015 的 A 端是明确拮抗，B 端是方向不能确认，因此整条观察不进入同型比较。PGPA_067、PGPA_068、PGPA_069 两端都只写有 FRET 活性，两端同为“方向不能确认”，没有负类，不计算 AUROC。

## 同型人数与 AUROC

AUROC 只在同一靶点对、同一种明确作用、并且两个方向的可评分人数都不少于 10 时计算。两个 PPAR 靶点对没有合并。分数使用已有 M0，没有重新对接。

PPARG/PPARA，两端都是激动或转录激活：dual 18、A_only 17、B_only 19、neither 5。dual 对 A_only 的点估计 0.5621，dual 对 B_only 的点估计 0.6067。正式 S0 的 `summary_min` 是 0.6020507812，较弱方向仍是 dual 对 B_only。这个子集人数更少，不能替换正式 110 条的结果。

PPARA/PPARD，两端都是激动或转录激活：dual 25、A_only 20、B_only 11、neither 8。dual 对 A_only 的点估计 0.7740，dual 对 B_only 的点估计 0.3236。正式 S0 的 dual 对 A_only 是 0.6787109375，dual 对 B_only 是 0.5107421875，`summary_min` 是 0.5107421875。子集里较弱方向的点估计更低，B_only 只有 11 人，没有自助区间。它只说明同型子集可以单独报告，不能据此改正式标签或正式 AUROC。

其余同型组合都没有同时形成不少于 10 人的正类和负类，只保留人数：

- PPARA/PPARD 方向不能确认：dual 1
- PPARA/PPARD 直接结合：dual 2、A_only 1、B_only 5、neither 4
- PPARA/PPARD 部分激动或其他调节：neither 1
- PPARG/PPARA 方向不能确认：dual 5
- PPARG/PPARA 直接结合：dual 3、A_only 8、B_only 2、neither 6

两端同为拮抗的观察数为 0，所以没有拮抗定义下的 AUROC。

## 对正式结果的影响

θ=6.0 的冻结类别没有改。正式 M0 与 S0–S4 没有重算。同型探索不能解释成“PPAR 的选择性已经在单一药理模式下得到正式证明”。PPARA/PPARD 激动子集的 dual 对 B_only 点估计低于正式较弱方向，这是新增的解释限制，需要在方法或讨论里说明子集范围。

## 尚未恢复

每条记录的原始论文实验方案、细胞与结合实验是否测量同一种受体构象，以及仓鼠 PPAR 与人源实验的逐条对应，都还没有从全文恢复。语义表里的“激动或转录激活”仍把激动和转录激活放在同一已有类别中，本次没有再拆开。
