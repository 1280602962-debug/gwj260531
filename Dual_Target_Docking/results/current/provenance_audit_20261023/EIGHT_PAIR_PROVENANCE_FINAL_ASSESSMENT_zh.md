# 八对靶点来源与选样总评估

当前正式人群有 869 行，其中 867 条合格观察：dual 248、A-only 256、B-only 252、neither 111。PIK3CA/mTOR 是 Q110 的 110 人，其余七对合计 757 人。这些人数在审计开始时已经与权威入口一致，审计没有为了凑人数修改正式文件。

867 条观察的唯一键是 pair 加当前 `panel_id`。来源由三部分拼成，没有按全局 InChIKey 删除跨靶点对的合法共有分子：

- 757 条其他七对：复用 2026-10-08 审计，`REUSED_20261008_AUDIT`。
- 9 条 Q110 与 PM48 的重叠：复用同一审计并核对当前 Q110 编号。
- 101 条 Q110 新增：本轮从 ChEMBL 37 提取，`NEW_Q110_20261023_AUDIT`。

合计 766 + 101 = 867。

## 来源覆盖

867 条观察两端都有 `STANDARD_OK` 原始活性。这表示可以追溯到 ChEMBL 37 的活性、assay 和文档记录，不是同一实验平台验证，也不是外部独立验证。

855 条满足指定文献标识的完整覆盖：两端每一个并列最大值都有 DOI、PubMed 或专利号之一。12 条缺少这种标识。其中 11 条是 2026-10-08 已经记录的 PubChem BioAssay 或 DrugMatrix 数据集文档，分布在 EGFR/HER2、JAK1/JAK2 和 JAK1/TYK2；新增的 1 条是 Q110_A21 的 mTOR 端，文档同样是 CHEMBL1201862。

81 条观察的两端最大值来自不同文档。123 条存在终点类型混合。两项至少占一项的是 168 条。35 条在同一端有跨越 θ=6.0 的重复测定。数值或四状态不一致的是 0 条。这些情况允许同时存在。最大值规则没有改成平均值或中位数。

## 选样

七对都是 `DOCUMENTED_NOT_REPLAYED`。Track B 的沉积 110 人与当前名单一致，但原脚本明确禁止用事后排序重跑替换面板，候选池缓存也不在库中。EGFR/HER2 的 109 人由沉积 110 人减去已记录的 EH120_059 缺失臂排除得到；完整重放依赖当时的 API，本轮未调用。AChE/BChE 做了一次只读重放，与沉积面板差 1 人；当前 98 人仍由沉积面板经过 AB_040 别名和 AB_087 排除闭合。

## 分级

没有 Critical 或 Major 的已确认错误。没有靶点身份错误，没有冻结类别错误，也没有必须因此更换评价集的成员差异。

Moderate：

- `METHOD_LIMITATION`：168 条观察存在跨文档或终点混合，35 条存在跨越阈值的重复测定。冻结标签仍按最大值和 θ=6.0 重现。
- `PROVENANCE_GAP`：Q110_A24 的 mTOR 审核排除只有 `source_adjudication`，实质理由未解决；类别仍是 A-only。
- `PROVENANCE_GAP`：七对选样程序不能独立完整重放。AChE 重放与沉积面板差 AB_056 对 CHEMBL173607。
- `METHOD_LIMITATION`：EGFR/HER2 沉积面板写明保留 panel40、TAK-285、lapatinib，以及 70 个 hard-negative 优先扩充名额。

Minor：

- `DOCUMENTATION_ISSUE`：AChE 脚本把 `gray_excluded=True` 写在全部输出行上，不能当成逐分子灰区排除。fold 文件仍有 AB_040 行；现行权威说明它不是独立观察。

审计结论是 PASS。12 项完整性检查通过。Track B 独立重抽和 EGFR 的 API 重放记为 NOT_VERIFIABLE，因为缺少可合法重跑的输入，没有写成 PASS。
