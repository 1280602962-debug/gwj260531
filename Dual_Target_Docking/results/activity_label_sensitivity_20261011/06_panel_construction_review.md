# 评价集构建核查

本文件只记录当前仓库能够复核的最终规则。不能从最终人数反推一套抽样程序。

## 可以复核的最终规则

合格观察是 `activity_eligible=1` 且类别属于 dual、A_only、B_only、neither。共 805 条，783 个不重复 `molecule_chembl_id`。人群表 807 行，mapping 808 行。

每个端点使用 ChEMBL 37 STANDARD_OK 记录的最大 pChEMBL。允许的类型是 IC50、Ki、Kd、EC50、Potency、IC50app、Ki app，且 pChEMBL 非空。没有追加 confidence、assay type 或物种过滤。θ=6.0，不低于 6 为有活性。两端都有活性为 dual，仅 A 为 A_only，仅 B 为 B_only，两端都低于 6 为 neither。恰好等于 6 算有活性。缺少任一端活性的记录没有被写成 neither。EH120_059 和 AB_087 的类别为空、`activity_eligible=0`。

AB_040 是 AB_046 的 alias，不单独占一条观察，也不单独准备和对接。身份键是去盐母体 InChIKey。冻结说明里的全局实体数是 785，pair 内独立母体单位是 807。

各对合格人数：JAK1/JAK2、JAK1/TYK2、F2/F10、PPARG/PPARA、PPARA/PPARD 各 110；EGFR/HER2 109；AChE/BChE 98；PIK3CA/mTOR 48。正式四类人数是 EGFR/HER2 28/37/32/12，五个 110 人的靶点对都是 32/32/32/14，AChE/BChE 28/27/28/15，PIK3CA/mTOR 18/14/12/4。顺序为 dual、A_only、B_only、neither。

严格 6.5/5.5 不是构建这 805 条的纳入规则。EGFR/HER2 有 51 条、PIK3CA/mTOR 有 18 条在严格规则下落入灰区，它们仍留在正式评价集里。灰区是未分类，不是 neither。

## 无法复现的环节

当前协议目录和正式脚本里没有候选池定义、配额算法、抽样随机种子，也没有可执行的 Bemis–Murcko 骨架上限。五个靶点对呈现 32/32/32/14，只能说明最终人数平衡，不能写成已经恢复了抽样过程。EGFR/HER2 少 1 条合格观察，是因为 EH120_059 不合格；这不能解释它的四类为什么不是 32/32/32/14，也不能解释 PIK3CA/mTOR 为什么是 48 条。

化学折文件含有 scaffold 字段，用于已冻结的折分配。补充材料写明本轮不建立骨架等价规则，也不用骨架共享率声称没有泄漏。这不是评价集抽样时的骨架删除记录。

因此，最终标签规则和最终成员可以核对，完整候选池和抽样细节不能复现。本次敏感性分析不重新抽样。
