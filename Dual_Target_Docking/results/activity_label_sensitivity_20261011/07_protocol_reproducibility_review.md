# 受体、配体和统计流程的可重复性

只读已有准备审计和正式配置。没有重新准备结构，没有重新对接。

## 已经有执行记录、可以写入方法的部分

M0 使用 Vina seed 42、mode 1。方向性分数是亲和力的相反数，较高值更有利。dual 对 A-only 使用 B 端分数，dual 对 B-only 使用 A 端分数。summary_min 是这两个 AUROC 的较小值，不是两个分数最小值的 AUROC。dual 必须两端 M0 都有效；A-only 需要 B 端有效，B-only 需要 A 端有效。缺失不填补。AChE/BChE 正式 M0 因此是 27 个有效 dual，而不是类别人数 28。

AUROC 是 Mann–Whitney，并列贡献 0.5。bootstrap 10000 次，种子 271828，每个靶点对重置 `numpy.random.default_rng`。每次重复按人群表原顺序依次抽取 dual、A-only、B-only，dual 抽样在两个方向共用。区间是 NumPy 线性插值的第 2.5 和 97.5 百分位数。点估计来自原始集合。S0 的人数、点估计和区间与 PRIMARY 的 M0 共 168 项一致。

一级受体是冻结的 14 个 PDB：3POZ、3RCD、6N7A、8BXH、3LXP、4L23、4JT6、4EY7、4BDS、4UDW、2JKH、9V8H、6LXA、5U3Q。4L23 去掉 p85，4JT6 去掉 mLST8，9V8H 保留肽，4UDW 保留 H/I/L，2JKH 保留 A/L。3LXP 的五个突变、4BDS 的 N→Q、2JKH 的 R372E 按沉积结构保留。6N7A、8BXH、3LXP 的磷酸化酪氨酸在 Meeko 前统一变成 TYR。这个决定没有使用对接分数或 AUROC。Meeko 对二硫键的覆盖限于晶体可核对的 CYX。6N7A A:1041、8BXH A:839、4UDW I:565 只保留晶体中观察到的重原子。

配体身份是去盐母体。AB_040 不单独生成三维结构。GNINA 审计确认 785 个冻结实体都有 SDF 和 PDBQT。29 个大环配体使用 CGn 真实原子类型；Gn 只是 glue 伪原子，不进入 SMILES 原子索引。

## 证据不足、不能写成已经核实的部分

V4.2 准备审计写明 785 个生产配体资产复用自 V4.1。当前仓库里 ETKDGv3、20 个构象、随机种子 42、MMFF94s、失败后改用 UFF，出现在共晶配体再对接脚本中。不能把这段参数写成 785 个生产配体在 V4.2 中重新生成构象的记录。Phase 0 冻结说明也写着生产用 ETKDGv3 当时尚未开始。

去盐、互变异构和手性标准化的逐步 RDKit 参数没有作为本次可执行协议留在正式目录。能复核的是冻结后的母体 InChIKey 和 alias 规则。

评价集候选池、配额和骨架抽样见 `06_panel_construction_review.md`，不能复现。

人源单蛋白靶点条目不等于每条测定都是人源野生型。EH40_27 的 A 端最大值来自 EGFR 突变体。PAPD_082 的 A 端最大值来自仓鼠 PPAR alpha，同端人源记录仍不低于 6。
