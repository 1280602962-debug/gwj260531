# 结果

八对靶点按固定顺序编号：1 EGFR/HER2，2 JAK1/JAK2，3 JAK1/TYK2，4 PIK3CA/mTOR，5 AChE/BChE，6 F2/F10，7 PPARG/PPARA，8 PPARA/PPARD。方向任务中，D_vs_A 使用 B 端分数区分 dual 与 A_only，D_vs_B 使用 A 端分数区分 dual 与 B_only。summary_min 是这两个方向 AUROC 的较小值。

## R1. 八对靶点挑战集的标签构成

本研究评价的是实验标签构成的方向任务，不是天然化合物库的命中率。Table 1 给出八对的 A/B 靶点、冻结 PDB 和四类标签人数。合格且有类别的 pair 内观察共 867 条，其中 dual 248、A_only 256、B_only 252、neither 111（Figure 1a）。dual 可以进入两个方向，因此标签方向记录为 1004 条，neither 不进入这两个任务。867 不是 867 个互不共享的小分子，跨靶点对可以共享配体实体。

869 行人群表在这 867 条之外，还保留 EGFR/HER2 的 EH120_059 和 AChE/BChE 的 AB_087 两条 activity_eligible=0、类别为空的记录；它们不是 neither。当前 mapping 为 870 行，多出的一行是 AChE/BChE 的 AB_040，它是 AB_046 的 alias。阈值、合格性和 alias 规则见补充方法第 1 节。三层用途不同：1004 条标签方向记录用于组成描述，1000 条正式 S1 方向记录用于 PRIMARY，994 条共同折外记录用于化学模型和 chemistry-matched S1。每臂人数见 Table S4A，评分可用性见 Table S4B 和 Figure S2，缺失名单见 Supplementary Data。

867 条合格观察的两端活性均可追溯到 ChEMBL 37 原始记录，冻结最大 pChEMBL 和 θ=6.0 类别均与来源核对一致。855 条观察的两端并列最大值均有 DOI、PubMed 或专利号，12 条缺少指定文献标识。81 条存在跨文档聚合，123 条存在终点类型混合，两者的并集为 168 条；35 条在同一端存在跨越活性阈值的重复测定。来源覆盖汇总见 Figure 1b 和 Table S1B，完整来源对照及逐条活性记录见 Supplementary Data。

## R2. 双方向表现与评分工作流差异

共晶配体重对接的 RMSD 结果见 Table S3 和 Figure S1。

Figure 2a,b 显示，S1 的两个方向并不对称。16 个方向点估计中，4 个低于 0.5，其余均高于 0.5。EGFR/HER2 的 D_vs_A 为 0.6274，D_vs_B 为 0.2467，summary_min 因此为 0.2467。F2/F10 的两个方向 AUROC 点估计均低于 0.5，分别为 0.4082 和 0.3623，summary_min 为 0.3623。PIK3CA/mTOR 的 D_vs_A 和 D_vs_B 分别为 0.7363 和 0.5801，summary_min 为 0.5801（95% 区间 0.4346 至 0.7129）。这些是点估计与已有区间的报告，不是按点值宣布显著失败。JAK1/JAK2 的 summary_min 区间来自 summary_min 自身行，上界为 0.7246，不等于较弱臂 D_vs_A 的上界 0.7275。

Figure 3 使用 PRIMARY 已有的四个配对差值，而不是把各方法自身人群的 AUROC 相减。S3−S1 的 summary_min 差值在 F2/F10 为 0.2471，95% 区间 0.0654 至 0.3896；在 PPARG/PPARA 为 −0.2915，区间 −0.4380 至 −0.1104。同一对照可以在不同靶点对上改变符号。S4−S2 是同一保存姿势集合上的重排对照，并带有已计算的配对区间。它不同于后文仅作 supporting 的 S4−S1 点估计和。S5−S1 比较的是独立对接工作流与 S1，不是单纯搜索算法效应。五方法的自身结果见 Table S5A 和 Figure S3，配对差值及完整区间见 Table S5B 和 Figure S11。

## R3. 负类组成改变部分方向的判别难度

Figure 4 在同一 dual 集合上比较单靶负类和 neither。差值为 neither−single。16 个臂中 12 个点为正、4 个为负，14 个配对区间包含 0。EGFR/HER2 的 D_vs_B 为 0.5004，区间 0.3025 至 0.6819；JAK1/TYK2 的 D_vs_B 为 0.3839，区间 0.1953 至 0.5693。这两个方向上，neither 比单靶负类更容易与 dual 分开。其余方向保留正负和宽区间，不能据此说所有单靶负类都更难，也不能把含 0 的区间写成等效。

PIK3CA/mTOR 的 neither 在两个臂都是 14 人；两个方向的 neither−single 差值分别为 0.0315 和 0.1164，区间均包含 0。完整人数和差值区间见 Table S6A，各条件 AUROC 区间见 Supplementary Data。

## R4. 化学标签可预测性与有限的 S1 增量

Figure 5 的三个化学模型和 chemistry-matched S1 使用同一第 3 层成员，共 994 条方向记录。ECFP4 的折外 AUROC 为 0.6633 至 0.9678，16 个臂都高于同成员 matched S1 的点估计。加入 S1 后的增量是 ECFP4+S1−ECFP4：11 个为负、5 个为正，没有零点，范围 −0.0164 至 0.0363。最大正点在 F2/F10 的 D_vs_B。PIK3CA/mTOR 两方向的增量分别为 0.0215 和 −0.0029。

化学模型在训练折上使用这些标签，S1 没有针对这些标签训练。二者的对照说明当前固定分折和模型下的可预测性与增量，不是同等训练条件的算法竞赛。Figure 5 没有区间，小的正增量不能写成稳定增益，小的负增量也不能写成显著损失。高 ECFP4 点估计也不说明已经找到来源偏差、数据泄漏或结合机制的原因。描述符和异类相似度只描述第 1 层组成，见 Figures S4–S6 及 Supplementary Data 中的 `module_D1_composition.csv`；模型表现与折外预测分别见 Table S7 和 `module_D_chemistry_oof_ligands.csv`。

## R5. 口袋对应提供的局部评分证据

Figure 6 以等比例坐标图比较同一组成员上的对应与非对应口袋 AUROC，完整差值区间见 Figure S12。D_vs_A 的对应口袋是 B，D_vs_B 的对应口袋是 A。差值为对应口袋−非对应口袋。16 个臂中 10 个点为正、5 个为负，JAK1/JAK2 的 D_vs_A 恰为 0；14 个区间包含 0。AChE/BChE 的 D_vs_A 为 0.1659，区间 0.0548 至 0.2830；PPARA/PPARD 的 D_vs_A 为 0.2354，区间 0.0693 至 0.4043。

16 个臂的非对应口袋额外减员都是 0。AChE/BChE D_vs_A 的 A_only 从标签层 27 人到正式 S1 的 25 人，发生在进入口袋对照之前，原因是 S1 超时，不是非对应口袋再删人。该对照只描述冻结口袋指定下的评分差异，不证明结合机制、姿势正确性或关键残基。结果汇总见 Table S6B，完整成员见 Supplementary Data。

## R6. 种子与受体条件的敏感性

Figure 7a 显示最小方向 AUROC 相对种子 42 的变化，行标签给出种子 42 的绝对值。种子为 17、29、42、71 和 101，24 个汇总都有 5 个有限点。图中的横线表示五种子相对变化的最小值到最大值，黑色短线标出中位数；观察范围不是置信区间。两个方向的完整结果在 Figure S7，汇总在 Table S8A，逐种子记录见 Supplementary Data。

PIK3CA/mTOR 的 summary_min 为 0.5596 至 0.5869，IQR 为 0.0190，五个种子的方向人数均为 dual 32、A_only 32、B_only 32。EGFR/HER2 的 summary_min 为 0.2305 至 0.2467，F2/F10 为 0.3594 至 0.3682；这两对的低表现没有因为更换种子而明显抬高。

除 AChE/BChE 外，七对的五种子有效人数不变。AChE/BChE 在 seed 17 的 A_only 为 24，其余种子为 25，因为 AB_046 在 BChE/4BDS 于 seed 17 超时。因此只有这一对的种子差异同时包含分数和成员变化。

当前替代受体结果覆盖 AChE/BChE、F2/F10 和 PPARA/PPARD 三对、六个条件，见 Figure 7b 和 Table S8B。未替换端的六行差值均为 0，这是只替换指定一端的构造核对。六个受影响端差值的区间均包含 0。PIK3CA/mTOR 的 Q110 未进行替代受体计算。

## R7. 方向判别变化与联合名单组成可一致或背离

联合名单来自 dual、A_only 和 B_only 构成的挑战集。名次使用两端 average midrank 的 worst rank，k 为 ceil(0.10×N)。Figure 8 的每个面板使用该方法与 S1 自己的共同成员。横轴是已有 summary_min 差值，纵轴是 Top-k 中单靶占比的变化，图中以百分点显示。预期同时有利的方向是横轴为正、纵轴为负。

S3 在 F2/F10 的 summary_min 差值为 0.2471，共同 Top-10 的 dual 数由 4 增至 7。PPARA/PPARD 为 0.1582，dual 由 3 增至 6。PPARG/PPARA 为 −0.2915，dual 由 6 降至 4。这三对的方向差值和名单变化同向。PIK3CA/mTOR 中，S3 的 summary_min 差值为 0.1240，共同 Top-10 的 dual 数由 4 增至 7；S2 的 summary_min 差值为 −0.0049，但 Top-10 的 dual 数由 4 增至 6。方向判别与名单组成可以背离。k 较小的名单里，一个分子就会带来较大的比例变化。

四态按舍入前点值计数，顺序为一致改善、一致变差、不一致、无变化：S2 为 2/2/3/1，S4 为 3/3/1/1，S3 为 5/3/0/0，S5 为 4/2/2/0。四种状态在全部比较中都出现。S3 的八个比较均为方向判别与名单变化一致，S5 没有无变化。无变化指两项差值至少有一项恰为 0，不是两项都不变。S4 的横轴是三组成员核对后已有 PRIMARY 点估计之和，不另带区间。各比较的 N、k 和 dual 数见 Table 2，共同名单的类别组成另见 Figure S9；方法自身的 S1 名单和四态计数见 Table S10A,C；方法间配对变化见 Table S10B。这些比例描述构造挑战集，不是天然库命中率或独立验证。

## R8. 限制同类骨架重复后的敏感性分析

在全部 867 条合格观察上，将同一靶点对、同一类别内的同一非空 Bemis–Murcko 骨架限制为最多 5 个分子后，保留 861 条、排除 6 条（Figure S10a）。排除发生在 JAK1/JAK2 的 3 条 A_only，以及 PPARA/PPARD 的 1 条 A_only 和 2 条 B_only。111 条 neither 均保留。CAP5 用于支持分析，主要评价仍采用未过滤的 867 条观察。

其余六对的 S1 summary_min 及区间不变。JAK1/JAK2 的点估计由 0.5908 变为 0.6315（CAP5 95% 区间 0.4871 至 0.7597）；PPARA/PPARD 由 0.5107 变为 0.4969（0.3427 至 0.6354）。两对的原始区间和 CAP5 区间均包含 0.5。16 个方向的 ECFP4 折外 AUROC 点变化绝对值均小于 0.02。

PPARA/PPARD 的 S5 在 dual 与 A_only 比较中的区间由包含 0.5 变为不包含 0.5，点估计由 0.6338 变为 0.6431。部分 Top-k 名单的成员和 k 值随排序池改变，不能直接沿用原 Top-10 的描述。人数、Vina 汇总及化学模型见 Table S9A,B,D，全部评分的方向性比较另见 Figure S8；完整方法结果、区间、排除成员及名单变化见 Supplementary Data 中的 CAP5 文件。PIK3CA/mTOR 的严格活性子集结果见 Figure S10b 和 Table S9C。
