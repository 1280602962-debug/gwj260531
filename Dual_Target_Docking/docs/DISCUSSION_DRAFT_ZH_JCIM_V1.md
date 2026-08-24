# Discussion（中文工作稿 · JCIM Articles）

## 4. 讨论

### 4.1 基准 formulation 改变了双靶对接的证据标准

本研究的首要发现不是某一种 docking scoring function 在双靶任务上取得了最高性能，而是 **benchmark formulation 可以改变“双靶对接成功”的表观含义**。传统 docking benchmark 通常将活性配体与 decoy 进行区分，而本文要求模型同时区分 dual-active 配体与两个方向上的 single-target selective hard negatives。后者在一个靶点上具有较强实验活性，因此不能被简单视为普通 decoy。同一套分数上的 Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是“the conventional dual-target benchmark”。供给审计显示，公开生物活性数据很难同时提供两个方向上足够数量的这类配体；49 个候选靶对中只有少数能够满足严格厚面板要求（Results 3.1；Figure 2）。

相对于 Zhou 等采用的 dual-target evaluation setting，[9] 本文进一步引入实验定义的 A-only/B-only 方向硬负，并直接比较不同 formulation 下的表观判别。EGFR/HER2 上，Dual versus neither 为 0.756，方向性 `summary_min` 为 0.430，混合库 Top-10 中 9/10 为选择性配体。**这是依赖靶对的清晰例子，不是四对定律或配对显著性检验：** AChE/BChE 与 PIK3CA/PIK3CB 只有小且区间重叠的增量，PIK3CA/mTOR comparator 则效能不足。

这一数据限制本身具有方法学意义。DUD、DUD-E 与 LIT-PCBA 已经表明，decoy construction、chemical bias 和真实 assay 标签会显著改变虚拟筛选的性能判断。^(5–7) 简单方法或不恰当的 unbiasing 也可以通过学习配体分布而高估 structure-based virtual screening。^(12) 近期基于 bioassay-derived data 的评价则进一步强调，与人工构造的 ligand/decoy 集相比，真实 assay-derived benchmarks 可以揭示模型在更接近实际筛选环境中的局限。^(13) DualFourClass-Bench 并未使用这些单靶数据集，也不评价 DiffDock-Pocket；它把同一关切延伸到双靶任务：评价结论取决于硬负样本如何被实验定义，而不是取决于候选靶对清单有多长。

因此，DualFourClass-Bench 的主要价值并不在于提供一个规模很大的数据集，而在于将双靶识别问题转化为一个实验标签驱动的 hard-negative discrimination task，并显式要求两个方向同时成立。该资源是 curated four-pair panel + evaluation protocol，不是 comprehensive dual-target suite。

错口袋对照属于同一证据标准下尚未解决的面板外失败模式。主面板上 pocket-matched 高于 wrong-pocket；unused-pool holdout 上点估计反转，但配对区间仍包含零（Results 3.5；Figure 6）。因此 wrong-pocket 不是在面板迁移下可靠的通用负对照；主面板 matched > wrong 也不作为口袋特异证明。

### 4.2 配体化学基线揭示了表观对接判别中的混淆来源

双靶点 docking 的困难并不能简单理解为两个单靶 docking 任务的性能相加。四个 `summary_min` bootstrap CI 均包含 0.5，因此在当前样本量下没有一个靶对获得排除随机水平的明确证据。在主分析冻结受体协议下，PIK3CA/mTOR 的点估计最高（0.692），但相对 best single-descriptor reference 的配对差值 CI 也包含 0（Results 3.2；Table S19）。因此该对是 **conditional directional signal**，不是可推广的成功案例。

AChE/BChE 上 TPSA 单独即可获得高于 docking 的 discrimination，而 ECFP4 scaffold-grouped baseline 在多个方向上进一步超过 docking（Results 3.3）。这一强 ligand-only baseline 表明，实验标签与无需受体信息即可利用的化学空间差异相关；因此 docking performance 必须相对于 ligand-only baselines 解释，而不能孤立阅读。在当前支架分组基准下，把 docking 加到 ECFP4 后只产生很小的 CV AUROC 增量改善（最大绝对变化 0.020），且若干方向为负。这一面板内结果不等于 docking 不提供任何额外信息或口袋特异信息。

化学型约束硬负给出同一方向的证据。T ≥ 0.7 的匹配子集为空。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only，0.691）降至 0.503（n_neg = 11），而远缘硬负（T < 0.3）升至 0.819。T ≥ 0.3 是 similarity-constrained subset，不是 chemically matched analogue set。

这一结果与近年来对虚拟筛选 benchmark 中化学偏倚的关注一致：简单模型或不恰当构造的 decoys 可以通过学习 ligand distribution 而获得看似优异的 performance。^(7,12) 如果不设置 A-only/B-only hard negatives 以及 ligand-property/chemical baselines，一个看似优秀的 dual-target docking result 可能实际上只是识别了与 dual label 相关的分子属性。

### 4.3 受体实现是评价条件的另一维度

PIK3CA/mTOR 提供了本研究中最值得进一步研究、但也最需要谨慎解释的案例。其主面板 summary_min 为 0.692，PM110 扩展面板为 0.648，unused-pool holdout 为 0.765。三者同向提示结果并非完全由 PM48 的具体成员构成所决定，但这些检查仍处于同一 ChEMBL 数据体系，只是稳定性分析而非独立外部验证。与此同时，该信号并非 receptor-invariant：替换 PIK3CA 后 summary_min 降至 0.486 和 0.505，而替换 mTOR 后为 0.639（Results 3.4；Figure 5A）。更准确的结论是在特定评价条件下观察到有限的计算判别信号，而不是可靠的生物学双靶识别。

同一套 PIK3CA 替换在 PIK3CA/PIK3CB 上方向相反：summary_min 由 0.500 升至 0.691 和 0.685，而 B 端受体保持冻结（Figure 5B）。因此受体选择不仅可以改变表观判别的幅度，也可以改变其方向。这一对比不支持把受体替换简单解释为 docking accuracy 的损失：受体实现是评价条件的另一独立维度，并可增强或削弱表观判别。两对、单口袋设计比单一 collapse 轶事更有解释力，但还不是普遍定律（K = 4；两对共享 PIK3CA）。不宣称具体分子机制。

耦合任务读法值得讨论，但仍是假说。`summary_min` 跟踪弱臂，因此替换 \(S_A\) 可以改变哪一臂成为瓶颈。PIK3CA/PIK3CB 原来的弱臂是 4L23 上的 D/B（0.500）；4JPS 后该臂升至 0.707，冻结 2WXF 臂（0.691）成为限制。其他非互斥假说包括局部侧链/口袋几何重排 dual 与选择性分数，以及两套面板面对同一 PIK3CA 晶体时不同的配体化学分布。由于没有面板范围的残基级 PLIF 分析，这些解释仍未解决。

5DXT 与 4L23 的局部口袋 Cα RMSD 仅为 0.343 Å，但 PIK3CA/mTOR 的 summary_min 仍降至 0.505。这说明“结构相似”与“screening discrimination 可迁移”并不是同一个问题。pose-generation QC 通过，也不等于 screening performance 不变。该结果与近期 cross-docking benchmark 中 receptor representation 被视为独立性能变量的观点相一致；那些工作使用的对接引擎与本文不同，不能当作同一协议的外推。^(14)

### 4.4 对双靶虚拟筛选与生成式设计的含义

这些结果对双靶点虚拟筛选与生成式设计具有直接意义。同时在两个口袋获得 favorable docking scores 并不能自动等同于 experimentally plausible dual-active ligands。若双靶生成模型把 docking 当作下游过滤器，就不应只报告 two-pocket score，而应在 selective hard-negative 和 ligand-only controls 下评价。即便在单靶超大规模对接之后，hit 的后处理与再打分也已被证明难以稳健地区分已知结合分子与无活性分子；^(15) 双靶场景额外要求同时压住两条实验硬负臂，因此更不能把两端有利分数或其简单平均当作充分证据。

本研究并不证明现有双靶生成模型无效，也没有直接评测 DualDiff、FuseDiff 或其他生成模型。^(10,11) DualDiff 的 Dual High Affinity 是生成分子在两个靶上均优于各自参考配体，不是均值池化；FuseDiff 的独立测试集为 DualDiff benchmark（DDF）。这些工作回答的是结构生成能否获得有利对接分数。DualFourClass-Bench 可以作为这类方法的 downstream evaluation layer，用于检验生成分子是否真正超过实验定义的 single-target hard negatives，而不是仅仅优化 docking score。

### 4.5 局限性

以下五条限制界定本文结论的解释范围。

第一，评价集仅含四对靶标，因为实验定义的双靶硬负样本稀缺。K = 4 是受数据供给约束的案例面板，不是 comprehensive suite。四个 `summary_min` 还混合了面板构建差异（严格 6.5/5.5 对 θ = 6.0；不等 n）与靶对生物学，不能读成 intrinsic docking performance 的总体排序。

第二，ground truth 来自 ChEMBL。unused-pool holdout 仍属同一抓取批次，不是跨数据库独立验证。BindingDB/PubChem 核对仅为计数。

第三，活性聚合对照之后，assay heterogeneity 仍然存在。主策展使用最大 pChEMBL；换成重复测定中位数后 pair-level 估计变化很小（Results 3.4；Table S29），因此 max 聚合是 controlled limitation，不是未关闭的 fatal ground-truth 风险。pChEMBL 仍非 assay-equivalent。confidence≥8 与 Homo sapiens 过滤未重建。

第四，受体实现可以提高或降低成对判别，但实验并未给出分子起源。两个受体敏感性例子均共享 PIK3CA，因此不能视为两个互不相关替换靶标上的独立证据。口袋局域 Cα RMSD 不能单独解释性能变化，残基级 PLIF/侧链分析未系统进行。PIK3CA/PIK3CB 的一次对接超时在原始与替代 PIK3CA 晶体上均一致报告，排除原因不是其标签（Table S30）。

第五，本研究评价的是计算判别，而不是对新预测双靶化合物做前瞻实验。benchmark 回答的是对接排序的可靠性，不是入选分子的前瞻生物学效力。本文不旨在证明 docking 对双靶发现普遍有效或无效；它问的是评价 formulation 本身是否改变双靶识别的表观证据。
