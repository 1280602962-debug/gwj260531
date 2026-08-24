# Discussion（中文工作稿 · JCIM Articles）

> 与 [`DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md`](DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md) 配套：本文件写解释、替代解释与使用边界；完整 13 条局限库存仍在 Limitations 稿，正文 4.6 只保留五条。收束段见 [`CONCLUSIONS_DRAFT_ZH_JCIM_V1.md`](CONCLUSIONS_DRAFT_ZH_JCIM_V1.md)。  
> 口径对照 [`DISCUSSION_RESULT_MAP_V1.md`](DISCUSSION_RESULT_MAP_V1.md)；引用核验 [`DISCUSSION_REFS_JCIM_V1.md`](DISCUSSION_REFS_JCIM_V1.md)。  
> 写法：Finding → Interpretation → Alternative explanation → Evidence → Implication。不把开放问题写成已解决。

---

## 4. Discussion

### 4.1 A Strict Dual-Target Benchmark Exposes a Task That Conventional Docking Evaluation Does Not Directly Test

本研究的首要发现不是某一种 docking scoring function 在双靶任务上取得了最高性能，而是 **benchmark formulation 可以改变“双靶对接成功”的表观含义**。传统 docking benchmark 通常将活性配体与 decoy 进行区分，而本文要求模型同时区分 dual-active 配体与两个方向上的 single-target selective hard negatives。后者在一个靶点上具有较强实验活性，因此不能被简单视为普通 decoy。同一套分数上的 Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是“the conventional dual-target benchmark”。供给审计显示，公开生物活性数据很难同时提供两个方向上足够数量的这类配体；49 个候选靶对中只有少数能够满足严格厚面板要求（Results 3.1；Figure 2）。

Zhou、Li 与 Hou 已经表明，对接用于双激酶筛选时，相对 noninhibitor 可以看起来有用，表现依赖结构，并且预测 dual 列表仍有较高 false-positive rate。^(9) DualFourClass-Bench 在同一套分数上追问更窄的问题：Dual-versus-neither（inactive）读出与方向性 Dual-versus-selective 读出是否一致。它们在 EGFR/HER2 上并不一致（Results 3.2，Table 3）：Dual versus neither 为 0.756，方向性 summary_min 为 0.430，混合库 Top-10 富集选择性配体（9/10）。该 gap 依赖靶对，不是四对 overestimation 定律，也是两种任务 formulation 的描述性对照，而不是配对显著性检验。相对 2013 年工作的增量是这一 formulation gap，而不是又一次四对对接普查。

这一数据限制本身具有方法学意义。DUD、DUD-E 与 LIT-PCBA 已经表明，decoy construction、chemical bias 和真实 assay 标签会显著改变虚拟筛选的性能判断。^(5–7) 简单方法或不恰当的 unbiasing 也可以通过学习配体分布而高估 structure-based virtual screening。^(12) 近期基于 bioassay-derived data 的评价则进一步强调，与人工构造的 ligand/decoy 集相比，真实 assay-derived benchmarks 可以揭示模型在更接近实际筛选环境中的局限。^(13) DualFourClass-Bench 并未使用这些单靶数据集，也不评价 DiffDock-Pocket；它把同一关切延伸到双靶任务：评价结论取决于硬负样本如何被实验定义，而不是取决于候选靶对清单有多长。

因此，DualFourClass-Bench 的主要价值并不在于提供一个规模很大的数据集，而在于将双靶识别问题转化为一个实验标签驱动的 hard-negative discrimination task，并显式要求两个方向同时成立。该资源是 curated four-pair panel + evaluation protocol，不是 comprehensive dual-target suite。

### 4.2 Why Dual-Target Docking Is Harder Than Two Single-Target Docking Tasks

双靶点 docking 的困难并不能简单理解为两个单靶 docking 任务的性能相加。若一个配体在靶点 A 上具有较强活性，则它在 A 口袋中获得有利 docking score 并不能说明其同时具有靶点 B 的活性；真正严格的 dual recognition 要求模型在另一靶点中进一步压制 A-only 配体，并同时在相反方向压制 B-only 配体。正因如此，单一 pooled score 容易掩盖较弱的一臂，而 pocket-matched directional AUROC 更直接地测试了“非选择性靶点上的额外识别”。

但即便采用这一任务对齐的评价方式，四对靶标仍表现出明显异质性。三个靶对的 summary_min 位于随机水平附近或低于 best single-descriptor reference，只有 PIK3CA/mTOR 表现出较高的点估计，且其置信区间仍与随机相容（Results 3.2）。该对是 **conditional directional signal**（配体面板可持续，受体实现可塌掉），不是可推广的成功案例。更重要的是，AChE/BChE 的结果显示，TPSA 单独即可获得高于 docking 的 discrimination，而 ECFP4 scaffold-grouped baseline 在多个方向上进一步超过 docking（Results 3.3）。由此可见，dual/selective 标签本身携带的 ligand-level information 可以产生强烈的 apparent signal。在当前支架分组基准下，把对接分数加到 ECFP4 后增量 AUROC 很小；这是对本面板的陈述，不是 docking 没有口袋特异信息的证明。

这一结果与近年来对虚拟筛选 benchmark 中化学偏倚的关注是一致的：简单模型或不恰当构造的 decoys 可以通过学习 ligand distribution 而获得看似优异的 virtual screening performance，因此 benchmark 必须区分 target-specific signal 与 chemical composition signal。^(7,12) 本研究进一步将这一问题扩展到双靶任务：如果不设置 A-only/B-only hard negatives 以及 ligand-property/chemical baselines，一个看似优秀的 dual-target docking result 可能实际上只是识别了与 dual label 相关的分子属性。

### 4.3 The PIK3CA/mTOR Case: Limited Directional Signal Rather Than a Generalizable Rule

PIK3CA/mTOR 提供了本研究中最值得进一步研究、但也最需要谨慎解释的案例。其主面板 summary_min 为 0.692，PM110 扩展面板为 0.648，unused-pool holdout 为 0.765，说明该方向性信号并非完全由 PM48 中少数配体驱动。与此同时，该信号并未表现出 receptor invariance：替换 PIK3CA 受体结构后，summary_min 降至 0.486 和 0.505，而替换 mTOR 结构后仍为 0.639（Results 3.4）。更准确的说法不是“PIK3CA/mTOR docking 可以可靠识别双靶配体”，而是该靶对在特定 receptor realization 下存在有限的 directional signal：配体面板替换后可以持续，受体实现替换后不能假定不变。

PIK3CA 和 mTOR 均具有可被 ATP-site chemotypes 访问的结合模式，因此某些 dual ligands 可以在两个口袋中形成合理的 hinge-oriented poses；但同样的 ATP-site compatibility 也可能使选择性硬负在另一靶点获得几何上合理的 pose，从而产生 false dual recognition。姿态级观察与这种可能性相符（Results 3.6 的 T2 / T5），但尚未进行完整的 residue-level PLIF analysis，因此不能将该解释提升为确定的结构机制。

更值得注意的是，5DXT 与 4L23 的局部口袋 Cα RMSD 仅为 0.343 Å，但对应的 summary_min 仍降至 0.505。这说明“结构相似”与“screening discrimination 可迁移”并不是同一个问题。pose-generation QC 通过，也不等于 screening performance robustness。该结果与近期 cross-docking benchmark 中 receptor representation 被视为独立性能变量的观点相一致；那些工作使用的对接引擎与本文不同，不能当作同一协议的外推。^(14)

### 4.4 Implications for Dual-Target Virtual Screening and Generative Design

这些结果对双靶点虚拟筛选与生成式设计具有直接意义。对于生成模型而言，同时在两个口袋获得 favorable docking scores 并不能自动等同于生成了 experimentally plausible dual-active ligands。若 scoring function 本身受到 ligand size、polarity 或 chemotype distribution 的影响，那么生成模型可能通过优化这些易被 scoring function 奖励的属性而获得较高的 dual docking score，而未真正获得两个靶点上的独立结合优势。即便在单靶超大规模对接之后，hit 的后处理与再打分也已被证明难以稳健地区分已知结合分子与无活性分子；^(15) 双靶场景额外要求同时压住两条实验硬负臂，因此更不能把两端有利分数或其简单平均当作充分证据。

因此，dual-target generative design 的 downstream evaluation 应至少包含三个层次：首先，dual-active 与 A-only/B-only selective hard negatives 的实验标签驱动 discrimination；其次，ligand-property 和 ligand-only chemical baselines；第三，receptor-structure sensitivity。现有 benchmark 的结果表明，只报告两个 pocket 的 docking scores 或其简单平均值，无法充分回答这些问题。

本研究并不证明现有双靶生成模型无效，也没有直接评测 DualDiff、FuseDiff 或其他生成模型。^(10,11) DualDiff 的 Dual High Affinity 是生成分子在两个靶上均优于各自参考配体，不是均值池化；FuseDiff 的独立测试集为 DualDiff benchmark（DDF）。这些工作回答的是结构生成能否获得有利对接分数。DualFourClass-Bench 可以作为这类方法的 downstream evaluation layer，用于检验生成分子是否真正超过实验定义的 single-target hard negatives，而不是仅仅优化 docking score。

### 4.5 Wrong-Pocket Reversal Is an Unresolved Benchmark Failure Mode

本研究还暴露出一个目前尚未解决的 benchmark failure mode：wrong-pocket control 在主面板与 unused-pool holdout 中表现出相反关系。在主面板中，pocket-matched score 均高于 wrong-pocket control；但在 holdout 中，三对靶标均表现为 wrong-pocket score 不低于 matched-pocket score。效价和尺寸匹配不能消除这一反转，而 scoring-independent contact_count 仅解释了其中部分 B-arm signal（Results 3.5）。

这一现象不应被简单归因于某一 docking engine 的 scoring artifact，因为它在同一 Vina protocol 下出现，并且粗粒度几何指标也显示了 ligand-size/burial-related signal。然而，contact_count 的效应幅度不足以重现 Vina wrong-pocket discrimination，因此当前数据尚不能确定其唯一来源。可能因素包括 ligand distribution shift、pose selection、receptor-specific interaction patterns 以及 scoring-function 的非线性 size dependence，但本研究未对这些机制进行系统拆分。主面板 matched > wrong 不作为口袋特异证明。wrong-pocket **不是在面板迁移下可靠的通用负对照**。

对 benchmark 而言，这一“未解决结果”本身具有价值。它说明一个在固定 benchmark 上表现合理的 pocket-specificity control，并不能保证在未见配体池中仍然成立。因此，未来双靶 docking benchmark 不应仅报告 matched-pocket performance，而应同时报告 wrong-pocket、chemical-property controls 和 panel-external holdout。

### 4.6 Limitations

正文只强调五条最高优先级限制；完整清单见 Limitations 稿。收束主张见 Conclusions，不在此重复。

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. K = 4 是 data-constrained case panel，不是 comprehensive suite。四个 `summary_min` 还混合了面板构建差异（严格 6.5/5.5 对 θ = 6.0；不等 n）与靶对生物学，不能读成 intrinsic docking performance 的总体排序。Existing ligand-level bootstrap (B = 2000) describes uncertainty inside a fixed panel; leftover strict hard negatives after the main panels and holdout cannot support anything close to 1000 non-overlapping balanced panels.

Second, activity labels were aggregated using the maximum available pChEMBL value. This may inflate apparent activity when measurements vary across assays. The current frozen data package does not retain sufficient assay-level metadata to perform a complete max-versus-median or assay-confidence sensitivity analysis.

Third, the unused-pool holdout remains within the ChEMBL-derived data ecosystem and therefore should not be considered an independent cross-database validation.

Fourth, receptor-swap experiments demonstrate structure dependence but do not identify its molecular origin. In particular, pocket-local Cα RMSD alone could not explain the observed performance change, and residue-level PLIF/side-chain conformational analyses were not systematically performed.

Finally, this study evaluates computational discrimination rather than experimentally testing newly predicted dual-target compounds. The benchmark therefore addresses the reliability of docking-based ranking, not the prospective biological efficacy of the selected molecules. 本文不旨在证明 docking 对双靶发现普遍有效或无效；它问的是评价 formulation 本身是否改变双靶识别的表观证据。
