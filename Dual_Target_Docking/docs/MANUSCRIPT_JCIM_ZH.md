# 基于对接的双靶识别：四靶对的评价设定审计

## 摘要

两端有利的对接分数能否作为双靶活性的证据，取决于评价时采用的负类。我们用 ChEMBL 衍生的 dual、A-only、B-only 与 neither 操作性状态，对四个靶对进行评价设定审计。两条口袋匹配方向 AUROC 分别将 dual 与对应单靶选择性配体比较；二者最小值 `summary_min` 仅作保守的描述性摘要。在 EGFR/HER2 上，Vina 的 Dual versus neither AUROC 为 0.756，而方向性 `summary_min` 为 0.430；独立 GNINA 姿态生成得到的对应数值为 0.783 与 0.220。其他靶对没有出现同样的设定差距，且 PIK3CA/mTOR 的 neither 对照只有 4 个分子。在支架分组模型中，把 docking 加入 ECFP4 后交叉验证 AUROC 的最大绝对变化为 0.020。替代受体使 PIK3CA/mTOR 的 `summary_min` 从 0.692 降至 0.486/0.505，但使 PIK3CA/PIK3CB 升高，说明结果依赖受体实现而不是具有结构稳健性。四个主 `summary_min` 的 95% CI 均包含 0.5。同一套冻结 EGFR/HER2 分数上，以 Dual 中位 `vina_worst` 做 AND 式双口袋过滤时，通过者多数仍是实验选择性配体（precision 0.298；硬负比例 0.702）。在四对完整 ChEMBL 图上，仅用配体 ECFP4 仍比 Dual versus 选择性更容易分开 Dual versus neither，说明设定问题不是 n ≈ 28 对接面板的抽样伪影。按文献 `document_id` 阻断后，EGFR/HER2 的弱方向臂仍为 0.430；预先冻结的 2018 文献年份分割没有两个可评估靶对，因此不作为外部验证。BindingDB 202608 原生归档在文献、结构与 ECFP4 < 0.70 过滤后，也没有两对达到预先冻结的主外部门槛，故未对接。因此，这一受数据供给约束的案例面板支持把选择性硬负与混淆感知对照作为评价要求，但不能建立靶标通用的对接性能或生物学识别结论。

**关键词：** 双靶对接；基准设定；选择性硬负样本；化学混淆；受体实现；虚拟筛选

## 1. 引言

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多药理学（polypharmacology）研究的重要方向。[1] 近年来，多靶点小分子的理性设计逐渐由经验性筛选转向结合结构生物学、计算化学与生成式模型的结构导向设计。[2] 分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的构象，再用打分函数对配体–受体互补性排序。[3][4] 因此，在双靶点药物发现中，一个自然策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在双靶结合能力。对接结果的解释高度依赖数据集构建。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观富集退化为粗粒度配体性质分离。[5][6] LIT-PCBA 采用实验 assay 标签并控制已知 decoy 与化学偏倚。[7] CASF-2016 评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

一个严格的双靶评价需要区分 **dual-active**、**A-selective**、**B-selective** 与 **neither** 四种实验状态（Figure 1A）。A-only 与 B-only 是该任务的**选择性硬负样本（selectivity hard negatives）**：它们在一个靶点上已有较强活性。计算终点因而检验对接能否在两个方向上将 dual-active 与对应单靶选择性配体区分开。Zhou、Li 与 Hou 曾在四对激酶上评价相对非抑制剂的 dual-target docking。[9] 本文在该设定上引入实验定义的方向性硬负，并在同一套分数上比较不同基准设定。Dual versus neither 按实验 inactive 计分，作为基准设定对照。平衡的四状态面板还受两端可比较测量与双向选择性硬负供给的约束。

近期双靶生成方法仍使用相对参考配体的对接成功指标。[10][11] 这些指标评价相对参考配体的计算双靶设计，而本文基准检验的是相对实验定义选择性硬负样本的判别。

本文要问的是：基准设定是否改变双靶识别的表观证据。我们构建 DualFourClass-Bench，作为具有两条方向主任务的四状态面板：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B），并以保守的最差方向判别摘要（`summary_min`）汇总。我们进一步考察该判别是否能够在不同配体、活性聚合方式及受体结构条件下保持。

## 2. 方法

### 2.1 数据与实验状态定义

双靶评价所需的配体活性通过 ChEMBL Web API 的公开 activity 端点获取，作为实验衍生标签。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干经标准化的定量效力或亲和力测量（如 IC50、EC50、Ki、Kd 和 Potency）转换为近似 −log10 活性尺度。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似。

同一配体–靶标若有多条可用 pChEMBL 记录，主策展采用**最大 pChEMBL** 作为一对一代表值。活性聚合敏感性分析从 ChEMBL activity 端点重拉 assay 级记录，并在同一 θ = 6.0 规则下用重复测定的**中位数**替换最大值，不改变面板成员、对接参数或 Vina 分数。API 重拉的最大值对中位数估计作为标签聚合敏感性，与 Table 2 并列报告（Table S29）。作为 post-hoc 标签稳健性，已打分配体–靶标记录于 2026-08-26 从当前 ChEMBL API 重拉，仅保留 Homo sapiens `SINGLE PROTEIN`、assay confidence ≥8、等式关系、IC50/Ki/Kd/EC50/Potency、无 validity 注释且无 `potential_duplicate` 的记录，再按 θ = 6.0 重标而不改变面板成员或对接分数（Table S36）。该日期化视图检验这些显式过滤器是否改变冻结标签，不是 2026-07-23 数据库状态的重建，也不统一 assay 条件、蛋白构建体或突变背景。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。ChEMBL 结构按连通片段拆分，保留重原子数最多的有机片段。

完整病例选择按冻结可用 pChEMBL 映射统计只在 A、只在 B 或两端均有值的结构（Table S37）。缺失不解释为无活性。来源文献集中度按高置信保留记录、分实验状态统计独立 `document_id` 与最大单篇文献份额。

对每一对靶标 A/B，配体被定义为四种实验状态：**dual**（两端较强）、**A-only**（仅 A 端较强）、**B-only**（仅 B 端较强）和 **neither**（两端均不足）。A-only 与 B-only 是选择性硬负样本。

**严格 6.5/5.5 规则**为：dual，两端 pChEMBL ≥ 6.5；A-only，A ≥ 6.5 且 B ≤ 5.5；B-only 对称；neither，两端 ≤ 5.5。5.5–6.5 灰区不进入该审计。金属依赖体系（如 HDAC）预先排除。严格 6.5/5.5 规则仅用于靶对供给资格审定。全部主基准标签随后统一按 θ = 6.0 定义（dual，两端 ≥ θ；A-only，A ≥ θ 且 B < θ；B-only 对称；neither，两端 < θ），并在查看对接结果之前冻结。两条阈值因此服务于预先规定的不同目的。建造规则在抽样前按供给审计冻结（Table 1）。θ ∈ {5.5, 6.5} 与严格 6.5/5.5 重标作为敏感性报告（Table S4）。样本量过小的格子在 Results 中标记效能不足。

为核对 ChEMBL 供给门槛，对冻结靶对另做 BindingDB / PubChem 计数核对（零对接、不重建面板；Table S12）。类型限于 IC50/Ki/Kd/EC50；配体身份分别用 BindingDB monomerid 与 PubChem CID，不做跨库结构合并。主比较采用等式测定。

随后从 BindingDB 202608 文章与专利 TSV 归档重建原生切片，而不是从 Table S12 的 REST pmax JSON 反推。[16] 规则预先写在 `external_slice_contract.yaml`：BindingDB 策展的文章或专利；人源野生型单链 UniProt；等式 IC50/Ki/Kd；两端均有测定；配体–靶–endpoint 内取中位数；θ = 6.0 四状态；去掉与开发面板共享的 PMID/DOI/专利；去掉已打分面板、未使用池留出集、PM110 或该对 ChEMBL 图中的 InChIKey/ChEMBL ID；对开发分子的最大 ECFP4 Tanimoto < 0.70。主外部门槛为 dual/A-only/B-only 各 n ≥ 20、每类至少 3 个来源、最大单文献配体份额 ≤ 50%。本会话未对接。ChEMBL 文献解析为 519/680，因此剩余计数是完全文献独立集的上界。

### 2.2 基准构建

DualFourClass-Bench 保留四种实验状态。主分析由两条方向性成对任务组成。neither 类保留给描述性的基准设定对照。

候选靶对按 2.1 的严格供给审计筛选。冻结评价集包含 PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB 与 EGFR/HER2。EGFR/HER2 保留为供给受限案例。配体按预先冻结的类别配额抽样，随机种子为 20260729。抽样时结构可用的面板施加 Bemis–Murcko 支架封顶：PIK3CA/mTOR（PM48）同一类别内同一支架最多 2 个分子，EGFR/HER2 最多 5 个。AChE/BChE 与 PIK3CA/PIK3CB 只采用类别配额和确定性随机顺序。观察对接分数后不再重抽面板。

AChE/BChE 与 PIK3CA/PIK3CB 按严格供给门槛抽样（目标 28 / 28 / 28 / 16；n_panel = 100）。EGFR/HER2（n_panel = 110）与 PIK3CA/mTOR PM48（n_panel = 48；建造 18 / 14 / 12 / 4）按主分析 θ = 6.0 标签构建（Table 1）。因此跨对 AUROC 同时混合靶对生物学与面板构建差异。对接失败的配体–受体组合被剔除，n_scored 可低于 n_panel（Table 1；Table S27）。扩面面板 PM110 保留 PM48 全部 48 个配体，用于检查面板规模增加后点估计是否同向。

**Table 1.** DualFourClass-Bench 评价集组成与对接设置。建造标签记录各靶对的供给/建面规则；Tables 2–3 的全部主 AUROC 均使用统一 θ = 6.0 实验状态标签。n_panel 为冻结面板成员数（含 neither）；n_scored 为两端均有有效 Vina 分数、进入方向性主 AUROC 的 dual / A-only / B-only 计数。

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | n_panel | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.3 受体准备与对接协议

受体取自含小分子共晶配体的 PDB 条目：PIK3CA/mTOR，4L23 / 4JT6（X6K / PI-103）；AChE/BChE，4EY7 / 4BDS（E20 / THA）；PIK3CA/PIK3CB，4L23 / 2WXF（X6K / 039）；EGFR/HER2，3POZ / 3RCD（03P / TAK-285）。结合位点由共晶配体定义。以共晶配体重原子计算轴对齐包围盒，三方向各外扩 5 Å；任一边若小于 20 Å，则设为至少 20 Å（Table S2）。去除水分子与共晶配体后，用 Meeko 生成 PDBQT。PIK3CA、mTOR、EGFR 与 HER2 使用冻结目录中已含氢的蛋白坐标（`mk_prepare_receptor.py --read_pdb`）。AChE、BChE 与 PIK3CB 从沉积 ATOM/TER 记录提取，并以 `mk_prepare_receptor`（默认 alternate location A）转换。主分析均为非共价小分子对接。

正式对接前对每个冻结受体做共晶配体重对接。生成 9 个姿态，计算与实验共晶构象的重原子 RMSD。预先通过标准为 \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å，即九个保留姿态中是否存在与共晶配体重原子 RMSD 小于 2.0 Å 的构象。若默认 exhaustiveness 未通过门槛，则提高至预先规定的备用水平。主分析因此采用 PIK3CA/mTOR exhaustiveness = 16、其余主面板为 8（Table S3）。EGFR/HER2 原始九姿态生产 PDBQT 未能找回，已按冻结协议重对接并标为 reconstructed QC，而非历史生产文件。拓扑核对后的 ranked RMSD：EGFR 3POZ top-1 9.505 Å、top-3 6.227 Å、best-of-9 0.760 Å（top-1/top-3 未过，搜索覆盖通过）；HER2 3RCD top-1 1.855 Å、top-3 1.394 Å（通过）。这些数值替换原先 NA 的 top-3 单元格，不改变预先规定的 best-of-nine 生产门槛。

配体从冻结 ChEMBL SMILES 统一准备：去盐并保留最大有机片段，RDKit 加显式氢，ETKDGv3 生成三维构象（种子 20260727），MMFF 局部优化最多 200 步，再经 Meeko 转为 PDBQT。不进行系统性质子化、互变异构或构象枚举。对接采用 AutoDock Vina 1.2.7 默认 `vina` 打分函数，保留 9 个姿态，`energy_range = 3` kcal mol\(^{-1}\)，随机种子 20260727（Table S1）。为检验打分函数依赖性，同一组 Vina 姿态另用 RTMScore（`rtmscore_model1`，取九姿态最高分）与 GNINA 1.3.2 CNN（`--cnn_scoring rescore --minimize`，Open Babel 转 SDF 后取九姿态最高分）重打分。Vina 主读出是 mode-1 能量；RTM 与 GNINA CNN 是 best-of-9 重打分。主终点始终由 Vina 定义。

另在 EGFR/HER2 与 PIK3CA/mTOR 上以 GNINA 1.3.2 对接搜索模式独立生成姿态（不是对 Vina 姿态重打分），复用冻结 Meeko 配体 PDBQT、受体坐标、对接盒、exhaustiveness（分别为 8 与 16）、九个保留姿态和种子 20260727。读出为 mode-1 `minimizedAffinity`。两个配体在两端口袋均失败（EGFR/HER2 的 neither 配体 EH120_109；PIK3CA/mTOR 的 A-only 配体 PM48_19），从需要完整分数的分析中剔除。该协议检验设定效应在更换姿态生成引擎后是否仍在，不是多引擎比赛（Table S32）。

### 2.4 主终点与统计分析

全文中的“双靶识别”指这一计算判别任务。对每个靶对计算两条二分类 AUROC。dual 对 A-only 使用口袋 B 分数，\( \mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B) \)；dual 对 B-only 使用口袋 A 分数，\( \mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A) \)。dual 始终为正类。Vina 输出结合能 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，越负表示预测结合越强）；\(S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}\)。

最差方向判别摘要定义为 \( \mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}) \)。它把两条方向 AUROC 保守地汇总为单值。算术平均、几何平均与调和平均作为聚合敏感性报告（Table S26）。全文唯一主终点是统一 θ = 6.0 下的口袋匹配 Vina `summary_min`（Table 2；PIK3CA/mTOR 主面板为 PM48）。预先指定的 RDKit 描述符面板（重原子数、分子量、cLogP、TPSA）按同一方向流程评价；其中 AUROC 最高者记为最佳单一描述符参考（Tables 2、S28、S19）。Dual versus neither（实验 inactive；`vina_mean`）与 Dual versus all non-duals 为同一套冻结分数上的基准设定对照（Table 3；Table S22）。PIK3CA/mTOR 的 neither n = 4 标记效能不足。

AUROC 与 summary_min 的不确定度用配体层 bootstrap：在保持类别结构的条件下对配体有放回重采样（\(B = 2000\)，种子 20260729，百分位数 95% CI）。配对比较在同一次重采样上计算（Tables S17、S19）。置信区间作描述性不确定度。可分辨效应模拟使用观察得的类别样本量、同一 bootstrap、双正态分数模型和一组真实 AUROC，报告 95% CI 排除 0.5 的概率，而不是观察后功效（Table S31；Figure S6）。

### 2.5 混淆、留出集与受体敏感性分析

将靶点 A 与 B 的分数对调作为证伪对照，配体、受体与其余设置不变。另在配体效率归一（\(S_{\mathrm{dock}}/N_{\mathrm{heavy}}\)）、效价约束（\(|\Delta\mathrm{pChEMBL}| \leq 0.5\)）和尺寸约束（\(|\Delta N_{\mathrm{heavy}}| \leq 2\)）后重算方向 AUROC。逻辑回归比较 docking alone 与 docking + 重原子数 + TPSA。Morgan/ECFP4（半径 2，2048 bit）加逻辑回归在 Bemis–Murcko 支架 `GroupKFold` 下提供配体化学基线（Tables S5、S20、S23、S24）。最近邻 Tanimoto 子集只作诊断。接触计数与全链序列一致性仅为探索性对照（Tables S7、S11）。

为检验结论是否依赖于冻结面板的具体成员，排除已用于主面板与 PM110 的 ChEMBL 条目后，在剩余未使用配体池中构建留出集（holdout）。配体仍来自同一 ChEMBL 抓取批次、同一靶对与同一标签规则。该留出集在 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 上构建（各 20 dual / 20 A-only / 20 B-only；`HOLDOUT_SEED = 20260731`）；EGFR/HER2 不具备同等抽样条件。受体、盒子、配体准备、exhaustiveness、打分与统计与主基准相同（Tables S8、S13）。

文献阻断分析使用同一套冻结分数。共享任一保留高置信 `document_id` 的配体连成一组，使同一篇文献的化合物不能同时进入训练与测试。`GroupKFold` 按该组划分；ECFP4、物化描述符与 docking logistic 使用相同折（Tables S39、S40）。缺少正负两类的折被丢弃；若有效折少于 2，则报告无法稳定估计，而不在看到 AUROC 后更换分组规则。document-cluster bootstrap 重采样的是这些文献连通组，而不是把同一系列化合物当作独立观察。

文献年份分割在计算 AUROC 之前冻结（`docs/TIME_SPLIT_PROTOCOL_FREEZE.md`）。配体年份取其保留高置信记录中最早的 `document.year`。主截止年为 2018（训练：first year < 2018；测试：first year ≥ 2018）；2015 与 2020 为预先指定的敏感性。晚期文献中的化合物不参与阈值、受体或指标选择。仅当测试集 dual、A-only、B-only 每类 n ≥ 10 时报告方向 AUROC；更小格子只报计数（Table S41）。主截止年至少两个靶对通过该门槛，才包装为外部验证。

对 352 个已打分配体中的 186 个优先分子提取了 assay-context 字段，包括 EGFR/HER2 全部方向类、PIK3CA/mTOR 的 4 个 neither、混合端点记录，以及对主 AUROC 影响最大的分子（`assay_context_audit.csv`）。随后一次元数据审核为全部 186 个优先分子填写纳入/排除（179 include / 7 uncertain / 0 exclude）；ChEMBL assay 自由文本不可用，故蛋白构建体与突变状态仍为 unknown。冻结 DualFourClass 标签未改，因此未重算 Table 2。该步骤不是文献级 assay 条件统一。

受体结构敏感性分析另选满足以下预先声明条件的替代晶体：（i）polymer entity 与目标蛋白真实对应；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.3 相同的共晶重对接 QC。实际对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。替换采用单口袋设计：在 PIK3CA/mTOR（PM48）上，4JPS/5DXT 替换口袋 A、口袋 B 仍用冻结 4JT6 分数，4JSX 替换口袋 B、口袋 A 仍用冻结 4L23 分数（exhaustiveness = 16）。在 PIK3CA/PIK3CB 上，同一套 4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 2WXF 分数（exhaustiveness = 8）。刚体 Cα 叠合作为探索性几何对照（Table S10）。在 Table S30 所用同一套 PM48 配体与 PIK3CA 晶体上，另做探索性接触快照：占有率定义为与 20 个冻结口袋残基的重原子距离 ≤ 4.5 Å（Table S33）。占有率变化只作为结构假说，不是残基层因果解释。

事后 θ = 6.0 四状态普查复用冻结的 J0 候选靶对名单与缓存 pChEMBL 图，去掉顺序别名后剩 49 对。dual/A-only/B-only 均 n ≥ 10 记为方向可评估，neither 也 n ≥ 10 记为设定可评估（Table S44）。这些计数只诊断标签供给，不对额外靶对做对接，也不把对接评价集扩到 K = 4 以外。冻结已打分面板上的多元物化匹配按 z 标准化 MW/cLogP/TPSA/重原子做 1:1 贪心配对，欧氏 caliper 为 0.5 与 1.0 SD（Table S45）；n_matched < 8 标记效能不足。AND 式双口袋过滤在 Dual+A-only+B-only 库上按 `vina_worst` 或 `vina_mean` 的 Dual 百分位截断（Table S46）。配体层 ECFP4 与四描述符逻辑回归则在四对完整 θ = 6.0 ChEMBL 图上拟合，每类最多抽 120 个分子（种子 20260729），支架 `GroupKFold`（Table S47）。该分析只用实验标签与二维结构，不是对接结果，也不替换 Table 2。

MCL1/Bcl-xL 冻结为 PPI/BH3 槽域外推候选，不是异质折叠对，也不是“首次非激酶对”（AChE/BChE 已是非激酶）。ChEMBL θ = 6.0 图为 dual/A-only/B-only/neither 82/77/24/122；按种子 20260729 抽 24/24/24/24，B-only 24 个全部纳入。主受体为 MCL1 3WIY 与 Bcl-xL 3WIZ（LC6 / Tanaka compound 10），替代 holo 为 6UDV 与 3SP7，选择依据是分辨率与野生型占位，不是 AUROC。[17] LC6 pose-gold gate 使用与生产相同的 Vina 协议（种子 20260727，exhaustiveness 8，九个姿态）。3WIY 通过 best-of-top3（1.689 Å）；3WIZ 失败（best-of-top3 2.011 Å；top-1 4.17 Å；Table S51）。因此该对只作为预先声明的 applicability stress-test 对接，不是第五个 Table 2 靶对（Table S53）。与 Zhou 2013、DUD-E、LIT-PCBA、CASF-2016、DOCKSTRING 的对照见表 S52。[5–9,18]

计算在 Python 3 环境下完成，主要软件为 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore。评价面板、对接分数、分析脚本与参数表见 Data and Software Availability。评价合约见 `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`。

## 3. 结果

### 3.1 实验数据供给限制了严格双靶基准的构建

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例（Table 1）。BindingDB 与 PubChem 的零对接计数核对支持同一供给稀缺结论（Table S12）。

因此，最终基准的规模主要由公开实验数据中方向性选择性硬负样本的可获得性所约束。严格 6.5/5.5 规则用于量化供给并记录面板构建，而 θ = 6.0 定义全部主 AUROC 的实验状态标签（Methods 2.1）。对同一 49 对在主规则 θ = 6.0 下重计（不再排除 5.5–6.5 灰区）后，有 17 对 dual/A-only/B-only 均 n ≥ 10，且这 17 对 neither 也均 n ≥ 10（Table S44；Figure S7A）。其中 16 对非金属酶。对接评价仍是原来的四对；普查是标签供给结果，不是 17 对对接基准。

### 3.2 基准设定改变了表观双靶判别

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina 对接分数进行评价（Figure 1B；Methods 2.4）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。预先指定的最差方向判别摘要为 `summary_min`。算术、几何与调和平均下四对排序不变（Table S26）。AChE/BChE 与 PIK3CA/PIK3CB 上，θ = 6.0 与更严格建造门槛给出相同分类与 AUROC（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 在严格规则下 B_only 过少。整张阈值网格内排序趋势保持一致（Figure S1A）。

EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的方向性 summary_min 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 4A）。EGFR/HER2 的 dual-versus-B-only AUROC 为 0.430，PIK3CA/PIK3CB 为 0.500；PIK3CA/mTOR 两个方向分别为 0.714 和 0.692。相对池化协议，口袋匹配抬高了点估计但未改变排序（Table S6）。

同一套冻结分数再按 Dual versus neither 以及 Dual versus all non-duals 计分（Table 3；Figure 3）。Dual versus neither 是使用实验 inactive 的描述性对照（`vina_mean`）。EGFR/HER2 提供了最清晰的设定对照。Dual versus neither 的 AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性 summary_min 仍为 0.430 [0.282, 0.578]。Dual versus all non-duals 降至 0.551 [0.443, 0.666]，表明引入选择性配体后判别任务明显变难。在 110 个 EGFR/HER2 配体的混合库中按 `vina_mean` 取 Top-10：1 个 dual、5 个 A-only、4 个 B-only、0 个 neither（EF10 = 0.393；硬负比例 = 0.90）；EF5 也低于随机（Table S25）。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小（0.649 与 0.559），区间与方向性臂重叠。PIK3CA/mTOR Dual versus neither 因 neither n = 4 而效能不足；该对 Dual versus all non-duals 为 0.674，接近 summary_min 0.692。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0），并列出四个预先指定描述符的 `summary_min`。表中类别样本量为 n_scored（dual / A-only / B-only），即对应方向分析中具有有效 Vina 分数的配体数，因而可以低于冻结面板定额；neither 不进入主终点。n_panel 与两端对接覆盖见 Table 1 和 Table S27。最高描述符是最佳单一描述符参考。错口袋与配体效率见 Table S6；描述符双臂见 Table S28。

| 靶对 | n_scored (dual / A-only / B-only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** 同一套 Vina 分数在 Dual-versus-neither 与方向性设定下的 AUROC（统一 θ = 6.0）。Dual-versus-neither 使用实验 inactive（`vina_mean`）；Dual versus all non-duals 把 A-only、B-only 与 neither 都计为负类。方向性 CI 来自 Table 2。负样本集合不同。PIK3CA/mTOR Dual versus neither 效能不足（n_neg = 4）。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.350, 0.650] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

四个靶对的 `summary_min` 95% bootstrap CI 均包含 0.5；因此在本研究的样本量下，没有一个靶对获得排除随机水平的明确证据。按观察得的类别样本量做可分辨效应模拟表明，当前样本更容易分辨较大的方向性效应，而对中等效应则较弱。当两臂真实 AUROC 均为 0.70 时，`summary_min` CI 排除 0.5 的概率在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 上分别为 0.62、0.50、0.56 和 0.22；真实 AUROC 为 0.60 时，相应概率为 0.03–0.07（Table S31；Figure S6）。因此，CI 未能排除 0.5 并不能建立与随机等价。

主面板两端均得分：EGFR/HER2 110/110，AChE/BChE 95/100，PIK3CA/PIK3CB 99/100，PIK3CA/mTOR 48/48（Table S27）。一个 A-only 配体因计算超时而持续无法完成 PIK3CA 对接，因此从需要该分数的分析中剔除（Tables S27、S30）。AUROC 因此以 AutoDock Vina 能够处理的化合物为条件。同一组 Vina 姿态上的替代打分器未改变总体排序（Tables S14–S15；Figure S1B）。独立姿态生成见 Results 3.7。

### 3.3 配体性质与化学型解释了相当一部分表观信号

口袋匹配对接首先与四种预先定义的物化性质进行比较（Figure 4B；Table 2）。相对于每个靶对的最佳单一描述符参考，docking summary_min 的配对差值在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 中分别为 −0.052、−0.128、−0.122 和 +0.229；四个 95% 置信区间均包含 0（Table S19；Figure S3C）。即使 PIK3CA/mTOR 的点估计表现出最大的正向差异，现有样本仍不足以将其与配体性质参考明确区分。

AChE/BChE 提供了一个较为直接的混淆案例。dual 配体平均 TPSA 约为 75，而选择性硬负配体约为 51（Figure 4C）；TPSA 单独获得约 0.769 的 AUROC，高于相同比较下的 Vina（约 0.56）。进一步加入重原子数和 TPSA 后，dual-versus-B-only AUROC 从 0.606 增至 0.807，而对接分数的优势比仅约为 1.18（Figure 7C）。该结果表明，该方向上的部分表观对接判别可由配体物化信息解释。

PIK3CA/mTOR 的情况有所不同。加入重原子数和 TPSA 后，AUROC 的变化约为 +0.07 至 +0.11，对接分数的优势比约为 2.19 和 3.08。与描述符的配对差值置信区间仍包含 0。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

二维化学结构基线进一步说明了这一问题（Figure 7A）。ECFP4 逻辑回归在 Bemis–Murcko 支架 GroupKFold 下多个方向获得约 0.78–0.91 的 fold AUROC，明显高于部分对应对接对照，例如 EGFR/HER2 dual-versus-B-only 中 ECFP4 AUROC 为 0.89，而对接 AUROC 仅为 0.43。实验标签因此包含可由配体二维结构表征捕获的信息，即使不使用受体信息也可以获得较强判别。PIK3CA/mTOR 上 \(n_{\mathrm{scaffolds}} \approx n\)，该折接近 leave-one-scaffold。同一设定下随机 `StratifiedKFold` 相对支架折的平均差为 +0.026（Table S20；Figure S3D）。

在当前支架分组基准下，加入口袋匹配对接分数后的最大绝对变化为 0.020（未四舍五入值：PIK3CA/mTOR dual versus A-only 的 −0.0198），且若干方向为负（Table S24）。T ≥ 0.3 时，PIK3CA/PIK3CB dual versus A-only 从 0.691 降至 0.503（n_neg = 11），远缘硬负（T < 0.3）为 0.819；T ≥ 0.4/0.5 的格子常为 n_neg ≤ 7，T ≥ 0.7 为空（Table S23）。效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52），各臂 n 常低于 15（Table S5；Figure 7D）。

### 3.4 活性聚合、配体面板与受体实现影响评价结果

主标签使用可用 pChEMBL 的最大值。将最大值换成重复测定中位数后，改变了 7/110、1/95、1/99 和 0/48 个配体状态分配，标签一致率分别为 93.6%、98.9%、99.0% 和 100%。API 重拉标签上的靶对 `summary_min` 仅适度变化（0.417→0.424、0.606→0.629、0.500→0.500 和 0.692→0.692；Table S29）。

测量频次仍呈类别不平衡（Table S35）。PIK3CA/mTOR 上 dual 配体两端 API 记录数中位数为 22，A-only、B-only 与 neither 分别为 3、3 和 2。EGFR/HER2 有一个 dual 配体含 318 条记录。这些计数不证明测量频次导致了观察得的 AUROC，但说明 max-versus-median 标签稳定并不能消除类别间的 profiling 强度差异。

post-hoc 当前 ChEMBL 高置信视图审计了 2748 条 activity records，在显式靶标、confidence、relation、endpoint、validity 与 duplicate 过滤后保留 1546 条；513 条因 potential duplicate 排除（Table S36）。全部 352 个已打分配体在两端均至少保留一条合格记录，且 352 个高置信 θ = 6.0 状态与冻结状态一致，因此四条方向性点估计不变。该结果针对当前数据库快照上的指定记录过滤，而不是 assay 条件等价或进入双端测定子集的选择机制。

覆盖审计量化了这一选择（Table S37）。在至少一端有可用 pChEMBL 的结构中，EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 与 PIK3CA/mTOR 两端均有值的比例分别为 14.5%、34.0%、23.3% 和 26.5%。缺失不解释为无活性。来源文献集中度也随类别变化。最突出的是 PIK3CA/mTOR 的 4 个 neither 配体及其 8 条保留记录全部来自同一篇 ChEMBL 文献。

在同一套冻结 Vina 分数上，按文献阻断交叉验证后，EGFR/HER2 的弱方向臂仍为 0.430（document-cluster bootstrap 95% CI [0.321, 0.617]；5 个有效折；23 个组；113 篇文献；Table S39）。Dual versus A-only 仍为 0.666。同一折上 ECFP4 logistic 的 OOF AUROC 为 0.623，低于支架分组的 0.89，说明同篇文献系列贡献了二维信号。八个方向臂中七个可估计；PIK3CA/mTOR Dual versus B-only 在 9 个文献连通组中只有 1 个折同时含两类（最大组 19/30），按原规则报告为无法稳定估计（Table S40）。

预先冻结的 2018 时间分割在任何靶对上都未达到样本量门槛（测试集 dual/A-only/B-only/neither：EGFR/HER2 6/3/14/2；AChE/BChE 8/5/15/6；PIK3CA/PIK3CB 12/11/0/3；PIK3CA/mTOR 2/0/1/0；Table S41）。2015 敏感性中仅 AChE/BChE 三类均 ≥10（11/11/24）。主截止年可评估靶对少于两个，因此不包装为外部验证。BindingDB REST 两端等式测定在去掉已打分面板 InChIKey 后，四对的 dual/A-only/B-only 仍均 ≥10（EGFR/HER2 1589/161/62；AChE/BChE 966/450/230；PIK3CA/PIK3CB 1341/261/257；PIK3CA/mTOR 2008/238/266；Table S43）。Table S43 仍是该历史 REST 供给计数，不是文献独立切片。随后从 BindingDB 202608 原生归档重建、并经文献/结构/ECFP4 过滤的结果见 Results 3.9（Tables S48–S49）。

机器 assay-context 审计标记了 186 个优先配体（1163 条保留 activity 行），包括 98 个 EGFR/HER2 方向类配体、40 个混合端点配体、22 个同时有生化与功能实验的配体，以及 4 个 PIK3CA/mTOR neither 配体（Table S42）。元数据审核随后给出 179 include / 7 uncertain / 0 exclude；7 个 uncertain 配体保留冻结类别。蛋白构建体与突变状态因 assay 自由文本不可用仍为 unknown。冻结标签未改。

将 exhaustiveness 从 16 降至 8 后，PIK3CA/mTOR 的 summary_min 从 0.692 降至 0.660（Figure S1D）。在 PM110 面板中（n_scored = 115；dual / A_only / B_only 各 30），Vina summary_min 为 0.648 [0.51, 0.76]，相比 PM48 下降约 0.04，排序趋势保持一致（Figure S1C）。在未使用配体池留出集中（每对 20 / 20 / 20，种子 20260731；EGFR/HER2 不具备同等配额），PIK3CA/mTOR 的 summary_min 为 0.765 [0.603, 0.891]，AChE/BChE 为 0.618 [0.422, 0.759]，PIK3CA/PIK3CB 则下降至 0.425 [0.241, 0.618]（Tables S8、S16）。该留出集共享同一 ChEMBL 抓取批次。

我们进一步测试受体实现：一端受体冻结，只替换另一端（Figure 5；Tables S9、S30）。三个替代晶体均通过共晶重对接 QC（best-of-9 RMSD 分别为 4JPS 0.607 Å、5DXT 0.624 Å、4JSX 0.515 Å）。在 PIK3CA/mTOR 上，当 PIK3CA 4L23 替换为 4JPS 或 5DXT、mTOR 4JT6 保持不变时，PM48 的 summary_min 分别由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]（Figure 5A）。变化主要发生在依赖替代 PIK3CA 结构的 D/B 方向，D/A 保持 0.714。将 mTOR 4JT6 替换为 4JSX 后 summary_min 为 0.639 [0.418, 0.776]。

同一套 PIK3CA 晶体再用于 PIK3CA/PIK3CB 面板、2WXF 保持冻结时，summary_min **上升**：由 0.500 至 0.691 [0.516, 0.779]（4JPS）和 0.685 [0.506, 0.768]（5DXT）（Figure 5B）。dual versus A-only 保持 0.691；dual versus B-only 由 0.500 升至 0.707 和 0.685。三种 PIK3CA 条件使用相同的 99 个已打分配体（Table S30）。因此，替换同一 PIK3CA 受体在两个靶对上使表观判别向相反方向变化。

### 3.5 错口袋对照揭示面板外尚未解决的失败模式

在主面板中，口袋匹配 summary_min 均高于错口袋对照，四对的 matched-minus-wrong 差值分别为 0.170、0.161、0.151 和 0.090；其中 EGFR/HER2 和 AChE/BChE 的差异置信区间排除 0，PIK3CA/PIK3CB 与 PIK3CA/mTOR 的区间包含 0（Tables S6、S17；Figure 6A；Figure S3A）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602。

未使用配体池留出集中的点估计关系发生反转（Figure 6B）。PIK3CA/mTOR、AChE/BChE 和 PIK3CA/PIK3CB 的错口袋 summary_min 分别为 0.788、0.643 和 0.520，而匹配口袋分别为 0.765、0.618 和 0.425。相应的 matched-minus-wrong 点差值均为负（−0.023 / −0.025 / −0.095），但其 95% 置信区间均包含 0（Table S17；Figure S3B）。效价或尺寸匹配后，错口袋仍 ≥ 匹配口袋（Table S13）。不依赖打分的接触计数在 B 方向获得 0.698–0.714 的 AUROC，但这一粗粒度替代不能解释观察到的 Vina 错口袋判别幅度（Figure 6D；Table S11）。因此，留出集中的错口袋反转是当前基准暴露出的尚未解决的失败模式。

### 3.6 探索性结构背景

结构相似性并不能预测筛选判别：5DXT 与 4L23 的局部口袋 Cα RMSD 仅 0.343 Å，但二者的基准表现仍有明显差异（Table S10）。在 PM48 姿态上，以重原子距离 ≤ 4.5 Å 对 20 个冻结 PIK3CA 口袋残基做几何占有率快照，占有率变化最大的残基为 Met772、Leu807、Gln859、Thr856、Cys838、Glu849、Phe930 和 Asp933（Table S33）。表观判别位移与这些接触模式变化同时出现，只提供受体敏感性的结构假说。占有率变化不能证明某一残基导致了 AUROC 改变，也不能解释 PIK3CA/PIK3CB 上的相反位移。姿态类型学保留在 Supporting Information（Table S7；Note S1）。

### 3.7 独立姿态生成仍保留 EGFR/HER2 的设定差距

Results 3.2 中的 GNINA CNN 与 RTMScore 是对冻结 Vina 姿态的重打分。为检验设定效应是否依赖于该姿态生成引擎，在同一套冻结 EGFR/HER2 与 PIK3CA/mTOR 配体、受体、对接盒和 exhaustiveness 上以 GNINA 1.3.2 对接搜索模式独立生成姿态（Methods 2.3；Table S32）。分析仅使用两端均有分数的配体（EGFR/HER2 在 neither 配体 EH120_109 两端失败后 n = 109；PIK3CA/mTOR 在 A-only 配体 PM48_19 两端失败后 n = 47）。方向性 n_scored 在 EGFR/HER2 上仍为 28 / 38 / 32，在 PIK3CA/mTOR 上变为 18 / 13 / 12。

EGFR/HER2 上 Dual versus neither 仍高（AUROC 0.783 [0.610, 0.922]；n_neg = 11），而方向性 `summary_min` 为 0.220 [0.109, 0.343]。按口袋均分取混合库 Top-10 仍含 1 个 dual 与 9 个实验选择性配体（4 个 A-only、5 个 B-only；EF10 = 0.389）。因此，在该独立姿态生成协议下，设定效应不是 Vina 特有的。PIK3CA/mTOR 的方向性 `summary_min` 符号不变（0.633 对 Vina 0.692）。该比较只检验设定差距是否仍在，不是主张 GNINA 优于或劣于 Vina。

### 3.8 标签供给规模、物化匹配与 AND 过滤

θ = 6.0 普查（Results 3.1；Table S44）表明四状态实验标签并不只存在于已对接的四对，但不表明这些额外靶对已被对接评价。

在冻结已打分面板上，1.0 SD caliper 的 1:1 物化匹配在有足够样本的 Dual-versus-B-only 格子上仍接近随机：EGFR/HER2 0.566 [0.356, 0.781]（n = 16），AChE/BChE 0.462 [0.243, 0.692]（n = 13），PIK3CA/PIK3CB 0.556 [0.284, 0.827]（n = 9）。0.5 SD caliper 以及 PIK3CA/mTOR 1.0 SD 的 Dual-versus-B-only（n = 5）效能不足（Table S45）。该匹配是用实验选择性配体做的 DUDE-Z 式对照，不是因果调整。

AND 式双口袋过滤在 Dual+A-only+B-only 库上把 Table S25 的 Top-10 观察写成明确工作点（Table S46；Figure S7B）。EGFR/HER2 上，Dual 中位 `vina_worst` 截断保留 14/28 个 Dual，但同时留下 33 个选择性配体（precision 0.298；硬负比例 0.702）。收到 Dual 第 90 百分位后 Dual precision 进一步降到 0.130，因为极端 AND 分数由 B-only 主导。AChE/BChE 在该尾部 precision 为 0.600（n_pass = 5）。这些截断是对虚拟筛选和生成式设计中双口袋过滤器的诊断，不是对某一生成模型的重打分。

在完整 ChEMBL 图而不是 n ≈ 28 对接面板上，配体层模型在每一对上都比 Dual versus 选择性更容易回收 Dual versus neither（Table S47；Figure S7C）。ECFP4 GroupKFold 的 Dual-versus-neither AUROC 在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 与 PIK3CA/mTOR 上分别为 0.921、0.851、0.835 与 0.956，而方向性 ECFP4 `summary_min` 为 0.801、0.744、0.720 与 0.887。EGFR/HER2 Dual versus B-only 在二维化学上仍为 0.864，而对接下采样面板为 0.430。因此实验标签在全图尺度上是化学可分的；这并不证明对接在硬臂上回收了该结构。Table S47 不替换 Table 2。

### 3.9 BindingDB 原生切片：供给冻结，不是外部验证

BindingDB 202608 文章与专利归档按对接前冻结的合约重建（`external_slice_contract.yaml`；Methods 2.1）。[16] 原生配对 θ = 6.0 的 dual/A-only/B-only/neither 为 EGFR/HER2 371/30/54/13（n = 468）、AChE/BChE 159/46/37/76（n = 318）、PIK3CA/PIK3CB 114/29/4/149（n = 296）、PIK3CA/mTOR 1000/115/30/48（n = 1193）、MCL1/Bcl-xL 32/22/6/30（n = 90；Table S48）。去掉共享文献、共享结构并要求对开发分子最大 ECFP4 Tanimoto < 0.70 后，剩余为 180/10/20/6、4/8/14/59、9/0/3/100、91/4/1/2 与 1/0/2/0（Table S49；Figure S8）。ChEMBL 文献解析为 519/680，因此这些剩余 n 是完全文献独立集的上界。即便作为上界，也没有任何一对达到预先冻结的主外部门槛（方向类各 n ≥ 20、每类至少 3 个来源、最大单文献份额 ≤ 50%），也没有一对达到 EGFR 式薄复制门槛。0.50–0.70 Tanimoto 敏感性仍不足（AChE/BChE 39/13/20；EGFR/HER2 A-only 11；PIK3CA/mTOR B-only 6）。因此没有一对达到预先冻结主外部门槛的配体被对接，该切片不包装为外部评价。停止规则使本稿保持为四靶对评价设定审计。

MCL1/Bcl-xL 仍是同源 BCL-2 折叠上的 PPI/BH3 槽候选，不是异质折叠对，也不是首次非激酶对（Table S50）。LC6 pose-gold gate 在 3WIZ 失败（best-of-top3 2.011 Å；3WIY 1.689 Å 通过；Table S51）后，冻结的 24/24/24/24 面板只作为预先声明的 applicability stress-test 对接。93/96 个配体获得两端分数（186/192 个作业；3 个嵌入失败）。Vina Dual versus neither 为 0.628 [0.462, 0.786]，Dual versus A-only 为 0.793 [0.655, 0.915]，Dual versus B-only 为 0.609 [0.439, 0.776]，`summary_min` 为 0.609（Table S53）。[17] Dual-versus-B-only 区间包含 0.5。该对不加入 Table 2，也不包装为筛选性能证据或外部验证。Table S52 把 DualFourClass 与 Zhou 2013、DUD-E、LIT-PCBA、CASF-2016、DOCKSTRING 并置为文献对照，不是 bake-off。[5–9,18]

## 4. 讨论

### 4.1 基准设定改变了双靶对接的证据标准

该基准将双靶识别具体化为两个方向上的实验选择性硬负判别任务。同一套冻结 EGFR/HER2 分数上，Dual versus neither 的 AUROC 为 0.756，而方向性 `summary_min` 为 0.430（Results 3.2；Table 3）。因此，相对 inactive 的设定可以给出比方向性硬负任务更有利的表观证据。在独立 GNINA 姿态生成协议下，Dual versus neither 仍为 0.783，方向性 `summary_min` 为 0.220，因此此处的设定效应不是 Vina 特有的（Results 3.7；Table S32）。

相对于 Zhou 等的双靶评价设定，[9] 本文比较的是对接能否把 dual-active 与实验定义的选择性配体分开，而不仅是与 inactive 分开。已有对接基准表明，decoy 构建、化学偏倚和真实 assay 标签会改变虚拟筛选解释；[5–7,12,13] 同一关切适用于双靶任务：评价结论取决于负类如何被实验定义。物化匹配的实验选择性配体（Table S45）与 Dual 百分位 AND 过滤（Table S46）是这条 decoy-bias 课的双靶版本：EGFR/HER2 上要求两端都像 Dual，优先留下的是 B-only。

四状态标签问题并不限于已对接面板。冻结 49 对名单中有 17 对在 θ = 6.0 下 dual/A-only/B-only 均 n ≥ 10（Table S44）。四对完整图上的配体层 ECFP4 仍显示 Dual versus neither 比 Dual versus 选择性更容易（Table S47）。这是化学标签结果，不能据此写成 17 对对接评价。

### 4.2 表观对接信号不能直接归因于对接分数

物化描述符与化学型已经携带了相当一部分实验标签信息。AChE/BChE 上 TPSA 单独即可获得高于对接的判别；加入重原子数和 TPSA 后 dual-versus-B-only AUROC 上升，而对接分数的优势比仍接近 1（Results 3.3）。支架分组 ECFP4 在多个方向上超过对接，例如 EGFR/HER2 dual-versus-B-only 中 0.89 对 0.43，说明标签包含可不依赖受体信息使用的配体结构信息。把口袋匹配对接分数加到 ECFP4 后，CV AUROC 的最大绝对变化仅为 0.020。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only）由 0.691 降至 0.503，而远缘硬负（T < 0.3）升至 0.819。若缺少这些配体层对照，一个看似优秀的双靶对接结果可能只是识别了与 dual 标签相关的分子属性。[7,12]

### 4.3 受体实现构成评价条件的另一重要维度，并可提高或降低表观判别

表观判别还取决于评价条件如何指定。将最大 pChEMBL 换成重复测定中位数后，靶对层结论基本不变（Results 3.4）。同一 ChEMBL 抓取批次中的未使用配体池留出集保持了部分排序趋势，也改变了另一部分。一端受体冻结、只替换另一端时，同一套 PIK3CA 替换在一个相关靶对上提高表观判别、在另一个靶对上降低表观判别（Figure 5）。这些受体实现变化与把受体表示视为性能变量的激酶交叉对接工作相一致。[14] PM48 的 PIK3CA 姿态几何占有率快照与包括 Met772、Leu807 在内的接触模式变化同时出现（Table S33）。该模式只是结构假说，不能证明某一残基导致了 AUROC 改变，也不能解释 PIK3CA/PIK3CB 上的相反位移。

### 4.4 对双靶虚拟筛选与生成式设计的含义

这些结果对双靶点虚拟筛选具有直接的方法学启示，并可能适用于将 docking 作为下游筛选环节的生成式设计流程。同时在两个口袋获得有利分数，并不能自动建立实验定义的双靶活性。EGFR/HER2 上以 Dual 中位 `vina_worst` 做 AND 截断时 Dual precision 仅为 0.298（Table S46）。这一关切与近期 JCIM 研究相一致：在有实验依据的筛选集上，docking 再打分表现可以有很大差异。[15] 在本文的双靶任务中，对接分数因此需要连同实验定义的选择性硬负样本和配体化学对照一起解读。

同一组四个检查可以作为一个实用诊断流程使用（Figure 8）。在双口袋分数看起来有利之后：（i）要求相对 A-only 与 B-only 硬负的方向性判别；（ii）询问配体层 ECFP 或物化模型在抗泄漏分割下是否回收相似信号；（iii）检验未使用配体池，并在有文献标识时做 document-blocked 分割；（iv）至少替换一端受体实现。任一步失败，就把该主张标为依赖于设定、化学、面板或受体的计算证据。

### 4.5 局限性

第一，评价集仅含四对靶标，因为严格建造门槛下实验定义的双靶硬负样本稀缺。主规则 θ = 6.0 下，49 对审计中有 17 对 Dual/A-only/B-only 均 n ≥ 10，但这些额外靶对未加入 Table 2。MCL1/Bcl-xL 仅在 LC6 pose-gold gate 失败后作为 applicability stress-test 对接，不是第五个评价靶对。K = 4 是受数据供给约束的对接案例面板，而不是全面的双靶基准套件。四个 `summary_min` 还混合了面板构建差异（严格 6.5/5.5 对 θ = 6.0；不等 n）与靶对生物学。当前类别样本量更容易分辨较大的方向性效应，而对中等效应较弱（Table S31）。

第二，实验标签来自 ChEMBL，并要求两端均有可用测定。完整病例只覆盖可用值并集的 14.5%–34.0%。未使用配体池留出集仍属同一抓取批次，因此不是独立外部验证。BindingDB/PubChem 的 Table S12 核对仅为计数。Table S43 仍是去掉已打分面板 InChIKey 后的历史 REST 供给计数。随后 BindingDB 202608 原生归档在文献、结构与 ECFP4 < 0.70 过滤后，没有一对达到预先冻结的主外部门槛（Tables S48–S49）。因 ChEMBL 文献解析不完整，剩余 n 是上界；即便作为上界，该切片也不包装为外部评价，且未对接。按文献阻断后 EGFR/HER2 弱臂仍为 0.430，且 PIK3CA/mTOR Dual versus B-only 无法稳定估计。预先冻结的 2018 文献年份分割没有两个可评估靶对，故不声称时间外验证。

第三，assay 异质性仍然存在。IC50、Ki、Kd、EC50 与 Potency 被映射到同一阈值。主策展使用最大 pChEMBL。2026-08-26 的当前 ChEMBL 高置信重建（人源单蛋白、confidence≥8、等式关系、允许端点、validity 与 duplicate 过滤）保留了全部 352 个已打分标签（Table S36），但不能等同于 assay 条件、蛋白构建体或突变背景的统一。186 个优先分子的元数据纳入/排除未改变任何冻结类别（179 include / 7 uncertain / 0 exclude），构建体与突变仍为 unknown。因此标签不应视为 assay-harmonized ground truth。

第四，受体替换可以提高或降低成对判别，但实验并未给出分子起源。两个受体敏感性例子均共享 PIK3CA。共晶 best-of-nine 只证明搜索覆盖。重建的 EGFR 3POZ QC 是显式例子：九个姿态中存在近晶构象（0.760 Å），但未排进 top-3（top-1 9.505 Å）。HER2 3RCD 重建 QC 通过 top-1（1.855 Å）。

第五，主协议为 AutoDock Vina；GNINA CNN 与 RTMScore 是对同一组 Vina 姿态的重打分。EGFR/HER2 与 PIK3CA/mTOR 上的独立 GNINA 对接搜索仍保留主要设定差距，不是多引擎比赛。本研究未对新预测双靶化合物做前瞻实验。

## 5. 结论

在所评价的四个案例靶对中，方向性对接判别的点估计依赖靶对、化学组成、面板成员与受体实现；全部四个主 `summary_min` 区间均包含 0.5。EGFR/HER2 的负类设定差距在独立 GNINA 姿态生成下仍然存在，但其他靶对没有支持同一普遍规律。

这些结果识别的是当前 ChEMBL 衍生面板与计算协议中的失败模式，不是靶标通用的可靠性边界。仅依据两个口袋中的有利 docking 分数不足以建立双靶活性的充分证据。对于双靶虚拟筛选及将 docking 用作下游筛选环节的生成式设计流程，可使用四步诊断：方向性硬负、配体层化学基线、未使用配体池或文献阻断分割，以及受体结构敏感性（Figure 8）。θ = 6.0 标签普查、物化 caliper 匹配与 AND 过滤工作点加强了该诊断，但并不把本研究变成 17 对对接基准。在这些面板上预先冻结的文献年份分割不能作为时间外验证。BindingDB 原生独立切片没有给出两对主外部靶对，故未对接。

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源表，包括独立 GNINA 姿态生成分数（Table S32）、PIK3CA 占有率位移（Table S33）、文献阻断交叉验证（Tables S39–S40）、冻结的文献年份分割（Table S41）、assay-context 元数据审核（Table S42）、BindingDB REST 独立性计数（Table S43）、θ = 6.0 靶对普查（Table S44）、物化 caliper 匹配（Table S45）、AND 过滤工作点（Table S46）、配体层全图 AUROC（Table S47）、BindingDB 原生候选流程与切片摘要（Tables S48–S49）、MCL1/Bcl-xL 面板、LC6 gate 与 applicability-stress-test AUROC（Tables S50–S51、S53）、文献对照（Table S52）、重建的 EGFR/HER2 共晶 QC（Table S3）、Figure S8、原生切片合约（`protocol/external_slice_contract.yaml`），以及评价合约（`DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`）。面向稿件的表 SHA-256 校验和见 `REVISION_CHECKSUM_MANIFEST_v1.csv`。ChEMBL 供给审计冻结于 2026-07-23；高置信 activity 视图抓取于 2026-08-26；BindingDB 原生归档锁定为 202608。GitHub Release 与 Zenodo DOI 将从打标签快照签发，而不是从当前仍可能变化的分支签发。分析环境与零新对接的复现命令见仓库 README。BindingDB TSV 归档本身不随仓库分发；CI 只核对已提交的 CSV。

## 参考文献

(1) Anighoro, A.; Bajorath, J.; Rastelli, G. Polypharmacology: Challenges and Opportunities in Drug Discovery. *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: 10.1021/jm5006463.

(2) Proschak, E.; Stark, H.; Merk, D. Polypharmacology by Design: A Medicinal Chemist’s Perspective on Multitargeting Compounds. *J. Med. Chem.* **2019**, *62*, 420–444. DOI: 10.1021/acs.jmedchem.8b00760.

(3) Kitchen, D. B.; Decornez, H.; Furr, J. R.; Bajorath, J. Docking and Scoring in Virtual Screening for Drug Discovery: Methods and Applications. *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: 10.1038/nrd1549.

(4) Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: 10.1021/acs.jcim.1c00203.

(5) Huang, N.; Shoichet, B. K.; Irwin, J. J. Benchmarking Sets for Molecular Docking. *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: 10.1021/jm0608356.

(6) Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: 10.1021/jm300687e.

(7) Tran-Nguyen, V.-K.; Jacquemard, C.; Rognan, D. LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening. *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: 10.1021/acs.jcim.0c00155.

(8) Su, M.; Yang, Q.; Du, Y.; Feng, G.; Liu, Z.; Li, Y.; Wang, R. Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: 10.1021/acs.jcim.8b00545.

(9) Zhou, S.; Li, Y.; Hou, T. Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors. *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: 10.1021/ci400065e.

(10) Zhou, X.; Guan, J.; Zhang, Y.; Peng, X.; Wang, L.; Ma, J. Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design. In *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; 2024. arXiv:2410.20688.

(11) Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. In *Proceedings of the 32nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining, Vol. 2*; ACM: New York, 2026; pp 12432–12443. DOI: 10.1145/3770855.3819050.

(12) Tran-Nguyen, V.-K.; Ballester, P. J. Beware of Simple Methods for Structure-Based Virtual Screening: The Critical Importance of Broader Comparisons. *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405. DOI: 10.1021/acs.jcim.3c00218.

(13) Ahmed, F.; Soellner, M. B.; Brooks, C. L., III. Real-World Assessment of Machine-Learned Docking Using Bioassay-Derived Benchmarks. *J. Chem. Inf. Model.* **2026**, *66*, 8752–8759. DOI: 10.1021/acs.jcim.5c03020.

(14) Schaller, D. A.; Christ, C. D.; Chodera, J. D.; Volkamer, A. Benchmarking Cross-Docking Strategies in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 8848–8858. DOI: 10.1021/acs.jcim.4c00905.

(15) Sindt, F.; Bret, G.; Rognan, D. On the Difficulty to Rescore Hits from Ultralarge Docking Screens. *J. Chem. Inf. Model.* **2025**, *65*, 5553–5566. DOI: 10.1021/acs.jcim.5c00730.

(16) Liu, T.; Hwang, L.; Burley, S. K.; Nitsche, C. I.; Southan, C.; Walters, W. P.; Gilson, M. K. BindingDB in 2024: a FAIR Knowledgebase of Protein-Small Molecule Binding Data. *Nucleic Acids Res.* **2025**, *53*, D1633–D1644. DOI: 10.1093/nar/gkae1075.

(17) Tanaka, Y.; Aikawa, K.; Nishida, G.; Homma, M.; Sogabe, S.; Igaki, S.; Hayano, Y.; Sameshima, T.; Miyahisa, I.; Kawamoto, T.; Tawada, M.; Imai, Y.; Inazuka, M.; Cho, N.; Imaeda, Y.; Ishikawa, T. Discovery of Potent Mcl-1/Bcl-xL Dual Inhibitors by Using a Hybridization Strategy Based on Structural Analysis of Target Proteins. *J. Med. Chem.* **2013**, *56*, 9635–9645. DOI: 10.1021/jm401170c.

(18) García-Ortegón, M.; Simm, G. N. C.; Tripp, A. J.; Hernández-Lobato, J. M.; Bender, A.; Bacallado, S. DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design. *J. Chem. Inf. Model.* **2022**, *62*, 3486–3502. DOI: 10.1021/acs.jcim.1c01334.
