# Methods（中文工作稿 · JCIM Articles）

## 2. 方法

### 2.1 数据与实验状态定义

双靶评价所需的配体活性通过 ChEMBL Web API 的公开 activity 端点获取，作为实验衍生标签。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干经标准化的定量效力或亲和力测量（如 IC50、EC50、Ki、Kd 和 Potency）转换为近似 −log10 活性尺度。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似。

同一配体–靶标若有多条可用 pChEMBL 记录，主策展采用**最大 pChEMBL** 作为一对一代表值。活性聚合敏感性分析从 ChEMBL activity 端点重拉 assay 级记录，并在同一 θ = 6.0 规则下用重复测定的**中位数**替换最大值，不改变面板成员、对接参数或 Vina 分数。API 重拉的最大值对中位数估计作为标签聚合敏感性，与 Table 2 并列报告（Table S29）。作为 post-hoc 标签稳健性，已打分配体–靶标记录于 2026-08-26 从当前 ChEMBL API 重拉，仅保留 Homo sapiens `SINGLE PROTEIN`、assay confidence ≥8、等式关系、IC50/Ki/Kd/EC50/Potency、无 validity 注释且无 `potential_duplicate` 的记录，再按 θ = 6.0 重标而不改变面板成员或对接分数（Table S36）。该日期化视图检验这些显式过滤器是否改变冻结标签，不是 2026-07-23 数据库状态的重建，也不统一 assay 条件、蛋白构建体或突变背景。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。ChEMBL 结构按连通片段拆分，保留重原子数最多的有机片段。

完整病例选择按冻结可用 pChEMBL 映射统计只在 A、只在 B 或两端均有值的结构（Table S37）。缺失不解释为无活性。来源文献集中度按高置信保留记录、分实验状态统计独立 `document_id` 与最大单篇文献份额。

对每一对靶标 A/B，配体被定义为四种实验状态：**dual**（两端较强）、**A-only**（仅 A 端较强）、**B-only**（仅 B 端较强）和 **neither**（两端均不足）。A-only 与 B-only 是选择性硬负样本。

**严格 6.5/5.5 规则**为：dual，两端 pChEMBL ≥ 6.5；A-only，A ≥ 6.5 且 B ≤ 5.5；B-only 对称；neither，两端 ≤ 5.5。5.5–6.5 灰区不进入该审计。金属依赖体系（如 HDAC）预先排除。严格 6.5/5.5 规则仅用于靶对供给资格审定。全部主基准标签随后统一按 θ = 6.0 定义（dual，两端 ≥ θ；A-only，A ≥ θ 且 B < θ；B-only 对称；neither，两端 < θ），并在查看对接结果之前冻结。两条阈值因此服务于预先规定的不同目的。建造规则在抽样前按供给审计冻结（Table 1）。θ ∈ {5.5, 6.5} 与严格 6.5/5.5 重标作为敏感性报告（Table S4）。样本量过小的格子在 Results 中标记效能不足。

为核对 ChEMBL 供给门槛，对冻结靶对另做 BindingDB / PubChem 计数核对（零对接、不重建面板；Table S12）。类型限于 IC50/Ki/Kd/EC50；配体身份分别用 BindingDB monomerid 与 PubChem CID，不做跨库结构合并。主比较采用等式测定。

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

正式对接前对每个冻结受体做共晶配体重对接。生成 9 个姿态，计算与实验共晶构象的重原子 RMSD。预先通过标准为 \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å，即九个保留姿态中是否存在与共晶配体重原子 RMSD 小于 2.0 Å 的构象。若默认 exhaustiveness 未通过门槛，则提高至预先规定的备用水平。主分析因此采用 PIK3CA/mTOR exhaustiveness = 16、其余主面板为 8（Table S3）。

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

对 352 个已打分配体中的 186 个优先分子提取了 assay-context 字段，包括 EGFR/HER2 全部方向类、PIK3CA/mTOR 的 4 个 neither、混合端点记录，以及对主 AUROC 影响最大的分子（`assay_context_audit.csv`）。蛋白构建体与突变状态需要阅读原文，本机器提取不改写冻结标签。

受体结构敏感性分析另选满足以下预先声明条件的替代晶体：（i）polymer entity 与目标蛋白真实对应；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.3 相同的共晶重对接 QC。实际对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。替换采用单口袋设计：在 PIK3CA/mTOR（PM48）上，4JPS/5DXT 替换口袋 A、口袋 B 仍用冻结 4JT6 分数，4JSX 替换口袋 B、口袋 A 仍用冻结 4L23 分数（exhaustiveness = 16）。在 PIK3CA/PIK3CB 上，同一套 4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 2WXF 分数（exhaustiveness = 8）。刚体 Cα 叠合作为探索性几何对照（Table S10）。在 Table S30 所用同一套 PM48 配体与 PIK3CA 晶体上，另做探索性接触快照：占有率定义为与 20 个冻结口袋残基的重原子距离 ≤ 4.5 Å（Table S33）。占有率变化只作为结构假说，不是残基层因果解释。

计算在 Python 3 环境下完成，主要软件为 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore。评价面板、对接分数、分析脚本与参数表见 Data and Software Availability。
