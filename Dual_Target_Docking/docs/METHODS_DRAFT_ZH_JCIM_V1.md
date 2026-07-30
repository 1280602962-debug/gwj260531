# Methods（中文工作稿 · JCIM Articles）

> 结构与语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843  
> （Dataset → Docking → Scoring → Evaluation；正文写清关键参数，其余见 SI Table S1）  
> 与 [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。投稿以英文为准；本稿供中文审改。

---

## 2. Methods

### 2.1 任务定义

本文将双靶对接评测建成四类配体判别任务。对每一对靶标 A/B，依据两端实测活性将配体分为：dual（两端均强）、A_only（仅 A 强）、B_only（仅 B 强）和 neither（两端均弱）。评价目标不是比较某一口袋上谁的对接分更高，而是判断对接分数能否把 dual 配体同时与两类单靶选择性配体（A_only 与 B_only）区分开。

为此，我们分别计算 dual 对 A_only、dual 对 B_only 的 AUROC，并以二者的最小值作为该靶对的汇总指标（记为 summary_min）。采用最小值是为了避免只报告较好一侧：若仅有一端可区分而另一端接近随机，则不足以支持“双靶判别”的主张。

### 2.2 数据策展与靶对选择

活性数据来自 ChEMBL。pChEMBL 是 ChEMBL 提供的近似统一活性标度，将若干摩尔浓度–响应型测量（如 IC50、Ki、Kd、EC50 等）转换为 −log10 形式，以便在公开数据整合中比较。不同 assay 类型与实验条件并不等价；本文将 pChEMBL 视为大规模策展中的常用近似，而非同一条件下的绝对亲和力。同一配体–靶标若有多条记录，取可用的最大 pChEMBL 作为代表值（最强报道活性）；该选择与公开 compound–target 数据集中同时报告 mean / median / maximum 聚合的做法一致，但偏向最强报道，且未按 assay 中位数或 ChEMBL 置信度过滤。任一端缺少可用 pChEMBL 的配体不进入四类主分析。上述局限在 Limitations 中进一步讨论。ChEMBL 结构条目常含盐形式；对接前将结构去盐并保留最大有机片段，以得到单一可对接分子。

四类标签首先按**严格规则**定义，以便在公开数据中识别真正“一端强、一端弱”的单靶选择性配体：dual 为两端 pChEMBL ≥ 6.5；A_only 为 A ≥ 6.5 且 B ≤ 5.5；B_only 为对称定义；neither 为两端 ≤ 5.5。介于 5.5 与 6.5 之间的灰区配体不进入严格面板。我们以该规则为供给审计标准：在 49 对可审计靶对上统计 A_only 与 B_only 的可用配体数。两端均 ≥ 50 的仅 4 对；去掉金属依赖、不适合作为常规对接主对象的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB。这三对构成评价集的主体。

EGFR/HER2 虽在同一严格规则下 B_only 仅 7 个，无法建成规模均衡的严格四类面板，但仍被纳入评价集，原因有二。第一，该对在双靶文献中出现频率高，若完全排除，读者无法看到“热门靶对在严格四类设定下究竟缺什么”。第二，它提供一个供给受限的对照案例：后续对接结果应结合标签供给上限解读，而不能默认再加大对接预算即可补齐硬负。EGFR/HER2 因此沿用既有统一 RDKit 准备面板，作为案例对与上述三对并列报告，而不是作为严格规则下的合格厚面板。

在实际抽样建面板时，并非所有靶对都能在严格规则下凑齐配额。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则完成配额抽样；EGFR/HER2 与 PIK3CA/mTOR 主面板（PM48）因严格规则下单端选择性配体过少，改用**单阈值规则**建成：选定 θ = 6.0 后，两端均 ≥ θ 为 dual，一端 ≥ θ 且对端 < θ 为对应单靶选择性类，两端均 < θ 为 neither。换言之，严格规则用于供给审计与优先建面板；单阈值规则是供给不足时的建面板折中，不是与严格规则并列的第二套主标准。为检查标签定义对结论的影响，我们在既有面板上按 θ ∈ {5.5, 6.0, 6.5} 重标四类，并用两端 Vina 分数均值重算方向 AUROC（Supporting Information Table S4）；正文主指标仍为口袋匹配 summary_min（Table 2）。

### 2.3 评价面板的构建

对选定的每一对靶标，我们按四类标签从 ChEMBL 候选集中抽样组成评价面板，并限制同一类别内重复骨架的过多富集。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则以 dual / A_only / B_only / neither ≈ 28 / 28 / 28 / 16 定额建成（面板 n = 100）。EGFR/HER2 沿用既有 n = 110 面板。PIK3CA/mTOR 则先建成较小的试点面板（PM48，n = 48；θ = 6.0），用于冻结受体、对接与配体准备协议；该面板中 dual / A_only / B_only 仅 18 / 14 / 12，方向性评价功效有限。对接失败的配体–受体组合从相应受体分数中剔除。主方向性评价仅使用 dual、A_only 与 B_only；neither 保留在面板中供描述，但不进入 dual 对 A_only / B_only 的 AUROC 计算。各面板的组成与对接设置汇总于 Table 1。

**Table 1.** 评价集组成与对接设置

| 靶对 | 标签规则 | 受体 PDB (A / B) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|----------|------------------|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 110 | 28 / 38 / 32 | 8 |

由于 PM48 样本量偏小、summary_min 的置信区间较宽，我们在 PIK3CA/mTOR 上另构建扩面面板：保留 PM48 的全部 48 个配体，并按严格规则追加新分子（目标配额 30 / 30 / 30 / 25），实际 n = 115。该扩面面板是 PM48 的超集，用于提高功效并收窄置信区间，而不是独立重复实验。主文以 PM48 报告与其他三对可比的协议冻结结果，并以扩面面板检验样本量增大后结论是否保持同向。

### 2.4 蛋白结构与结合位点定义

受体结构取自 Protein Data Bank 中含小分子共晶配体的条目。PIK3CA 与 mTOR 分别使用 4L23 与 4JT6；AChE 与 BChE 使用 4EY7 与 4BDS；PIK3CB 使用 2WXF；EGFR 与 HER2 使用 3POZ 与 3RCD。对每个受体，结合位点定义为以共晶配体为中心的盒子：先取共晶配体的轴对齐包围盒，再向外扩展 5 Å，并将每边长度下限设为 20 Å。受体准备为对接所用的 PDBQT 格式。完整盒子坐标见 Supporting Information Table S2。

八个冻结受体均以共晶配体重对接作协议检查：门槛为重原子 RMSD 在输出的 9 个姿态中的最小值（best_of_9）&lt; 2.0 Å。在 exhaustiveness = 8 时，4L23、4EY7、4BDS、2WXF、3POZ 与 3RCD 的 best_of_9 均 &lt; 2 Å；mTOR（4JT6）上 PI-103 的 best_of_9 为 5.003 Å，未过门槛。将 exhaustiveness 提高到 16 后，4JT6 的 best_of_9 降至 0.445 Å（近晶姿态为 mode 3；Vina mode1 仍约 7.1 Å）。因此 PIK3CA/mTOR 主面板及其扩面、单靶对照均采用 exhaustiveness = 16；其余靶对采用 8，并在 PIK3CA/mTOR 上另行报告 exhaustiveness = 8 的对照结果。另需说明：EGFR（3POZ）在 E = 8 时 best_of_9 已过关（约 0.96 Å），但 Vina 排序第一的姿态 RMSD 约 9.5 Å；升高 exhaustiveness 不能把 mode1 翻成近晶。完整 cognate 表见 Supporting Information Table S3。

### 2.5 对接与打分

配体由 ChEMBL 结构去盐后保留最大有机片段，用 RDKit ETKDGv3 生成三维构象，再经 meeko 转为 PDBQT。对接使用 AutoDock Vina 1.2.7：每个配体保留 9 个姿态，`energy_range = 3`，exhaustiveness 按 Table 1（PIK3CA/mTOR 为 16，其余靶对为 8）。构象生成、面板抽样与对接均使用固定随机种子；完整参数见 Supporting Information Table S1。

在同一组 Vina 姿态上另用两种函数重打分，以检查结论是否依赖单一打分通道。RTMScore（公开权重 `rtmscore_model1`）对 9 个姿态取最高分。GNINA 1.3.2 在 CPU 模式下，将 Vina 排序第一的姿态转为 SDF 后，以 `--cnn_scoring rescore --minimize` 重打分。正文以 Vina 亲和力（取负号，使越高越好）报告口袋匹配方向 AUROC；RTMScore 与 GNINA 仅作通道对照，不据此更换主指标。

### 2.6 评价指标与基线

对每个配体，分别记录其在口袋 A 与口袋 B 上的分数（方向约定见 2.5）。

主指标为口袋匹配的方向 AUROC：比较 dual 与 A_only 时使用口袋 B 的分数；比较 dual 与 B_only 时使用口袋 A 的分数。summary_min 为这两条 AUROC 的最小值。两端分数的平均值（池化分数）仅作对照，不作为主指标。

作为平凡基线，我们用同一方向 AUROC 流程，但以配体描述符代替对接分数，包括重原子数、分子量、cLogP 与拓扑极性表面积（TPSA）。对接相对最优平凡基线的差值及其置信区间用于判断对接是否提供超出简单理化性质的信息。

不确定度用配体层 bootstrap 估计：对每个靶对重采样 2000 次，报告 summary_min 的 95% 百分位区间。按 Murcko 支架重采样的结果另作对照，主文以配体层区间为准。

### 2.7 对照与敏感性分析

为检查方向性信号是否可能来自分子属性而非口袋匹配，我们设置下列对照。错口袋对照将对调后的口袋分数用于同一方向比较。配体效率对照将分数除以重原子数后再计算口袋匹配 AUROC。匹配子集对照分别在 \|ΔpChEMBL\| ≤ 0.5 或 \|Δheavy atoms\| ≤ 2 的近邻匹配子集上重算单侧 AUROC。协变量对照在逻辑回归中同时纳入对接分数、重原子数与 TPSA，比较仅分数模型与加入协变量后的判别 AUROC。此外，以 ECFP4 指纹（Morgan 半径 2，2048 bit）与逻辑回归建立仅依赖二维结构的配体基线；交叉验证按 Murcko 支架分组，避免同一骨架跨训练/测试折。

作为单靶虚拟筛选参照，我们在 4L23（PIK3CA）与 4JT6（mTOR）上分别构建活性–decoy 集合：活性分子 pChEMBL ≥ 6.5；decoy 为同靶已测定但 pChEMBL ≤ 5.5 的弱效分子，并按分子量（±50）、logP（±1.5）与 TPSA（±25）与活性分子属性匹配。目标规模约为 50 个活性分子与 150 个 decoy，对接设置与 PIK3CA/mTOR 主面板一致（Vina，exhaustiveness = 16）。该对照回答单靶富集是否存在，不替代四类方向评价。

标签阈值（Table S4）与匹配子集（Table S5）的敏感性分析见 Supporting Information。

### 2.8 软件与数据可用性

计算使用 Python 3，以及 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7 与 GNINA 1.3.2。RTMScore 使用公开预训练权重。评价面板、对接分数、分析脚本与完整参数表将随公开数据包发布；DOI 见 Data and Software Availability。

---

## 写法说明（不进正文）

对照 Vu et al. 2025 Methods 的取舍：

| Vu 怎么写 | 本稿怎么对齐 |
|-----------|--------------|
| Dataset 先给规模与规则，再列输入文件类型 | 2.2–2.3：ChEMBL 规则 → 选对 → 面板规模（Table 1） |
| Docking 里写清盒子、配体准备、引擎参数；细节进 Table S1 | 2.4 受体/盒子/cognate；2.5 合并配体准备 + Vina + 重打分 |
| Scoring 可并入对接节，或单独一小节 | 本稿并入 2.5 第二段（通道对照，不换主指标） |
| Evaluation 单独一节写指标 | 2.6 |
| 不写仓库路径、内部角色名、未做实验的 holdout | 已去掉 protocol.yaml、boxes/、主开发对、NLRP3 等；SI 只汇编已有数字 |
| 过去时 / 陈述句，少口号 | 全文按此收紧 |

先前稿问题主要来自把项目笔记（角色名、黑话、路径、决策口吻）直接写进 Methods，而不是按已发表评测文把“读者要复现实验所需的信息”写清楚。本版从 2.3 起按该标准重写；2.1–2.2 保留上一轮已改口径；2.4 已与真实 cognate 表对齐；2.5–2.7 原三分节过碎，已合并为 2.5。
