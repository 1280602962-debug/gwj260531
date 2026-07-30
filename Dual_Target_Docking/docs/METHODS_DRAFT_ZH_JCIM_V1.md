# Methods（中文工作稿 · JCIM Articles）

> 结构与语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843  
> （Dataset → Task → Docking → Scoring → Evaluation；正文写清关键参数，其余见 SI Table S1）  
> 与 [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。投稿以英文为准；本稿供中文审改。

---

## 2. Methods

### 2.1 数据策展与靶对选择

活性数据来自 ChEMBL Web API（公开 activity 端点；靶对供给审计锁定于 2026-07-23）。pChEMBL 是 ChEMBL 提供的近似统一活性标度，将若干摩尔浓度–响应型测量（如 IC50、Ki、Kd、EC50 等）转换为 −log10 形式，以便在公开数据整合中比较。不同 assay 类型与实验条件并不等价；本文将 pChEMBL 视为大规模策展中的常用近似，而非同一条件下的绝对亲和力。同一配体–靶标若有多条记录，取可用的最大 pChEMBL 作为代表值；未按 assay 中位数或 ChEMBL 置信度过滤。任一端缺少可用 pChEMBL 的配体不进入四类主分析。上述局限在 Limitations 中进一步讨论。ChEMBL 结构条目常含盐形式；对接前将结构按连通片段拆分并保留重原子数最多的有机片段，以得到单一可对接分子。

四类标签首先按**严格规则**定义，以便在公开数据中识别真正“一端强、一端弱”的单靶选择性配体：dual 为两端 pChEMBL ≥ 6.5；A_only 为 A ≥ 6.5 且 B ≤ 5.5；B_only 为对称定义；neither 为两端 ≤ 5.5。介于 5.5 与 6.5 之间的灰区配体不进入严格面板。我们以该规则为供给审计标准：在 49 对可审计靶对上统计 A_only 与 B_only 的可用配体数。两端均 ≥ 50 的仅 4 对；去掉金属依赖、不适合作为常规对接主对象的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB。这三对构成评价集的主体。

EGFR/HER2 在同一严格规则下 B_only 仅 7 个，无法建成规模均衡的严格四类面板，但仍纳入评价集：该对在双靶文献中常见，且可作为供给受限案例，使对接结果与标签供给上限一并解读。EGFR/HER2 因此沿用既有统一 RDKit 准备面板，作为案例对与上述三对并列报告，而不作为严格规则下的合格厚面板。

在实际抽样建面板时，并非所有靶对都能在严格规则下凑齐配额。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则完成配额抽样；EGFR/HER2 与 PIK3CA/mTOR 主面板（PM48）因严格规则下单端选择性配体过少，改用**单阈值规则**建成：选定 θ = 6.0 后，两端均 ≥ θ 为 dual，一端 ≥ θ 且对端 < θ 为对应单靶选择性类，两端均 < θ 为 neither。换言之，严格规则用于供给审计与优先建面板；单阈值规则是供给不足时的建面板折中，不是与严格规则并列的第二套主标准。为检查标签定义对结论的影响，我们在既有面板上按 θ ∈ {5.5, 6.0, 6.5} 重标四类，并用两端 Vina 分数均值重算方向 AUROC（Supporting Information Table S4）；正文主指标仍为口袋匹配 summary_min（Table 2）。

### 2.2 任务定义

本文将双靶对接评测建成四类配体判别任务。对每一对靶标 A/B，依据两端实测活性将配体分为：dual（两端均强）、A_only（仅 A 强）、B_only（仅 B 强）和 neither（两端均弱）。评价目标不是比较某一口袋上谁的对接分更高，而是判断对接分数能否把 dual 配体同时与两类单靶选择性配体（A_only 与 B_only）区分开。

为此，我们分别计算 dual 对 A_only、dual 对 B_only 的 AUROC，并以二者的最小值作为该靶对的汇总指标（记为 summary_min）。采用最小值是为了避免只报告较好一侧：若仅有一端可区分而另一端接近随机，则不足以支持“双靶判别”的主张。

### 2.3 评价面板的构建

对选定的每一对靶标，我们按 2.1 所述标签规则从 ChEMBL 候选集中抽样组成评价面板。PIK3CA/mTOR（PM48）在抽样时限制同一类别内同一 Murcko 支架最多 2 个分子；EGFR/HER2（n = 110）相应上限为 5；AChE/BChE 与 PIK3CA/PIK3CB 在建面时以类别配额控制组成，并以标识符前缀作弱多样性约束。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则定额抽取 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2 沿用既有 θ = 6.0 面板（n = 110）。PIK3CA/mTOR 因严格规则下单端选择性配体不足，按 θ = 6.0 建成 n = 48 的主比较面板（PM48；dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。对接失败的配体–受体组合从相应受体分数中剔除，故分析用类别计数可低于建面定额（Table 1）。主评价仅使用 dual、A_only 与 B_only 计算 AUROC；neither 保留在面板中供描述，但不进入上述对比。各面板组成与对接设置汇总于 Table 1。

**Table 1.** 评价集组成与对接设置

| 靶对 | 标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|----------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

由于 PM48 的 dual / A_only / B_only 样本量较小、summary_min 置信区间较宽，我们在 PIK3CA/mTOR 上另构建扩面面板：保留 PM48 的全部 48 个配体，并按严格规则追加新分子（目标配额 30 / 30 / 30 / 25），实际 n = 115。该扩面面板是 PM48 的超集，用于收窄区间并检验结论是否同向，而不是独立重复实验。主文以 PM48 报告与其他三对在同一协议下的比较；扩面结果另报。

### 2.4 蛋白结构与结合位点定义

受体结构取自 Protein Data Bank 中含小分子共晶配体的条目。PIK3CA 与 mTOR 分别使用 4L23 与 4JT6（共晶配体 X6K / PI-103）；AChE 与 BChE 使用 4EY7 与 4BDS（E20 / THA）；PIK3CB 使用 2WXF（039）；EGFR 与 HER2 使用 3POZ 与 3RCD（03P / TAK-285）。分辨率见表 1。对接用受体均为去除水分子与共晶配体后的蛋白坐标，并转为 PDBQT。PIK3CA、mTOR、EGFR 与 HER2 采用冻结目录中已含氢的蛋白坐标，经 meeko `mk_prepare_receptor.py --read_pdb` 生成 PDBQT；AChE、BChE 与 PIK3CB 由沉积 PDB 的 ATOM/TER 记录提取蛋白（去除水与异源原子），以 meeko `mk_prepare_receptor`（默认 alternate location A）生成 PDBQT。对每个受体，结合位点定义为以共晶配体为中心的盒子：先取共晶配体的轴对齐包围盒，再向外扩展 5 Å，并将每边长度下限设为 20 Å。完整盒子坐标见 Supporting Information Table S2。

八个冻结受体均以共晶配体重对接作协议检查。通过门槛为输出 9 个姿态中重原子 RMSD 的最小值（best_of_9）&lt; 2.0 Å；RMSD 在对接坐标系中计算，不做蛋白叠合。PIK3CA/mTOR 与 EGFR/HER2 的主 QC 采用 meeko `REMARK SMILES IDX` 原子映射，并在模板图自同构上取最小 CalcRMS；AChE/BChE 与 PIK3CB 冻结 QC 采用重原子坐标的匈牙利匹配。在 exhaustiveness = 8 时，4L23、4EY7、4BDS、2WXF、3POZ 与 3RCD 的 best_of_9 均 &lt; 2 Å；mTOR（4JT6）上 PI-103 的 best_of_9 为 5.003 Å，未过门槛。将 exhaustiveness 提高到 16 后，4JT6 的 best_of_9 降至 0.445 Å（近晶姿态为 mode 3；Vina mode1 仍约 7.1 Å）。因此 PIK3CA/mTOR 主面板及其扩面、单靶对照均采用 exhaustiveness = 16；其余靶对采用 8，并在 PIK3CA/mTOR 上另行报告 exhaustiveness = 8 的对照结果。另需说明：EGFR（3POZ）在 E = 8 时 best_of_9 已过关（约 0.96 Å），但 Vina 排序第一的姿态 RMSD 约 9.5 Å；升高 exhaustiveness 不能把 mode1 翻成近晶。完整 cognate 表见 Supporting Information Table S3。

### 2.5 对接与打分

配体由 ChEMBL SMILES 去盐（保留最大有机片段）后，用 RDKit 加氢，以 ETKDGv3 生成三维构象（随机种子 20260727），并用 MMFF 力场局部优化（最多 200 步），再经 meeko 默认参数转为 PDBQT。对接使用 AutoDock Vina 1.2.7，打分函数为默认 `vina`：每个配体保留 9 个姿态，`energy_range = 3`，exhaustiveness 按 Table 1（PIK3CA/mTOR 为 16，其余靶对为 8）。构象生成、面板抽样与对接均使用固定随机种子；完整参数见 Supporting Information Table S1。

在同一组 Vina 姿态上另用两种函数重打分，以检查结论是否依赖单一打分通道。RTMScore（公开权重 `rtmscore_model1`）对 9 个姿态取最高分。GNINA 1.3.2 在 CPU 模式下，将 Vina 排序第一的姿态经 Open Babel 转为 SDF 后，以 `--cnn_scoring rescore --minimize` 重打分。正文以 Vina 亲和力（取负号，使越高越好）按 2.6 所述计算 summary_min；RTMScore 与 GNINA 仅作通道对照，不据此更换主指标。

### 2.6 评价指标与基线

每个配体在口袋 A 与口袋 B 上各有一个对接分数。主评价计算两条二分类 AUROC：dual 对 A_only 使用口袋 B 的分数；dual 对 B_only 使用口袋 A 的分数。该设定对应口袋匹配的判别任务：单靶选择性配体在对端口袋上预期较弱，故以对端分数检验其是否可与 dual 区分。两条 AUROC 的较小值记为 summary_min，作为该靶对的汇总指标，以避免仅报告较好一侧。两端分数取平均后再计算 AUROC（池化）仅作对照。

为判断对接信号是否可由简单分子属性解释，我们采用同一 AUROC 流程，但以配体描述符替代对接分数，包括重原子数、分子量、cLogP 与 TPSA。每个靶对取其中 AUROC 最高的描述符作为对照基线，并报告对接 summary_min 与该基线之差 Δ 及其置信区间。若 Δ 的区间包含 0，则不足以支持对接提供了超出上述描述符的额外信息。

summary_min 的不确定度以 bootstrap 估计：在每个靶对内对配体有放回重采样 2000 次，报告 2.5%–97.5% 百分位区间（随机种子 20260729）。重采样以配体为单元。另报告按 Murcko 支架重采样的区间作为对照；正文以配体层区间为准。本文以置信区间作描述性不确定度报告，不对多靶对、多对照的全部对比作多重比较校正或正式假设检验。

### 2.7 对照与敏感性分析

为检验 2.6 所述 summary_min 是否可能主要由分子属性而非正确口袋上的对接分数驱动，我们设置下列对照：

1. **错口袋对照**：将口袋 A 与口袋 B 的分数对调后，仍按 dual 对 A_only / dual 对 B_only 计算 AUROC 并取 summary_min。  
2. **配体效率对照**：将各口袋分数除以重原子数后，再计算 summary_min。  
3. **匹配子集对照**：在 \|ΔpChEMBL\| ≤ 0.5（效价匹配）或 \|Δheavy atoms\| ≤ 2（尺寸匹配）的子集上，分别重算 dual 对 A_only 与 dual 对 B_only 的 AUROC。  
4. **协变量对照**：以逻辑回归比较“仅对接分数”与“对接分数 + 重原子数 + TPSA”两类模型的判别 AUROC。  
5. **二维结构基线**：以 ECFP4 指纹（Morgan 半径 2，2048 bit）与逻辑回归建立仅依赖二维结构的基线；交叉验证按 Murcko 支架分组（GroupKFold），使同一骨架不跨训练/测试折。

作为单靶参照，我们在 4L23（PIK3CA）与 4JT6（mTOR）上分别构建活性–decoy 集合：活性分子 pChEMBL ≥ 6.5；decoy 为同靶已测定且 pChEMBL ≤ 5.5 的弱效分子，并按分子量（±50）、logP（±1.5）与 TPSA（±25）与活性分子匹配。目标规模约为 50 个活性分子与 150 个 decoy；对接参数与 PIK3CA/mTOR 主面板一致（Vina，exhaustiveness = 16）。该分析用于估计单靶富集水平，不替代四类面板上的 summary_min 评价。

标签阈值敏感性（Table S4）与匹配子集完整结果（Table S5）见 Supporting Information。

### 2.8 软件与数据可用性

计算在 Python 3 环境下完成，主要软件包括 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore（公开预训练权重 `rtmscore_model1`）；Vina 姿态转为 SDF 时使用 Open Babel。AUROC、逻辑回归与交叉验证等分析使用常规 Python 科学计算栈（NumPy、SciPy、scikit-learn、pandas；版本见公开复现环境说明）。评价面板、对接分数、分析脚本与完整参数表将随公开数据包提供，详见 Data and Software Availability。

---

## 写法说明（不进正文）

对照 Vu et al. 2025 Methods 的取舍：

| Vu 怎么写 | 本稿怎么对齐 |
|-----------|--------------|
| Dataset 先给规模与规则 | 2.1 数据 → 2.2 任务 → 2.3 面板（Table 1） |
| Docking 写清盒子、准备、引擎；细节进 Table S1 | 2.4 受体/盒子/cognate；2.5 配体 + Vina + 重打分 |
| Evaluation 单独写指标 | 2.6；对照列表化 2.7 |
| 不写未做实验、不写仓库路径 | SI 只汇编已有数字；受体准备按冻结文件如实区分含氢坐标 vs ATOM+meeko |
| 过去时 / 陈述句 | 全文按此收紧 |

**本轮据仓库补全（有出处才写）：** ChEMBL API 锁定日；Murcko/多样性上限分面板；分辨率；meeko 受体路径差异；cognate RMSD 两种实现；Vina `scoring_function: vina`；ETKDG+MMFF200；bootstrap 描述性 CI、不做多重比较校正；不捏造 numpy/sklearn 版本号。
