# Methods（中文工作稿 · JCIM Articles）

> 结构与语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843  
> （Dataset → Task → Docking → Scoring → Evaluation；正文写清关键参数，其余见 SI）  
> 与 [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。投稿以英文为准；本稿供中文审改。  
> 本轮已吸收审稿式 Methods 批评（见 [`METHODS_REVIEWER_CRITIQUE_RESPONSE_V1.md`](METHODS_REVIEWER_CRITIQUE_RESPONSE_V1.md)）：统一标签作跨对主稳健分析、分数方向写死、GNINA 不对称披露、decoy 窗口引文献、cognate 细节压 SI；**不编造** median pChEMBL 表与 1000 次 panel 重抽结果。

---

## 2. Methods

### 2.1 数据策展与靶对选择

活性数据来自 ChEMBL Web API（公开 activity 端点；靶对供给审计锁定于 2026-07-23）。pChEMBL 是 ChEMBL 提供的近似统一活性标度，将若干摩尔浓度–响应型测量（如 IC50、Ki、Kd、EC50 等）转换为 −log10 形式，以便在公开数据整合中比较。不同 assay 类型与实验条件并不等价；本文将 pChEMBL 视为大规模策展中的常用近似，而非同一条件下的绝对亲和力。

同一配体–靶标若有多条记录，**代表值取可用的最大 pChEMBL**（activity inflation 风险见 Limitations）。未按 assay 中位数、ChEMBL 置信度阈值或物种字段过滤：冻结数据包（`mols_*.json`）仅存每个配体–靶标一对一的最大 pChEMBL 浮点数，不具备在本地重算 median / confidence≥8 / Homo sapiens 过滤的条件。因此本稿**不报告** max 对 median 的数值敏感性表；该分析需重新拉取逐条 assay 记录后方可进行，列入后续工作而非现有结果。任一端缺少可用 pChEMBL 的配体不进入四类主分析。ChEMBL 结构条目常含盐形式；对接前将结构按连通片段拆分并保留重原子数最多的有机片段，以得到单一可对接分子。

**供给审计用严格规则（construction gate，非跨对主比较的唯一标签）。** 严格规则：dual 为两端 pChEMBL ≥ 6.5；A_only 为 A ≥ 6.5 且 B ≤ 5.5；B_only 为对称定义；neither 为两端 ≤ 5.5；灰区（介于 5.5 与 6.5）不进入严格面板。在 49 对可审计靶对上，两端严格硬负均 ≥ 50 的仅 4 对；排除金属依赖的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 作为厚面板候选。EGFR/HER2 在同一严格规则下 B_only 仅 7 个，无法建成规模均衡的严格四类面板，但仍纳入评价集作为**供给受限案例**（文献常见双靶对；不宣称其为严格厚面板）。

**面板建造（panel construction）允许在供给不足时使用单阈值规则，但这是建造协议，不是为抬高 AUROC 而事后改阈值。** AChE/BChE 与 PIK3CA/PIK3CB 按严格规则完成配额抽样；EGFR/HER2 与 PIK3CA/mTOR（PM48）因严格规则下单端选择性配体过少，按 θ = 6.0 建成（两端 ≥ θ 为 dual；一端 ≥ θ 且对端 < θ 为对应单靶类）。建造规则在抽样前按供给审计结果冻结，并写入 Table 1；阈值选择的动机是**凑齐可分析配额**，不是在观察对接分数后回改标签。

**跨对主稳健分析（primary cross-pair robustness）：统一标签重标。** 为消除“不同靶对用不同阈值”的质疑，我们在既有面板配体与既有 Vina 分数上，对四对统一施加 θ ∈ {5.5, 6.0, 6.5} 与严格 6.5/5.5 规则重标四类，并重算口袋匹配 summary_min（Supporting Information Table S4）。正文跨对排序与稳健性结论以该统一重标为准（Results 3.8）；Table 2 并列报告各面板**建造时**标签下的点估计，作为 construction readout，二者不得混称为两套互相竞争的“主标准”。统一重标下若某类 n 过小，标记 underpowered，不作功效充足主张。

### 2.2 任务定义

本文将双靶对接评测建成四类配体判别任务。对每一对靶标 A/B，依据两端实测活性将配体分为：dual（两端均强）、A_only（仅 A 强）、B_only（仅 B 强）和 neither（两端均弱）。评价目标不是比较某一口袋上谁的对接分更高，而是判断对接分数能否把 dual 配体同时与两类单靶选择性配体（A_only 与 B_only）区分开。

为此，我们分别计算 dual 对 A_only、dual 对 B_only 的 AUROC，并以二者的最小值作为该靶对的汇总指标（记为 summary_min）。采用最小值是为了避免只报告较好一侧：若仅有一端可区分而另一端接近随机，则不足以支持“双靶判别”的主张。

### 2.3 评价面板的构建

对选定的每一对靶标，我们按 2.1 所述**建造规则**从 ChEMBL 候选集中抽样组成评价面板。PIK3CA/mTOR（PM48）在抽样时限制同一类别内同一 Murcko 支架最多 2 个分子；EGFR/HER2（n = 110）相应上限为 5；AChE/BChE 与 PIK3CA/PIK3CB 在建面时以类别配额控制组成，并以标识符前缀作弱多样性约束。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则定额抽取 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2 沿用既有 θ = 6.0 面板（n = 110）。PIK3CA/mTOR 按 θ = 6.0 建成 n = 48 的主比较面板（PM48；dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。对接失败的配体–受体组合从相应受体分数中剔除，故分析用类别计数可低于建面定额（Table 1）。主评价仅使用 dual、A_only 与 B_only 计算 AUROC；neither 保留在面板中供描述，但不进入上述对比。

**Table 1.** 评价集组成与对接设置（建造规则）

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

由于 PM48 的 dual / A_only / B_only 样本量较小、summary_min 置信区间较宽，我们在 PIK3CA/mTOR 上另构建扩面面板：保留 PM48 的全部 48 个配体，并按严格规则追加新分子（目标配额 30 / 30 / 30 / 25），实际 n = 115。该扩面面板是 PM48 的超集，用于收窄区间并检验结论是否同向，而不是独立重复实验。主文以 PM48 报告与其他三对在同一协议下的比较；扩面结果另报。

面板组成对结论的影响：现报告配体层与 Murcko 支架层 bootstrap CI（2.6）。**从严格供给池按同一配额重新抽取 1000 个独立 panel 并重算 summary_min 的分布**，需要尚未对接的池内分子分数，超出当前冻结分数包范围；该分析列入后续对接工作（可与 unused-pool holdout 合并执行），本稿不虚构其数值。

### 2.4 蛋白结构与结合位点定义

受体结构取自 Protein Data Bank 中含小分子共晶配体的条目。PIK3CA 与 mTOR 分别使用 4L23 与 4JT6（共晶配体 X6K / PI-103）；AChE 与 BChE 使用 4EY7 与 4BDS（E20 / THA）；PIK3CB 使用 2WXF（039）；EGFR 与 HER2 使用 3POZ 与 3RCD（03P / TAK-285）。分辨率见表 1。结合位点定义为以共晶配体为中心的盒子：共晶配体轴对齐包围盒外扩 5 Å，每边长度下限 20 Å（坐标见 Supporting Information Table S2）。

**受体准备（按实际冻结流程如实写，不补充未执行步骤）。** 对接用受体均去除水分子与共晶配体后转为 PDBQT。PIK3CA、mTOR、EGFR 与 HER2：采用冻结目录中已含氢的蛋白坐标，经 meeko `mk_prepare_receptor.py --read_pdb` 生成 PDBQT。AChE、BChE 与 PIK3CB：由沉积 PDB 的 ATOM/TER 记录提取蛋白（去除水与异源原子），以 meeko `mk_prepare_receptor`（默认 alternate location A）生成 PDBQT。未另行运行 PDBFixer 补全缺失原子，亦未用 Reduce 做独立的 pH 7.4 质子化枚举；组氨酸互变异构与局部质子态采用 meeko 默认处理。本课题受体均为非共价小分子共晶体系，对接网格中不保留金属辅因子作为可对接位点组分。上述省略构成复现边界，完整命令与输入文件随公开数据包提供。

**Cognate 协议检查（细节表见 SI）。** 八个冻结受体均以共晶配体重对接：通过门槛为 9 个姿态中重原子 RMSD 的最小值（best_of_9）&lt; 2.0 Å（对接坐标系，不做蛋白叠合）。PIK3CA/mTOR 与 EGFR/HER2 用 meeko `REMARK SMILES IDX` + 图自同构最小 CalcRMS；AChE/BChE 与 PIK3CB 用重原子匈牙利匹配。除 mTOR（4JT6）在 exhaustiveness = 8 时未过门槛、升至 16 后通过外，其余受体在 E = 8 时 best_of_9 &lt; 2 Å；EGFR（3POZ）存在 best_of_9 过关但 Vina mode1 远离共晶的情形。完整数值、mode1 与 best_of_9 对照见 Supporting Information Table S3。因此 PIK3CA/mTOR 主分析采用 exhaustiveness = 16，其余靶对采用 8，并在 PIK3CA/mTOR 上报告 E = 8 对照。

### 2.5 对接与打分

配体由 ChEMBL SMILES 去盐（保留最大有机片段）后，用 RDKit 加氢，以 ETKDGv3 生成三维构象（随机种子 20260727），并用 MMFF 力场局部优化（最多 200 步），再经 meeko 默认参数转为 PDBQT。对接使用 AutoDock Vina 1.2.7，打分函数为默认 `vina`：每个配体保留 9 个姿态，`energy_range = 3`，exhaustiveness 按 Table 1。构象生成、面板抽样与对接均使用固定随机种子；完整参数见 Supporting Information Table S1。

在同一组 Vina 姿态上另用两种函数重打分，以检查结论是否依赖单一打分通道。RTMScore（公开权重 `rtmscore_model1`）对 **全部 9 个姿态**取最高分。GNINA 1.3.2 在 CPU 模式下，将 **Vina 排序第一的姿态**（mode 1）经 Open Babel 转为 SDF 后，以 `--cnn_scoring rescore --minimize` 重打分。**因此 RTM 与 GNINA 的姿态覆盖不对称**：前者为 best-of-9，后者为 mode-1 rescore；GNINA 通道可能系统性吃亏，不作“三引擎已公平对齐”的主张，正文仅把 RTM/GNINA 作定性通道对照。若需公平比较，应对全部 9 个 Vina 姿态分别 GNINA rescore 后取最优（后续工作）。

### 2.6 评价指标与基线

**分数方向（固定，便于复现）。** 记 Vina 输出的结合能为 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，通常为负且越负表示预测结合越强）。定义对接分数  
\[
S = -E_{\mathrm{Vina}},
\]  
使 \(S\) **越大越好**。RTMScore 与 GNINA CNN 分数本身已是“越高越好”，不再取负。所有 AUROC 均以 \(S\)（或对应重打分）为连续预测值；在 dual 对 A_only 与 dual 对 B_only 的二分类中，**dual 始终为正类（positive class）**，A_only 或 B_only 为负类。

每个配体在口袋 A 与口袋 B 上各有一个分数 \(S\)。主评价计算两条二分类 AUROC：dual 对 A_only 使用口袋 B 的分数；dual 对 B_only 使用口袋 A 的分数。两条 AUROC 的较小值记为 summary_min。作为聚合对照，另报告池化与 worst-pocket summary_min（不作主指标）。

为判断对接信号是否可由简单分子属性解释，我们采用同一 AUROC 流程，但以配体描述符替代对接分数，包括重原子数、分子量、cLogP 与 TPSA。每个靶对取其中 AUROC 最高的描述符作为对照基线，并报告对接 summary_min 与该基线之差 Δ 及其置信区间。若 Δ 的区间包含 0，则不足以支持对接提供了超出上述描述符的额外信息。

summary_min 的不确定度以 bootstrap 估计：在每个靶对内对配体有放回重采样 2000 次，报告 2.5%–97.5% 百分位区间（随机种子 20260729）。重采样以配体为单元。另报告按 Murcko 支架重采样的区间作为对照；正文以配体层区间为准。本文以置信区间作描述性不确定度报告，不对多靶对、多对照的全部对比作多重比较校正或正式假设检验。

### 2.7 对照与敏感性分析

为检验 2.6 所述 summary_min 是否可能主要由分子属性而非正确口袋上的对接分数驱动，我们设置下列对照：

1. **错口袋对照**：将口袋 A 与口袋 B 的分数对调后，仍按 dual 对 A_only / dual 对 B_only 计算 AUROC 并取 summary_min。  
2. **配体效率对照**：将各口袋分数除以重原子数后，再计算 summary_min。  
3. **匹配子集对照**：在 \|ΔpChEMBL\| ≤ 0.5（效价匹配）或 \|Δheavy atoms\| ≤ 2（尺寸匹配）的子集上，分别重算 dual 对 A_only 与 dual 对 B_only 的 AUROC。  
4. **协变量对照**：以逻辑回归比较“仅对接分数”与“对接分数 + 重原子数 + TPSA”两类模型的判别 AUROC，并报告对接分数的回归系数与优势比（OR）。  
5. **二维结构基线**：以 ECFP4 指纹（Morgan 半径 2，2048 bit）与逻辑回归建立仅依赖二维结构的基线；交叉验证按 Murcko 支架分组（GroupKFold），使同一骨架不跨训练/测试折。  
6. **统一标签重标**（2.1；Table S4）：跨对主稳健分析。

作为单靶参照（非主指标），我们在 4L23（PIK3CA）与 4JT6（mTOR）上分别构建活性–decoy 集合：活性分子 pChEMBL ≥ 6.5；decoy 为同靶已测定且 pChEMBL ≤ 5.5 的弱效分子，并按分子量（±50 Da）、logP（±1.5）与 TPSA（±25 Å²）与活性分子做性质匹配。上述窗口沿用 property-matched decoy 的常见设定（如 DUD-E 对分子量、logP 等物理化学性质的匹配思想；Mysinger et al., *J. Med. Chem.* 2012, 55, 6582–6594），TPSA 窗口为在同一思想下增加的极性匹配项，用于降低“用完全无关分子当 decoy”造成的虚高富集。目标规模约为 50 个活性分子与 150 个 decoy；对接参数与 PIK3CA/mTOR 主面板一致（Vina，exhaustiveness = 16）。该分析用于估计单靶富集水平，不替代四类面板上的 summary_min 评价。

匹配子集完整结果（Table S5）与分数聚合对照（Table S6）见 Supporting Information。

### 2.8 软件与数据可用性

计算在 Python 3 环境下完成，主要软件包括 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore（公开预训练权重 `rtmscore_model1`）；Vina 姿态转为 SDF 时使用 Open Babel。AUROC、逻辑回归与交叉验证等分析使用常规 Python 科学计算栈（NumPy、SciPy、scikit-learn、pandas；版本见公开复现环境说明）。评价面板、对接分数、分析脚本与完整参数表将随公开数据包提供，详见 Data and Software Availability。

---

## 写法说明（不进正文）

对照 Vu et al. 2025 Methods 的取舍：

| Vu 怎么写 | 本稿怎么对齐 |
|-----------|--------------|
| Dataset 先给规模与规则 | 2.1 数据 → 2.2 任务 → 2.3 面板（Table 1） |
| Docking 写清盒子、准备、引擎；细节进 Table S1 | 2.4 受体/盒子；cognate 数字进 SI S3 |
| Evaluation 单独写指标 | 2.6（含 \(S=-E\)）；对照列表化 2.7 |
| 不写未做实验、不写仓库路径 | 不编造 PDBFixer/Reduce/median 表/1000-panel 分布 |

**明确不做/未做（防审稿追问时撒谎）：** max vs median 数值表；GNINA 九姿态公平重打；从供给池重抽 1000 panel 的 summary_min 分布；PDBFixer+Reduce 质子化流程。
