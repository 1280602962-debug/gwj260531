# Methods（中文工作稿 · JCIM Articles）

> 结构与语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843  
> （Dataset → Task → Docking → Scoring → Evaluation；正文写清关键参数，其余见 SI）  
> 与 [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套。投稿以英文为准；本稿供中文审改。  
> 本轮已吸收审稿式 Methods 批评（见 [`METHODS_REVIEWER_CRITIQUE_RESPONSE_V1.md`](METHODS_REVIEWER_CRITIQUE_RESPONSE_V1.md)）：统一标签作跨对主稳健分析、分数方向写死、GNINA 不对称披露、decoy 窗口引文献、cognate 细节压 SI；并补入已做实验的可复现边界（holdout 三对范围、换晶单口袋替换、叠合与接触计数参数）。**不编造** median pChEMBL 表与 1000 次 panel 重抽结果。

---

## 2. Methods

### 2.1 数据策展与靶对选择

活性数据来自 ChEMBL Web API（公开 activity 端点；靶对供给审计锁定于 2026-07-23）。pChEMBL 是 ChEMBL 提供的近似统一活性标度，将若干摩尔浓度–响应型测量（如 IC50、Ki、Kd、EC50 等）转换为 −log10 形式，以便在公开数据整合中比较。不同 assay 类型与实验条件并不等价；本文将 pChEMBL 视为大规模策展中的常用近似，而非同一条件下的绝对亲和力。

同一配体–靶标若有多条记录，**代表值取可用的最大 pChEMBL**（activity inflation 风险见 Limitations）。未按 assay 中位数、ChEMBL 置信度阈值或物种字段过滤：冻结数据包（`mols_*.json`）仅存每个配体–靶标一对一的最大 pChEMBL 浮点数，不具备在本地重算 median / confidence≥8 / Homo sapiens 过滤的条件。因此本稿**不报告** max 对 median 的数值敏感性表；该分析需重新拉取逐条 assay 记录后方可进行，列入后续工作而非现有结果。任一端缺少可用 pChEMBL 的配体不进入四类主分析。ChEMBL 结构条目常含盐形式；对接前将结构按连通片段拆分并保留重原子数最多的有机片段，以得到单一可对接分子。

**供给审计用严格规则（construction gate，非跨对主比较的唯一标签）。** 严格规则：dual 为两端 pChEMBL ≥ 6.5；A_only 为 A ≥ 6.5 且 B ≤ 5.5；B_only 为对称定义；neither 为两端 ≤ 5.5；灰区（介于 5.5 与 6.5）不进入严格面板。在 49 对可审计靶对上，两端严格硬负均 ≥ 50 的仅 4 对；排除金属依赖的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 作为厚面板候选。EGFR/HER2 在同一严格规则下 B_only 仅 7 个，无法建成规模均衡的严格四类面板，但仍纳入评价集作为**供给受限案例**（文献常见双靶对；不宣称其为严格厚面板）。

**跨库计数核对（零对接，不重建面板）。** 为检验上述 ≥50 双侧硬负门槛是否只是 ChEMBL 覆盖假象，我们对冻结的 K = 4 四对另查 BindingDB REST（`getLigandsByUniprots`，cutoff = 1 mM，以免把弱端测定截掉）与 PubChem PUG REST（`protein/accession/…/concise`）。类型限于 IC50/Ki/Kd/EC50；代表值取最大转换 p 活性；分类规则与 J0 严格门槛相同。配体身份分别用 BindingDB monomerid 与 PubChem CID，**不做**跨库结构合并。主比较采用与 pChEMBL 更接近的**等式测定**（去掉 `>`/`<` 截尾）；将不等式数值当作点估计计入只作敏感性。该核对只报告计数（Supporting Information Table S12），不进入对接或改写冻结面板。

**面板建造（panel construction）允许在供给不足时使用单阈值规则，但这是建造协议，不是为抬高 AUROC 而事后改阈值。** AChE/BChE 与 PIK3CA/PIK3CB 按严格规则完成配额抽样；EGFR/HER2 与 PIK3CA/mTOR（PM48）因严格规则下单端选择性配体过少，按 θ = 6.0 建成（两端 ≥ θ 为 dual；一端 ≥ θ 且对端 < θ 为对应单靶类）。建造规则在抽样前按供给审计结果冻结，并写入 Table 1；阈值选择的动机是**凑齐可分析配额**，不是在观察对接分数后回改标签。

**主表述采用单一统一标签规则（θ = 6.0）。** 为消除“不同靶对用不同阈值”的质疑，正文 Table 2 对全部四对统一采用 θ = 6.0 规则报告主结果。对 EGFR/HER2 与 PIK3CA/mTOR，这与建造时直接采用的规则相同；对 AChE/BChE 与 PIK3CA/PIK3CB，建造时按更严格的 6.5/5.5 规则完成供给配额抽样，但在本数据上 θ = 6.0 给出与该严格规则**完全相同**的配体分类与 AUROC（Supporting Information Table S4），即标签在阈值网格内对这两对不敏感。作为支持性稳健性分析，我们进一步在 θ ∈ {5.5, 6.5} 与严格 6.5/5.5 规则下重标四类并重算口袋匹配 summary_min（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 在严格规则下 B_only 样本量过小，标记 underpowered，不作功效充足主张。该敏感性网格不是与 Table 2 竞争的第二套主标准，只用于证明排序不随阈值网格翻转（Results 3.2）。

### 2.2 任务定义

本文将双靶对接评测建成四类配体判别任务。对每一对靶标 A/B，依据两端实测活性将配体分为：dual（两端均强）、A_only（仅 A 强）、B_only（仅 B 强）和 neither（两端均弱）。评价目标不是比较某一口袋上谁的对接分更高，而是判断对接分数能否把 dual 配体同时与两类单靶选择性配体（A_only 与 B_only）区分开。

为此，我们分别计算 dual 对 A_only、dual 对 B_only 的 AUROC，并以二者的最小值作为该靶对的汇总指标（记为 summary_min）。采用最小值是为了避免只报告较好一侧：若仅有一端可区分而另一端接近随机，则不足以支持“双靶判别”的主张。

### 2.3 评价面板的构建

对选定的每一对靶标，我们按 2.1 所述**建造规则**从 ChEMBL 候选集中抽样组成评价面板。PIK3CA/mTOR（PM48）在抽样时限制同一类别内同一 Murcko 支架最多 2 个分子；EGFR/HER2（n = 110）相应上限为 5；AChE/BChE 与 PIK3CA/PIK3CB 在建面时以类别配额控制组成，并以标识符前缀作弱多样性约束。AChE/BChE 与 PIK3CA/PIK3CB 按严格规则定额抽取 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2 沿用既有 θ = 6.0 面板（n = 110）。PIK3CA/mTOR 按 θ = 6.0 建成 n = 48 的主比较面板（PM48；dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。对接失败的配体–受体组合从相应受体分数中剔除；任一端缺少可用分数的配体不进入需要两端分数的口袋匹配 AUROC，故分析用类别计数可低于建面定额（Table 1）。主评价仅使用 dual、A_only 与 B_only 计算 AUROC；neither 保留在面板中供描述，但不进入上述对比。

**Table 1.** 评价集组成与对接设置（建造规则）

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

由于 PM48 的 dual / A_only / B_only 样本量较小、summary_min 置信区间较宽，我们在 PIK3CA/mTOR 上另构建扩面面板（记为 PM110，为历史命名）。该面板保留 PM48 的全部 48 个配体，并按严格规则追加新分子，目标配额 dual / A_only / B_only / neither = 30 / 30 / 30 / 25，实际面板 n = 115。口袋匹配主对比只用前三类（各 30 个）。扩面面板是 PM48 的超集，用于收窄区间并检验点估计是否同向，而不是独立重复实验。主文以 PM48 报告与其他三对在同一协议下的比较；扩面结果另报（Results 3.5）。

面板组成对结论的影响：现报告配体层与 Murcko 支架层 bootstrap CI（2.6）。从严格供给池按同一配额重复抽取**互不重叠**的 30/30/30 panel 在供给上做不到接近 1000 次：主面板与 holdout 用掉硬负后，剩余严格 A_only/B_only 为 PIK3CA/mTOR 37/39、AChE/BChE 141/30、PIK3CA/PIK3CB 8/19、EGFR/HER2 22/0（`C_CLASS_EXPERIMENT_NECESSITY_VERDICT_V1.md`）。硬负总量本身只有数十至一百余个，不是“把未用池对接完就能得到独立 1000-panel 分布”。本稿以一次 unused-pool holdout（Results 3.9）作为面板外检验，不虚构 1000 次重抽数值，也不把已对接配体的有放回 bootstrap 称作供给池重抽。

### 2.4 蛋白结构与结合位点定义

受体结构取自 Protein Data Bank 中含小分子共晶配体的条目。PIK3CA 与 mTOR 分别使用 4L23 与 4JT6（共晶配体 X6K / PI-103）；AChE 与 BChE 使用 4EY7 与 4BDS（E20 / THA）；PIK3CB 使用 2WXF（039）；EGFR 与 HER2 使用 3POZ 与 3RCD（03P / TAK-285）。分辨率见表 1。结合位点定义为以共晶配体为中心的盒子：共晶配体轴对齐包围盒外扩 5 Å，每边长度下限 20 Å（坐标见 Supporting Information Table S2）。

**受体准备（按实际冻结流程如实写，不补充未执行步骤）。** 对接用受体均去除水分子与共晶配体后转为 PDBQT。PIK3CA、mTOR、EGFR 与 HER2：采用冻结目录中已含氢的蛋白坐标，经 meeko `mk_prepare_receptor.py --read_pdb` 生成 PDBQT。AChE、BChE 与 PIK3CB：由沉积 PDB 的 ATOM/TER 记录提取蛋白（去除水与异源原子），以 meeko `mk_prepare_receptor`（默认 alternate location A）生成 PDBQT。未另行运行 PDBFixer 补全缺失原子，亦未用 Reduce 做独立的 pH 7.4 质子化枚举；组氨酸互变异构与局部质子态采用 meeko 默认处理。本课题受体均为非共价小分子共晶体系，对接网格中不保留金属辅因子作为可对接位点组分。上述省略构成复现边界，完整命令与输入文件随公开数据包提供。

**Cognate 协议检查（细节表见 SI）。** 八个冻结受体均以共晶配体重对接：通过门槛为 9 个姿态中重原子 RMSD 的最小值（best_of_9）&lt; 2.0 Å（对接坐标系，不做蛋白叠合）。PIK3CA/mTOR 与 EGFR/HER2 用 meeko `REMARK SMILES IDX` + 图自同构最小 CalcRMS；AChE/BChE 与 PIK3CB 用重原子匈牙利匹配。除 mTOR（4JT6）在 exhaustiveness = 8 时未过门槛、升至 16 后通过外，其余受体在 E = 8 时 best_of_9 &lt; 2 Å；EGFR（3POZ）存在 best_of_9 过关但 Vina mode1 远离共晶的情形。完整数值、mode1 与 best_of_9 对照见 Supporting Information Table S3。因此 PIK3CA/mTOR 主分析采用 exhaustiveness = 16，其余靶对采用 8，并在 PIK3CA/mTOR 上报告 E = 8 对照。

### 2.5 对接与打分

配体由 ChEMBL SMILES 去盐（保留最大有机片段）后，用 RDKit 加氢，以 ETKDGv3 生成三维构象（随机种子 20260727），并用 MMFF 力场局部优化（最多 200 步），再经 meeko 默认参数转为 PDBQT。对接使用 AutoDock Vina 1.2.7，打分函数为默认 `vina`：每个配体保留 9 个姿态，`energy_range = 3`，exhaustiveness 按 Table 1，随机种子 20260727（与 ETKDG 相同）。构象生成、面板抽样与对接均使用固定随机种子；完整参数见 Supporting Information Table S1。

在同一组 Vina 姿态上另用两种函数重打分，以检查结论是否依赖单一打分通道。RTMScore（公开权重 `rtmscore_model1`）对 **全部 9 个姿态**取最高分。GNINA 1.3.2 在 CPU 模式下最初仅将 **Vina 排序第一的姿态**（mode 1）经 Open Babel 转为 SDF 后以 `--cnn_scoring rescore --minimize` 重打分，与 RTM 的 best-of-9 覆盖不对称。**2026-08-24 已补做全 9 姿态 GNINA 公平重打**：对每个配体–靶标的全部 9 个 Vina 姿态分别转 SDF 并重打分，取每端最高 CNNscore，与 RTM 姿态覆盖对齐；mode-1 结果保留为历史对照（`scores_gnina_*_mode01_backup.csv`）。全 9 姿态重打后，GNINA 口袋匹配 summary_min 相对 mode-1 的变化很小且方向不一致（K=4：AChE/BChE −0.03、PIK3CA/PIK3CB −0.02、PIK3CA/mTOR +0.08、EGFR/HER2 −0.04；`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`），姿态覆盖不对称因而**不是**GNINA 表现偏弱的主要原因。三对（EGFR/HER2、AChE/BChE、PIK3CA/mTOR）上 GNINA 不超过同面板 Vina 口袋匹配；PIK3CA/PIK3CB 上 GNINA 略高于 Vina（0.533 对 0.500），但该关系在 mode-1 时已存在（0.554 对 0.500），并非 best-of-9 新产生，且二者均接近随机、区间重叠，不构成“GNINA 优于 Vina”的主张。正文仍把 RTM/GNINA 作定性通道对照，不作“三引擎已公平对齐即验证同一决策臂”的主张。

### 2.6 评价指标与基线

**分数方向（固定，便于复现）。** 记 Vina 输出的结合能为 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，通常为负且越负表示预测结合越强）。定义对接分数  
\[
S = -E_{\mathrm{Vina}},
\]  
使 \(S\) **越大越好**。RTMScore 与 GNINA CNN 分数本身已是“越高越好”，不再取负。所有 AUROC 均以 \(S\)（或对应重打分）为连续预测值；在 dual 对 A_only 与 dual 对 B_only 的二分类中，**dual 始终为正类（positive class）**，A_only 或 B_only 为负类。

每个配体在口袋 A 与口袋 B 上各有一个分数 \(S\)。主评价计算两条二分类 AUROC：dual 对 A_only 使用口袋 B 的分数；dual 对 B_only 使用口袋 A 的分数。两条 AUROC 的较小值记为 summary_min。作为聚合对照，另报告池化与 worst-pocket summary_min（不作主指标）。

为判断对接信号是否可由简单分子属性解释，我们采用同一 AUROC 流程，但以 RDKit 计算的配体描述符替代对接分数，包括重原子数（GetNumHeavyAtoms）、分子量（MolWt）、cLogP（MolLogP）与 TPSA。每个靶对取其中 AUROC 最高的描述符作为对照基线，并报告对接 summary_min 与该基线之差 Δ 及其置信区间。若 Δ 的区间包含 0，则不足以支持对接提供了超出上述描述符的额外信息。

summary_min 的不确定度以 bootstrap 估计：在每个靶对内对配体有放回重采样 2000 次，报告 2.5%–97.5% 百分位区间（随机种子 20260729）。重采样以配体为单元。另报告按 Murcko 支架重采样的区间作为对照；正文以配体层区间为准。本文以置信区间作描述性不确定度报告，不对多靶对、多对照的全部对比作多重比较校正或正式假设检验。

### 2.7 对照与敏感性分析

为检验 2.6 所述 summary_min 是否可能主要由分子属性而非正确口袋上的对接分数驱动，我们设置下列对照：

1. **错口袋对照**：将口袋 A 与口袋 B 的分数对调后，仍按 dual 对 A_only / dual 对 B_only 计算 AUROC 并取 summary_min。  
2. **配体效率对照**：将各口袋分数除以重原子数后，再计算 summary_min。  
3. **匹配子集对照**：在 \|ΔpChEMBL\| ≤ 0.5（效价匹配）或 \|Δheavy atoms\| ≤ 2（尺寸匹配）的子集上，分别重算 dual 对 A_only 与 dual 对 B_only 的 AUROC。  
4. **协变量对照**：以逻辑回归（scikit-learn `LogisticRegression`，C = 1.0，max_iter = 2000）比较“仅对接分数”与“对接分数 + 重原子数 + TPSA”两类模型的判别 AUROC，并报告对接分数的回归系数与优势比（OR）。  
5. **二维结构基线**：以 ECFP4 指纹（Morgan 半径 2，2048 bit）与同一设定的逻辑回归建立仅依赖二维结构的基线；交叉验证按 Murcko 支架分组（`GroupKFold`），折数取 min(5, 正类数, 负类数, 支架数) 且不少于 2，使同一骨架不跨训练/测试折。  
6. **统一标签重标**（2.1；Table S4）：阈值敏感性支持性分析。
7. **跨对结构决定因素（探索性）**：从各冻结受体 `*_protein.pdb` 中，用 Biopython `PDBParser` 提取最长蛋白链的一级序列，仅计入标准氨基酸 ATOM 记录；以 `Bio.Align.PairwiseAligner`（BLOSUM62 替换矩阵，全局比对，gap open = −11、extend = −1）对每对靶标内的两条链做两两比对，报告全链序列一致性（分别以比对长度与较短链长度归一，Supporting Information Table S7）。该指标为整体结构相似度的粗粒度代理，不涉及口袋残基对应或结构叠合，不用于口袋 RMSD 或 PLIF 主张。

作为单靶参照（非主指标），我们在 4L23（PIK3CA）与 4JT6（mTOR）上分别构建活性–decoy 集合：活性分子 pChEMBL ≥ 6.5；decoy 为同靶已测定且 pChEMBL ≤ 5.5 的弱效分子，并按分子量（±50 Da）、logP（±1.5）与 TPSA（±25 Å²）与活性分子做性质匹配。上述窗口沿用 property-matched decoy 的常见设定（如 DUD-E 对分子量、logP 等物理化学性质的匹配思想；Mysinger et al., *J. Med. Chem.* 2012, 55, 6582–6594），TPSA 窗口为在同一思想下增加的极性匹配项，用于降低“用完全无关分子当 decoy”造成的虚高富集。目标规模约为 50 个活性分子与 150 个 decoy；对接参数与 PIK3CA/mTOR 主面板一致（Vina，exhaustiveness = 16）。该分析用于估计单靶富集水平，不替代四类面板上的 summary_min 评价。

匹配子集完整结果（Table S5）与分数聚合对照（Table S6）见 Supporting Information。

作为**面板外冻结验证**（post-hoc unused-pool holdout），holdout **只覆盖有足够 unused-pool 配额的三对**：PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB。EGFR/HER2 为供给受限案例，严格规则下 B_only 仅 7 个，不具备与另外三对同等的未用池抽样条件，故不进入 holdout。对上述三对，我们在建面时已用严格规则筛出、但未进入对应冻结面板的 ChEMBL 候选池中抽样：以 ChEMBL ID 精确排除已用条目（PIK3CA/mTOR 排除的是 PM110 超集，从而覆盖 PM48），再以新种子 `HOLDOUT_SEED=20260731` 按 dual / A_only / B_only = 20 / 20 / 20 定额抽取（Murcko 支架封顶 3 个/类）。现查 SMILES 后按与主面板完全相同的 RDKit/meeko 制备与 Vina 协议对接（受体、盒子、exhaustiveness、种子均不重调）。评价仍用 2.6 的口袋匹配 summary_min 与配体层 bootstrap；平凡基线在同一 holdout 配体上并列计算。含硼配体 HOAP_028 因 AutoDock 原子类型 `B` 不支持、两端均未得分，已从 AUROC 装配中剔除（59/60 配体进入该对分析）。该 holdout 用于检验“分数规则与协议是否只在建面板时凑效”，**不是**跨数据库的独立外部验证集；抽样清单先冻结、后看分数（Supporting Information Table S8）。

结构稳健性方面，替代晶体的入选要求为：polymer entity 核实为真 PIK3CA α 或真 mTOR（非嵌合体）、含 ATP 位点小分子共晶，并先做与 Methods 2.4 相同的共晶配体重对接 QC（best_of_9 &lt; 2 Å）。通过者为 PIK3CA **4JPS**（共晶 1LT）、**5DXT**（5H5）与 mTOR **4JSX**（Torin2 / 17G；与 4JT6 同属 mTORΔN–mLST8 截短构建体家族，但共晶配体与晶型不同）。嵌合体结构（如曾误用的 3T8M / PIK3CG 骨架）一律排除。换晶时**一次只替换一个口袋**：4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结的 4JT6 分数；4JSX 替换口袋 B，口袋 A 仍用冻结的 4L23 分数。新盒子按该替代晶体自身共晶配体、以 2.4 的同一 AABB + 5 Å / 边长下限 20 Å 规则生成；受体准备与 2.4 相同，exhaustiveness = 16、随机种子 20260727。仅在冻结 PM48 配体上重对接并重算 summary_min。

**受体依赖的探索性结构对照（零新对接）。** 为对照换晶体后 PIK3CA 端崩溃、mTOR 端仅小幅下降的不对称，我们直接在已冻结晶体坐标上做刚体叠合：以 Biopython `PDBParser` 提取各受体最长蛋白链的 Cα 坐标，按残基编号与残基名精确匹配（不匹配即剔除），用 `Superimposer` 对全部匹配 Cα 做一次 Kabsch 拟合得到全域 RMSD；口袋残基由参考结构自身共晶配体的重原子 ≤5 Å 界定，在**同一变换**下计算口袋局域 RMSD，不做二次局部拟合。再将替代结构自身的共晶配体坐标按同一变换投影，计算其质心与参考结构共晶配体质心的距离，检验二者是否落在同一大类口袋。5DXT 匹配的 Cα 少于 4JPS，全域 RMSD 不是等覆盖比较；本对照仅含 PIK3CA 替代结构 n = 2、mTOR n = 1，不预设 Cα RMSD 能够定量解释 AUROC 变化（Supporting Information Table S10）。

**Holdout 错口袋对照的几何对照（探索性，零新对接）。** 为检验 holdout 上 `wrong_pocket_control_vina` 不低于 `pocket_matched_vina` 是否仅为打分函数伪象，我们直接在已冻结的 **mode-1** 姿态上计算一个不依赖打分函数的几何量：配体重原子中与受体重原子距离 ≤4.0 Å 的原子数（`contact_count`；4.0 Å 为粗粒度接触阈值，非经验证 PLIF）。用该量在口袋 A 上比较 dual 对 A_only、在口袋 B 上比较 dual 对 B_only，与错口袋对照的同口袋比较同构，作为与 Vina 结果并列的几何对照，**不预设**其幅度与 Vina 错口袋一致（Supporting Information Table S11）。

**Holdout 错口袋的效价/尺寸匹配诊断（零新对接）。** 为检验该悖论是否来自 unused-pool 抽样相对主面板的效价或尺寸偏移，我们用与 Table S5 相同的最近邻匹配（效价：共享活性端 \|ΔpChEMBL\| ≤ 0.5；尺寸：\|Δheavy\| ≤ 2）在 holdout 上重算口袋匹配与错口袋 AUROC（Supporting Information Table S13）。该诊断不改写 Table S8 的主 holdout 数字。

### 2.8 软件与数据可用性

计算在 Python 3 环境下完成，主要软件包括 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore（公开预训练权重 `rtmscore_model1`）；Vina 姿态转为 SDF 时使用 Open Babel。刚体叠合与全链序列比对使用 Biopython（`PDBParser`、`Superimposer`、`PairwiseAligner`）。AUROC、逻辑回归与交叉验证等分析使用常规 Python 科学计算栈（NumPy、SciPy、scikit-learn、pandas；版本见公开复现环境说明）。评价面板、对接分数、分析脚本与完整参数表将随公开数据包提供，详见 Data and Software Availability。

---

## 写法说明（不进正文）

对照 Vu et al. 2025 Methods 的取舍：

| Vu 怎么写 | 本稿怎么对齐 |
|-----------|--------------|
| Dataset 先给规模与规则 | 2.1 数据 → 2.2 任务 → 2.3 面板（Table 1） |
| Docking 写清盒子、准备、引擎；细节进 Table S1 | 2.4 受体/盒子；cognate 数字进 SI S3 |
| Evaluation 单独写指标 | 2.6（含 \(S=-E\)）；对照列表化 2.7 |
| 不写未做实验、不写仓库路径 | 不编造 PDBFixer/Reduce/median 表/1000-panel 分布 |

**明确不做/未做（防审稿追问时撒谎）：** 全面板 max vs median 数值表（仅有 27 配体 API 诊断样本，不得升格为 SI 表）；从供给池重抽 1000 个**互不重叠独立** panel 的 summary_min 分布（硬负供给不够，不是“未对接所以还没做”）；PDBFixer+Reduce 质子化流程；主面板残基级 PLIF（git 无 K=4 production poses）。GNINA 九姿态公平重打已于 2026-08-24 补做完成（见 2.7；`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`），不再列入未做项。本轮补写的 holdout / 换晶 / 叠合 / 接触计数参数均来自已冻结实验，不是新对接。必要性复核见 `data/jcim_strengthen_t0t1_v0/analysis/C_CLASS_EXPERIMENT_NECESSITY_VERDICT_V1.md`。
