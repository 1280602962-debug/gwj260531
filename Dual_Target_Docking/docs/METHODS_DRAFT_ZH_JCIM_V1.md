# Methods（中文工作稿 · JCIM Articles）

> 结构：预定义评价协议 + robustness / sensitivity / falsification，而不是实验记录压缩版。  
> 语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843（写清做什么、参数与软件；数字进 Results / SI）。  
> 配套：[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md)、[`INTRODUCTION_DRAFT_ZH_JCIM_V1.md`](INTRODUCTION_DRAFT_ZH_JCIM_V1.md)、英文稿 [`METHODS_SECTION_JCIM_EN_V1.md`](METHODS_SECTION_JCIM_EN_V1.md)。  
> **Methods 只写协议**；供给计数、cognate RMSD、AUROC、holdout 点估计等一律在 Results / SI。  
> **不编造** median pChEMBL 表、1000 次互不重叠 panel 重抽，也不重建已冻结面板。  
> DualFourClass-Bench 是 **four-state curated benchmark**；primary endpoint 是两条方向 pairwise AUROC，不是四分类器。

---

## 2. Methods

### 2.1 数据来源与活性数据整理

双靶评价所需的配体活性作为 **experimentally derived activity labels**，通过 ChEMBL Web API 的公开 activity 端点获取。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干摩尔浓度–响应型测定（如 IC50、Ki、Kd、EC50）转换为近似 −log10 活性尺度，便于大规模公开数据整合。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似，而不解释为同一条件下可直接比较的绝对结合亲和力。

同一配体–靶标若有多条可用 pChEMBL 记录，冻结数据包采用**最大 pChEMBL** 作为一对一代表值。该聚合保持覆盖率，但可能引入 activity inflation，其作为数据来源限制在 Discussion 中讨论。冻结文件（`mols_*.json`）仅保存该代表浮点数，不具备在本地重算 assay 中位数、ChEMBL 置信度阈值或物种过滤的条件；本文因此不在冻结包上重建 median / confidence≥8 / *Homo sapiens* 活性表，以免在面板冻结后改写标签定义。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。

ChEMBL 结构常含盐、溶剂化物或多组分形式。对接前按连通片段拆分，并保留重原子数最多的有机片段作为计算母体。

为检验 ChEMBL 供给门槛是否仅为单一库的覆盖假象，对进入冻结评价集的靶对另做 **BindingDB / PubChem 计数核对**（零对接、不重建面板）。BindingDB 使用 REST `getLigandsByUniprots`（cutoff = 1 mM，以免截掉弱端测定）；PubChem 使用 PUG REST `protein/accession/…/concise`。类型限于 IC50/Ki/Kd/EC50；代表值取最大转换 p 活性；分类规则与 2.2 的严格供给门槛相同。配体身份分别用 BindingDB monomerid 与 PubChem CID，**不做**跨库结构合并。主比较采用去掉 `>`/`<` 截尾的**等式测定**；将不等式当作点估计只作敏感性。该核对只报告计数（Supporting Information Table S12）。

### 2.2 靶对供给审计与实验配体状态定义

为判断公开数据能否支持严格的双靶 benchmark，先对候选靶对做数据供给审计。对每一对靶标 A/B，按两端实验活性将配体定义为四种**实验状态**：

- **dual**：两端均有较强活性；
- **A-only**：仅 A 端较强，B 端较弱；
- **B-only**：仅 B 端较强，A 端较弱；
- **neither**：两端均缺乏足够活性。

A-only 与 B-only 是选择性硬负样本，不是 DUD/DUD-E 式假定 decoy。

**严格供给审计规则（construction gate，不是全部最终比较的唯一标签）。** Dual：两端 pChEMBL ≥ 6.5。A-only：A ≥ 6.5 且 B ≤ 5.5。B-only 对称。Neither：两端 ≤ 5.5。介于 5.5 与 6.5 的灰区不进入严格审计。该规则用于判断某一靶对在两个方向上是否具有足够的选择性硬负，以支持规模较均衡的面板。审计通过门槛与最终进入评价集的靶对名单见 Results 3.1；金属依赖体系（如 HDAC）按预先声明排除，不作为常规非共价对接主对象。

**正文主比较采用预先统一的 θ = 6.0 标签。** Dual：两端 ≥ θ；A-only：A ≥ θ 且 B < θ；B-only 对称；neither：两端 < θ。建造阶段允许在严格规则下单端选择性过少时改用该单阈值规则凑齐配额；建造规则在抽样前按供给审计冻结，并写入 Table 1。阈值选择服务于可分析配额，不是在观察对接分数后回改标签。作为支持性敏感性分析，在 θ ∈ {5.5, 6.5} 与严格 6.5/5.5 规则下重标四种状态并重算口袋匹配 summary_min（Table S4）。该网格不是与 Table 2 竞争的第二套主标准。样本量过小的格子在 Results 中标记 underpowered，Methods 不预判其数值。

### 2.3 DualFourClass-Bench 面板构建

**The benchmark preserves four experimentally defined ligand states—dual, A-only, B-only, and neither—whereas the prespecified primary endpoint focuses on two directional pairwise discrimination tasks, dual versus A-only and dual versus B-only.** neither 保留以描述完整四状态实验空间及后续辅助分析，不进入 primary directional AUROC。

候选靶对按 2.2 的严格供给审计筛选。最终冻结评价集包含 PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB 与 EGFR/HER2。EGFR/HER2 按预先批准的角色保留为**供给受限案例**（`PAIR_ROLES_APPROVED_JCIM.yaml`），其组成不与其余靶对按同一厚面板供给条件等价。

每个靶对从符合相应实验标签的候选池中按预先冻结的类别配额抽样。面板抽样使用固定随机种子 20260729。在能够计算 Bemis–Murcko 支架的面板上施加支架封顶，以降低同一化学系列过度代表：PIK3CA/mTOR（PM48）同一类别内同一支架最多 2 个分子；EGFR/HER2 最多 5 个。AChE/BChE 与 PIK3CA/PIK3CB 在建面时 SMILES 尚未并入抽样表，无法施加 Murcko 封顶；实际抽样按类别配额与确定性随机顺序进行，并以 ChEMBL identifier 前缀作为当时的占位过滤。该前缀**没有化学意义，不解释为化学多样性约束**。各面板的最终成员、状态标签、ChEMBL identifier、SMILES、Murcko 支架（若事后可算）与抽样脚本随冻结数据包提供；本文不在观察对接分数后重抽面板。

配额与建造标签如下。AChE/BChE 与 PIK3CA/PIK3CB：严格 6.5/5.5，目标 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2：沿用既有 θ = 6.0 面板（n = 110）。PIK3CA/mTOR：θ = 6.0，主比较面板 PM48（n = 48；建造 dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。

对接失败的配体–受体组合从该受体分数中剔除；任一端缺少可用分数的配体不进入需要两端分数的口袋匹配 AUROC，故分析用计数可低于建造定额（Table 1）。

PIK3CA/mTOR 另构建扩面面板（历史名 PM110）：保留 PM48 全部 48 个配体，并按严格规则追加分子，目标配额 dual / A_only / B_only / neither = 30 / 30 / 30 / 25。PM110 是 PM48 的超集，用于评价面板规模增加后点估计是否同向，不是与其他靶对独立等价的 primary benchmark，也不是独立重复实验。主文跨对比较以 PM48 为准。

本文不以从供给池重复抽取互不重叠平衡面板的分布作为稳健性读出。该路径受硬负供给限制（定量见 Results）；正式的配体侧外推是一次 unused-pool holdout（2.11）。配体层有放回 bootstrap（2.8）描述固定面板内的不确定度，不称作供给池重抽。

**Table 1.** DualFourClass-Bench 评价集组成与对接设置（建造规则）

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.4 蛋白结构与结合位点定义

受体取自 Protein Data Bank 中含实验确定结构与小分子共晶配体的条目。冻结主分析使用：PIK3CA/mTOR，4L23 / 4JT6（共晶配体 X6K / PI-103）；AChE/BChE，4EY7 / 4BDS（E20 / THA）；PIK3CA/PIK3CB，4L23 / 2WXF（X6K / 039）；EGFR/HER2，3POZ / 3RCD（03P / TAK-285）。分辨率见表 1。

结合位点由各结构的共晶配体定义。以共晶配体重原子坐标计算轴对齐包围盒（AABB），三方向各外扩 5 Å；任一边若小于 20 Å，则将该边设为至少 20 Å。盒子中心与边长冻结于 JSON，并汇总于 Supporting Information Table S2。

受体准备时去除水分子与共晶配体，再用 Meeko 生成 PDBQT。PIK3CA、mTOR、EGFR 与 HER2：使用冻结目录中已含氢的蛋白坐标，经 `mk_prepare_receptor.py --read_pdb` 转换。AChE、BChE 与 PIK3CB：从沉积 PDB 的 ATOM/TER 记录提取蛋白（去除水与异源原子），以 `mk_prepare_receptor`（默认 alternate location A）转换。未额外运行 PDBFixer 补全缺失原子，也未用 Reduce 做独立的 pH 依赖质子化或组氨酸互变异构枚举；质子化属于冻结准备协议的一部分。主分析均为非共价小分子对接，不在盒子中把金属离子或其他辅因子当作额外可对接组分。完整命令与输入文件随公开数据包提供。

### 2.5 共晶配体重对接质量控制

正式对接前，对每个冻结受体做共晶配体重对接，以检验对接盒子、受体准备与搜索参数能否在**保留的姿态集合**中生成近似共晶构象。

每个共晶配体生成 9 个 docking poses，计算其与实验共晶构象的重原子 RMSD（对接坐标系，不做蛋白叠合）。PIK3CA/mTOR 与 EGFR/HER2 使用 meeko `REMARK SMILES IDX` 映射后，在图自同构上取最小 CalcRMS；AChE/BChE 与 PIK3CB 使用重原子匈牙利匹配（`linear_sum_assignment`）。定义

\[
\mathrm{RMSD}_{\mathrm{best9}} = \min_{i=1,\ldots,9} \mathrm{RMSD}_i.
\]

预先通过标准为 \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å。

该 QC 检验的是 **pose-generation capability**：协议能否在保留的 pose ensemble 中产生近晶构象。它**不等于**要求 Vina 排名第一的 pose（mode 1）必须为近晶构象。best-of-9 QC 与 mode-1 scoring 是不同层面的评价。

若默认 exhaustiveness 下未满足预设 QC，则在不改变盒子、受体与随机种子的条件下，将搜索强度提高至预先规定的备用水平并重新 QC。主分析因此采用受体特异的冻结 exhaustiveness：PIK3CA/mTOR 为 16，其余主面板为 8。各受体的 QC 数值、mode-1 与 best-of-9 对照见 Supporting Information Table S3。

### 2.6 配体准备与分子对接

配体从冻结 ChEMBL SMILES 统一准备：去盐并保留最大有机片段，RDKit 加显式氢，ETKDGv3 生成三维构象（随机种子 20260727），MMFF 局部优化最多 200 步，再经 Meeko 默认参数转为 PDBQT。不进行系统性的质子化状态、互变异构体或构象枚举；各靶对使用同一 ligand-preparation protocol。不使用 Schrödinger LigPrep。

分子对接采用 AutoDock Vina 1.2.7，默认 `vina` 打分函数。每个配体–受体组合生成 9 个 poses，`energy_range = 3` kcal mol\(^{-1}\)，随机种子 20260727（与 ETKDG 相同）。exhaustiveness 按 Table 1 的受体特异冻结值。配体准备、盒子生成规则与打分函数在各主面板上相同；仅受体坐标、盒子数值与预先定义的 exhaustiveness 按靶标变化。完整参数见 Supporting Information Table S1。

### 2.7 Alternative scoring channels

为检验主观察是否依赖单一打分函数，在**同一组 Vina-generated poses** 上另用 RTMScore 与 GNINA CNN 重打分。

RTMScore 使用公开权重 `rtmscore_model1`，对每个配体–受体组合的 9 个 Vina poses 分别打分，取该口袋最高 RTMScore。

GNINA 1.3.2 在 CPU 模式下做 CNN rescoring（`--cnn_scoring rescore --minimize`）。最终协议对全部 9 个 Vina poses 分别转 SDF（Open Babel）并重打分，取每端最高 CNNscore，与 RTM 的姿态覆盖对齐。仅使用 Vina mode 1 的 GNINA 结果保留为历史 sensitivity control，不是最终通道读出。

Vina 主读出是 mode 1 能量；RTM 与 GNINA 是 best-of-9 重打分。三者对 9 个姿态的聚合并不相同，因此 **不作为 head-to-head docking-engine competition**，而作为 scoring-channel sensitivity analysis。Primary endpoint 始终由 Vina 定义。

### 2.8 Primary endpoint 与统计分析

#### 2.8.1 口袋匹配方向 AUROC

对每个靶对 A/B 计算两条二分类 AUROC。dual 对 A-only 使用**口袋 B** 的分数：

\[
\mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B),
\]

以检验对接能否利用非选择性靶点 B 的结构信息，把 dual-active 与已在 A 端强效的 A-only 分开。dual 对 B-only 使用口袋 A 的分数：

\[
\mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A).
\]

dual 始终为正类。neither 不进入上述对比。

Vina 输出结合能 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，通常越负表示预测结合越强）。定义

\[
S_{\mathrm{Vina}} = -E_{\mathrm{Vina}},
\]

使所有 primary scores 遵循“越大表示预测结合越强”。RTMScore 与 GNINA CNN 分数本身已是越高越好，不再取负。

#### 2.8.2 summary_min

靶对汇总为较弱一臂：

\[
\mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}).
\]

该规则是与双靶任务同构的 **worst-arm aggregation**，不是新的 scoring function。选择最小值是为了避免一端较强的 discrimination 掩盖另一端失败。全文只有一个主终点：统一 θ = 6.0 下的口袋匹配 Vina `summary_min`（Table 2；PIK3CA/mTOR 主面板为 PM48）。预指定次级终点为两条方向臂、RTMScore 口袋匹配、GNINA CNN best-of-9 口袋匹配，以及 2.8.3 的描述符面板。稳健性 / 证伪终点为 θ 网格、PM110、E = 8、unused-pool holdout、受体替换、错口袋对照（含配对 Δ）。探索性终点为 ECFP4、contact_count（非 PLIF）以及 pooled `vina_mean` 的 Top-10 硬负计数。完整层级见 Supporting Information Table S16。`vina_mean` 池化方向 AUROC **不是** Table 2。

#### 2.8.3 物化描述符对照

用 RDKit 计算预先指定的描述符面板：重原子数（GetNumHeavyAtoms）、分子量（MolWt）、cLogP（MolLogP）与 TPSA。每个描述符按与对接分数相同的方向 AUROC 流程评价，正文与 SI 报告全部四个。其中 AUROC 最高者记为 **best-performing physicochemical descriptor**，只作为描述性强基线，**不是**预先指定的独立假设检验基准。为避免先选最优描述符再做正式比较的选择偏倚，docking 与描述符的配对 Δ 不以“击败 best descriptor”作为 confirmatory test，而用于描述对接是否明显超出简单配体属性信号（Table S19）。

#### 2.8.4 分数聚合对照

作为辅助分析，同时计算两端分数的 pooled mean、wrong-pocket assignment（定义见 2.9.1）以及 worst-pocket aggregation。它们不是 primary endpoint，只用于判断不同聚合是否改变双靶识别结论（Table S6）。

#### 2.8.5 Bootstrap 不确定度

AUROC 与 summary_min 的不确定度用配体层 bootstrap：在保持类别标签结构的条件下对配体有放回重采样，每次重算两条方向 AUROC 与 summary_min。\(B = 2000\)，随机种子 20260729，百分位数 95% CI 为 \([P_{2.5}, P_{97.5}]\)。错口袋与描述符等配对比较在**同一次**重采样上计算 \(\Delta = \mathrm{Metric}_1 - \mathrm{Metric}_2\)，得到 paired bootstrap 区间（Table S17、S19）。另报 Murcko 支架重采样区间作为对照；正文以配体层为准。置信区间作描述性不确定度；除预先定义的主终点外，不对多靶对、多对照做多重比较意义上的 confirmatory testing，也不把“CI 是否跨越 0.5”单独等同于正式显著性。

### 2.9 Confounder、falsification 与化学对照

#### 2.9.1 Wrong-pocket falsification control

将靶点 A 与 B 的分数对调，配体、受体与其余分析设置不变，重算方向 AUROC 与 summary_min。该分析是 **falsification control**：检验观察到的 discrimination 在错误口袋赋值下是否仍然保持。错口袋明显低于匹配口袋，可提供 pocket-specificity 的支持性证据；错口袋接近或高于匹配口袋，则视为对 pocket-specific interpretation 的反证。它不是用来“证明口袋特异”的阳性对照。

#### 2.9.2 配体效率归一

各口袋分数除以重原子数，\(S_{\mathrm{LE}} = S_{\mathrm{dock}} / N_{\mathrm{heavy}}\)，再按 primary 流程计算方向 AUROC 与 summary_min，以检验对接分是否主要反映分子大小。

#### 2.9.3 效价与尺寸匹配子集

分别构建 \(|\Delta\mathrm{pChEMBL}| \leq 0.5\) 的 potency-matched 子集与 \(|\Delta N_{\mathrm{heavy}}| \leq 2\) 的 size-matched 子集，在子集上重算方向 AUROC。匹配会减小样本量；该分析只判断方向是否明显改变，不把低样本量子集的点估计当作独立强证据（Table S5）。

#### 2.9.4 Covariate-adjusted analysis

逻辑回归比较

\[
\mathrm{Model}_1:\ Y \sim S_{\mathrm{dock}}, \qquad
\mathrm{Model}_2:\ Y \sim S_{\mathrm{dock}} + N_{\mathrm{heavy}} + \mathrm{TPSA},
\]

其中 \(Y\) 为 dual 对相应选择性硬负的二分类标签。使用 scikit-learn `LogisticRegression`（\(C = 1.0\)，`max_iter = 2000`）。报告模型 AUROC、对接分数回归系数及其优势比（OR）。该分析询问对接分在控制分子大小与极性后是否仍有 residual discrimination，协变量模型不是 primary predictor。

#### 2.9.5 二维化学基线

Morgan/ECFP4（半径 2，2048 bit）加与 2.9.4 相同的逻辑回归，建立仅依赖二维结构的基线。评价采用 Bemis–Murcko scaffold `GroupKFold`，折数 \(\min(5, N_{+}, N_{-}, N_{\mathrm{scaffold}})\) 且至少两折，使同一骨架不跨训练/测试折。该分析识别 chemotype–label association，不是口袋物理证据。随机 `StratifiedKFold` 仅作泄漏核对（Table S20），不以寻找更大 gap 为目的。

#### 2.9.6 Scoring-independent contact count

在已冻结的 Vina **mode-1** 姿态上计算不依赖打分函数的几何量：配体重原子中与受体重原子距离 ≤ 4.0 Å 的原子数

\[
N_{\mathrm{contact}} = \#\{i:\ \min_j d_{ij} \le 4.0\,\text{Å}\}.
\]

该描述符不使用对接能量函数。用 \(N_{\mathrm{contact}}\) 在口袋 A 上比较 dual 对 A-only、在口袋 B 上比较 dual 对 B-only，与错口袋对照的同口袋比较同构，作为 scoring-independent geometric confounder control，检验错口袋判别是否可能只反映更大分子产生更多埋藏接触。4.0 Å 为粗粒度接触阈值，**不是** PLIF。不预设其幅度与 Vina 错口袋一致（Table S11）。

#### 2.9.7 跨对序列一致性（探索性）

从各冻结受体 `*_protein.pdb` 用 Biopython `PDBParser` 提取最长蛋白链一级序列（仅标准氨基酸 ATOM），以 `PairwiseAligner`（BLOSUM62，全局比对，gap open = −11、extend = −1）计算靶对内全链序列一致性（分别以比对长度与较短链归一；Table S7）。该指标是整体相似度的粗粒度代理，不涉及口袋残基对应或结构叠合，不用于口袋 RMSD 或 PLIF 主张。

### 2.10 单靶富集参照

在 PIK3CA 4L23 与 mTOR 4JT6 上分别构建单靶 active–weak-active 集合。活性分子：pChEMBL ≥ 6.5。弱效分子：同靶已测定且 pChEMBL ≤ 5.5，并按分子量 ±50 Da、cLogP ±1.5、TPSA ±25 Å² 与活性分子做性质匹配。分子量与 logP 窗口沿用 property-matched decoy 的常见设定（Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594）；TPSA 窗口为同一思想下增加的极性匹配。目标规模约 50 个活性分子与 150 个弱效分子。配体准备、受体、盒子与 Vina 协议与 PIK3CA/mTOR 主面板相同（exhaustiveness = 16）。报告 AUROC、EF1% 与 EF5%。该实验只提供单靶 docking enrichment 的背景参照，不替代 dual-target 的 summary_min。

### 2.11 Unused-pool holdout

为检验结论是否依赖于冻结面板的具体成员，从严格标签池中排除已用于主面板与 PM110 的 ChEMBL 条目，在剩余 unused pool 中构建 **unused-pool, panel-external holdout**。它不是跨数据库或跨实验体系的 external validation：配体仍来自同一 ChEMBL 抓取批次、同一靶对与同一标签规则。

Holdout 只在 unused-pool 配额足以按 dual / A-only / B-only 各抽 20 个配体的靶对上构建。预先冻结为 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB；EGFR/HER2 不具备同等未用池抽样条件，记为 not eligible，不补做不等价抽样。PIK3CA/mTOR 排除的是 PM110 超集，从而覆盖 PM48。抽样种子 `HOLDOUT_SEED = 20260731`（区别于建面种子），同一 Murcko 支架在每一状态类中最多 3 个成员。清单先冻结、后对接。

Holdout 不参与主面板构建、对接协议调整或 primary endpoint 选择。受体、盒子、配体准备、exhaustiveness、打分与统计与主 benchmark 相同，并使用同一 `summary_min` 与配体层 bootstrap。未能产生 Vina 分数的配体–受体组合按 2.3 从需要该分数的分析中剔除。同一 holdout 配体上并列计算描述符对照；错口袋、效价/尺寸匹配与 contact count 按 2.9 在 holdout 上重算（Table S8、S13）。效价/尺寸匹配诊断不改写 Table S8 的主 holdout 数字。

### 2.12 Receptor-structure sensitivity analysis

为评价 benchmark 结论对受体结构选择的敏感性，另选满足以下**预先声明**条件的替代晶体：（i）polymer entity 与目标蛋白真实对应，排除嵌合体或非目标同源骨架；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.5 相同的 cognate redocking QC。实际进入对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。该分析是 **receptor-structure sensitivity analysis**，不是用来证明某一晶体“更正确”，也不是把 PIK3CA/mTOR 预设为结构不变的 positive case。

替换采用**单口袋**设计：每次只替换靶对中的一个受体，另一端保持主 benchmark 的冻结结构与既有分数。4JPS/5DXT 替换口袋 A（PIK3CA），口袋 B 仍用冻结 4JT6 分数；4JSX 替换口袋 B（mTOR），口袋 A 仍用冻结 4L23 分数。新盒子按该替代晶体自身共晶配体、以 2.4 的同一 AABB 规则生成。配体准备、exhaustiveness（16）、随机种子（20260727）、打分函数与 primary endpoint 与 PM48 主分析相同。仅在冻结 PM48 配体上重对接并重算 summary_min。

作为探索性、零新对接的几何对照，在已冻结晶体坐标上做刚体叠合：Biopython `PDBParser` 提取最长链 Cα，按残基编号与残基名精确匹配，`Superimposer` 一次 Kabsch 拟合得全域 RMSD；口袋残基由参考结构共晶配体重原子 ≤5 Å 界定，在**同一变换**下计算口袋局域 RMSD，不做二次局部拟合。再将替代结构共晶配体按同一变换投影，计算与参考共晶配体质心的距离。不同结构匹配的 Cα 数目可以不同，全域 RMSD 因此不是等覆盖比较。本对照仅含有限数目的替代晶体，不预设 Cα RMSD 能够定量解释 AUROC 变化（Table S10）。

### 2.13 软件与数据可用性

计算在 Python 3 环境下完成。主要软件：RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2、RTMScore（`rtmscore_model1`）；Vina 姿态转 SDF 使用 Open Babel。刚体叠合与全链序列比对使用 Biopython（`PDBParser`、`Superimposer`、`PairwiseAligner`）。AUROC、逻辑回归与交叉验证使用 NumPy、SciPy、scikit-learn 与 pandas（版本见公开复现环境）。评价面板、对接分数、分析脚本与完整参数表随公开数据包提供，见 Data and Software Availability。

---

## 写法说明（不进正文）

| 原则 | 本稿处理 |
|------|----------|
| Methods 不含结果数字 | 49→4、cognate 谁失败、GNINA Δ、HOAP_028、holdout AUROC、换晶 0.486 等全部在 Results / SI |
| four-state ≠ 四分类 | 2.3 英文定位句；primary 只用 dual / A-only / B-only |
| summary_min | worst-arm 任务约束，不是 novel metric |
| best descriptor | 四描述符全报；max 只作 descriptive strong baseline |
| wrong-pocket | falsification control |
| holdout | unused-pool, panel-external；不是 external validation |
| 换晶 | receptor-structure sensitivity，不是“验证稳健” |
| identifier 前缀 | 如实写为 SMILES 缺失时的占位过滤；**不**称为多样性；**不**重抽已冻结面板 |
| cognate QC | pose-generation，不是 top-ranked pose recovery |
| RTM/GNINA | scoring-channel sensitivity，不是三引擎竞赛 |
| PIK3CA/mTOR | Methods 不预设 positive case |

**Methods → Results 对应**

| Methods | Results |
|---------|---------|
| 2.1 策展 + 跨库计数协议 | 3.1 供给计数（含 Table S12） |
| 2.2 状态定义与 θ 网格协议 | 3.2 标签稳健 |
| 2.3 面板建造 | 3.1 Table 1 |
| 2.4–2.5 受体与 cognate QC 协议 | 3.5 / Table S3 |
| 2.6–2.7 对接与通道 | 3.3 |
| 2.8 主终点 | 3.3 Table 2 |
| 2.9 混淆 / 证伪 / contact | 3.4、3.9 |
| 2.10 单靶富集 | 3.5 |
| 2.11 unused-pool holdout | 3.9 |
| 2.12 receptor-structure sensitivity | 3.10 / 3.11 |

**明确不做/未做：** 全面板 max vs median 数值表；1000 个互不重叠独立 panel；PDBFixer+Reduce；主面板残基级 PLIF。GNINA 九姿态公平重打已完成。必要性复核见 `C_CLASS_EXPERIMENT_NECESSITY_VERDICT_V1.md`。
