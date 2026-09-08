# 基于对接的双靶识别：多靶对的评价设定审计

## 摘要

对接常被用来解释双靶识别，但负类通常是非结合配体或未匹配 decoy，而不是实验测定的单靶选择性配体。本文将 DualFourClass 建成该选择的四状态配方审计。最初的主评价集是 2026-07-23 的四靶对冻结。事后 ChEMBL 37 宇宙普查之后，五个通过预先声明门槛的普通非共价靶对被加入同一分析栈——这是收集补全，不是预先注册的八靶对冻结。原先对接的第四个候选靶对 PIK3CA/PIK3CB 在受体身份审计中被撤回——其"PIK3CB"受体（PDB 2WXF）实为小鼠 PIK3CD，并非人源 PIK3CB，因此仅作为已记录的受体身份错误在 Supporting Information 中报告。主表因此为八行：PIK3CA/mTOR、AChE/BChE、供给受限的 EGFR/HER2，以及 F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA、PPARA/PPARD。在 AutoDock Vina 下，实验定义的双靶配体相对 A-only 与 B-only 选择性配体排序。最弱臂 AUROC 分别为 EGFR/HER2 0.430、AChE/BChE 0.606、PIK3CA/mTOR 0.692、F2/F10 0.345、JAK1/TYK2 0.365、JAK1/JAK2 0.588、PPARG/PPARA 0.649 与 PPARA/PPARD 0.446。仅 PPARG/PPARA 的 95% 置信区间完全高于 0.5，且该区间对面板成员与打分公式敏感。同一套 EGFR/HER2 分数在双靶对 neither 时 AUROC 为 0.756，定向评价则降至 0.430；JAK1/TYK2 再现同一配方差距（0.770 对 0.365）。独立 GNINA 姿态生成再现了这两处差距。支架分组的纯配体模型已捕获大部分表观排序，替换备选 PIK3CA 晶体结构则使 PIK3CA/mTOR 对比降至接近随机水平。双口袋分数过滤因此需要选择性感知的负类；两个口袋同时给出有利分数，本身并不是双靶识别的证据。

**关键词：** 双靶对接；选择性；硬负例；AutoDock Vina；GNINA；虚拟筛选

## 1. 引言

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多药理学（polypharmacology）研究的重要方向。[1] 近年来，多靶点小分子的理性设计逐渐由经验性筛选转向结合结构生物学、计算化学与生成式模型的结构导向设计。[2] 分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的构象，再用打分函数对配体–受体互补性排序。[3][4] 因此，在双靶点药物发现中，一个自然策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在双靶结合能力。对接结果的解释高度依赖数据集构建。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观富集退化为粗粒度配体性质分离。[5][6] LIT-PCBA 采用实验 assay 标签并控制已知 decoy 与化学偏倚。[7] CASF-2016 评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

一个严格的双靶评价需要区分 **dual-active**、**A-selective**、**B-selective** 与 **neither** 四种实验状态（Figure 1A）。A-only 与 B-only 是该任务的**选择性硬负样本（selectivity hard negatives）**：它们在一个靶点上已有较强活性。计算终点因而检验对接能否在两个方向上将 dual-active 与对应单靶选择性配体区分开。Zhou、Li 与 Hou 曾在四对激酶上评价相对非抑制剂的 dual-target docking。[9] 本文在该设定上引入实验定义的方向性硬负，并在同一套分数上比较不同基准设定。Dual versus neither 按实验 inactive 计分，作为基准设定对照。平衡的四状态面板还受两端可比较测量与双向选择性硬负供给的约束。

近期工作已经把双靶问题做实，但回答的并不是同一评价问题。Wu 等表明，大规模对接可以前瞻性地得到选定靶对的联合结合配体，同时也报告后续优化仍然困难。[19] 该研究问的是对接能否找到双靶活性分子；本文问的是：当负类由实验选择性配体而不是非结合配体定义时，回顾性双靶识别证据是否改变。POLYGON 在十万级结合数据上生成双靶化学空间，并合成 32 个 MEK1/mTOR 化合物做实验验证，[20] 因此不宜把近期生成式多药理学一概写成对接成功指标。生成式双靶方法也使用相对参考配体的对接度量。[10][11] Kinase-Bench 汇集 6875 个选择性配体、75 个激酶与 422,799 个 decoy，检验相对激酶特异性 decoy 的选择性富集；[22] DualFourClass 则在每一对上直接构造实验测定的 dual、A-only、B-only 与 neither 四状态。一项 147 靶的 AI 对接基准进一步表明，方法排序取决于负类是实验低活性 TrueDecoy 还是商业库随机 decoy，[21] 说明配方本身就是科学主张的一部分。

本文要问的是：基准设定是否改变双靶识别的表观证据。我们构建 DualFourClass-Bench，作为具有两条方向主任务的四状态面板：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B），并以二者较弱一臂汇总为最弱臂 AUROC（`summary_min`）。最初的主评价集是 2026-07-23 的四靶对冻结；普查合格的五个普通靶对随后加入同一分析栈，PIK3CA/PIK3CB 被撤回，主表为八行。我们进一步考察该判别是否能够在不同配体、活性聚合方式及受体结构条件下保持。这是收集补全，不是“八对在看分之前已冻结”的主张。

## 2. 方法

### 2.1 数据与实验状态定义

双靶评价所需的配体活性通过 ChEMBL Web API 的公开 activity 端点获取，作为实验衍生标签。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干经标准化的定量效力或亲和力测量（如 IC50、EC50、Ki、Kd 和 Potency）转换为近似 −log10 活性尺度。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似。

同一配体–靶标若有多条可用 pChEMBL 记录，主策展采用**最大 pChEMBL** 作为一对一代表值。活性聚合敏感性分析从 ChEMBL activity 端点重拉 assay 级记录，并在同一 θ = 6.0 规则下用重复测定的**中位数**替换最大值，不改变面板成员、对接参数或 Vina 分数。API 重拉的最大值对中位数估计作为标签聚合敏感性，与 Table 2 并列报告（Table S29）。作为 post-hoc 标签稳健性，已打分配体–靶标记录于 2026-08-26 从当前 ChEMBL API 重拉，仅保留 Homo sapiens `SINGLE PROTEIN`、assay confidence ≥8、等式关系、IC50/Ki/Kd/EC50/Potency、无 validity 注释且无 `potential_duplicate` 的记录，再按 θ = 6.0 重标而不改变面板成员或对接分数（Table S36）。该日期化视图检验这些显式过滤器是否改变冻结标签，不是 2026-07-23 数据库状态的重建，也不统一 assay 条件、蛋白构建体或突变背景。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。ChEMBL 结构按连通片段拆分，保留重原子数最多的有机片段。

完整病例选择按冻结可用 pChEMBL 映射统计只在 A、只在 B 或两端均有值的结构（Table S37）。缺失不解释为无活性。来源文献集中度按高置信保留记录、分实验状态统计独立 `document_id` 与最大单篇文献份额。

对每一对靶标 A/B，配体被定义为四种实验状态：**dual**（两端较强）、**A-only**（仅 A 端较强）、**B-only**（仅 B 端较强）和 **neither**（两端均不足）。A-only 与 B-only 是选择性硬负样本。

**严格 6.5/5.5 规则**为：dual，两端 pChEMBL ≥ 6.5；A-only，A ≥ 6.5 且 B ≤ 5.5；B-only 对称；neither，两端 ≤ 5.5。5.5–6.5 灰区不进入该审计。金属依赖体系（如 HDAC）预先排除。严格 6.5/5.5 规则仅用于靶对供给资格审定。全部主基准标签随后统一按 θ = 6.0 定义（dual，两端 ≥ θ；A-only，A ≥ θ 且 B < θ；B-only 对称；neither，两端 < θ），并在查看对接结果之前冻结。两条阈值因此服务于预先规定的不同目的。建造规则在抽样前按供给审计冻结（Table 1）。θ ∈ {5.5, 6.5} 与严格 6.5/5.5 重标作为敏感性报告（Table S4）。样本量过小的格子在 Results 中标记效能不足。

为核对 ChEMBL 供给门槛，对冻结靶对另做 BindingDB / PubChem 计数核对（零对接、不重建面板；Table S12）。类型限于 IC50/Ki/Kd/EC50；配体身份分别用 BindingDB monomerid 与 PubChem CID，不做跨库结构合并。主比较采用等式测定。

随后从 BindingDB 202608 文章与专利 TSV 归档重建原生切片，而不是从 Table S12 的 REST pmax JSON 反推。[16] 规则预先写在 `external_slice_contract.yaml`：BindingDB 策展的文章或专利；人源野生型单链 UniProt；等式 IC50/Ki/Kd；两端均有测定；配体–靶–endpoint 内取中位数；θ = 6.0 四状态；去掉与开发面板共享的 PMID/DOI/专利；去掉已打分面板、未使用池留出集、PM110 或该对 ChEMBL 图中的 InChIKey/ChEMBL ID；对开发分子的最大 ECFP4 Tanimoto < 0.70。主外部门槛为 dual/A-only/B-only 各 n ≥ 20、每类至少 3 个来源、最大单文献配体份额 ≤ 50%。因没有任何靶对满足预先规定的外部切片门槛，故未进行外部切片对接。ChEMBL 文献解析为 519/680，因此剩余计数是完全文献独立集的上界。该切片不作为外部验证。

### 2.2 基准构建

DualFourClass-Bench 保留四种实验状态。主分析由两条方向性成对任务组成。neither 类保留给描述性的基准设定对照。

候选靶对按 2.1 的严格供给审计筛选。2026-07-23 原始冻结包含 PIK3CA/mTOR、AChE/BChE、EGFR/HER2 与 PIK3CA/PIK3CB。EGFR/HER2 保留为供给受限案例。PIK3CA/PIK3CB 同样通过供给筛选并完成对接，但事后受体身份审计发现其"PIK3CB"口袋所用 PDB 条目（2WXF）实际对应小鼠 PIK3CD（UniProt O35904），并非人源 PIK3CB（P42338，尚无沉积晶体结构）；共晶重对接 QC 无法识别该错误，因为共晶配体本身就属于这个（错误的）蛋白。该靶对因此从主评价集撤回，仅作为已记录的受体身份错误在 Supporting Information 中报告，而不作为经修正的常规结果。

事后 ChEMBL 37 人源 SINGLE PROTEIN 宇宙普查之后，五个通过预先声明 G1–G5 门槛的普通非共价靶对——F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA、PPARA/PPARD——被加入同一分析栈。这是收集补全，不是预先注册的八靶对冻结。CTSK/CTSS 在普查上可评，但因可逆共价晶体不进入普通 Vina。配体按预先冻结的类别配额抽样，随机种子为 20260729。抽样时结构可用的面板施加 Bemis–Murcko 支架封顶：PIK3CA/mTOR（PM48）同一类别内同一支架最多 2 个分子，EGFR/HER2 最多 5 个，普查后五对最多 3 个。AChE/BChE 只采用类别配额和确定性随机顺序。观察对接分数后不再重抽面板。

AChE/BChE 按严格供给门槛抽样（目标 28 / 28 / 28 / 16；n_panel = 100）。EGFR/HER2（n_panel = 110）与 PIK3CA/mTOR PM48（n_panel = 48；建造 18 / 14 / 12 / 4）按主分析 θ = 6.0 标签构建。普查后五对均用 θ = 6.0 配额 32 / 32 / 32 / 14（n_panel = 110；Table 1）。因此跨对 AUROC 同时混合靶对生物学与面板构建差异。对接失败的配体–受体组合被剔除，n_scored 可低于 n_panel（Table 1；Table S27）。扩面面板 PM110 保留 PM48 全部 48 个配体，用于检查面板规模增加后点估计是否同向。

**Table 1.** DualFourClass-Bench 评价集组成与对接设置。建造标签记录各靶对的供给/建面规则；Tables 2–3 的全部主 AUROC 均使用统一 θ = 6.0 实验状态标签。n_panel 为冻结面板成员数（含 neither）；n_scored 为两端均有有效 Vina 分数、进入方向性主 AUROC 的 dual / A-only / B-only 计数。已撤回的 PIK3CA/PIK3CB（受体身份错误；见 2.2 节）仅见 Supporting Information。普查后五行在 2026-07-23 冻结之后加入。

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | n_panel | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |
| F2/F10 | θ = 6.0；普查后 | 4UDW / 2JKH | 1.16 / 1.25 | 110 | 31 / 32 / 32 | 8 |
| JAK1/TYK2 | θ = 6.0；普查后 | 6N7A / 3LXP | 1.33 / 1.65 | 110 | 31 / 32 / 32 | 8 |
| JAK1/JAK2 | θ = 6.0；普查后 | 6N7A / 8BXH | 1.33 / 1.30 | 110 | 32 / 32 / 32 | 8 |
| PPARG/PPARA | θ = 6.0；普查后 | 9V8H / 6LXA | 1.39 / 1.23 | 110 | 32 / 31 / 32 | 8 |
| PPARA/PPARD | θ = 6.0；普查后 | 6LXA / 5U3Q | 1.23 / 1.50 | 110 | 32 / 32 / 32 | 8 |

### 2.3 受体准备与对接协议

受体取自含小分子共晶配体的 PDB 条目：PIK3CA/mTOR，4L23 / 4JT6（X6K / PI-103）；AChE/BChE，4EY7 / 4BDS（E20 / THA）；EGFR/HER2，3POZ / 3RCD（03P / TAK-285）；F2/F10，4UDW / 2JKH（N6L / BI7）；JAK1/TYK2，6N7A / 3LXP（KEV / IZA，TYK2 JH1）；JAK1/JAK2，6N7A / 8BXH（KEV / C87）；PPARG/PPARA，9V8H / 6LXA（BRL / EPA）；PPARA/PPARD，6LXA / 5U3Q（EPA / 7UJ）。结合位点由共晶配体定义。以共晶配体重原子计算轴对齐包围盒，三方向各外扩 5 Å；任一边若小于 20 Å，则设为至少 20 Å（Table S2）。去除水分子与共晶配体后，用 Meeko 生成 PDBQT。PIK3CA、mTOR、EGFR 与 HER2 使用冻结目录中已含氢的蛋白坐标（`mk_prepare_receptor.py --read_pdb`）。AChE 与 BChE 从沉积 ATOM/TER 记录提取，并以 `mk_prepare_receptor`（默认 alternate location A）转换。主分析均为非共价小分子对接。

正式对接前对每个冻结受体做共晶配体重对接。生成 9 个姿态，计算与实验共晶构象的重原子 RMSD。预先通过标准为 \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å，即九个保留姿态中是否存在与共晶配体重原子 RMSD 小于 2.0 Å 的构象。若默认 exhaustiveness 未通过门槛，则提高至预先规定的备用水平。主分析因此采用 PIK3CA/mTOR exhaustiveness = 16、其余主面板为 8（Table S3）。EGFR/HER2 原始九姿态生产 PDBQT 未能找回，已按冻结协议重对接并标为 reconstructed QC，而非历史生产文件。拓扑核对后的 ranked RMSD：EGFR 3POZ top-1 9.505 Å、top-3 6.227 Å、best-of-9 0.760 Å（top-1/top-3 未过，搜索覆盖通过）；HER2 3RCD top-1 1.855 Å、top-3 1.394 Å（通过）。这些数值替换原先 NA 的 top-3 单元格，不改变预先规定的 best-of-nine 生产门槛。

配体从冻结 ChEMBL SMILES 统一准备：去盐并保留最大有机片段，RDKit 加显式氢，ETKDGv3 生成三维构象（种子 20260727），MMFF 局部优化最多 200 步，再经 Meeko 转为 PDBQT。不进行系统性质子化、互变异构或构象枚举。对接采用 AutoDock Vina 1.2.7 默认 `vina` 打分函数，保留 9 个姿态，`energy_range = 3` kcal mol\(^{-1}\)，随机种子 20260727（Table S1）。为检验打分函数依赖性，同一组 Vina 姿态另用 RTMScore（`rtmscore_model1`，取九姿态最高分）与 GNINA 1.3.2 CNN（`--cnn_scoring rescore --minimize`，Open Babel 转 SDF 后取九姿态最高分）重打分。Vina 主读出是 mode-1 能量；RTM 与 GNINA CNN 是 best-of-9 重打分。主终点始终由 Vina 定义。

另在配方差距靶对 EGFR/HER2、PIK3CA/mTOR 与 JAK1/TYK2 上以 GNINA 1.3.2 对接搜索模式独立生成姿态（不是对 Vina 姿态重打分），复用冻结 Meeko 配体 PDBQT、受体坐标、对接盒、exhaustiveness（分别为 8、16 与 8）、九个保留姿态和种子 20260727。读出为 mode-1 `minimizedAffinity`。两端均失败的配体从需要完整分数的分析中剔除（EGFR/HER2 的 neither 配体 EH120_109；PIK3CA/mTOR 的 A-only 配体 PM48_19；JAK1/TYK2 四个 CG0 解析失败与一次超时）。该协议检验设定效应在更换姿态生成引擎后是否仍在，不是多引擎比赛，也不扩到 F2/F10、JAK1/JAK2 或两个 PPAR 对（Table S32）。预先规定的五种子 Vina 敏感性复用生产种子 20260727，并增加 20260811–20260814 四个冻结种子，配体、受体、对接盒与搜索设置保持不变；Dual versus neither 使用与 Table 3 相同的逐配体 `vina_mean` 估计量（Table S54）。

### 2.4 主终点与统计分析

全文中的“双靶识别”指这一计算判别任务。对每个靶对计算两条二分类 AUROC。dual 对 A-only 使用口袋 B 分数，\( \mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B) \)；dual 对 B-only 使用口袋 A 分数，\( \mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A) \)。dual 始终为正类。Vina 输出结合能 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，越负表示预测结合越强）；\(S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}\)。

最弱臂 AUROC 定义为 \( \mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}) \)。后文使用“最弱臂 AUROC”；表格保留 `summary_min` 列名。它把两条方向 AUROC 保守地汇总为单值。算术平均、几何平均与调和平均作为聚合敏感性报告（Table S26）。全文唯一主终点是统一 θ = 6.0 下的口袋匹配 Vina 最弱臂 AUROC（Table 2；PIK3CA/mTOR 主面板为 PM48）。预先指定的 RDKit 描述符面板（重原子数、分子量、cLogP、TPSA）按同一方向流程评价；其中 AUROC 最高者记为最佳单一描述符参考（Tables 2、S28、S19）。Dual versus neither（实验 inactive；`vina_mean`）与 Dual versus all non-duals 为同一套冻结分数上的基准设定对照（Table 3；Table S22）。PIK3CA/mTOR 的 neither n = 4 标记效能不足。

Table 2 的不确定度取自 `unified_threshold_sensitivity_v2.csv` 中的配体层非分层百分位 bootstrap（\(B = 2000\)；由靶对与 \(\theta\) 经 SHA-256 派生的确定性子种子）。配对比较在同一次重采样上计算（Tables S17、S19）。置信区间作描述性不确定度。可分辨效应模拟使用观察得的类别样本量、同一 bootstrap、双正态分数模型和一组真实 AUROC，报告 95% CI 排除 0.5 的概率，而不是观察后功效（Table S31；Figure S6）。

### 2.5 混淆、留出集与受体敏感性分析

将靶点 A 与 B 的分数对调作为证伪对照，配体、受体与其余设置不变。另在配体效率归一（\(S_{\mathrm{dock}}/N_{\mathrm{heavy}}\)）、效价约束（\(|\Delta\mathrm{pChEMBL}| \leq 0.5\)）和尺寸约束（\(|\Delta N_{\mathrm{heavy}}| \leq 2\)）后重算方向 AUROC。逻辑回归比较 docking alone 与 docking + 重原子数 + TPSA。Morgan/ECFP4（半径 2，2048 bit）加逻辑回归在 Bemis–Murcko 支架 `GroupKFold` 下提供配体化学基线（Tables S5、S20、S23、S24）。最近邻 Tanimoto 子集只作诊断。接触计数与全链序列一致性仅为探索性对照（Tables S7、S11）。

为检验结论是否依赖于冻结面板的具体成员，排除已用于主面板与 PM110 的 ChEMBL 条目后，在剩余未使用配体池中构建留出集（holdout）。配体仍来自同一 ChEMBL 抓取批次、同一靶对与同一标签规则。该留出集在 PIK3CA/mTOR、AChE/BChE 与普查后五对上构建（各 20 dual / 20 A-only / 20 B-only；`HOLDOUT_SEED = 20260731`；JAK1/JAK2 leftover B-only 偏薄，实抽 20 / 20 / 18）；EGFR/HER2 不具备同等抽样条件，已撤回的 PIK3CA/PIK3CB 也未纳入本项或后续任何敏感性分析。受体、盒子、配体准备、exhaustiveness、打分与统计与主基准相同（Tables S8、S13）。

文献阻断分析使用同一套冻结分数。共享任一保留高置信 `document_id` 的配体连成一组，使同一篇文献的化合物不能同时进入训练与测试。`GroupKFold` 按该组划分；ECFP4、物化描述符与 docking logistic 使用相同折（Tables S39、S40）。缺少正负两类的折被丢弃；若有效折少于 2，则报告无法稳定估计，而不在看到 AUROC 后更换分组规则。document-cluster bootstrap 重采样的是这些文献连通组，而不是把同一系列化合物当作独立观察。

文献年份分割在计算 AUROC 之前冻结（`docs/TIME_SPLIT_PROTOCOL_FREEZE.md`）。配体年份取其保留高置信记录中最早的 `document.year`。主截止年为 2018（训练：first year < 2018；测试：first year ≥ 2018）；2015 与 2020 为预先指定的敏感性。晚期文献中的化合物不参与阈值、受体或指标选择。仅当测试集 dual、A-only、B-only 每类 n ≥ 10 时报告方向 AUROC；更小格子只报计数（Table S41）。主截止年至少两个靶对通过该门槛，才包装为外部验证。

对 352 个已打分配体中的 186 个优先分子提取了 assay-context 字段，包括 EGFR/HER2 全部方向类、PIK3CA/mTOR 的 4 个 neither、混合端点记录，以及对主 AUROC 影响最大的分子（`assay_context_audit.csv`）。随后一次元数据审核为全部 186 个优先分子填写纳入/排除（179 include / 7 uncertain / 0 exclude）；ChEMBL assay 自由文本不可用，故蛋白构建体与突变状态仍为 unknown。冻结 DualFourClass 标签未改，因此未重算 Table 2。该步骤不是文献级 assay 条件统一。

受体结构敏感性分析另选满足以下预先声明条件的替代晶体：（i）polymer entity 与目标蛋白真实对应；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.3 相同的共晶重对接 QC。实际对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。替换采用单口袋设计：在 PIK3CA/mTOR（PM48）上，4JPS/5DXT 替换口袋 A、口袋 B 仍用冻结 4JT6 分数，4JSX 替换口袋 B、口袋 A 仍用冻结 4L23 分数（exhaustiveness = 16；Table S10）。在受体身份审计发现问题之前，同一套 4JPS/5DXT 也曾用于替换已撤回的 PIK3CA/PIK3CB 口袋 A（口袋 B 仍用冻结 2WXF 分数，exhaustiveness = 8）；该结果仅作为受体身份错误记录的一部分列入 Supporting Information，不作为主要受体敏感性结果。刚体 Cα 叠合作为探索性几何对照（Table S10）。在 Table S30 所用同一套 PM48 配体与 PIK3CA 晶体上，另做探索性接触快照：占有率定义为与 20 个冻结口袋残基的重原子距离 ≤ 4.5 Å（Table S33）。占有率变化只作为结构假说，不是残基层因果解释。

事后 θ = 6.0 四状态普查复用冻结的 J0 候选靶对名单与缓存 pChEMBL 图，去掉顺序别名后剩 49 对。dual/A-only/B-only 均 n ≥ 10 记为方向可评估，neither 也 n ≥ 10 记为设定可评估（Table S44）。这些 J0 名单计数只诊断标签供给，不是对接菜单；对接主表为撤回 PIK3CA/PIK3CB 后的八行（原冻结三对加普查后五对）。冻结已打分面板上的多元物化匹配按 z 标准化 MW/cLogP/TPSA/重原子做 1:1 贪心配对，欧氏 caliper 为 0.5 与 1.0 SD（Table S45）；n_matched < 8 标记效能不足。AND 式双口袋过滤在 Dual+A-only+B-only 库上按 `vina_worst` 或 `vina_mean` 的 Dual 百分位截断（Table S46）。配体层 ECFP4 与四描述符逻辑回归则在三个主评价靶对的完整 θ = 6.0 ChEMBL 图上拟合，每类最多抽 120 个分子（种子 20260729），支架 `GroupKFold`（Table S47）。该分析只用实验标签与二维结构。

MCL1/Bcl-xL 在 LC6 pose-gold 未建立后退出主评价，面板对接仅作为 Supporting Information 中的探索性归档（Tables S50–S53）。[17]

计算在 Python 3 环境下完成，主要软件为 RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2 与 RTMScore。评价面板、对接分数、分析脚本与参数表见 Data and Software Availability。评价合约见 `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`。

## 3. 结果

### 3.1 双端实验数据供给与四状态评价集构建
ChEMBL 中能够支持四状态评价的双端活性数据随样本要求提高而迅速减少。在所有至少存在双端测量的靶对中，2,164,618 对至少有 1 个双端测量配体，63,790 对具有不少于 10 个双端测量配体。进一步要求在 \(\theta=6.0\) 下 dual、A-only 和 B-only 均不少于 10 个后，仅剩 5,253 对；采用严格 6.5/5.5 双向选择性供给标准后，进一步减少至 86 对。

对满足数据供给条件的靶对进一步检查后，部分体系因共同高通量测定背景、靶点类型、结合位点性质或缺少合适的人源实验结构而不适合纳入统一非共价对接评价。结合双向选择性样本供给、靶点类型、结合位点性质以及人源共晶结构可用性，最终主要评价包括 PIK3CA/mTOR、AChE/BChE、EGFR/HER2、F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD 八个靶对（Table 1）。EGFR/HER2 的严格 B-only 样本供给仍然有限，因此保留为双向选择性数据受限的评价体系。四状态分类要求同一配体在两个靶点上均具有实验测量，双端实验测量覆盖和双向选择性样本数量共同限制了严格四状态评价集的规模。

### 3.2 实验状态定义在部分靶对中改变 docking 性能解释
在统一 \(\theta=6.0\) 标签下，八个靶对的方向性判别表现差异较大，\(\mathrm{summary}_{\min}\) 介于 0.345 至 0.692 之间。仅 PPARG/PPARA 的配体水平 bootstrap 95% 置信区间完全高于 0.5（0.649 [0.504, 0.751]），其余靶对的置信区间均包含或低于 0.5（Table 2）。PPARG/PPARA 是唯一在主要 Vina 评价中区间完全高于 0.5 的靶对。

为单独考察实验状态比较本身的影响，我们首先固定相同的靶点评分通道，仅改变负类组成（Figure 2C；Table S34）。EGFR/HER2 使用靶点 A（EGFR）评分时，dual-versus-B-only 的 AUROC 为 0.430，而将负类替换为双端低活性的 neither 配体后，AUROC 升至 0.808，差值达到 0.378 [0.205, 0.547]；另一方向上的负类差异明显较小。JAK1/TYK2 同样出现显著的负类依赖差异：使用靶点 A（JAK1）评分时，dual-versus-neither 与 dual-versus-B-only 的差值达到 0.444 [0.263, 0.620]，表明这一现象并非单一体系特有。这种负类依赖差异在不同靶对之间表现出明显的特异性：在 AChE/BChE、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD 中，相同评分通道下的负类差值均较小，95% 置信区间均包含 0（Figure 2C；Table S34）。

在进一步采用双口袋平均评分进行描述性比较时，EGFR/HER2 和 JAK1/TYK2 的 dual-versus-neither AUROC 分别达到 0.756 [0.562, 0.920] 和 0.770 [0.597, 0.906]，而其方向性 \(\mathrm{summary}_{\min}\) 分别仅为 0.430 和 0.365（Table 3；Figure 2B）。由于该比较同时改变了评分聚合形式和负类定义，负类组成的单独影响以上述固定评分通道分析为准。

评价集的实际排序直观反映了这一设定的筛选后果。在包含 110 个配体的 EGFR/HER2 评价集按双口袋平均评分排序时，Top-10 中仅有 1 个 dual，其余 9 个均为单靶选择性配体，而没有任何 neither 配体（Table S25）。该结果直观显示，较高的 dual-versus-neither AUROC 并未转化为对单靶选择性配体的有效排除。若以 dual 配体的中位 worst 评分进行双口袋联合过滤，在保留 14 个 dual 的同时仍会保留 33 个单靶选择性配体，precision 仅为 0.298（Table S46）。

**Table 2.** 八个主评价靶对上的口袋匹配方向 AUROC（Vina，统一 \(\theta=6.0\)）及四个物化描述符的 \(\mathrm{summary}_{\min}\) 参考。表中类别样本量为 n_scored（dual / A-only / B-only）。最高描述符代表最佳单一描述符基线。

| 靶对 | n_scored (dual / A-only / B-only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 [0.211, 0.477] | 0.432 | 0.482 | 0.515 | 0.345 |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 [0.231, 0.503] | 0.369 | 0.389 | 0.580 | 0.425 |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 [0.444, 0.725] | 0.578 | 0.565 | 0.480 | 0.570 |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 [0.504, 0.751] | 0.507 | 0.478 | 0.485 | 0.627 |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 [0.296, 0.584] | 0.490 | 0.436 | 0.564 | 0.351 |

**Table 3.** 同一套 Vina 对接分数在方向性与 Dual-versus-neither 设定下的 AUROC 比较（统一 \(\theta=6.0\)）。Dual-versus-neither 采用双口袋平均评分 vina_mean。PIK3CA/mTOR 的 neither 样本量（n = 4）较少。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (vina_mean) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 | 0.405 [0.273, 0.534] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 | 0.527 [0.411, 0.638] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 | 0.668 [0.561, 0.770] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 | 0.675 [0.571, 0.780] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 | 0.522 [0.406, 0.640] |

### 3.3 配体化学基线可解释部分方向性判别
部分靶对中，配体自身的物理化学特征即可对 dual 和单靶选择性配体产生一定区分。以 AChE/BChE 为例，仅使用 TPSA 对 dual-versus-A-only 和 dual-versus-B-only 的判别 AUROC 分别达 0.733 和 0.801（Figure 3C；Table 2）。这种单一物化性质的区分能力具有靶对特异性，PIK3CA/mTOR 中最佳单一描述符（重原子数）的 \(\mathrm{summary}_{\min}\) 仅为 0.463（Table 2；Table S28）。在八个靶对中，Vina 与各自最佳单一描述符之间的 \(\mathrm{summary}_{\min}\) 差值大多数接近或包含 0（Figure S4；Table S19）。

在 Bemis–Murcko 骨架分组交叉验证下，不使用受体结构的 ECFP4 分子指纹在多个方向上获得与 docking 相当甚至更高的 AUROC（Figure 3A）。将对应方向的 docking score 作为输入特征加入 ECFP4 逻辑回归模型后，全部 16 个方向上的 AUROC 变化介于 −0.023 至 +0.021 之间，最大绝对变化仅为 0.023（Figure 3B；Table S24）。在当前评价集和交叉验证设置下，未观察到 docking score 在配体二维化学信息基础上提供稳定的增量判别。

### 3.4 方向性判别的结构对应性与计算敏感性
共晶配体重对接显示，所有主受体的保存姿态中均存在重原子 RMSD < 2.0 Å 的近天然构象，但最低 RMSD 姿态并不总是排名第一。例如，EGFR 3POZ 的 top-1 RMSD 为 9.505 Å，而全部保存姿态中的最低 RMSD 为 0.760 Å；HER2 3RCD 的 top-1 RMSD 为 1.855 Å，最低为 1.394 Å（Table S3）。这表明当前搜索设置能够覆盖近天然构象，但打分函数在姿态排序上仍存在局限。

为检验方向性判别是否具有与活性差异靶点一致的结构特异性，我们比较了使用正确对应口袋评分与非对应口袋评分的判别结果（Figure 5A）。在主评价集中，EGFR/HER2 和 AChE/BChE 的匹配评分 \(\mathrm{summary}_{\min}\) 高于非对应评分，差值分别为 0.170 [0.060, 0.280] 和 0.161 [0.037, 0.269]，置信区间均排除 0；其余六个靶对的差值 95% 置信区间均包含 0。而在未使用池留出集（unused-pool holdout）中，所有七个可评靶对的对应与非对应口袋评分差值均包含 0（差值介于 −0.079 至 +0.150 之间；Figure 5B），表明对应口袋优势在独立样本中并未稳定重现。

在 PIK3CA/mTOR 中，将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.470, 0.813] 降至 0.486 [0.259, 0.692]；替换为 5DXT 后降至 0.505 [0.292, 0.696]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.418, 0.776]（Figure 4B；Table S30）。这表明方向性评价结果对受体晶体结构的选择具有较高敏感性。

采用 GNINA 1.3.2 独立生成姿态和评分后，EGFR/HER2 的 dual-versus-neither AUROC 为 0.783 [0.610, 0.922]，而方向性 \(\mathrm{summary}_{\min}\) 仅为 0.220 [0.109, 0.343]；JAK1/TYK2 同样呈现一致模式，其 dual-versus-neither 为 0.705，而方向性 \(\mathrm{summary}_{\min}\) 仅为 0.317 [0.183, 0.463]（Figure 4A；Table S32）。这表明在更换独立姿态搜索和基于深度学习的评分流程后，选择性负类与 neither 之间的表观判别差异仍然稳固存在。

PPARG/PPARA 在主要 Vina 评价中是唯一 \(\mathrm{summary}_{\min}\) 置信区间完全高于 0.5 的靶对（0.649 [0.504, 0.751]），但在同姿态 RTMScore 重评分下降至 0.369 [0.233, 0.475]，GNINA CNN 重评分降至 0.500，且在未使用池留出集中降至 0.535 [0.350, 0.717]，表明单一软件实现下的表观优势对评分方式和样本组成较为敏感。此外，不同 Vina 随机种子及 PIK3CA/mTOR 的 exhaustiveness 设置产生的数值波动相对有限（Figure 4C；Figure 6C；Table S54）。

### 3.5 方向性评价的标签与样本敏感性
改变活性阈值（\(\theta = 5.5, 6.0, 6.5\) 及严格 6.5/5.5）后，样本较充足靶对的 \(\mathrm{summary}_{\min}\) 总体平稳（如 AChE/BChE 稳定在 0.606，F2/F10 介于 0.345–0.353，JAK1/TYK2 介于 0.365–0.385）；但当单侧选择性样本量严重萎缩时，估计波动明显增大（如 EGFR/HER2 在严格条件下仅剩 7 个 B-only 配体，\(\mathrm{summary}_{\min}\) 降至 0.324；Figure 6A；Table S4）。将重复活性记录的汇总方式由最大 pChEMBL 改为中位数后，标签一致率保持在 93.6%–100%，方向性汇总指标的最大绝对变化仅为 0.023（Table S29）。

在排除主评价集分子后，基于剩余候选分子构建的未使用池留出集（unused-pool holdout）进一步评估了样本组成的影响。各靶对表现出明确的靶对依赖性：AChE/BChE（0.618 [0.422, 0.759]）、PIK3CA/mTOR（0.765 [0.603, 0.891]）和 JAK1/JAK2（0.619 [0.420, 0.749]）保持了相近水平；JAK1/TYK2（0.475 [0.282, 0.660]）略有上升；F2/F10（0.392 [0.214, 0.573]）与 PPARA/PPARD（0.445 [0.241, 0.559]）仍处于较低水平；而主评价中较高的 PPARG/PPARA 则由 0.649 降至 0.535 [0.350, 0.717]（Figure 5C；Table S13）。

在具备完整文献来源标识、能够进行来源簇重采样的评价集中，改以文献连通簇作为重采样单位后，部分方向的置信区间显著扩大（如 EGFR/HER2 dual-versus-B-only 的 95% 置信区间扩大至 [0.321, 0.617]，PIK3CA/mTOR 相应区间扩大至 [0.000, 0.818]；Tables S39–S40）。可检测效应模拟进一步表明，在当前类别样本量下，即便真实方向性 AUROC 达到 0.70，\(\mathrm{summary}_{\min}\) 的 95% 置信区间排除 0.5 的概率也仅为 0.219–0.621；当真实 AUROC 为 0.60 时该概率仅为 0.025–0.065（Table S31；Figure S6）。

### 3.6 独立外部评价受双端实验数据供给限制
对 BindingDB 和 PubChem 中的双端实验活性数据进行检索与核对后，经过人源单蛋白、定量测量、排除与主评价集共享来源和结构高相似分子（ECFP4 Tanimoto < 0.70）等多重准入过滤后，没有任何靶对同时满足 dual、A-only 和 B-only 样本量（各 n ≥ 20）及来源多样性（各 ≥ 3 个独立来源）要求，因此均未进入外部对接流程（Figure 6D；Tables S48–S49）。

按配体发表年份进行时间切分同样受到双向选择性样本供给的限制。以 2018 年为主要界点时，没有任何靶对在 2018 年之后的测试集中同时具有不少于 10 个 dual、A-only 和 B-only 配体，因此未能构建出可同时独立评价两个方向的时间切分测试集（Table S41）。综上，当前公开数据供给尚不足以构建满足严格准入标准的独立双端四状态外部评价集。

## 4. 讨论

### 4.1 基准设定改变了双靶对接的证据标准

该基准将双靶识别具体化为两个方向上的实验选择性硬负判别任务。冻结 EGFR/HER2 分数上，Dual versus neither 明显强于方向性硬负任务（Results 3.2；Table 3）。独立 GNINA 姿态生成再现了该差距，因此此处的设定效应不是 Vina 特有的。EGFR/HER2 的设定差距在五个 Vina 种子上均为正（Table S54）。JAK1/TYK2 是第二对配方差距（0.770 对 0.365；独立 GNINA 0.705 对 0.317）。其余靶对并不支持同一普遍规律，八行也不平均。相对于 Zhou 等，[9] 本文比较的是对接能否把 dual-active 与实验定义的选择性配体分开，而不仅是与 inactive 分开。已有对接基准表明，decoy 构建、化学偏倚和真实 assay 标签会改变虚拟筛选解释；[5–7,12,13,21] 同一关切适用于双靶任务。Wu 等表明大规模对接可以前瞻性得到联合结合配体；[19] DualFourClass 回答的是另一问题——当负类由实验选择性配体而不是非结合配体定义时，回顾性证据是否改变。Kinase-Bench 检验相对激酶 decoy 的选择性富集，[22] DualFourClass 则在每一对上直接构造四种实验测定状态。

### 4.2 配体化学与受体实现

物化描述符与化学型已经携带了相当一部分实验标签信息。对接相对支架分组纯配体模型几乎没有增量判别（Results 3.3；Table S24）。若缺少这些配体层对照，一个看似优秀的双靶对接结果可能只是识别了与 dual 标签相关的分子属性。[7,12]

表观判别还取决于受体实现。mTOR 端固定、只替换 PIK3CA 端时，PIK3CA/mTOR 的表观判别降低（Results 3.4；Figure 4B），与把受体表示视为性能变量的激酶交叉对接工作相一致。[14] 同一套替换用于受体身份存疑的 PIK3CA/PIK3CB 时方向相反，但由于该口袋后来被证实是非人源脱靶受体（见 Supporting Information），本文不以此作为"受体效应因靶对而方向相反"的普遍证据；在唯一经受体身份核实的靶对上，替换只降低了表观判别。

### 4.3 对双靶虚拟筛选的含义

同时在两个口袋获得有利分数，并不能自动建立实验定义的双靶活性。在双口袋分数看起来有利之后：（i）要求相对 A-only 与 B-only 硬负的方向性判别；（ii）询问配体层 ECFP 或物化模型在抗泄漏分割下是否回收相似信号；（iii）检验未使用配体池或文献阻断分割；（iv）至少替换一端受体实现（Figure 1；TOC graphic）。POLYGON 等生成式方法表明双靶设计可以被实验验证，[20] 但对接型双靶主张在把双口袋占有当作识别证据之前，仍需要选择性感知评价。这一关切与近期 JCIM 研究相一致：在有实验依据的筛选集上，docking 再打分表现可以有很大差异。[15]

### 4.4 局限性

分析覆盖八对靶标，而不是靶标通用基准：2026-07-23 冻结撤回 PIK3CA/PIK3CB 后的三对，加上 ChEMBL 37 普查后加入的五个普通靶对。这是收集补全，不是预先注册的八靶对冻结。PPARG/PPARA 是唯一主面板最弱臂 CI 完全高于 0.5 的一对，且该数字对面板成员（unused-pool holdout）与打分公式（RTM）敏感。PIK3CA/PIK3CB 已完成对接，但事后受体身份审计发现其预期的 PIK3CB 受体实为小鼠 PIK3CD，因此撤回并仅作为已记录的受体身份错误保留在 Supporting Information，提示共晶配体重对接 QC 本身无法验证受体身份，而不作为主结果。标签来自异质公开生物活性记录，并要求两端均有完整测定（完整病例 14.5%–34.0%）。受体依赖性只在选定结构上评价，且未对新预测双靶配体做前瞻实验。没有任何靶对满足预先规定的 BindingDB 外部切片门槛，因此该切片不作为外部验证。共晶 best-of-nine QC 证明搜索覆盖，不是 top-ranked pose 验证。MCL1/Bcl-xL 已退出主文，仅见 Supporting Information。

## 5. 结论

用实验定义的单靶选择性配体作为硬负，在部分设定下实质改变了双靶识别的表观证据。配体化学与受体实现进一步贡献了这一变异，表明两个口袋中的有利分数应当连同选择性感知对照一起解读。这些发现来自八个靶对案例——2026-07-23 冻结保留的三对，以及 ChEMBL 37 普查后加入的五对——旨在推动更广泛的验证，而不是建立靶标通用的对接性能。八对不是预先注册冻结，也不平均。原先对接的一个候选靶对因事后发现的受体身份错误而撤回，提示共晶配体重对接 QC 本身不能验证受体身份。失败的 BindingDB 原生重建不作为外部验证。

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源表，包括独立 GNINA 姿态生成分数（Table S32）、PIK3CA 占有率位移（Table S33）、文献阻断交叉验证（Tables S39–S40）、冻结的文献年份分割（Table S41）、assay-context 元数据审核（Table S42）、BindingDB REST 独立性计数（Table S43）、θ = 6.0 靶对普查（Table S44）、物化 caliper 匹配（Table S45）、AND 过滤工作点（Table S46）、配体层全图 AUROC（Table S47）、BindingDB 原生候选流程与切片摘要（Tables S48–S49）、MCL1/Bcl-xL 面板、LC6 gate 与 applicability-stress-test AUROC（Tables S50–S51、S53）、文献对照（Table S52）、与 Table 3 同一 Dual-versus-neither 估计量的五种子 Vina 敏感性（Table S54）、重建的 EGFR/HER2 共晶 QC（Table S3）、Figure S8、原生切片合约（`protocol/external_slice_contract.yaml`），以及评价合约（`DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`）。面向稿件的表 SHA-256 校验和见 `REVISION_CHECKSUM_MANIFEST_v1.csv`。ChEMBL 供给审计冻结于 2026-07-23；高置信 activity 视图抓取于 2026-08-26；BindingDB 原生归档锁定为 202608。GitHub Release 与 Zenodo DOI 将从打标签快照签发，而不是从当前仍可能变化的分支签发。分析环境与零新对接的复现命令见仓库 README。BindingDB TSV 归档本身不随仓库分发；CI 只核对已提交的 CSV。

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

(19) Wu, Y.; Vigneron, S.; Braz, J.; Srinivasan, K.; Fink, E. A.; Huang, X.-P.; Xu, X.; Huebner, H.; Kim, J. Y.; Wang, J.; Pfeiffer, T.; Sakamoto, K.; Moroz, Y. S.; Radchenko, D. S.; Rodriguiz, R. M.; Irwin, J. J.; Gmeiner, P.; Billesboelle, C.; Roth, B. L.; Basbaum, A. I.; Manglik, A.; Wetsel, W. C.; Shoichet, B. K. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.

(20) Munson, B. P.; Chen, M.; Bogosian, A.; Kreisberg, J. F.; Licon, K.; Kuenzi, B. M.; Ideker, T. De novo generation of multi-target compounds using deep generative chemistry. *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.

(21) Gu, S.; Shen, C.; Zhang, X.; Sun, H.; Cai, H.; Luo, H.; Zhao, H.; Liu, B.; Du, H.; Zhao, Y.; Fu, C.; Zhai, S.; Deng, Y.; Liu, H.; Hou, T.; Kang, Y. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.

(22) Wei, T.-H.; Zhou, S.-S.; Jing, X.-L.; Liu, J.-C.; Sun, M.; Zhao, Z.-H.; Li, Q.-Q.; Wang, Z.-X.; Yang, J.; Zhou, Y.; Wang, X.; Ling, C.-X.; Ding, N.; Xue, X.; Yu, Y.-C.; Wang, X.-L.; Yin, X.-Y.; Sun, S.-L.; Cao, P.; Li, N.-G.; Shi, Z.-H. Kinase-Bench: Comprehensive Benchmarking Tools and Guidance for Achieving Selectivity in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.
