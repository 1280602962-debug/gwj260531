# 方法初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> 引言：[`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。  
> 口径与引言一致：产出为 URAT1–NLRP3 **双节点候选分子**的计算假说，不是已验证的双靶抑制剂；对接评分仅作池内排序，不解释为结合亲和力。

## 2 材料与方法

### 2.1 研究设计

本研究的总体流程如图 1 所示。引言已说明，URAT1 与 NLRP3 两侧的公开活性证据在规模、测定条件和化学空间覆盖上并不对等，因此不能对两靶采用同一套数据驱动策略。据此，NLRP3 分类模型仅用于缩小临床阶段化合物库；URAT1 侧则以结构对接作为主要排序依据。为避免临床库中的候选身份反过来影响方法选择，对接协议的比较在解释临床库结果之前完成，且不使用临床库排名。比较在 URAT1 活性物–诱饵基准上进行，并同时报告性质匹配诱饵与随机诱饵上的富集，以降低单一负样本构造带来的表观优势（Gu 等，2025）。

随后，满足 NLRP3 模型阈值的分子分别对接至 URAT1 与 NLRP3。两靶结构分及 NLRP3 模型分先在同一成功对接集合内转换为百分位，再作 Pareto 非支配分析，该分析仅用于审计。跟进提名与 Pareto 成员身份分开：以双对接百分位门控（\(S_U\) 与 \(S_{N,\mathrm{dock}}\) 同时达到阈值）筛选较高优先级假说，再经结构警报、类药性、分子量窗口和骨架去冗余。URAT1 回归模型仅作对照，不参与临床库主排序。生产对接读出为事先锁定的 gnina CNNaffinity（记为 P2）。全文将上述名单视为可供后续实验验证的计算线索，而不将其表述为已发现的双靶抑制剂。

### 2.2 活性数据收集与临床化合物库

人源 URAT1（SLC22A12，ChEMBL 标识 CHEMBL6120）与 NLRP3（CHEMBL1741208）的活性记录取自 ChEMBL（Zdrazil 等，2024）。URAT1 仅保留标准关系符为“=”且具有有效 SMILES 的记录。优先采用 ChEMBL 提供的 pChEMBL 值；当该值缺失且标准值为正、单位为 nM 或 µM 时，先换算为摩尔浓度，再按式（1）计算：

\[
\mathrm{pActivity}=-\log_{10}C ,
\tag{1}
\]

其中 \(C\) 为以 mol/L 表示的标准活性浓度。活性值限定在 4–10 之间。同一来源 SMILES 上标准差超过 0.5 或极差超过 1.0 个对数单位的记录组视为测定冲突并剔除，其余取中位数；随后进行结构规范化并按规范 SMILES 去重。清洗后得到 822 个 URAT1 化合物，对应 218 个 Bemis–Murcko 骨架（Bemis 与 Murcko，1996）。

NLRP3 保留测定描述中含“IL-1”、ChEMBL 测定类型代码为 B、标准关系符为“=”且 pChEMBL 位于 4–10 的记录；该靶点直接使用 pChEMBL，不再由浓度单位另行换算。多数此类记录为细胞体系中的 IL-1 相关读出，因此模型学习的是测定条件下的相对活性，而不是与 NACHT 口袋的直接结合常数。仅纳入至少包含 5 个化合物的测定，并以 \(\mathrm{pActivity}\ge6\) 定义活性。该数据集包含 609 条测定–分子记录、513 个唯一化合物和 39 个测定；建模时保留全部测定标识。URAT1 与 NLRP3 清洗数据之间不存在相同的规范 SMILES，因此未构建共享训练集或多任务模型。

临床重定位库的纳入口径如下。化合物来自 ChEMBL 临床阶段导出集以及一级、二级 ATC 导出集，合并后按 InChIKey 第一段去重。结构去盐后保留重原子数最多的片段，并排除蛋白质、肽、寡核苷酸和抗体等非小分子实体。分子量过滤范围为 150–800 Da；若导出表提供分子量则优先采用，否则由标准化结构计算。经 RDKit 消毒后转为规范 SMILES。ChEMBL 最高临床阶段标注用于描述分子所处开发阶段，但不作为进入主分析库的硬性门槛：ATC 来源且尚无阶段标注的小分子仍可纳入。已知对照药若未出现在导出集中，则作为基准分子补入。最终临床库包含 8,319 个唯一小分子。该库用于回答引言提出的转化问题，即在已经进入人体研究或临床使用的小分子中寻找双节点计算线索；其规模数字只在方法中给出操作定义，不作为发现主张。

### 2.3 URAT1 活性物–诱饵基准

URAT1 对接协议的选择不依赖临床库。活性基准由清洗数据中 \(\mathrm{pActivity}\ge6\) 的分子组成（469 个）。TrueDecoy 负样本由两部分构成：（i）清洗数据中 \(\mathrm{pActivity}<5\) 的实验弱活性分子（80 个），作为有实测依据的阴性种子；（ii）自商业类药小分子库经水库抽样与物化包络预过滤后的候选池中抽取的未标注分子，按分子量、\(\log P\)、拓扑极性表面积、氢键供体数、氢键受体数和可旋转键数作性质匹配，使活性物–诱饵比例达到 1:10。对活性分子 \(a\) 与候选匹配诱饵 \(x\)，归一化性质距离定义为

\[
d(a,x)=
\sqrt{
\frac{1}{|\mathcal D|}
\sum_{k\in\mathcal D}
\left(\frac{a_k-x_k}{w_k}\right)^2
},
\tag{2}
\]

其中 \(\mathcal D\) 为上述六项描述符，严格匹配窗口 \(w_k\) 依次为 40 Da、1.0、25 Å\(^2\)、1.5、2.5 和 2.5。后三项为整数描述符，实际允许差值分别不超过 1、2 和 2。为降低近邻活性类似物被误作诱饵的风险，以半径 2、长度 2,048 位的 Morgan 指纹计算候选分子相对全部活性分子的最大 Tanimoto 相似度（Rogers 与 Hahn，2010），并要求 \(\mathrm{TC}_{\max}(x)\le0.50\)。性质匹配采用不重复的 round-robin 分配，随机种子固定为 42；在扩大后的候选池上，全部 4,610 条匹配均落在严格窗口内。最终 TrueDecoy 含 469 个活性分子与 4,690 个负样本（80 个实验弱活 + 4,610 个性质匹配）。

该构建参照 Gu 等（2025）的 TrueDecoy/RandomDecoy 评测思路，但负样本并非 BindingDB 多靶实验无活集，而是单靶 URAT1 上的“实验弱活种子 + 商业库性质匹配分子”。RandomDecoy 与 TrueDecoy 共用同一批活性分子；负样本从同一预过滤商业库池中随机抽取，数目同为 4,690，且禁止与 TrueDecoy 负样本的 SMILES 重叠。两套诱饵集合的 SMILES 交集为零。对接在活性物与两套诱饵的结构并集上进行（9,849 个独特结构），再分别汇总 TrueDecoy 与 RandomDecoy 上的富集指标。早期允许诱饵重叠的蒸馏子集构建已废弃，不以该旧集作为协议筛选依据。

### 2.4 NLRP3 分类模型与 URAT1 对照模型

分子表示由半径为 2 的 2,048 位 Morgan 指纹与 12 个 RDKit 描述符拼接而成（Landrum）。描述符包括分子量、\(\log P\)、拓扑极性表面积、氢键供体数、氢键受体数、可旋转键数、环数、芳香环数、sp\(^3\) 碳比例、Bertz 复杂度、杂原子数和定量类药性估计（QED；Bickerton 等，2012）。NLRP3 模型进一步加入出现频率最高的 25 个测定标识及一个“其他测定”类别的 one-hot 编码，并采用 XGBoost 二分类器（Chen 与 Guestrin，2016）。模型包含 400 棵树，最大深度为 5，学习率为 0.05，行采样率和列采样率均为 0.8，\(L_2\) 正则化系数为 1.0。每条训练记录的权重设为其所属测定样本数平方根的倒数，以降低大样本测定的支配作用。

模型性能采用 5 折 Bemis–Murcko 骨架分组交叉验证评估，同一分子的全部测定记录始终位于同一折。各折中拼接后的特征仅使用训练折拟合标准化器，再将相同变换应用于测试折。交叉验证中，同一测试分子在其已有测定记录上的预测值取最大值，得到分子级输出。报告指标包括 AUROC、AUPRC 和前 10% 富集因子。用于临床库筛选的最终模型在完整训练记录上拟合；等渗映射使用按骨架分组划出的校准子集拟合。对于每个临床库分子，模型分别在出现频率最高的 5 个测定条件下计算输出，取其中最大值作为 NLRP3 分类分数 \(q_N\)。由于交叉验证聚合方式与部署时的五测定推理并不完全相同，该分数只作为相对筛选分数，而不解释为严格校准的活性概率。满足 \(q_N\ge0.5\) 的分子进入双靶对接池（8,319 个临床库分子中 1,588 个）。

URAT1 对照模型采用相同的分子表示和 XGBoost 回归器，在未使用有机阴离子转运体辅助迁移的条件下训练，并通过 5 折骨架分组交叉验证计算均方根误差、决定系数和 Spearman 相关系数。另以 lesinurad、benzbromarone、verinurad 和 dotinurad 进行预设命名药物的回顾性检查。由于这些分子并非全部位于训练化学空间之外，该检查只用于观察已知药物的回收行为，不作为严格的外部验证。URAT1 回归分数不用于临床库主排序。

### 2.5 蛋白与配体准备

URAT1 采用与 lesinurad 复合的 inward-open 结构 9DKB（分辨率 2.55 Å；Suo 等，2025）。NLRP3 采用 NACHT 结构 7ALV；其共晶配体为磺酰脲类类似物 NP3-146，而不是 MCC950 本身（Dekker 等，2021）。两种结构均保留 A 链，并去除结晶水、共晶配体和其他异原子。蛋白结构经 Gemmi 提取后，用 Open Babel 在 pH 7.4 条件下加氢并转换为刚性受体 PDBQT（O’Boyle 等，2011）。配体由规范 SMILES 出发，使用 RDKit ETKDGv3 生成三维构象，随后进行最多 200 步 MMFF 几何优化，再经 Meeko 生成 PDBQT。构象嵌入使用固定随机种子，以保证配体准备可复现。

URAT1 搜索盒以 9DKB 共晶配体的几何中心为中心，坐标为 \((99.966,\,102.967,\,105.699)\) Å，三个方向的边长均为 22 Å。NLRP3 搜索盒以 7ALV 共晶配体为中心，坐标为 \((16.756,\,35.449,\,125.714)\) Å，三个方向的边长均为 20 Å。除上述两处生产口袋外，本文不对 URAT1 的其他构象状态做对接。

9DKB 中直接提取的共晶配体用于自对接，以分别评价构象生成和排序能力。报告原生排名第一构象相对晶体配体的重原子 RMSD、多构象集合中的最低 RMSD，以及 RTMScore 所选构象的 RMSD。几何通过阈值设为 2.0 Å。RMSD 根据完成原子对应后的配体重原子坐标计算。该检验只回答搜索盒与配体准备是否能够采样到近原生姿，以及某一读出能否把近原生姿排在第一；它不替代活性物–诱饵基准上的协议锁定，也不把生产读出的第一构象自动视为晶体结合模式。

### 2.6 分子对接、协议评价与生产读出

协议比较在 9DKB 上进行。AutoDock Vina 1.2.5 以 exhaustiveness 32 运行，每个配体最多输出 9 个构象，能量窗口设为 3.0 kcal/mol（Trott 与 Olson，2010；Eberhardt 等，2021）。gnina 使用相同受体和搜索盒，在 CPU 模式下以 exhaustiveness 32 运行，并启用卷积神经网络重打分（McNutt 等，2021）。RTMScore 用于重新评价 Vina 与 gnina 产生的构象集（Shen 等，2022）。比较的读出包括 gnina CNNscore、Vina affinity、gnina CNNaffinity、gnina affinity、Vina 构象的 RTMScore，以及 gnina 构象的 RTMScore。CNNscore 预注册为负对照，其余五项纳入生产协议选择。

为消除不同评分方向带来的混淆，全部评分先转换为“数值越高、排名越前”的统一方向。对 Vina affinity 和 gnina affinity 等低值优指标，取 \(s=-E\)；对 CNNscore、CNNaffinity 和 RTMScore 等高值优指标，直接取原始分数。对接评分始终作为同一基准或同一临床对接池内的相对排序，而不转换为实验 \(K_i\) 或结合自由能。

设基准集中共有 \(N\) 个分子，其中 \(y_i=1\) 表示活性分子，\(y_i=0\) 表示诱饵。排名前比例 \(f\) 所含分子数为 \(n_f=\max(1,\lfloor fN\rfloor)\)，富集因子定义为

\[
\mathrm{EF}_{f}=
\frac{\displaystyle
\frac{1}{n_f}\sum_{i\in\mathrm{Top}(f)}y_i}
{\displaystyle
\frac{1}{N}\sum_{i=1}^{N}y_i}.
\tag{3}
\]

各读出在 TrueDecoy 与 RandomDecoy 上分别计算 ROC-AUC、EF\(_{1\%}\) 和 EF\(_{5\%}\)。TrueDecoy 的 EF\(_{1\%}\) 为首要判据；结果相同时依次比较 EF\(_{5\%}\) 和 ROC-AUC。RandomDecoy 作为否决对照，用于识别仅在容易区分的随机负样本上出现的表观优势。若主要指标仍无法区分候选协议，则比较四种已知尿酸药的排名百分位。选定后，生产协议以相同的配体准备、搜索深度和评分设置应用于 9DKB 与 7ALV。临床库生产筛选采用 gnina CNNaffinity（P2），exhaustiveness 为 32，并启用 CNN rescoring；每个配体按其首选构象的 CNNaffinity 进入后续百分位排序。P2 在自对接中可以采样到近原生姿，但第一构象不一定通过 2.0 Å 门控，因此不以该读出的 Top-1 作为构象金标准。具体富集数值与自对接 RMSD 见结果部分。

### 2.7 双靶百分位排序与 Pareto 审计

NLRP3 分类筛选后的分子分别在 9DKB 和 7ALV 上按选定协议对接。由于不同靶点和评分指标的原始数值不具有直接可比性，所有有效分数均在同一对接成功集合内转换为平均秩百分位。令 \(\mathrm{PR}(x_i)\) 表示 \(x_i\) 的平均秩除以有效分子数，取值范围为 0–1，并列值采用平均秩。统一为高值优方向后，三个证据分量定义为

\[
\begin{aligned}
S_U(i)&=100\,\mathrm{PR}\!\left[s_U(i)\right],\\
S_{N,\mathrm{dock}}(i)&=100\,\mathrm{PR}\!\left[s_N(i)\right],\\
S_{N,\mathrm{ML}}(i)&=100\,\mathrm{PR}\!\left[q_N(i)\right].
\end{aligned}
\tag{4}
\]

其中 \(s_U\) 和 \(s_N\) 分别表示方向统一后的 URAT1 与 NLRP3 结构评分。Pareto 审计将 NLRP3 轴定义为

\[
S_N(i)=
\max\left\{
S_{N,\mathrm{ML}}(i),
S_{N,\mathrm{dock}}(i)
\right\}.
\tag{5}
\]

该定义保留分类证据或结构证据较强的分子，并在审计中标明证据来源。作为敏感性分析，另分别使用 \(S_{N,\mathrm{ML}}\) 和 \(S_{N,\mathrm{dock}}\) 构建 NLRP3 单证据轴。机器学习分数可以改变 Pareto 审计中的 NLRP3 轴，但不能单独把分子抬进下文的优选跟进名单。

在同时具有两靶有效分数的集合中，同时最大化 \(S_U\) 和 \(S_N\)。若不存在另一分子 \(j\) 满足 \(S_U(j)\ge S_U(i)\)、\(S_N(j)\ge S_N(i)\) 且至少一个不等式严格成立，则分子 \(i\) 被定义为 Pareto 非支配分子。除非另有说明，Pareto 前沿不施加额外分数阈值。排序敏感性通过两轴前 1%、2%、5% 和 10% 的交集，以及双轴阈值 85、90 和 95 下的候选集合进行比较。有放回重采样仅用于描述候选集合组成对抽样的敏感性，不解释为对接评分不确定性或重复对接稳定性。原始非支配集合完整保留，用于展示大环内酯等高分子量分子如何占据对接优势，而不是作为跟进名单。

### 2.8 化学适用域、结构警报与候选提名

URAT1 化学适用域以 Morgan 指纹的最近邻 Tanimoto 相似度定义。首先对 822 个 URAT1 清洗分子执行留一法最近邻计算，并将所得分布的第 5 百分位设为阈值 \(\theta_{\mathrm{AD}}\)。查询分子相对 URAT1 清洗集的最大相似度不低于该阈值时，标记为位于该参照化学空间内。另分别计算候选分子相对 URAT1 清洗集和 NLRP3 已知活性分子的最近邻相似度，用于描述化学新颖性；这些相似度不作为活性分数。

药物化学审计使用 RDKit FilterCatalog 标注 PAINS-A、PAINS-B、PAINS-C、Brenk 和 NIH 结构警报（Baell 与 Holloway，2010；Brenk 等，2008），并计算分子量、cLogP、拓扑极性表面积、氢键供体数、氢键受体数和 QED。Lipinski 条件包括分子量不超过 500 Da、cLogP 不超过 5、氢键供体数不超过 5 及氢键受体数不超过 10；违反项不多于 1 项时记为通过（Lipinski 等，2001）。Veber 条件定义为可旋转键数不超过 10 且拓扑极性表面积不超过 140 Å\(^2\)（Veber 等，2002）。该 Veber 门槛同时作为口服吸收替代标志。Ghose 条件用于化学排序中的软标记，不作为硬性剔除规则（Ghose 等，1999）。Egan 规则、NIH 警报和胶体聚集启发式（cLogP 不低于 3.5、芳香环数不少于 3 且拓扑极性表面积低于 75 Å\(^2\)）一并计算并保留在审计表中，但不作为排序键。结构警报用于候选降级和结果解释，不被视为实验假阳性的直接证明。

上述步骤属于类药性与结构警报过滤，而不是细胞色素、hERG 或药物性肝损伤等药代毒性预测。无 PAINS/Brenk、同时通过 Lipinski 与 Veber、以及分子量 200–550 Da，用于优选候选的硬性门控；QED、分子量居中性、Ghose 通过和羧酸软标记仅用于门控后的相对排序。这些规则用于把大环内酯等对接优势分子降级，不能替代临床库分子已有的人体药代或安全性信息，也不把预测 ADMET 写成活性证据。MM-GBSA、MM-PBSA 或其他结合自由能不参与临床库排序。

候选提名在 Pareto 审计之后独立进行，且明确与“仅按对接分取非支配前沿”区分。跟进分子不以对接绝对分或未经验证的线性加权总分为唯一排序依据。优选假说采用双结构门控：同时满足 \(S_U\ge90\) 和 \(S_{N,\mathrm{dock}}\ge90\)。另保留 \(S_U\ge90\) 且 \(S_N\ge90\) 的较宽集合供审计；其中仅由 NLRP3 模型分抬入、对接百分位未达 90 的分子标记为仅模型证据风险，不进入优选短名单。对较宽审计集合中的每个分子标注是否命中 PAINS 或 Brenk、是否同时通过 Lipinski 与 Veber、分子量是否落在口服小分子窗口 \(200\le\mathrm{MW}\le550\) Da。无 PAINS、无 Brenk，且同时通过 Lipinski 与 Veber 者标记为清洁候选；在此基础上进一步满足双结构门控、口服分子量窗口与 Veber/口服吸收替代标志者标记为优选候选。红霉素、epothilone 等大环骨架因分子量超出该窗口而被降级，即使其对接百分位较高也不作为优选假说。已知对照药另作标记，不与新的重定位候选混为“新命中”。

排序时依次优先考虑优选候选、清洁候选、分子量窗口、非已知对照药、NLRP3 结构证据支持、较高临床阶段，以及双结构平衡分

\[
B_{\mathrm{struct}}(i)=
\min\left\{
S_U(i),S_{N,\mathrm{dock}}(i)
\right\},
\tag{6}
\]

再辅以由 QED、分子量居中性、Ghose 通过与羧酸软标记构成的化学排序分；最后按 Murcko 骨架贪心去冗余生成跟进短名单，避免同一大环骨架占据多个名额。该规则保留每个候选的结构证据来源和降级原因，从而使数学非支配分子与较高优先级假说可被分开阅读。

### 2.9 分子动力学模拟方案

对药物化学提名后的代表性分子，拟以其 gnina P2 生产构象为起点进行分子动力学模拟，用以检验对接占位姿在有限轨迹内能否保持关键接触。该步骤属于跟进计算，不重新选择生产对接协议，也不把端点结合能解释为实验亲和力。URAT1 为跨膜转运体，复合物置于膜–脂双层环境中构建；NLRP3 NACHT 为可溶性结构域，复合物置于显式溶剂水盒子中构建。对照体系为 lesinurad 与 9DKB、以及 MCC950 与 7ALV（若存在可用生产姿）。配体质子化状态在生理 pH 下复核后，采用通用小分子力场参数化；体系经能量最小化与平衡后进入生产模拟。轨迹分析优先报告配体及口袋残基的 RMSD/RMSF、与晶体对照残基的接触维持，以及氢键占有率。MM-GBSA 或 MM-PBSA 仅在生产轨迹完成后、对同一受体和同一模拟方案下的提名分子与对照药作相对比较；该类端点能不参与 1,588 分子对接池的排序，也不解释为实验结合自由能或 \(K_i\)。因 P2 第一构象不一定近原生，本文不以对接单帧直接计算结合自由能。URAT1 为膜蛋白，隐式溶剂单点能尤其不能代替膜体系采样。自由能微扰或热力学积分不在本文范围内。在轨迹完成前，本文不报告数值结果，也不将尚未运行的模拟写成已完成的构象结论。当前计划跟进的代表性分子为 GSK-3008348（URAT1 膜体系）和 Vecabrutinib（NLRP3 水盒子）；二者均为计算假说，不是已验证的双节点活性分子。

### 2.10 软件与可复现性

分子标准化、指纹、描述符和结构警报计算使用 RDKit；机器学习使用 scikit-learn、XGBoost 与 NumPy；受体处理使用 Gemmi 和 Open Babel；配体准备使用 RDKit 与 Meeko；对接使用 AutoDock Vina 1.2.5 和 gnina 1.3.1；候选构象重打分使用 RTMScore。搜索盒尺寸、exhaustiveness、CNN 重打分开关及分析阈值均预先固定。随机种子在诱饵分配、构象嵌入和交叉验证中保持一致。分析脚本随研究资料一并版本化，以便按相同纳入口径复现临床库构建、协议比较和提名规则。

---

## 引用清单（定稿时改为期刊格式）

方法部分在引言已列文献之外，主要补充工具、数据库与结构来源。与引言重复的条目保留相同作者—年份写法。

1. Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res.* 2024;52(D1):D1180-D1192. doi:10.1093/nar/gkad1004  
2. Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. *J Med Chem.* 1996;39(15):2887-2893. doi:10.1021/jm9602928  
3. Gu S, Zhang X, Shen C, et al. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat Mach Intell.* 2025. doi:10.1038/s42256-025-00993-0  
4. Rogers D, Hahn M. Extended-connectivity fingerprints. *J Chem Inf Model.* 2010;50(5):742-754. doi:10.1021/ci100050t  
5. Landrum G. RDKit: Open-source cheminformatics. https://www.rdkit.org  
6. Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the chemical beauty of drugs. *Nat Chem.* 2012;4(2):90-98. doi:10.1038/nchem.1243  
7. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD.* 2016. doi:10.1145/2939672.2939785  
8. Suo Y, Fedor JG, Zhang H, et al. *Nat Commun.* 2025;16:5178. doi:10.1038/s41467-025-60480-3 （PDB 9DKB；Suo 与 Fedor 同等贡献）  
9. Dekker C, Mattes H, Wright M, et al. Crystal structure of NLRP3 NACHT domain with an inhibitor. *J Mol Biol.* 2021;433(24):167309. doi:10.1016/j.jmb.2021.167309 （PDB 7ALV；配体为 NP3-146）  
10. O’Boyle NM, Banck M, James CA, et al. Open Babel: an open chemical toolbox. *J Cheminform.* 2011;3:33. doi:10.1186/1758-2946-3-33  
11. Trott O, Olson AJ. AutoDock Vina. *J Comput Chem.* 2010;31(2):455-461. doi:10.1002/jcc.21334  
12. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0. *J Chem Inf Model.* 2021;61(8):3891-3898. doi:10.1021/acs.jcim.1c00203  
13. McNutt AT, Francoeur P, Aggarwal R, et al. GNINA 1.0: molecular docking with deep learning. *J Cheminform.* 2021;13:43. doi:10.1186/s13321-021-00522-2  
14. Shen C, Zhang X, Deng Y, et al. Boosting protein–ligand binding pose prediction and virtual screening based on residue–atom distance likelihood potential and graph transformer. *J Med Chem.* 2022;65(15):10691-10706. doi:10.1021/acs.jmedchem.2c00991 （RTMScore）  
15. Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability. *Adv Drug Deliv Rev.* 2001;46(1-3):3-26. doi:10.1016/s0169-409x(00)00129-0  
16. Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD. Molecular properties that influence the oral bioavailability of drug candidates. *J Med Chem.* 2002;45(12):2615-2623. doi:10.1021/jm020017n  
17. Baell JB, Holloway GA. New substructure filters for removal of pan assay interference compounds (PAINS). *J Med Chem.* 2010;53(7):2719-2740. doi:10.1021/jm901137j  
18. Brenk R, Schipani A, James D, et al. Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. *ChemMedChem.* 2008;3(3):435-444. doi:10.1002/cmdc.200700139  
19. Ghose AK, Viswanadhan VN, Wendoloski JJ. A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. *J Comb Chem.* 1999;1(1):55-68. doi:10.1021/cc9800071  
