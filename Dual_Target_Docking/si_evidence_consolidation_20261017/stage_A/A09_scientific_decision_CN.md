# 阶段 A 科学决策

判定：`STOP_NO_QUALIFYING_EXTERNAL_DATA`

这个判定只说明：按本轮预先写明的 E2 定义和每类至少 20 个合格分子的工作门槛，没有一套数据可以进入方案 C 的外部对接。它不授权阶段 B，也不授权任何对接或 AUROC。

## 基准

远程分支 `cursor/v4-2-final-audited-results-20260930` 的 HEAD 是 `7777e42e8f7f204d4ba275c7f1cbefc5f8f6bd5b`（2026-10-11）。默认 `git fetch origin` 的 refspec 只跟踪另一条分支，所以随后按该提交单独取了这条远程引用，没有修改 git 配置。本地对照提交与远程 HEAD 一致。工作树 `/mnt/d/CADD paper exercise/dual target docking/repo_jcim_plan_c` 从该提交新建，分支为 `cursor/jcim-plan-c-external-validation-20261012`。原有约 1904 个未提交改动的工作区没有被切换或清理。

从正式人群表复核：`activity_eligible=1` 且 `parent_independent=1` 的观察是 805 条。交叉表中的不重复 ChEMBL 分子是 783 个。方向性标签记录是 938 条（dual 对 A-only 472，dual 对 B-only 466）。其中 M0 两端或相应端点有效的记录是 934 条。化学模型共同折外方向记录是 928 条。正式活性规则是 pChEMBL ≥ 6.0。bootstrap 为 B = 10000、seed = 271828。八对受体与正式人群表一致。这些数字与冻结说明相符，正式评价集没有重建。

去盐母体的连接性用 InChIKey 第一段表示。正式 805 条的冻结实体编号，其连接性段与异立体 SMILES 的 InChIKey 连接性段一致（805/805）。完整立体身份用冻结实体编号的全长 InChIKey，它与异立体 SMILES 的 InChIKey 一致（805/805）。ECFP4 使用半径 2、2048 位、不用手性、不用特征。E3 没有用任何 Tanimoto 阈值自动赋值。

## 1. 实际检查过的数据源

Davis 等，*Nature Biotechnology*，2011，DOI `10.1038/nbt.1990`，PMID 22037378。读到的活性表是 Supplementary Table 4（MOESM5）：72 个化合物对 442 个激酶实验的 Kd，单位 nM。Supplementary Table 1（MOESM3）是激酶水平的 Skinase 汇总，不是配体活性矩阵，不参与分类。

Drewry 等，*PLOS ONE*，2017，DOI `10.1371/journal.pone.0181585`，PMID 28767711。S4 表有 646 个化合物的 SMILES，以及 1 µM、单次的百分抑制。表名是 `Table 4 - PKIS2 %Inh`。论文写明该读数由 DiscoverX 的 percent control 换算而来。

Klaeger 等，*Science*，2017，DOI `10.1126/science.aan4368`，PMID 29191878，PMCID PMC6542668。期刊页返回 403，Kd_app 工作簿没有拿到。正式活性长表里已有引用该 DOI 的 9 条 Kd 记录、8 个面板分子。

Karaman 等，2008，正式表中的 DOI 写作 `10.1038/nbt1358`，PMID 18183025。原始矩阵没有重新下载。正式活性长表里有 49 条 Kd 记录、8 个面板分子。

BindingDB REST `getLigandsByUniprots`，2026-10-12 下载，亲和力截止值 1,000,000 nM。十四个人源 UniProt 返回的记录，在按靶点名称过滤并去重后为 171001 条。字段只有 query、monomerid、SMILES、affinity type、affinity、PMID、DOI。

PubChem BioAssay 的原始工作表没有下载，没有进入四状态计数。

没有为了凑数据而新增靶点对。Davis 和 PKIS2 里的其他激酶对使用的是同一批已经判定为 E1 或不可按 θ = 6.0 分类的矩阵，不能据此扩展正式八对。

## 2. E0–E3

没有集合被评为 E2 或 E3。E3 需要在 E2 之上另有化学系列独立性和测定可靠性的逐项证据。由于没有 E2，E3 无从升级。相似度数字只作为化学轴的描述。

Davis 2011 是 E1。DOI 和 PMID 已经出现在正式 805 条所依据的 ChEMBL 37 活性行中：103 条 Kd 行、11 个面板分子。这是同一 KINOMEscan 实验系列，不是独立外部实验。其中被正式标签直接使用的那些测量属于 E0 式的重复收录；同一论文里其余化合物仍是 E1。MOESM4/MOESM5 只有化合物名称，没有 SMILES。本机对 PubChem 名称解析没有返回，因此 72 个化合物的结构重叠没有算完。不能把“结构重叠计数为 0”理解成已经证明分子独立。

Karaman 2008 是 E0。这 49 条 Kd 已经在正式活性提取里。

Klaeger 2017 是 E1。文献标识已经在正式活性行中。即使补充表后来被读到，同一论文的其他化合物也仍然是相关来源。该文的终点是 Kd_app，正式长表把引用行存成 Kd。本轮没有用这个存储名把 Kd_app 和生化 Kd 合并。

PKIS2 S4 的 DOI 不在正式 805 条的活性行中，实验来源这一轴与上述三篇不同。分子身份并不干净：去盐母体连接性与正式集重合的有 2 个，即 PFE-PKIS 29（CHEMBL1234354，PM48_08）和 UNC10225449A（CHEMBL29197，EH120_054）。EGFR/HER2 与 PIK3CA/mTOR 上的最大 ECFP4 Tanimoto 为 1.0，与这两个母体重合相符。共同 Bemis–Murcko 骨架在 EGFR/HER2 为 7，在 PIK3CA/mTOR 为 2，在 JAK1/TYK2 为 1，在 JAK1/JAK2 为 0。因为 1 µM 百分抑制不能按 θ = 6.0 分类，这套数据的总评是 `NOT_CLASSIFIABLE`，不是 E2。

BindingDB 中，DOI 或 PMID 已经出现在正式活性行里的文献队列评为 E1。只有 `10.7270` 数据集 DOI、没有 PMID 的队列，以及虽有 PMID 但接口没有突变体、结构域、物种回显和原始行位置的队列，评为 `UNRESOLVED_PROVENANCE`。这两种都不是 E2。接口还把非请求靶点夹在返回值里，例如 EGFR 文件中有 780 条名称为 ERBB2 的记录，PIK3CA 文件中有调节亚基 p85，PPARG 文件中有 GAL4。这些行按名称排除，没有当成请求靶点的活性。

## 3. 四状态人数

分类只在同一终点家族、同一 DOI 内进行。活性定义为 p ≥ 6，即数值型亲和力 ≤ 1000 nM。不等式只有在不可能跨过 1000 nM 时才分类。一端缺失记为测量不完整，不记为 neither。空白若不能唯一决定类别，记为 `INDETERMINATE`。

Davis，仅两端都是数值 Kd 时分类。Nature 补充表说明写明：空白表示已测试但结合弱（Kd > 10 µM），或者 10 µM 初筛未检出。这两种含义在单元格里没有分开，所以空白全部是 `INDETERMINATE`。

| 靶点对 | dual | A-only | B-only | neither | 无法分类 |
| --- | ---: | ---: | ---: | ---: | ---: |
| EGFR/HER2 | 7 | 7 | 1 | 1 | 56 |
| JAK1/JAK2 | 13 | 0 | 4 | 1 | 54 |
| JAK1/TYK2 | 12 | 1 | 3 | 2 | 54 |
| PIK3CA/mTOR | 4 | 0 | 0 | 0 | 68 |

AChE/BChE、F2/F10 和两个 PPAR 对不在该激酶面板中。mTOR 使用 Entrez 符号 FRAP1、激酶名 MTOR 的野生型行；PIK3CA 只用名称恰好为 PIK3CA 的行，不用突变体行。四个激酶对的三个方向性类别都低于 20。

PKIS2 不进入上表。四个激酶对各有 645 个化合物两端都有百分抑制读数，1 个化合物缺一端。这些读数保持为单浓度数据。

BindingDB 里，dual、A-only、B-only 都不少于 20 的队列只有两个，而且是同一篇论文的两种靶点对切分：

- PPARG/PPARA，IC50，DOI `10.1016/j.bmc.2012.04.005`，PMID 22564380：dual 45，A-only 50，B-only 36，neither 35。
- PPARA/PPARD，同一 DOI 和终点：dual 46，A-only 35，B-only 73，neither 12。

这篇 *Bioorganic & Medicinal Chemistry* 2012 年文章的标题是 Activity landscape modeling of PPAR ligands with dual-activity difference maps。正式活性长表里已有它的 98 条 IC50、39 个面板分子。BindingDB 队列与正式母体的重合是 33。总评 E1。PPAR 的正式药理核查已经说明激动、结合和方向不明不能混用；BindingDB 的 IC50 字段没有给出作用方式。因此人数达标不改变 E1。

其余带 PMID、且不是 `10.7270` 数据 DOI 的队列，三个方向性类别的最小值最高只有 6。例如 AChE/BChE 的 `10.1021/jm050578p` 是 11/11/6，而且已是 E1。没有 PMID 的 `10.7270` 队列可以出现数百个 dual，但 B-only 经常为 0，三个方向性类别不能同时达到 20。这些队列本身也缺少原始论文行位置。

Klaeger 的 Kd_app 矩阵没有在手，四状态人数不编造。PIK3CA 或 mTOR 上“没有观察到结合”没有被当成无活性。

## 4. 分子、来源和骨架重叠

实验来源重叠、分子身份重叠和化学骨架重叠是三列，不合成一个相似度门槛。

来源重叠的直接证据是正式活性长表中的 DOI 和 PMID。Davis、Karaman、Klaeger 和 2012 年 PPAR 景观论文都在其中。PKIS2 的 DOI 不在这 3555 行里。

分子身份重叠以去盐母体连接性为准。PKIS2 有 2 个母体与正式集重合，见上。Davis 的 72 个名称没有完成结构对照。BindingDB 的重合数写在 A05 每个队列的 `n_parent_overlap_any_formal`。2012 年 PPAR 队列是 33。

化学轴只对有结构的来源计算。PKIS2 的最大 Tanimoto 和共同骨架数见 A04。没有 E2/E3 集合，所以没有做分子量、MolLogP 和 TPSA 的匹配可行性分析，也没有做一对一匹配。

样本量是否达到每类 20 人，与上述三轴分开写在 A05 的 `gate_20_each_directional_class`。

## 5. 结构与软件

八对的冻结受体和盒子仍在正式人群表中：EGFR/HER2 为 3POZ/3RCD，JAK1/JAK2 为 6N7A/8BXH，JAK1/TYK2 为 6N7A/3LXP，PIK3CA/mTOR 为 4L23/4JT6，AChE/BChE 为 4EY7/4BDS，F2/F10 为 4UDW/2JKH，PPARG/PPARA 为 9V8H/6LXA，PPARA/PPARD 为 6LXA/5U3Q。

`receptor_component_policy.csv` 里读到的特殊组分是：4EY7 二聚体保留、配体 E20 的 A 链用于盒子；4L23 的 p85 删除；4JT6 的第二拷贝和 mLST8 删除；4UDW 的水蛭素保留；2JKH 的 FXa 轻链保留；9V8H 的肽段保留；5U3Q 的 PPARD B 链为正式拷贝。本轮没有新选受体，也没有改盒子。

当前环境有 AutoDock Vina 1.2.7 和 Meeko 0.7.1，以及 `mk_prepare_ligand.py`。RDKit 为 2025.09.6。没有找到 GNINA 和 Open Babel。M0 若在以后得到授权，Vina 和 Meeko 是已有工具；GNINA 不在本机。大环规则沿用已冻结的 M3 增补：Meeko 0.7.1 的 `--rigid_macrocycles`，真实原子保留，伪原子只是闭环胶水，删掉伪原子后的开环结构不是原大环分子。本轮没有配体需要按这条规则处理，因为没有准备新配体，也没有运行对接。

## 6. 最接近阶段 B 的候选，以及没有通过的原因

人数上最接近预先门槛的是 2012 年 PPAR IC50 队列。它没有通过，因为文献和 39 个正式面板分子已经在开发数据的活性提取里，等级是 E1，而且接口没有药理作用方式。

Davis 的数值 Kd 子集方向性类别最多只有 13 人，空白又不能分类，文献本身是 E1。

PKIS2 覆盖面大，但终点是 1 µM 百分抑制，不能换算成 pChEMBL，也不能为了达到 READY 而另设百分抑制阈值。

Klaeger 的原始 Kd_app 表不在手，且 DOI 已在正式提取中。

BindingDB 里来源尚未核完的期刊队列，限制性类别不超过 6 人。即使以后读到原文并把等级从 `UNRESOLVED_PROVENANCE` 升到 E2，也过不了每类 20 人的门槛。`10.7270` 的大队列缺少 PMID，而且 B-only 侧为空，不能当成双方向供给。

因此没有一个集合同时具备 E2 和三个方向性类别各不少于 20 人。

## 7. 未通过准入的原因

未通过的共同原因是下面几条里的至少一条：实验记录已经在正式 ChEMBL 37 提取中；终点不是可按 θ = 6.0 分类的双端摩尔亲和力；空白或缺失不能唯一分类；接口没有原始工作表行、物种回显和突变体信息；或者三个方向性类别没有同时达到 20。人数大的队列没有因此被降级标准放行。单方向人数较多的队列，例如若干 JAK 的 `10.7270` IC50 队列，只说明 dual 很多，不能称为完整双方向验证。

## 8. 方案 C 是否具备继续的科学条件

不具备把现有发现放到独立外部对接里检验的条件。可以写进论文的阶段性事实是：在本轮能够追溯到原始表或正式活性行的公开数据中，平衡的双端摩尔活性要么已经进入正式提取，要么在同一终点、同一文献内达不到预先的人数门槛，要么根本不是摩尔亲和力。这个阴性可行性结果本身有方法学意义，但它不是外部验证已经完成，也不是方向性评价结论已经被独立数据重复。

`READY_FOR_STAGE_B_DESIGN` 不成立。`EXPLORATORY_ONLY` 要求已经存在 E2 或 E3 而只是人数不足，本轮没有这样的集合。`BLOCKED_UNRESOLVED_PROVENANCE` 也不作为总判定：尚未核完的 BindingDB 队列，其已有计数离三个方向性类别各 20 人还有明显距离，或者根本没有 PMID。把总判定写成 BLOCKED，会暗示只要补齐来源就可能进入外部对接。当前计数不支持这种预期。

总判定保持 `STOP_NO_QUALIFYING_EXTERNAL_DATA`。阶段 B 至 E 不启动。
