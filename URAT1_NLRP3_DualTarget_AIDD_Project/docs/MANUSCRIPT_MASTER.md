# 临床酸性化学空间中的URAT1–NLRP3双节点结构候选

状态：基于已存数据的统一工作稿，尚未投稿定稿。唯一数字源为[事实锁](../MANUSCRIPT_FACTS.yaml)及[论文数据层](../data/manuscript/README.md)。英文润色、文献逐项核验和尚未执行的MD应在证据完成后补入。

## 摘要工作稿

URAT1介导的尿酸转运与NLRP3相关炎症提供了双节点候选研究的背景。本研究使用NLRP3相关活性分类器将8,319个临床库分子缩减至1,588个，再通过酸性药效团及药化规则获得303个酸性等价物和156个对接分子。URAT1与NLRP3结构评估分别采用9DKB和7ALV。URAT1自由自对接的Top-1姿势恢复不足，故采用晶体相对几何门；A1在228个羧酸活性分子与64个羧酸诱饵上获得LR+ 2.20。NLRP3结构门在seed42下保留9/9个家族集中的阳性与11/40个匹配诱饵；排除5个对接失败诱饵后为11/35。多种子交集及记录明确的药化审计形成12个冻结primary，另列21个reserve。PF-03882845被选为后续验证优先对象。目前结果支持可检验的结构假设，没有直接双靶活性或动态模拟结果。

## 研究设计

主线包含六个证据包：数据与NLRP3模型、酸性化学空间、URAT1结构验证、NLRP3结构验证、双臂交集及候选冻结、PF-03882845后续验证计划。P2分数迁移分析和URAT1回归模型的阴性结果用于解释方法选择，放在SI。

## Methods

### 数据与模型

NLRP3训练数据包含609条记录、513个唯一标准化分子及39个assay。原始描述、系统和细胞信息见[assay provenance](../data/manuscript/si/nlrp3_assay_provenance.csv)。该数据支持NLRP3相关活性分类，不能直接作为NACHT结合标签。模型使用Morgan2048与RDKit描述符，为频率最高的25个assay分别编码，另加一个other维度（共26维）。按分子骨架进行五折交叉验证，同一分子的记录留在同一折。分类活性阈值为pActivity≥6.0。临床库推断在5个冻结assay上下文中取最大输出，模型分数≥0.5的已存池为生产成员依据。

39个assay、25个编码类别和5个推断上下文是不同口径。存档CV的AUROC为0.8934、AUPRC为0.9143、EF@10%为1.5706。重新训练的5.5/6.0/6.5阈值敏感性单列SI；其中6.0重训池为1,377，不能替换1,588个生产池。模型文件与1,588分数表于2026-07-01一并提交至本仓库；2026-08未重打分。screening summary未写入模型哈希或软件环境锁，故不以当前环境比特级重现该池。已存哈希见[模型登记](../data/frozen/model_manifest.json)。检索记录见[环境检索](../data/manuscript/audit/nlrp3_scoring_env_search.json)。

### 酸性空间与配体准备

冻结1,588池中使用COOH/carboxylate、tetrazole及acylsulfonamide SMARTS得到303个分子，类别允许重叠。MW 200–550、TPSA≤140、可旋转键≤10、HBD≤5及HBA≤10得到156个。实际脚本没有logP硬阈值，不声称完整Ro5通过。配体按pH 7.4经Dimorphite-DL及Meeko准备，使用清单选定微观状态，并非合理状态穷举。描述符和Murcko骨架在整理时重算为注释，不改动冻结名单。

酸等价物入池规则与姿势规则覆盖范围不同：已执行酸基识别仅匹配COOH/carboxylate氧，不识别没有COOH的四唑或酰基磺酰胺。这些分子保留在分母并单独报告，其不通过不能解释为普遍无URAT1活性。

### 受体、对接与姿势选择

生产受体为URAT1 9DKB和NLRP3 7ALV，GNINA使用exhaustiveness 32、最多9个pose及CNN rescore；seed42为初始筛选，43/44为重复。盒子及已准备受体见[docking_final](../config/docking_final.yaml)，受体用途见[receptor manifest](../data/manuscript/tables/receptor_manifest.csv)。CPU/GPU实际执行信息以原始逐任务日志为准，统一配置是整理后的描述，不替代历史日志。

其他URAT1结构用于晶体锚定、自对接审计或有限刚体姿势转移；8ETR用于NLRP3刚体转移敏感性。它们不构成全库ensemble docking。接触采用已准备受体实际编号，[residue mapping](../data/manuscript/tables/residue_mapping.csv)记录9DKB Arg477→A/476及7ALV重新分段链/编号。URAT1预期12个关键残基中仅11个匹配；Q437在9DKB为LEU，不强行计为GLN接触。

### 结构门及验证

重对接RMSD在受体固定坐标系计算，使用重原子匹配，禁止独立配体叠合降低pose RMSD。分别报告CNNscore Top-1与best-of-9；best-of-9仅用于解释采样与选择步骤差别。

A1先选CNNscore最高姿势，再检验COOH/carboxylate存在、酸氧至Arg NE/NH1/NH2最短距离≤7.7027 Å及重原子质心距离≤6 Å。7.7027来自共晶A1AIL距离6.7027 Å加1 Å裕量，仅表示晶体相对邻近性。A1本身没有IFP或显式全局clash门。A2先取几何通过姿势，再选CNNscore最高者，用于可容纳性。W2采用沉积晶体坐标锚定的IFP、接触与clash检查，是冻结名单结构注释；不能把W2条件追溯加到全部12个primary上。[规则](../config/urat1_gate_definition.yaml)

NLRP3以CNNscore最高pose为代表，要求CNNscore≥0.5、质心偏移≤6 Å、参考配体重原子邻域覆盖比例≥0.5、IFP Jaccard≥0.5、7个关键残基至少5个接触且无小于2.2 Å的重原子冲突。覆盖率指候选重原子中距参考重原子小于2.5 Å的比例，不是实际体积重叠。氢键位为距离代理，不包含完整角度验证。[规则](../config/nlrp3_gate_definition.yaml)

回顾性验证报告TP/FN/FP/TN、Wilson区间、LR+、MLE OR、Haldane校正OR和Fisher检验。保留历史bootstrap区间并标明方法，避免混合不同OR口径。W4的40个诱饵含5个三seed均失败分子，完整面板以计算不通过计，另列35个成功对接分子的complete-case结果。三个seed分别报告，不把重复pose视为独立阳性。9个阳性集中于磺酰脲家族，不能推断跨家族灵敏度。

W4诱饵来自URAT1 true-decoy数据中的性质匹配抽样，不是经实验证实的NLRP3阴性。因此这里的特异度是所构造回顾性面板的区分指标，不能直接解释为临床化学空间中排除NLRP3非活性物的概率。

### 候选形成与后续选择

156个分子中，54个在至少2/3 seed满足URAT1 A2及NLRP3结构门；与历史40个药化合格集合取交集得37个。A1 seed42通过者13个，去除结构对照得12个primary；其余24个去除3个β-lactam得21个reserve。历史40个集合包含PAINS/Brenk排除及按药名排除，NIH为注释。β-lactam采用SMARTS末级排除。见[决策台账](../data/manuscript/tables/screening_decision_ledger.csv)及[药化规则](../config/chemistry_final.yaml)。本次整理复现冻结成员，不重新赋予前瞻性预注册身份。

PF-03882845优先级综合多seed结构一致性、NLRP3关键接触与药化注释，详见[提名说明](NOMINATION_RATIONALE.md)。优先级不是效力排名。MD目前仅有[计划](../config/md_protocol.yaml)，未报告轨迹、占有率或MM/GBSA数值。

## Results

8319→1588→303→156经唯一ID和嵌套检查。303酸池含272个carboxylate、9个tetrazole、27个acylsulfonamide标记，允许重叠；156池相应为135、6、16。按当前RDKit对沉积SMILES重算，303、156、12集合分别有273、148、12种Murcko骨架；303中含3个无环分子，空骨架计一个类别。

Lesinurad@9DKB Top-1 RMSD为4.298/4.312/4.878 Å；benzbromarone@9DKA为3.585/3.603/3.595 Å，其best-of-9约1.11 Å，提示近晶体姿势可被采样但Top-1选择不足。NP3-146@7ALV为0.821/0.668/0.683 Å，均通过2 Å检验。[URAT1表](../data/manuscript/tables/urat1_redocking_summary.csv)、[NLRP3表](../data/manuscript/tables/nlrp3_selfdocking_summary.csv)

A1的TP/FN/FP/TN为102/126/13/51，LR+ 2.20；A2为217/11/61/3，LR+约1；W2为96/132/12/52，属A1更严格子集。W4 seed42结构门为9/0/11/29，完整面板特异度0.725；complete-case为9/0/11/24，特异度0.686。9/9的Wilson 95%下限约0.701。该门对20个临床酸背景未比宽松门增加逐行区分能力。[URAT1验证](../data/manuscript/tables/urat1_gate_validation.csv)、[NLRP3验证](../data/manuscript/tables/nlrp3_gate_validation.csv)

冻结12个primary可由已存集合运算复现，8个在至少2/3 seed通过W2注释门。21个reserve用于后续探索，其A2依据不支持“次优活性”解释。Lintitript对更严格Arg裕量及9DKA刚体转移敏感。PF-03882845在其他URAT1结构刚体转移后保持酸基邻近性；NLRP3 8ETR转移存在2.14 Å近接触，不能称跨受体无冲突结合。刚体转移不等于重对接或动力学稳定性。

## Discussion

价值在于明确证据适用范围并形成可实验检验候选：NLRP3模型提供缩库依据，URAT1 A1提供有限回顾富集，NLRP3结构门提供家族内相容性证据。URAT1回归与P2评分失败保留为SI的方法选择依据，不作为现行筛选分数。

限制包括：A1漏掉已知URAT1药物Top-1姿势，不能把入选者解释为优于已知降尿酸药；W4阳性少且家族集中，并含失败诱饵；酸等价物定义宽于几何匹配器；历史药化集合含人工药名排除；1,588生产池以冻结成员为准，缺少同时生成的软件环境锁，不能比特级重算；微观状态为单选；构象转移只提供静态敏感性。PF-03882845的MR拮抗作用提供上游炎症抑制解释，MD不能排除该混杂或确认直接NACHT结合。

当前支持12个冻结结构候选及明确的后续选择，不支持已确认双靶抑制剂、新靶点注释或临床治疗结论。后续需区分URAT1尿酸转运功能、NLRP3直接结合/ATPase及上游通路效应。最终摘要和结论应随真实实验/MD证据更新。

## 主图表与SI

主图1：最终漏斗和交集；主图2：酸池与primary化学空间；主表1：双靶结构验证；主表2：12个primary的分seed指标。SI包括assay provenance、阈值敏感性、URAT1 ML、P2 negative-transfer分析、reserve、泄漏/相似性、受体映射与刚体转移。旧稿的51/7集合只能作为legacy audit set出现在归档或历史SI。
