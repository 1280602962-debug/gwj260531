# 面向 Molecular Diversity 的 URAT1–NLRP3 双靶筛选重构与全文计划

## 一、总体判断

当前项目不必改成一篇纯方法学论文，也不必放弃“双靶候选筛选”主线。需要修改的是证据表述和筛选架构：

- 不再声称单一对接分数能够预测或排序实验活性；
- 仍可用经过靶点特异回顾验证的结构计算来缩小候选范围；
- 双靶候选由多个相互独立的证据共同支持，而不是由一个对接分决定；
- 文章结论是“筛选并优先排序了可实验验证的双靶候选”，不是“发现了已验证双靶抑制剂”。

建议文章类型定为：**带有方法学严谨性的计算药物重定位与候选发现论文**。方法学是可信度来源，但候选发现仍是主结果。

## 二、新的一句话故事线

> 针对 URAT1 与 NLRP3 数据和结构证据不对称的问题，本研究建立了靶点特异、偏倚审计和多证据一致的临床阶段化合物筛选流程，在不把对接分解释为活性值的前提下，从 8,319 个分子中优先提出具有双靶结构兼容性、可接受药物化学性质和互补生物学证据的可验证候选。

## 三、对现有五种读出的重新解释

现有 P0–P5 不应统称为五个完整“对接协议”，而应拆成姿态来源和评分读出：

| 原编号 | 实际角色 | 新处理 |
|---|---|---|
| P0 GNINA CNNscore | 更接近姿态/复合物合理性判别信号 | 从“负对照”改为姿态选择和质控指标，不直接当亲和力 |
| P1 Vina affinity | 物理经验评分基线 | 保留作正交基线和共识输入，不单独决定候选 |
| P2 GNINA CNNaffinity | 学习型库内排序信号 | 保留为主排序分量，但不称实验活性或亲和力 |
| P3 minimizedAffinity | 局部优化后经验能量 | 当前富集差，可降为补充材料 |
| P4 RTMScore@Vina pose | Vina姿态的学习型重打分 | 覆盖不足且富集弱，降为补充材料 |
| P5 RTMScore@GNINA pose | GNINA姿态的学习型重打分 | 保留作敏感性分析；RandomDecoy失败不能简单用“0命中”一票否决 |

核心改变：**用 CNNscore 选姿态，用 CNNaffinity/正交评分排分子，用靶点相互作用和偏倚校正做最终确认。**

## 四、推荐主方案：GNINA 保留，改为靶点特异的分层共识筛选

这是最推荐、与现有资产兼容度最高、最容易在 Molecular Diversity 完成全文的一条路线。

### 4.1 模块 A：统一并冻结真正的生产设置

基准和生产必须使用完全相同的参数：

- 相同受体文件；
- 相同搜索盒中心和尺寸；
- 相同 GNINA 版本与 CNN ensemble；
- 相同 exhaustiveness；
- 相同质子化/互变异构体规则；
- 每个分子保留 9–20 个构象，而不是生产阶段只保留 1 个；
- 固定至少 3 个随机种子进行采样稳定性检查。

推荐：`exhaustiveness=32` 可保留；`num_modes` 改为 9 或 20。每个分子先由 CNNscore 选择合理姿态，再读取该姿态的 CNNaffinity，而不是让“只输出一个构象”同时承担搜索和排序两种任务。

### 4.2 模块 B：两靶分别验证，不再把 URAT1 结果直接等同于 NLRP3

#### URAT1

- 保留 9DKB 主结构；
- 使用 lesinurad 做自对接；
- 报告 Top-1、Top-3、best-of-N RMSD；
- 使用现有 469 活性物和两类诱饵做早期富集；
- 增加按电荷、分子量和骨架的偏倚审计；
- 可把其他 URAT1 构象作为候选级敏感性分析，而非重新扩大生产主线。

#### NLRP3

- 必须完成 NP3-146@7ALV 自对接；
- MCC950只能作相似药理对照，不能代替共晶配体自对接；
- 从文献/ChEMBL中另建“直接或位点相关 NLRP3 抑制剂集合”；
- 细胞 IL-1 数据和直接结合/位点数据分开，前者用于生物学缩库，后者用于结构协议验证；
- 如果可靠的直接阴性分子不足，使用性质匹配诱饵并明确其为推定阴性，同时把 NLRP3 结构验证定为探索性。

### 4.3 模块 C：新的协议评价指标

不再用“TrueDecoy EF@1%最高且 RandomDecoy 非零”作为唯一选择规则。推荐指标矩阵：

#### 排序能力

- BEDROC，预先固定 α=80.5 与 160.9；
- adjusted LogAUC；
- EF@1%、EF@5%；
- hits@1%、hits@5%；
- ROC-AUC与PR-AUC作为全局辅助指标。

#### 姿态能力

- Top-1 RMSD≤2 Å成功率；
- Top-3和best-of-N成功率；
- 重原子冲突；
- 口袋占位与关键相互作用恢复率。

#### 稳健性

- 3个随机种子的Spearman相关；
- 前1%、5%、10%集合Jaccard；
- 质子化/互变异构体改变后的排名稳定性；
- 按Bemis–Murcko骨架进行配对cluster bootstrap。

#### 工程可用性

- 覆盖率；
- 失败率；
- 单分子计算成本。

若多个读出差异不显著，选择覆盖最完整、计算最简单、跨两类诱饵最稳定的读出，而不是声称它在统计上“最优”。

### 4.4 模块 D：显式处理大分子偏倚

裸 Pareto 被大环和高分子量分子占据，说明对接分存在明显尺寸偏倚。建议同时生成三种排序：

1. 原始 CNNaffinity 百分位；
2. 配体效率型分数，如分数/重原子数，只作敏感性分析；
3. 物化残差分数：在完整临床池内用分子量、重原子数、cLogP、形式电荷、可旋转键数解释对接分，使用残差表示“超过同类物化性质预期的结构分”。

最终候选要求在原始排名与至少一种偏倚校正排名中均保持较高位置。不能仅因校正后某候选上升而重新定义规则。

### 4.5 模块 E：构建可解释的双靶结构证据

每个靶点不再只有一个分数，而有三个分量：

1. **Rank evidence**：CNNaffinity/BEDROC验证后的库内百分位；
2. **Pose evidence**：CNNscore、无冲突、口袋占位；
3. **Interaction evidence**：与共晶配体或已知药效团的一致性。

URAT1可重点检查：

- 芳香笼 Phe241/Phe360/Phe364/Phe365/Phe449；
- Arg477与酸性基团；
- 与lesinurad相似的口袋空间占位。

NLRP3可重点检查：

- 7ALV中央口袋占位；
- 与NP3-146已知关键相互作用的恢复；
- Ala227/Ala228、Arg351、Met408、Tyr443、Phe575、Arg578附近的接触模式。

相互作用不应简单按数量求和。建议使用二元必要特征＋指纹相似度，并通过已知配体回顾验证其区分能力。

### 4.6 模块 F：双靶提名规则

不建立不可解释的总分。使用分层门控：

#### 第一门：生物学相关性

- NLRP3细胞活性模型处于可靠适用范围；或已有明确NLRP3/炎症药理证据；
- 对模型域外分子不直接删除，但标记低置信度。

#### 第二门：双靶结构一致性

- 两靶原始排名均进入前10%；
- 两靶至少在一种偏倚校正排名中进入前15%；
- 两靶姿态均通过无冲突和口袋占位质控；
- 关键相互作用不与已知结构机制明显矛盾。

#### 第三门：跨设置稳定性

- 三随机种子中至少2次通过；
- 质子化/互变异构体敏感性不导致候选跌出前20%；
- 对候选级结构扰动仍保留主要口袋接触。

#### 第四门：药物化学与转化

- PAINS、Brenk、聚集和反应性风险；
- 溶解性、可购性、临床暴露；
- 已知靶点是否可能解释NLRP3细胞读出；
- 骨架去冗余。

### 4.7 候选分层

- **Tier 1：双证据一致候选**——Vecabrutinib优先；要求重算后仍通过两靶结构和NLRP3生物学证据。
- **Tier 2：结构驱动检验候选**——GSK-3008348；明确其NLRP3 ML证据弱，用于检验结构筛选能否发现模型未支持的候选。
- **Tier 3：备选化学型**——Praliciguat、MLN-0415、Zelenirstat等，按重算稳定性重新确定，不保证原名单全部保留。
- **Bias controls**——大环Pareto成员，用来展示未校正对接分的尺寸偏倚。
- **Reference controls**——lesinurad、verinurad、MCC950、NP3-146。

## 五、如果 GNINA 确认性验证仍失败，替换什么

### 路线 B：物理采样＋AI重打分

推荐组合：

- 姿态生成：AutoDock Vina/Vinardo、rDock、LeDock或DOCK系列；
- 姿态过滤：PoseBusters式物理合理性检查＋关键相互作用；
- 重打分：GNINA CNNaffinity、RTMScore、EquiScore中选择经过本靶点回顾验证者；
- 最终排序：至少一个物理经验分＋一个学习型分形成rank consensus。

优点是评分来源更正交；缺点是工作量增加。不能把多个失败评分简单平均，所有进入共识的分量都应先通过最低验证门。

### 路线 C：KarmaDock或CarsiDock生成姿态＋正交重打分

适用于GNINA采样明显失败时：

- KarmaDock、CarsiDock生成候选姿态；
- 不直接把其预测分当实验亲和力；
- 对姿态做物理合理性过滤；
- 用RTMScore、GNINA score-only或物理评分重排；
- 必须在URAT1和NLRP3各自完成自对接和活性物–诱饵验证。

这两个方法可以作为候选级确认工具，但不建议未经本靶点验证直接替换P2并重筛全库。近期综合基准显示，AI方法可能有较好的构象准确性，但物理合理性并不总是优于传统方法，因此仍需正交检查。

### 路线 D：结构药效团/Open-ComBind式相互作用共识

如果任一打分函数都没有稳定富集，可降低对绝对评分函数的依赖：

- 以lesinurad@9DKB和NP3-146@7ALV建立结构药效团；
- 用已知同靶配体形成相互作用指纹共识；
- 多姿态对接后按药效团满足度和相互作用相似性筛选；
- 对接能量只设宽松门槛，不负责精细排序。

这条路线特别适合已有明确共晶配体但通用评分函数表现差的靶点。Open-ComBind思路也可作为参考：利用一组已知配体间共享相互作用信息选择姿态，而不是完全依赖单个通用打分函数。

### 不推荐作为直接替代

- DiffDock等只用于姿态生成，不能单独承担活性排序；
- 单次MM-GBSA不能取代经验证的虚拟筛选协议；
- 一条短MD轨迹不能证明结合稳定或亲和力；
- 把P1–P5简单求平均不会自动得到可靠共识；
- 在看到最终候选后反复更换协议，容易造成选择偏倚。

## 六、MD与后处理如何放入文章

Molecular Diversity并不要求每篇计算筛选文章必须有MD，但加入规范的候选级动态验证可增强化学与结构解释。

建议仅对最终2个主候选和两靶参考配体执行：

- 体系：2候选×2靶点＋每靶1个参考，共6个复合物；
- 每体系3条独立重复；
- 若资源充足，每条100 ns；资源有限可做3×50 ns并明确探索性质；
- 报告配体RMSD、口袋保留率、关键接触占据、重复间一致性；
- MM/PBSA或MM/GBSA仅作同体系内补充描述，不写成实验结合自由能；
- 预先设定分析指标，避免只展示最稳定的一条轨迹。

MD不是候选入选的首要门控，而是对冻结候选的结构解释。

## 七、面向 Molecular Diversity 的完整论文结构

### 推荐题目

**A target-asymmetric and bias-audited virtual screening framework prioritizes clinical-stage candidates for coordinated URAT1–NLRP3 modulation in gout**

更偏候选发现的备选题目：

**Multi-evidence virtual screening of clinical-stage compounds identifies testable URAT1–NLRP3 dual-node candidates for gout**

### 摘要逻辑

1. 痛风同时涉及尿酸重吸收和炎性小体激活；
2. URAT1和NLRP3证据不对称，使传统同权双靶评分不可靠；
3. 建立靶点特异、偏倚审计、早期富集验证和多证据门控流程；
4. 8,319→生物学缩库→双靶结构筛选→药化/稳定性筛选；
5. 得到2个主候选及若干备选化学型；
6. 结论为可实验验证的双节点候选，而非已证实双靶抑制剂。

### 1 Introduction

- 痛风的代谢–炎症双过程；
- URAT1与NLRP3的互补临床意义；
- 单靶药物和组合治疗的不足；
- 计算双靶筛选的评分偏倚、证据不对称和假阳性问题；
- 本研究目的：从临床阶段分子中提出可验证双靶候选，同时建立可信的筛选证据链。

### 2 Materials and methods

1. Study design and preregistered information flow
2. Target evidence and compound data curation
3. Clinical-stage library and structure standardization
4. URAT1 and NLRP3 benchmark construction
5. Protein and ligand preparation
6. Multi-pose docking and target-specific pose validation
7. Early-recognition metrics and scaffold bootstrap
8. Physicochemical-bias audit and corrected ranks
9. NLRP3 biological relevance model
10. Dual-target evidence gates and candidate nomination
11. Chemical diversity and scaffold analysis
12. Candidate-level dynamic/interaction validation
13. Statistics and reproducibility

### 3 Results

#### 3.1 Curated evidence reveals target asymmetry

展示两靶数据规模、测定层级、骨架和化学空间差异。

#### 3.2 Target-specific validation distinguishes pose and ranking performance

分别报告两靶自对接、早期富集、BEDROC、LogAUC、EF和覆盖率；明确没有一个通用分能预测实验活性，但可筛出高优先级区域。

#### 3.3 Bias audit supports a corrected multi-evidence workflow

展示分数与MW、重原子数、电荷的关系；说明为什么裸Pareto富集大环；给出偏倚校正后的稳定性。

#### 3.4 Biological triage and dual-target screening reduce the clinical library

展示8319→缩库→完整案例→双结构门控的漏斗。

#### 3.5 Consensus structural evidence prioritizes diverse candidate chemotypes

展示原始分、校正分、姿态、相互作用和跨种子稳定性；给出Tier 1–3候选。

#### 3.6 Vecabrutinib and GSK-3008348 represent complementary testable hypotheses

- Vecabrutinib：结构＋NLRP3生物学一致；
- GSK-3008348：结构驱动、NLRP3模型不支持；
- 二者代表两种不同可证伪机制，而不是为了凑出“双靶hit”。

#### 3.7 Candidate-level structural dynamics/orthogonal validation

若完成MD或第二姿态方法，在此报告；未完成时可改为多姿态和相互作用指纹稳定性。

### 4 Discussion

- 双靶候选筛选仍然是主贡献；
- 贡献不是“某个分数预测活性”，而是对不可靠分数进行靶点特异校准与证据汇合；
- Vecabrutinib和GSK-3008348的不同证据类型；
- 通用打分迁移、推定阴性诱饵、NLRP3细胞终点和缺乏湿实验的局限；
- 清晰提出URAT1转运和NLRP3炎症实验验证顺序。

### 5 Conclusion

得到的是具有双靶结构兼容性和互补生物学证据的实验优先候选。不能写“确认双靶活性”，可以写“prioritized testable dual-node candidates”。

## 八、推荐图表

### 主图

1. 图1：研究流程与单向信息流；
2. 图2：两靶数据、骨架与证据不对称；
3. 图3：协议姿态能力、BEDROC/LogAUC/EF及bootstrap；
4. 图4：分数–MW/电荷偏倚与校正效果；
5. 图5：8319到最终候选的筛选漏斗和双靶二维景观；
6. 图6：最终候选结构、骨架多样性与证据雷达/矩阵；
7. 图7：两个主候选在两靶口袋中的相互作用或动态稳定性。

### 主表

1. 表1：数据集与用途；
2. 表2：两靶结构和验证结果；
3. 表3：对接/评分读出比较；
4. 表4：候选提名规则与每层数量；
5. 表5：最终候选的结构、生物学、药化和转化证据。

## 九、投稿前最低补算包

### 必须完成

1. 基准和生产参数完全统一；
2. URAT1重新计算多姿态版本；
3. NLRP3 NP3-146自对接；
4. BEDROC、LogAUC、EF@1/5%、PR-AUC；
5. 配对骨架bootstrap；
6. 电荷、MW、重原子数和cLogP偏倚审计；
7. 原始排名、配体效率和物化残差排名敏感性；
8. 7个候选的多种子、多质子化状态稳定性；
9. 两靶结构相互作用指纹；
10. 按新门控重新冻结候选，不保证原7个名单不变。

### 强烈建议

1. NLRP3直接/位点相关活性物基准；
2. 第二种姿态生成或重打分方法用于候选级确认；
3. 最终2个候选与参考配体的重复MD；
4. 临床暴露、可购性和实验干扰审计。

### 可放补充材料

- P3、P4等表现较差读出；
- 所有阈值敏感性；
- 全部51个结构门控分子；
- 大环偏倚案例；
- 对接失败清单和所有姿态QC。

## 十、执行顺序与停止规则

### 第1阶段：两周内完成协议确认

- 重跑精确生产配置基准；
- 计算新增指标；
- 做NLRP3自对接；
- 判断GNINA主方案是否继续。

**继续GNINA的最低条件：** 至少在URAT1上稳定早期富集，两靶均能产生物理合理姿态，候选跨种子和偏倚校正后不是完全重排。

### 第2阶段：两至三周完成候选重筛

- 对1,580完整案例做多姿态/偏倚校正；
- 完成相互作用指纹；
- 冻结2主＋3备候选。

### 第3阶段：三至六周完成候选级确认与全文图表

- 第二方法或MD；
- 化学空间和药化分析；
- 自动生成全部主图和表格。

### 停止或换路线条件

如果精确生产配置下URAT1的BEDROC/LogAUC、EF@1/5%都与随机不可区分，或候选在三个随机种子和质子化状态下完全不稳定，则停止把GNINA分数作为排序核心，转入“结构药效团/Open-ComBind式相互作用共识”或“新姿态生成器＋正交重打分”路线。

## 十一、最终推荐

首选不是立刻替换所有软件，而是：

1. 保留GNINA作为多姿态生成器；
2. CNNscore负责姿态选择和质控；
3. CNNaffinity只负责库内排序分量；
4. 增加BEDROC/LogAUC、骨架bootstrap和物化偏倚校正；
5. URAT1和NLRP3分别完成验证；
6. 最终以结构排名、姿态、相互作用、生物学与药化一致性筛选双靶候选。

这条路线既能保留现有8319→1588→1580的主要工作，也能保留Vecabrutinib/GSK-3008348作为待确认候选，同时把文章从“对接分挑药”升级为符合Molecular Diversity范围的“分子多样性、靶点不对称和偏倚审计驱动的双靶候选发现”。

---

## 十二、执行层：哪些必须本地（C1 战役）

产品目标仍是**可测的双靶候选**。本节不改第四–十一节的科学意图，只把它们接到一台有 gnina/GPU 的机器上，并写死两条轨道，避免 4.5（酸根–Arg477）与 4.7（Vecabrutinib Tier 1）同时当真。

**作战文件（读这个再开对接）：** [`LOCAL_C1_CANDIDATE_CAMPAIGN.md`](LOCAL_C1_CANDIDATE_CAMPAIGN.md)  
**科学锁：** [`config/campaign_c1.yaml`](../config/campaign_c1.yaml)  
**gnina 引擎配置：** [`config/docking_c1.yaml`](../config/docking_c1.yaml)

| 必须本地（无 gnina 做不了） | 云端可做、不算过关 | 现在不要做 |
|-----------------------------|---------------------|------------|
| L2 自对接（lesinurad 羧酸根@9DKB；NP3-146@7ALV） | 羧酸根准备脚本、SDF 读出解析器 | 覆盖冻结 `data/repurposing/p2/` |
| L3 全诱饵 URAT1（GPU；~9,849，`num_modes=9`） | 冻结 SDF 上的 CNN_VS / C1_P2star **frozen-prep** 对照 | 未过 L2 就开全诱饵 |
| L4 NLRP3 位点对接 | ChEMBL 位点阳性策展 | 与 L3 并行换 KarmaDock/DiffDock |
| L5 临床库（仅 `pass_fail.json` 之后） | 旧 1,580 行偏倚图、BEDROC SI | 把旧 7 个名字预承诺为 C1 hit |
| L7 MD（短名单冻结后；URAT1 必须膜） | Methods 空数字草稿 | 未过 Rank 轨就为 Vecabrutinib 开 URAT1 膜 MD |

**Rank 轨不过关就关闭。** Acid 轨仍可出羧酸姿态假说，但文章不能写“对接排出了活性”。Vecabrutinib 仅当 Rank 过关才可能留在 URAT1 臂；GSK-3008348 按羧酸规则重判，不因旧百分位自动入选。

冻结 P2 是九构象里 **CNNaffinity 最大**；规划书 4.1 要的是 **CNNscore 选姿后再读 CNNaffinity**。两者不是同一读出。本地若只改 `num_modes` 跑现有 `run_gnina_batch.py`，得到的仍是旧定义。细节与停止规则见 C1 执行书第 12–16 节。

本地第一次只做 L0–L2。L3 需要单独授权。`run_funnel_p2.sh` 是旧战役（`num_modes=1`），不要用来跑 C1。


