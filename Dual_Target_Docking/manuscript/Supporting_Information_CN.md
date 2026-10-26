# 补充信息

本文件说明八对双靶点结果材料的成员、统计定义、扩展结果和文件位置。它不是审计日志。结果数值来自已验收的 PRIMARY 与第二阶段表。G、J、K 没有执行。

## 1. 标签、实体和来源范围

primary 阈值为 6.0。activity_eligible=0 不是 neither。AB_040 是 AB_046 的 alias，不作为独立观察进入 869 行当前人群表。pair 内观察可以跨靶点对共享 global ligand entity，所以 867 不是 867 个独立分子。

当前标签与身份使用 results/current/pair_ligand_mapping.csv，共 870 行。PIK3CA/mTOR 使用冻结 Q110 名单，共 110 个分子，dual、A_only、B_only、neither 分别为 32、32、32、14；其他七对保留既有 V4.2 标签。

活性来源于 ChEMBL 37，使用人源、SINGLE PROTEIN、单一蛋白组分、pChEMBL 非空且标准类型为 IC50、Ki、Kd、EC50、Potency、IC50app 或 Ki app 的记录（STANDARD_OK）。同一分子和靶点取未取整的最高 pChEMBL，按 θ=6.0 确定四状态。

全部 867 条合格观察已与两端原始活性建立对应。855 条的两端全部并列最大值均带有 DOI、PubMed 或专利号；12 条缺少指定标识，包含 Q110_A21 的 mTOR 数据集文档。81 条存在跨文档聚合，123 条存在终点类型混合，并集为 168 条；35 条存在同一端跨越 6.0 的重复测定。冻结最大值及类别不一致为 0 条。上述标记可重叠，不能将各项人数相加作为独立观察数。

当前来源总表见 Table S20。101 个 Q110 分子的 238 条 ChEMBL 37 原始活性见 Table S21；其余观察通过来源总表关联已有原始记录。Q110_A24 的 mTOR 冻结最大值保持 5.83；后续审核表曾以 source_adjudication 排除该记录，具体理由仍未解决，但两种汇总下的类别均为 A_only。活性来源追溯不代表同一测定条件，也不属于外部独立验证。

## 2. 成员和缺失

第 1 层 1004 条标签方向记录供 D1 使用。第 2 层 1000 条是正式 M0 方向记录，少掉的 4 条来自 AChE/BChE 的 M0 超时。第 3 层 994 条是化学共同 OOF，在 1000 上再排除 6 条缺少冻结 fold 的方向记录。标签层缺 fold 为 10 条、7 个配体：AB_001、AB_053、AB_054、EH40_31、F2F10_080、J1TYK2_092、PGPA_030。M0 之后仍缺 fold 的是后四个配体的 6 条记录。AB_001、AB_053 和 AB_054 仍在 D1 中；它们离开第 2 层是因为超时，不是因为缺 fold。缺 fold 是没有冻结交叉验证分折，不是活性不可信。

fold_id 使用 0、1、2、3、4。0 是一个实际折号，不是缺失。AChE/BChE D_vs_A 的 A_only 由 27 到 25 发生在 M0 筛选，错口袋额外减员为 0。Table S2 给出每臂人数和名单。Table S13 中 expected = valid + missing。pair 展开记录数不等于物理 docking job 数。

## 3. PRIMARY 与配对统计

AUROC 按 higher-better 分数计算，ties=0.5。summary_min 是两个方向 AUROC 的较小值，不是对逐分子最低分再算 AUROC。D_vs_A 使用 score_B，D_vs_B 使用 score_A。已有区间来自 B=10000、seed=271828、每个分析单元重置的 default_rng，以及 linear 2.5 和 97.5 百分位。负类对照共享 dual 索引、两种负类分别抽取，差值在每次重采样内计算。口袋对照和替代受体对照使用同成员及共享索引。正文报告保留含 0 的区间，不把区间重叠写成等效，也不另作多重比较校正。

## 4. 化学模型与 D1

化学模型使用八个描述符、ECFP4，或 ECFP4 加相应方向的 M0。ECFP4 为 radius 2、2048 bit、useChirality=false、useFeatures=false。逻辑回归为 L2、C=1、liblinear、class_weight=None、max_iter=1000、random_state=271828。标准化只在训练折拟合；加入的 M0 列不使用测试折的均值和标准差。其他七对沿用原五个冻结折。Q110 使用读取分数前写定的分折：规范 Murcko 骨架去重排序后，按 0–4 循环分配，同一骨架进入同一折，共 105 个骨架。两种分折分别记录来源，均拼接折外预测得到 OOF，没有在本轮重新拟合模型。OOF 不是外部独立验证。正式输出不含 scaler_probe，训练折缩放的实现证据留在既有人工验收，不在本轮补造。

D1 在第 1 层、同一 canonical Mol 上记录八个描述符和 Murcko 字符串。空 Murcko 表示无环，不是解析失败。最大异类 ECFP4 Tanimoto 只在同一 pair×arm 内相对另一标签类计算。CAP5 的骨架定义与一致性检查见第 11 节。

## 5. 负类与口袋的完整结果

Table S5 和 Figure 2 给出 16 个臂的单靶负类、neither、差值和区间。Table S6 和 Figure 4 给出正确口袋、对侧口袋、差值和额外减员。反向点、零点和含 0 的区间都保留。

## 6. 五种子与替代受体

Table S10 为 40 个 pair×seed 原点，Table S11 为 24 个汇总。Figure S5 同时显示 auc_B、auc_A 和 summary_min。Figure 5a 只显示 summary_min。IQR 按各自 metric 的 linear P75−P25 读取。七对的五种子人数不变。AChE/BChE 仅 seed 17 因 AB_046 在 4BDS 超时而少 1 个 A_only。

Table S12 和 Figure S6 覆盖 AChE/BChE、F2/F10 和 PPARA/PPARD 三对、六个替代条件。未替换端差值为 0 是构造核对。Q110 没有替代受体计算，因此当前材料不采用 PM48 的替代受体统计。

## 7. 联合排序

方法自身描述和与 M0 的共同比较分开。共同比较中，新方法使用 pairwise_common_vs_M0，M0 使用 pairwise_common_vs_该方法，不用 method_specific 的 M0 行充当全部方法的基线。名次为 average midrank，联合名次为两端名次的 worst rank。并列先按 rank_A+rank_B，再按 global ID。k=ceil(0.10×N)。in_topk 取 YES 或 NO。Table S16 保留全部 9734 行，不在本 PDF 中逐行排印。M1、M2、M3 的方向差值读取 PRIMARY。M1b−M0 只在三组 dual、A_only、B_only 成员分别一致后，把已有的 M1b−M1 与 M1−M0 点估计相加，不加新区间。四态定义使用舍入前点值：x>0 且 y<0 为一致改善，x<0 且 y>0 为一致变差，两项同号且都非零为不一致，至少一项恰为 0 为无变化。Table S17 的四行与 Table S15 的 32 行计数一致。M2 和 M3 的无变化计数均为 0。

## 8. 复现、环境和数据可用性

主评价分数与统计的来源提交为 dc028b3350e804cbf8f40acea08261ed4f89c4fa，当前八对的主结果由 5d992b95 汇总。来源追溯见提交 196639e9，CAP5 支持分析见提交 7ef3b25f。本次整合读取上述已验收文件，更新文稿、补充表及数据包，不重新对接、拟合模型或计算 bootstrap。
原 V4.2 正式统计环境为 Python 3.11.15、NumPy 2.4.6、RDKit 2025.09.5 和 scikit-learn 1.9.0。Q110 配体准备记录为 RDKit 2025.09.6、Meeko 0.7.1、Vina 1.2.7，使用 ETKDGv3、randomSeed 20260727、单构象和最多 200 步 MMFF。Q110 统计程序在正式计算前通过 5649 项 PM48 回归；回归报告保留，旧 PM48 派生统计表通过 Git 历史查询。

Supplementary_Source_Data.zip 收录当前八对的主结果、人群、分折、来源对照、逐条活性证据，以及 CAP5 的保留名单、统计、折外预测和排名。原 805 条来源附件按来源记录保留，由当前 867 条总表映射到本次评价成员。包内 SHA256 清单只用于校验文件字节，包内不含 ChEMBL 数据库原库。来源与支持分析的逐表对应见 Source_Register.csv。仓库为 https://github.com/1280602962-debug/gwj260531，分支 cursor/v4-2-final-audited-results-20260930。

## 9. 严格阈值支持分析

严格 6.5/5.5 同类子集共 56 个分子：dual 25、A_only 12、B_only 14、neither 5。M0 的 summary_min 为 0.4571，M2 为 0.6286。完整结果见 Table S18。这是 Q110 内的支持分析，不替换主评价人群。

## 10. 选样证据及其覆盖范围

七对靶点的沉积名单、脚本记录及当前身份核对见 Table S22。它们均记为 DOCUMENTED_NOT_REPLAYED，不能写成原始抽样程序已完整独立重现。F2/F10、JAK1/JAK2、JAK1/TYK2、PPARG/PPARA 和 PPARA/PPARD 的构建记录使用严格 6.5/5.5 条件，正式标签按 θ=6.0。EGFR/HER2 的构建记录包含具名参照分子和负类优先选择。AChE/BChE 的旧脚本使用 ChEMBL ID 前缀作多样性键，不能称为 Bemis–Murcko 骨架限制。其一次只读重放与沉积名单有 1 个成员不同；当前 98 条合格观察仍可从沉积名单经别名和不合格记录处理闭合。PIK3CA/mTOR 采用冻结 Q110 名单。

CAP5 不改变这些选样记录，也不把八对原有的构建条件改写为同一种抽样规则。

## 11. CAP5 化学多样性敏感性分析

CAP5 是在正式评价集冻结后开展的事后敏感性分析。结构取自当前 mapping 的 canonical_smiles，使用 jcim_stage2_compute._murcko，即 MurckoScaffold.GetScaffoldForMol 后转换为 SMILES。在同一靶点对、同一实验类别内，同一非空骨架最多保留 5 个分子。无环分子不占非空骨架名额，同一骨架在不同类别中分别计数，不回补被排除的名额。

超额组按固定字段生成的 SHA-256 升序保留前 5 个，哈希相同再按 panel_id 排序；字段、分隔符和描述性带宽在敏感性统计前写入 CAP5_SELECTION_FREEZE.json。保留规则不读取对接分数。全部 867 条中保留 861 条、排除 6 条，4 条无环记录及 111 条 neither 均保留，见 Tables S23–S24。

方向性 AUROC、方法配对差、负类对照、口袋对照及联合排名使用原统计定义；化学模型在 CAP5 成员上重新拟合，排除成员既不训练也不测试，沿用原冻结 fold。保留成员中的原有缺 fold 情况仍按原规则处理，没有新的 fold 连接失败。完整统计和与主分析的对照见 Tables S25–S31。折外逐分子预测、全排名及分折检查随 CSV 数据包提供。

点变化的 0.02 带宽只用于描述，不是显著性或等效性界值。完整点估计、区间和人数同时保留。PPARA/PPARD 的 M3 在 dual 对 A_only 方向有区间边界变化：主分析为 0.4922–0.7676，CAP5 为 0.5030–0.7772。该边界变化单独列出，不以总体稳定概括掩盖。

CAP5 没有重评五随机种子和替代受体，也没有重新对接或重评分。CAP5 排名使用自身的共同成员和 k=ceil(0.10×N)，不能把其中 k=9 的名单称为原 Top-10。主分析继续采用 867 条合格观察，PIK3CA/mTOR 的 Q110 及其主结果不变。
