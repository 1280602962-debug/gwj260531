# 补充信息

本文件说明八对双靶点结果材料的成员、统计定义、扩展结果和文件位置。它不是审计日志。结果数值来自已验收的 PRIMARY 与第二阶段表。G、J、K 没有执行。

## 1. 标签、实体和来源范围

primary 阈值为 6.0。activity_eligible=0 不是 neither。AB_040 是 AB_046 的 alias，不作为独立观察进入 807 行人群表。pair 内观察可以跨靶点对共享 global ligand entity，所以 805 不是 805 个独立分子。

原 mapping 有标签、结构和身份字段，没有逐条连续 pA/pB、assay 端点或完整文献出处。ACTIVITY_LABEL_SEMANTICS_AUDIT.csv 只检查少量特殊行。因此，本文结果数值可以回溯到正式表，但不能写成 805 条标签都已逐条核对原始实验记录，也不能写成所有标签来自同一测定平台。已核对的来源和缺口列在 Source_Register.csv。

## 2. 成员和缺失

第 1 层 938 条标签方向记录供 D1 使用。第 2 层 934 条是正式 M0 方向记录，少掉的 4 条来自 AChE/BChE 的 M0 超时。第 3 层 928 条是化学共同 OOF，在 934 上再排除 6 条缺少冻结 fold 的方向记录。标签层缺 fold 为 10 条、7 个配体：AB_001、AB_053、AB_054、EH40_31、F2F10_080、J1TYK2_092、PGPA_030。M0 之后仍缺 fold 的是后四个配体的 6 条记录。AB_001、AB_053 和 AB_054 仍在 D1 中；它们离开第 2 层是因为超时，不是因为缺 fold。缺 fold 是没有冻结交叉验证分折，不是活性不可信。

fold_id 使用 0、1、2、3、4。0 是一个实际折号，不是缺失。AChE/BChE D_vs_A 的 A_only 由 27 到 25 发生在 M0 筛选，错口袋额外减员为 0。Table S2 给出每臂人数和名单。Table S13 中 expected = valid + missing。pair 展开记录数不等于物理 docking job 数。

## 3. PRIMARY 与配对统计

AUROC 按 higher-better 分数计算，ties=0.5。summary_min 是两个方向 AUROC 的较小值，不是对逐分子最低分再算 AUROC。D_vs_A 使用 score_B，D_vs_B 使用 score_A。已有区间来自 B=10000、seed=271828、每个分析单元重置的 default_rng，以及 linear 2.5 和 97.5 百分位。负类对照共享 dual 索引、两种负类分别抽取，差值在每次重采样内计算。口袋对照和替代受体对照使用同成员及共享索引。正文报告保留含 0 的区间，不把区间重叠写成等效，也不另作多重比较校正。

## 4. 化学模型与 D1

化学模型使用八个描述符、ECFP4，或 ECFP4 加相应方向的 M0。ECFP4 为 radius 2、2048 bit、useChirality=false、useFeatures=false。逻辑回归为 L2、C=1、liblinear、class_weight=None、max_iter=1000、random_state=271828。标准化只在训练折拟合；加入的 M0 列不使用测试折的均值和标准差。五个冻结折的折外预测拼接成 OOF，没有调参或重新分折。OOF 不是外部独立验证。正式输出不含 scaler_probe，训练折缩放的实现证据留在既有人工验收，不在本轮补造。

D1 在第 1 层、同一 canonical Mol 上记录八个描述符和 Murcko 字符串。空 Murcko 表示无环，不是解析失败。最大异类 ECFP4 Tanimoto 只在同一 pair×arm 内相对另一标签类计算。本轮不建立骨架等价规则，也不用新的骨架共享率声称没有泄漏。

## 5. 负类与口袋的完整结果

Table S5 和 Figure 2 给出 16 个臂的单靶负类、neither、差值和区间。Table S6 和 Figure 4 给出正确口袋、对侧口袋、差值和额外减员。反向点、零点和含 0 的区间都保留。

## 6. 五种子与替代受体

Table S10 为 40 个 pair×seed 原点，Table S11 为 24 个汇总。Figure S5 同时显示 auc_B、auc_A 和 summary_min。Figure 5a 只显示 summary_min。IQR 按各自 metric 的 linear P75−P25 读取。七对的五种子人数不变。AChE/BChE 仅 seed 17 因 AB_046 在 4BDS 超时而少 1 个 A_only。

Table S12 和 Figure S6 覆盖四对、八个替代条件。未替换端差值为 0 是构造核对。PIK3CA/mTOR 的 ALT_4JT5 提升了被替换的 B 端，但 summary_min 点差值仍为 0，因为未替换的 A 端仍是较小值。其 summary_min 区间不是零宽度。

## 7. 联合排序

方法自身描述和与 M0 的共同比较分开。共同比较中，新方法使用 pairwise_common_vs_M0，M0 使用 pairwise_common_vs_该方法，不用 method_specific 的 M0 行充当全部方法的基线。名次为 average midrank，联合名次为两端名次的 worst rank。并列先按 rank_A+rank_B，再按 global ID。k=ceil(0.10×N)。in_topk 取 YES 或 NO。Table S16 保留全部 9058 行，不在本 PDF 中逐行排印。M1、M2、M3 的方向差值读取 PRIMARY。M1b−M0 只在三组 dual、A_only、B_only 成员分别一致后，把已有的 M1b−M1 与 M1−M0 点估计相加，不加新区间。四态定义使用舍入前点值：x>0 且 y<0 为一致改善，x<0 且 y>0 为一致变差，两项同号且都非零为不一致，至少一项恰为 0 为无变化。Table S17 的四行与 Table S15 的 32 行计数一致。M2 的无变化计数为 0。

## 8. 复现、环境和数据可用性

科学实现与正式计算基准为 3f3e0d26052178c38ccb43dec637683a60f97a9d。已验收结果参考为 319375f10cb468d3a36e9cd0cf3e420bb5aa9632。正式统计环境为 Python 3.11.15、NumPy 2.4.6、RDKit 2025.09.5 和 scikit-learn 1.9.0。本轮绘图环境另行记录在整合记录中，二者不混用。版面几何和碰撞的专用检查工具本轮没有运行，该项为 NOT_VERIFIED；替代检查是正式图 PDF 页宽测量，以及实际查看最终图和两份审阅 PDF。该替代不能代替专用碰撞检查。运行日志中的 RDKit 与 scikit-learn 弃用警告被记录后继续执行；ConvergenceWarning 为 0。这些警告不是性能结论，本轮也不为消除警告改接口、solver 或 max_iter。

数据包 Supplementary_Source_Data.zip 按原字节收录结果 CSV、manifest、人群、alias、fold、四份权威 YAML 和三份指定评分主表，并附 SHA256 清单。SHA256 只核对文件字节，不证明标签已经逐条实验复核，也不表示稿件可以投稿。仓库为 https://github.com/1280602962-debug/gwj260531 ，分支 cursor/v4-2-final-audited-results-20260930。本材料没有 DOI。其他方法的评分主表由 FORMAL_AUTHORITY_PATHS.yaml 指向固定版本，不在包内另造统计来源。

## 可采用的 Methods 对齐

若下一轮写入全文 Methods，应保持下列已执行事实。M0 为 Vina seed 42、mode 1，higher-better 分数是亲和力的相反数。M1 为同一 mode 1 姿势的 GNINA CNNscore。M1b 为已保存姿势的最大 CNNscore。M2 为同一姿势的 RTMScore model1，列名 M2_RTMScore，不使用 max_saved_RTMScore。M3 为独立 GNINA docking 的正式 mode 1 CNNscore。AUROC、summary_min、bootstrap、化学模型、五种子 IQR 和联合排名的定义见上文第 3、4、7 节。没有证据支持把全部标签写成同一实验平台，也没有证据支持把冻结五折写成严格无骨架泄漏的外部测试。四态规则是计算前补记，不是历史预注册。
