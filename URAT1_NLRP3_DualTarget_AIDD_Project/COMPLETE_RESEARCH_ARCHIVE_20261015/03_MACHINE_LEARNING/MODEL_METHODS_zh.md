# 模型方法与已核实性能

数字抄自 `data/models/training_report.json` 和 `results/nlrp3_model_v2/performance_summary.csv`。两种评价不能写成同一个模型。

## URAT1 回归

原始表 `data/raw/URAT1_CHEMBL_cf12.csv` 有 1264 行，靶点 CHEMBL6120，UniProt Q96S37，标准类型全部是 IC50。ChEMBL 版本号无法从现有文件恢复。清洗保留等号关系、可解析 SMILES 和 pActivity 在 4 到 10 之间的记录；同一 SMILES 的标准差大于 0.5 或极差大于 1.0 时整分子丢弃，否则取中位数。得到 `data/processed/urat1_curated.csv` 的 822 个化合物。

特征是 Morgan 半径 2、2048 位，加 12 个 RDKit 描述符。模型是 XGBoost 回归，`n_estimators=400`、`max_depth=5`、`learning_rate=0.05`、`subsample=0.8`、`colsample_bytree=0.8`、`reg_lambda=1.0`、`random_state=42`。划分是五折 Murcko 骨架 GroupKFold。

合并折外结果：RMSE 0.663，MAE 0.498，R² 0.508，Spearman 0.726。按 p≥6 的 EF@10% 为 1.71，按 p≥7 的 ROC-AUC 为 0.852。Y 扰乱后 Spearman 接近 0。这些数字说明模型学到了训练集里的趋势，不说明它能给临床库排序。

回测：lesinurad 预测 5.64，未过 6；benzbromarone 6.56，过线；verinurad 在训练集中，预测 7.08；dotinurad 5.07，未过。probenecid 没有存档预测。

临床库没有这列预测。`scripts/screen_repurposing_library.py` 只有在显式加上 `--with-urat1-ml` 时才计算它，正式导出没有这个开关的产物。

## NLRP3 v1，正式缩库模型

原始表 `data/raw/NLRP3_CHEMBL_4807.csv` 有 1334 行，靶点 CHEMBL1741208，标准类型全部是 IC50。清洗后 `nlrp3_records.csv` 为 609 条记录、513 个分子、39 个测定。标签是记录 pActivity ≥ 6.0，分子标签取最大值。分子水平活性比例约 0.604，所以 EF@10% 的理论上限大约是 1.65。报告的 1.57 已经接近这个上限，不能单独当成强富集。

特征是 Morgan 2048 加 12 个描述符，再加测定频次最高的 25 个测定的 one-hot。五折按分子骨架划分，同一分子的记录不拆开。`training_report.json` 里的 AUROC 0.8934017784193573、AUPRC 0.9143469018214805、EF@10% 1.5706359287488578 是五个折的指标的平均值。灵敏度、特异度和精确率没有写进这份 v1 报告。

1588 个成员由这个模型在阈值 0.5 上产生，成员表已存档。同一规则的补充材料重训得到 1377，原因没有查实。因此 1588 是冻结的成员表，不是现在可以逐分子重现的概率。

适用域按训练分子留一最近邻 Morgan Tanimoto 的第 5 百分位，阈值约 0.567。8319 个临床化合物里只有 5 个在域内，1588 里 3 个，12 个主候选里 0 个，32 个新分子里 0 个。概率是模型输出，不是 NLRP3 抑制的实验概率，而且几乎全部是外推。

## NLRP3 v2，没有替换缩库

v2 是外层五折骨架、内层校正的嵌套验证。结构特征、不做概率校正的合并外层测试 AUROC 是 0.8896，AUPRC 0.9131，召回 0.8742，特异度 0.8128，Brier 0.1261。测定条件变体的 AUROC 是 0.8807。泄漏检查记录拟合、校正和测试之间的 SMILES 与骨架重叠为 0。

用全部 513 个分子拟合、不做校正的模型重打了临床库，p ≥ 0.5 得到 1525 个，与 1588 的交集 1249、新进入 276、离开 339，Jaccard 0.670。lesinurad 的历史概率 0.597 在 v2 上是 0.407，离开新池；dotinurad 从 0 变为 0.537，进入新池。这些是模型差异，不是活性变化。v2 没有改写 1588，也没有改写 12/21。
