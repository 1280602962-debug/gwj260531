# URAT1 数据审计

| 项目 | 已核实内容 | 文件 |
| --- | --- | --- |
| 来源 | ChEMBL 活性导出，文件名 `cf12`，版本号无法恢复 | `data/raw/URAT1_CHEMBL_cf12.csv` |
| 靶点 | CHEMBL6120，UniProt Q96S37，SLC22A12 | 同上，1264 行全部为此靶点 |
| 活性类型 | 存档导出只有 IC50。配置里写过 Ki 和 EC50，这份表里没有 | 原始表 `Standard Type` |
| 关系符 | 等号 951，大于 182，小于 46。只保留等号 | `scripts/utils_ml.py` |
| 最终样本 | 822 个化合物，218 个骨架，pActivity 均值 6.231 | `data/processed/urat1_curated.csv` |
| 划分 | 五折 Murcko 骨架 | `data/splits/scaffold_fold_*.csv` |
| 是否用于缩库 | 否 | 临床库评分表无此列 |

lesinurad 与 benzbromarone 因多测定冲突没有进入 822。这解释了它们出现在回测而不在训练集。dotinurad 的预测失败因此不是“训练集记忆”。
