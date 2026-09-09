# 数据与软件可用性

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。

排入 Supporting Information 的核实表为 **Tables S1–S13**（软件与种子、对接盒子与共晶 RMSD、阈值与 pChEMBL 聚合、固定评分负类、ECFP4 增量、对应/非对应口袋、未使用池留出集、PIK3CA/mTOR 晶体替换、独立 GNINA 与五种子 Vina、文献簇与可检测效应、BindingDB 外部准入、2018 年切分计数、EGFR/HER2 操作点）。历史工作表 S1–S54（含原 Table S54 多种子表）的合并对照见 `data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv`。

更长的逐配体分数、holdout 成员、多种子长表、物化 caliper、chemotype 硬负、聚合均值、完整病例与 assay-context 底表、J0 候选对普查、BindingDB REST 历史计数、leave-cognate-out、PIK3CA 占有率快照、接触计数、全链序列一致性，以及 MCL1/Bcl-xL 适用性压力测试（不进入 Table 2），见 **Note S14**，随代码与 SHA-256 清单归档至 GitHub Release；Zenodo DOI 将从打标签快照签发，而不是从当前仍可能变化的分支签发。

`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源 CSV。面向稿件的表 SHA-256 校验和见 `REVISION_CHECKSUM_MANIFEST_v1.csv`。原生切片合约为 `protocol/external_slice_contract.yaml`；评价合约为 `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`。ChEMBL 供给审计冻结于 2026-07-23；高置信 activity 视图抓取于 2026-08-26；BindingDB 原生归档锁定为 202609。分析环境与零新对接的复现命令见仓库 README。BindingDB TSV 归档本身不随仓库分发；CI 只核对已提交的 CSV。
