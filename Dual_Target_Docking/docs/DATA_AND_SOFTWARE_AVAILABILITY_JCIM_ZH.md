# 数据与软件可用性

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源表，包括独立 GNINA 姿态生成分数（Table S32）、PIK3CA 占有率位移（Table S33）、文献阻断交叉验证（Tables S39–S40）、冻结的文献年份分割（Table S41）、assay-context 元数据审核（Table S42）、BindingDB REST 独立性计数（Table S43）、θ = 6.0 靶对普查（Table S44）、物化 caliper 匹配（Table S45）、AND 过滤工作点（Table S46）、配体层全图 AUROC（Table S47）、BindingDB 原生候选流程与切片摘要（Tables S48–S49）、冻结的 MCL1/Bcl-xL 面板与受体（Tables S50–S51）、文献对照（Table S52）、重建的 EGFR/HER2 共晶 QC（Table S3）、Figure S8、原生切片合约（`protocol/external_slice_contract.yaml`），以及评价合约（`DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`）。面向稿件的表 SHA-256 校验和见 `REVISION_CHECKSUM_MANIFEST_v1.csv`。ChEMBL 供给审计冻结于 2026-07-23；高置信 activity 视图抓取于 2026-08-26；BindingDB 原生归档锁定为 202608。GitHub Release 与 Zenodo DOI 将从打标签快照签发，而不是从当前仍可能变化的分支签发。分析环境与零新对接的复现命令见仓库 README。BindingDB TSV 归档本身不随仓库分发；CI 只核对已提交的 CSV。
