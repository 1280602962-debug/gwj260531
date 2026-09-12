# 审核分层（不要混用“已通过”）

三种检查回答的问题不同，不能互相替代。

| 层级 | 实际做了什么 | 能证明什么 | 不能证明什么 |
|---|---|---|---|
| 自动字段筛查 | 脚本按物种、测定类型等规则写 `human_include_exclude` 等列；审核者名可能是 `local_agent_metadata_pass` | 记录是否触发预设规则 | 不是逐篇阅读原文 |
| 数字一致性 | 稿面数字、锁定 CSV、图注是否对得上 | 表与表、表与稿是否一致 | 上游对接或 ChEMBL 原始记录是否真实 |
| 产物回放 | 从保存姿态读 `REMARK VINA RESULT`，或按化学对应重算 RMSD | 已存文件能否重建该条分数 / 该条 RMSD | 未归档的运行；也不是实验活性本身 |

缺失项不得标 PASS。旧审计报告里预填的 “yes / complete” 只说明当时脚本这样写，不作为原始真实性证明。

与稿件的对应：

- “高置信人源单蛋白复核”属于自动字段筛查，见 Methods 2.2 与 Table S3。
- 主 Vina 分数回放见 `local_track_b_v0/analysis/PRODUCTION_POSE_INVENTORY_V1.md`。
- 14 个主受体共晶 RMSD 统一表见 `primary_cognate_rmsd_calcrrms_v1.csv`。新增八受体化学对应复核见 `LAYER3_COGNATE_RMSD_REAUDIT_V1.md`。
