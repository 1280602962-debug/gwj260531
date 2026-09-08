# 当前冻结资产

- [primary_candidates_frozen.csv](primary_candidates_frozen.csv)：12个primary，唯一当前候选成员表。
- [backup_candidates_frozen.csv](backup_candidates_frozen.csv)：21个reserve/admissible，不能称为次优活性候选。
- [model_manifest.json](model_manifest.json)：已存模型与生产池的SHA-256及历史生成链限制。

这些是原C5冻结成员的丰富视图，由生成器检查集合一致性后输出；原C5文件只作不可改写审计源。原数据中的n_seed_pass是A2宽松双靶通过数，新增structural_seed_passes及W2分seed字段避免口径混淆。nlrp3_percentile不进入最终表；模型分数保留其缩库含义。

描述符重算字段带recomputed后缀，原QED/MW保留；InChIKey原值与重算值分别命名。临床阶段是沉积数据库的max_phase，不声称是整理日最新临床状态。变更名单需版本化决策与新溯源，不手工编辑CSV。
