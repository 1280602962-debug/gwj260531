# 中英文 SI 压缩摘要（2026-09-13）

## 2026-09-14 ChEMBL 37 统一

正文、方法、讨论和正式 SI 不再区分 2026-07-23 REST 收获与 ChEMBL 37 转储。八对已评分配体的最大 pChEMBL 与本地 ChEMBL 37 一致（0 missing / 0 mismatch）。敏感性只保留同一套 ChEMBL 37 记录上的最大值/中位数重标（Table S3）。2026-08-26 API 快照和三对高置信字段筛查改为仓库归档，不排版。

---

# 中英文 SI 压缩摘要（2026-09-13）

正式投稿 SI 从 Tables S1–S14 / Figures S1–S5 压到 **Tables S1–S9 / Figures S1–S4**。从正式 SI 删除不等于从仓库删除。上一轮 PR37 的点估计方向性区间、类别分层 `summary_min` 敏感性 CSV、以及可执行分子身份规则均予保留；分层敏感性不再排成正式表。

## 判断

当前 SI 的主要问题不是“数字太多”，而是混入了历史对照、代码审计、post hoc 稳健性检查和重复证明。正式 SI 只保留四类问题：数据和 docking 怎么做；核心 fixed-channel 结果从哪里来；ligand chemistry / pocket correspondence 是否支持归因；主要结论对标签、样本和计算实现是否明显不稳定。

## 本轮结构

1. **S1** 计算设置与统计定义（压缩失败案例解释）。
2. **S2** 受体、对接盒与统一 CalcRMS（删除历史 Hungarian S2b）。
3. **S3** 标签/聚合敏感性（只展开组成变化的 EGFR/HER2、PIK3CA/mTOR）。
4. **S4** 固定评分通道 + 两个最大 Δ 的骨架/文献簇重采样。
5. **S5** 16 个 ECFP4 方向 + 每对最佳描述符（删除四描述符完整矩阵）。
6. **S6** 对应/非对应口袋 + 未使用池留出（删除单向拆分表；增加 weaker-arm-switched）。
7. **S7** 计算实现：受体替换 + 独立 GNINA + PPARG 重评分。
8. **S8** 独立性过滤后的外部准入（删除原始供给清点）。
9. **S9** 靶对审计（七个纳入对一句话概括；逐项列排除对和 EGFR/HER2 例外）。

## 从正式 SI 移出、留仓库

S2b Hungarian RMSD；S4 dual-versus-all-nonduals 小表；S6b 单向拆分；S9c–S9e 五种子/固定成员表；S10 分层 bootstrap 与样本量情景；Figure S5 可检测效应模拟；S11a 原始供给清点；整个 S12 发表年份表；S13 操作点表。

正文 Results 3.6 删除 2018 子集句，把 `not external validation` / `不作为外部验证` 留在 BindingDB 准入结论上。五种子结论改为一句正文 + 仓库 CSV。
