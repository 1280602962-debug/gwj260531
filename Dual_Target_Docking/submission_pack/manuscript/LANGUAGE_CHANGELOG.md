# 中英文润色修改摘要（2026-09-12，叙事回收轮）

上一轮并行改中英文时，把项目时间线和内部审计语又写回了正文。本轮按数据源重写，并同时改中英文，不再把英文当作单独“再翻译”一遍。

## 必须改的六项

1. **删除 original / added / early candidate 时间线。** 聚合敏感性按 ChEMBL API 快照与 ChEMBL 37 转储分源描述，并点名靶对。主文 Methods 删除 PIK3CA/PIK3CB（2WXF）撤出史。
2. **Figure 1 caption。** (C) 只写普查是供给可用性；八对按 Table 1 纳入。不再写 after which / 86→8。
3. **matched/mismatched。** Methods 写明交换已计算评分通道、不重新对接；删掉 `isolated effect` / “单独影响”。
4. **Redocking RMSD。** 本轮未统一重算 14 个受体。正文和 SI 改为定性搜索覆盖检查，不把 Figure S4 写成统一定量基准。
5. **Data Availability。** 压成正式的两句：仓库位置、分数表可复现、不要求生产姿态档案。内部 checksum / Table S54 / DOI 未发 全部移出主文。
6. **Discussion 4.3 / 4.5。** 4.3 只留对应口袋主结论；Figure S5 覆盖率细节压成一句并指向 SI。

## 其他

- 标题改为 `Evaluating Dual-Target Molecular Docking: Experimental-State Comparisons across Multiple Target Pairs and Sources of Discrimination`。
- 摘要：ECFP4 为 comparable or higher；0.023 限定在骨架分组逻辑回归设定下。
- 引言第三问改为 observed docking discrimination，不再写 remaining discrimination。
- 2.3 准入改为可核对条件：人源结构、非共价位点、共晶配体可定义对接盒、双向供给足够。
- 4.2 补回 dual 配体连接/融合设计两句，并立即限定未作架构标注。
- 参考文献 (16) 后补上 (17) DOCKSTRING，其后顺延至 (22)。
- 数字未改。

## 投稿面残留（2026-09-12 续）

上一轮已改主文 15 项。本轮只清投稿还会看到的残留：

- Figure 1C 图面：普查与八对评价之间加分隔线和 “Independent of the census”，去掉紧贴 86 的漏斗观感；图注写明八对不是从 86 对继续筛出。
- Figure S4 图注：删除 “eight added / original six”。
- Table S2c 脚注：不再写 “all eight”。
- 作图脚本与 figure audit 同步禁止这些时间线措辞。
