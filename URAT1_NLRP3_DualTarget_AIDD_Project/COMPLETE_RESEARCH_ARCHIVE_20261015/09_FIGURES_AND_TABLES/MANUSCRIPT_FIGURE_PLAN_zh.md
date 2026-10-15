# 图的计划

不在本次生成新图。已有图在 `data/manuscript/figures/` 和 `figures/generated/`。

| 图 | 应表达的内容 | 已有文件 | 不要画进去的东西 |
| --- | --- | --- | --- |
| 1 | 8319→1588→303→156→54→37→12/21，旁边单列第三阶段 B 的 39→32 | `figure_1_screening_funnel` | 不要把 13 个新分子加进 12 |
| 2 | 六种对接读出的 EF@1% 及区间，标出被拒绝的方法 | 需由 `protocol_ef_ci.csv` 重绘 | 不要只画 P2 |
| 3 | lesinurad 与 NP3-146 的 top-1 RMSD 和 best-of-9 RMSD | 表已有，图需新绘 | 不要用最低 RMSD 代替选中姿势 |
| 4 | A1 与 A2 的四格表 | `urat1_gate_validation.csv` | 不要把 LR+≈1 画成成功富集 |
| 5 | 12、21 和 13 的化学空间，三组分开着色 | `figure_2_chemical_space` 只覆盖旧集合 | 不要合并成一个候选云 |
| 6 | 适用域：8319 的最近邻 Tanimoto 分布，标出 0.567 | 需由 `clinical_library_scores.csv` 重绘 | 不要省略域外这一事实 |

分子动力学图只有在轨迹入库并且对照规则被遵守之后才进入正文。现在不要画。
