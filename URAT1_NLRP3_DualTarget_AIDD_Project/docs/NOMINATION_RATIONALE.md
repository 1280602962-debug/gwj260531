# PF-03882845后续验证优先级

身份：PF-03882845，REP_07580，CHEMBL1215331。成员来自[冻结primary](../data/frozen/primary_candidates_frozen.csv)，结构指标来自[分seed表](../data/manuscript/tables/candidate_seed_metrics.csv)。本文件集中说明历史选择，不是本次重新筛选或新的前瞻性预登记。

该分子在三个seed均通过A2双靶结构门及W2注释；NLRP3关键接触为7/7、7/7、7/7，IFP Jaccard均0.85，覆盖比例均约0.867。NLRP3模型分数为1.0，应解释为模型输出而非直接结合概率。URAT1 seed42酸基—Arg距离3.238 Å、IFP Jaccard 0.786、10/11关键接触；QED约0.736。历史优先规则偏重双臂一致性和完整接触，这支持把它优先安排跟进，不支持效力排名。

受体刚体转移分析中，PF在9DKB/9DKA/9DKC/9B1I的酸基邻近距离约3.24/3.23/3.20/3.18 Å。该分析没有重对接，也没有考虑完整松弛。在NLRP3的8ETR转移中最短关键距离约2.14 Å，存在近接触，不能称为跨构象无冲突结合。

已知药理来自原始论文，不是新靶点发现。Meyers 等（2010，DOI 10.1021/jm100505n）将 PF-3882845 报告为非甾体盐皮质激素受体拮抗剂（MR 结合 IC50 2.7 nM），并在 Dahl 盐敏感大鼠肾病模型中显示降压与肾脏保护。Orena 等（2013，PMC3796291）在醛固酮肾损伤模型中观察到肾脏 *Il-6* 等炎症基因表达下降。因此，即使将来 IL-1β 下降，优先解释仍是 MR 上游通路，而不是直接结合 NACHT。糖尿病肾病适应症或肾脏保护动物数据**不能**推断近端小管顶膜游离暴露。对 `PF-03882845 AND (URAT1 OR NLRP3 OR gout)` 的公开检索无命中，只说明文献检索未把该药连到这些词，**不是**已确认的新靶点注释。

HNW005（Sun 等，2025，DOI 10.1016/j.ejmech.2025.117644）与 Nature Communications compound 32（Zhang 等，2025，DOI 10.1038/s41467-025-62645-6；SI 化合物名为 1-((4-bromonaphthalen-1-yl)methyl)-N-(thiophen-2-ylsulfonyl)-1H-indole-2-carboxamide）已经是实验报道的 URAT1/NLRP3（或多转运体）抗痛风化学型。本工作不得写成首次双靶策略。结构来源见[文献核对](../data/manuscript/si/literature_primary_sources.md)。

优先级具有规则依赖：更偏URAT1接触会偏向Lanifibranor；不同通路混杂考量可改变后续选择。固定选择PF的目的在于使后续验证可追溯，不产生“所有指标最优”的结论。

下一阶段见[MD计划](../config/md_protocol.yaml)与[实验清单](FOLLOWUP_EXPERIMENTS.md)。MD需两个晶体配体对照，四个体系各3×200 ns；当前尚未授权、尚未执行，不得写入轨迹数字。MM/GBSA不能回头改变12个primary或PF选择。直接URAT1转运抑制、NLRP3结合/ATPase及MR干预条件下的通路实验，是不同问题的验证方向。
