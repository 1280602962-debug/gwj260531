# Confirmed Contribution（V2）

> PaperSpine V4 Contribution-First gate · 2026-07-21 更新
> 本页取代早期 Glide XP / 8973 主验证叙事。

## Core Contribution

| Field | Content |
|---|---|
| Main contribution statement | We establish a **protocol-selected, asymmetric computational funnel** for URAT1–NLRP3 clinical-library repurposing. The workflow first evaluates open docking/rescoring combinations for URAT1 9DKB against paired TrueDecoy and RandomDecoy benchmarks, then couples the selected ranking protocol to NLRP3 model-based library reduction, dual-structure percentile ranking, Pareto non-domination, and medicinal-chemistry-aware nomination. |
| Contribution type | **new system / new analysis-or-benchmark**：靶点专属双诱饵协议评价 + 不对称双节点漏斗 + Pareto 与审计提名分离；不是新打分函数，也不是经实验验证的新药发现。 |
| One-sentence reviewer payoff | The paper shows how a gout dual-node screen can choose rather than assume its URAT1 docking protocol, assign different evidential roles to transporter and inflammasome data, and retain the reason why a mathematically strong hit is promoted or downgraded. |

## Why This Contribution Is Needed

| Field | Content |
|---|---|
| Field problem | 痛风治疗同时涉及 URAT1 介导的尿酸重吸收与 NLRP3/IL-1β 炎症，但多数计算筛选将两靶视为同质任务，并直接采用单一默认对接分数。 |
| Specific gap | 缺少一条同时满足以下条件的公开流程：用困难诱饵检验 URAT1 排序协议；显式处理 URAT1/NLRP3 数据不对称；将 Pareto 非支配与药物化学审计后的提名分开。 |
| Concrete challenge | URAT1 命名药物回收不足以支持 ML 主排；NLRP3 标签跨 assay 异质；随机诱饵可能夸大富集；两个靶点的原始分数不可直接比较。 |
| Why prior work leaves it unresolved | 湿法双靶研究解决的是合成与活性验证，不是临床库复用；PLK1/NLRP3 类不对称流程面对的是激酶而非膜转运体；单靶 NLRP3 ML/docking 研究不包含 URAT1 协议选择及 Pareto–审计分离。 |

## How This Paper Responds

| Field | Content |
|---|---|
| Design response | URAT1 TrueDecoy/RandomDecoy 协议评价 → NLRP3 分类分数缩库 → 9DKB/7ALV 双结构百分位 → Pareto 非支配 → PAINS/Brenk、类药性、证据来源和化学空间审计后的透明提名。 |
| Evidence required | （1）双诱饵上各协议的 EF\(_{1\%}\)、EF\(_{5\%}\) 与 AUC；（2）自对接的 top-1、集合最优和重打分选姿 RMSD；（3）NLRP3 骨架交叉验证与缩库；（4）生产协议下的双靶合并、对照药与 Pareto；（5）审计前后候选变化；（6）完整版本和命令。 |
| Evidence available | TrueDecoy/RandomDecoy 数据集及构建记录；NLRP3 模型与临床库分数；URAT1 ML 命名药物回顾；结构警报、类药性、化学空间、Pareto 敏感性和提名脚本。 |
| Evidence missing | 开源协议 P0–P5 的完整服务器结果与最终协议锁定；生产协议下的临床库双靶重算；RTMScore 运行归档；若保留构象讨论，则还需实际 MD 参数和结果。 |

## Claim Boundary

| Field | Content |
|---|---|
| Strong claims allowed after completion | 双诱饵基准可以区分 URAT1 排序协议；数据不对称支持不同证据角色；Pareto 前沿与审计提名是两个不同决策层；最终名单为可追溯的计算假说。 |
| Claims to soften or avoid | 双口袋直接抑制、药效协同、临床推荐、对接分数等同亲和力、富集证明结合位点、某候选已是 URAT1/NLRP3 双抑制剂。 |
| Novelty risk | “只是 ML→dock→Pareto”。回答必须依赖：**先选协议而非默认协议 + 双诱饵困难度对照 + 不对称证据角色 + Pareto≠提名**。 |
| Significance risk | 无湿实验。正文只能将产出定义为可证伪、供摄取实验与 MSU–IL-1β 实验检验的计算假说。 |
