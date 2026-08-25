# 投稿文稿指南（唯一写作入口）

> 2026-08-25 起：原目的“对接排出活性 → 双靶候选”已不成立。改构见 [`PROJECT_RECONCEPTION.md`](PROJECT_RECONCEPTION.md)。  
> **主投 *Journal of Computer-Aided Molecular Design***；不再按 *Molecular Diversity* 的 hit 文模板写。

## 一句话贡献

在预注册的 TrueDecoy/RandomDecoy 上，开源 URAT1 对接读出**不能作为活性检索器**（P2 按规则锁定，但 RandomDecoy 不优于随机；阳性集还排除了 lesinurad 等教科书药）。把该读出迁到临床库做双百分位 AND，得到的 7 个名字是**失败排序器的审计产出**，不是双节点候选。贡献是协议比较、标签病理与一次阴性迁移，不是命中发现。

## 写作文件

| 用途 | 文件 |
|------|------|
| 引言正文 | [`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md) |
| Methods 正文 | [`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md) |
| Results 正文 | [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md) |
| Discussion 正文 | [`DISCUSSION_DRAFT_CN.md`](DISCUSSION_DRAFT_CN.md) |
| 结论正文 | [`CONCLUSIONS_DRAFT_CN.md`](CONCLUSIONS_DRAFT_CN.md) |
| 协议锁定表 | [`PROTOCOL_SELECTION_RESULT.md`](PROTOCOL_SELECTION_RESULT.md) |
| 补充分析（assay 重叠 / EF 区间 / MCC950@7ALV） | [`SI_SUPPLEMENT_ANALYSES.md`](SI_SUPPLEMENT_ANALYSES.md) |
| 自对接烟雾（pose vs 排序） | [`REDOCK_SMOKE_ANALYSIS.md`](REDOCK_SMOKE_ANALYSIS.md) |
| 与 PLK1/NLRP3 差异 | [`DIFFERENTIATION_VS_PLK1_NLRP3.md`](DIFFERENTIATION_VS_PLK1_NLRP3.md) |
| 贡献 / 动机 / 期刊 | [`paper_spine_ars_analysis/`](paper_spine_ars_analysis/) |
| 漏斗复现命令 | [`LOCAL_AGENT_TASKS.md`](LOCAL_AGENT_TASKS.md)、[`WORKFLOW.md`](WORKFLOW.md) |
| 对照 Mol Divers 的修改计划（不重锁 P2） | [`MOL_DIVERS_REVISION_PLAN.md`](MOL_DIVERS_REVISION_PLAN.md)（期刊策略已改为 JCAMD 优先，此文件只作历史对照） |
| **课题改构（现行目的）** | [`PROJECT_RECONCEPTION.md`](PROJECT_RECONCEPTION.md) |
| 若仍要候选：新实验 vs 不能做的事后抬分 | [`IF_STILL_WANT_CANDIDATES.md`](IF_STILL_WANT_CANDIDATES.md) |
| MD 应跑哪些体系 | [`MD_RUN_PLAN.md`](MD_RUN_PLAN.md) |

## 正文结构（JCAMD：方法学主文；名单降为审计）

1. Introduction — 已有中文稿  
2. Methods — [`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)；协议筛选在临床库之前完成；生产对接为 P2（`num_modes=1`）；含 2.12 统计分析小节（配对 bootstrap、超几何检验）；MD 六体系清单（2.11，轨迹未报数值）  
3. Results — [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)  
   - R1 诱饵相似性泄漏审计：RandomDecoy 无近邻泄漏，弱活分子的重叠符合设计意图  
   - R2 双诱饵 P0–P5，锁定 P2（不选 P5）；配对 bootstrap 显示 P2 与 P5 在 TrueDecoy 上不可区分，P5 被否决靠 RandomDecoy 零命中；自对接说明 P2 Top-1 不是构象金标准  
   - R3 **阳性集构成**：469 均为 URAT1 IC50，但 lesinurad/benzbromarone/dotinurad 不在阳性里；类似物偏置  
   - R4 8319 → 1588 → 归档 1,580 行（有效双靶分数 **1,579**）；已知痛风药对照表（未进门控 = 阴性迁移）  
   - R5 裸 Pareto 4 vs 百分位门控 51 / 名单 7：**审计表**，不是候选短名单；Vecabrutinib 不再当 URAT1 跟进；GSK 仅作羧酸姿态假说（讨论）  
   - R6 姿态 QC（lesinurad 生产姿 Arg477 ≈ 14 Å）  
   - R7 MD：对照体系必须；GSK 仅当酸根假说讨论时保留；轨迹未报数值  
4. Discussion — 原目的不成立；标签失败与打分失败分开写；百分位 AND 不能提名双节点；酸根/Arg477 只作前瞻，不回写主表  
5. Conclusions — 对接不是 URAT1 活性检索器；7 个名字是审计产出；不鉴定双靶候选

## 主张边界（全文禁止）

- 发现 / 鉴定双靶抑制剂、双口袋、协同、1+1>2、临床推荐  
- **双节点计算候选、待双通路验证的跟进分子**（改构后 7 个名字只做审计）  
- 对接分或 MM-GBSA = \(K_i\) / 亲和力；P2 已能排序活性所以筛库成立  
- 富集 = 结合位点证明；P2 Top-1 = 近原生姿  
- 看过临床名单后重锁 Π\*，或把 lesinurad 加回阳性再重选协议  
- EGCG 或 canagliflozin 作为当前主推荐 lead  
- Glide XP / 默认 Vina 作为生产读出  
- 未实现模块：三态 \(S_{\mathrm{trap}}\)、生成式路径、SLC22 迁移主创新、Teacher 蒸馏

## 当前分子角色（不再当双节点候选）

| 角色 | 分子 | 说明 |
|------|------|------|
| 审计名单（漏斗产出） | 原 7 个百分位优选 | 证明弱排序器 + 化学规则会吐出什么；不是候选 |
| 讨论中的羧酸姿态假说 | GSK-3008348 | 仅因酸根化学，不因 S_U=97.5；不升格为双节点候选 |
| 不再跟进 URAT1 | Vecabrutinib | 入选只靠 CNN 百分位 |
| 对照 | lesinurad @ 9DKB；NP3-146/MCC950 @ 7ALV | 校准；lesinurad 还说明阳性标签未覆盖口袋配体 |
| 方法学负例 | EGCG、红霉素类大环 | 裸 Pareto 审计 |

审计漏斗门控（冻结，不据此提名候选）：\(S_U\ge 90\) 且 \(S_{N,\mathrm{dock}}\ge 90\)；Veber + Ro5 的 HBD/HBA/logP；MW 200–550；显式降级红霉素/epothilone 骨架。

## 数据三套（禁止混用）

| 集合 | 用途 |
|------|------|
| 临床库 8319 / 对接池 1588 / P2 完整案例 **1580** | 阴性迁移实验（百分位名单仅审计） |
| TrueDecoy / RandomDecoy | 只选 Π\*，不参与临床库排名 |
| 8973 distill | 仅 URAT1 回顾（可选 SI），不作双靶 Pareto |

## 已退役（不要再打开或引用）

已删除：旧全文稿 / 多份大纲、JMM 双稿策略、TAPE-GATE / MASFL / Teacher 蒸馏设计、历史商业对接（Glide XP）生产分数/短名单/图、三态对接计划。  
`data/repurposing/p2/` 为 gnina P2 生产漏斗归档（完整案例 1,580）。`data/repurposing/pareto/` 不再存放对接分数。

