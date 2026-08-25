# 投稿文稿指南（唯一写作入口）

> 2026-08-18 起取代所有旧大纲 / 旧全文稿。  
> 先投 ***Molecular Diversity***；若拒，转投 ***Journal of Computer-Aided Molecular Design***（均为 Springer hybrid，可走非 OA）。

## 一句话贡献

在 ChEMBL 证据不对称条件下，痛风双节点临床库重定位应：**先用 TrueDecoy/RandomDecoy 选定 URAT1 对接协议（Π\* = P2，gnina CNNaffinity）**，再用 NLRP3 分类缩库、双靶结构百分位排序，并把 **Pareto 非支配与药物化学提名分开**。产出可证伪假说，不是双靶抑制剂发现。

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
| 对照 Mol Divers 的修改计划（不重锁 P2） | [`MOL_DIVERS_REVISION_PLAN.md`](MOL_DIVERS_REVISION_PLAN.md) |
| MD 应跑哪些体系 | [`MD_RUN_PLAN.md`](MD_RUN_PLAN.md) |

## 正文结构（Mol Divers；转 JCAMD 时把 R1 提前加重）

1. Introduction — 已有中文稿  
2. Methods — [`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)；协议筛选在临床库之前完成；生产对接为 P2（`num_modes=1`）；含 2.12 统计分析小节（配对 bootstrap、超几何检验）；MD 六体系清单（2.11，轨迹未报数值）  
3. Results — [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)  
   - R1 诱饵相似性泄漏审计：RandomDecoy 无近邻泄漏，弱活分子的重叠符合设计意图  
   - R2 双诱饵 P0–P5，锁定 P2（不选 P5）；配对 bootstrap 显示 P2 与 P5 在 TrueDecoy 上不可区分，P5 被否决靠 RandomDecoy 零命中；自对接说明 P2 Top-1 不是构象金标准  
   - R3 8319 → 1588 → 归档 1,580 行（有效双靶分数 **1,579**）；已知痛风药对照表；NLRP3 阈值敏感性；Tegoprazan/Levotofisopam 未进门控  
   - R4 裸 Pareto 4（大环，审计）vs 双结构门控 51 / 优选 7；τ 敏感性不替换生产短名单；跟进 GSK-3008348 + Vecabrutinib 为待双通路验证的双节点计算候选  
   - R5 姿态 QC（7 个优选均在口袋内；lesinurad 生产姿 Arg477 ≈ 14 Å）  
   - R6 MD 六体系已指定、轨迹数值未报  
4. Discussion — [`DISCUSSION_DRAFT_CN.md`](DISCUSSION_DRAFT_CN.md)：假说边界；对照未回收的读法；P2 ≠ 构象金标准；MD 只做口袋压力测试；Unlike 湿法双靶与计算筛选模板  
5. Conclusions — [`CONCLUSIONS_DRAFT_CN.md`](CONCLUSIONS_DRAFT_CN.md)：7 个双结构门控计算候选；GSK-3008348 与 Vecabrutinib 为待双通路验证的双节点计算候选（不是单靶假说，也不是已验证双靶药）

## 主张边界（全文禁止）

- 发现 / 鉴定双靶抑制剂、双口袋、协同、1+1>2、临床推荐  
- 对接分或 MM-GBSA = \(K_i\) / 亲和力  
- 富集 = 结合位点证明；P2 Top-1 = 近原生姿  
- EGCG 或 canagliflozin 作为当前主推荐 lead  
- Glide XP / 默认 Vina 作为生产读出  
- 未实现模块：三态 \(S_{\mathrm{trap}}\)、生成式路径、SLC22 迁移主创新、Teacher 蒸馏

## 当前跟进分子（双节点计算候选，待双通路验证；非已验证 hit）

| 角色 | 分子 | 说明 |
|------|------|------|
| 双节点计算候选 | GSK-3008348 | 羧酸；双对接均衡；NLRP3 模型分弱；I 期 αvβ6 已停。两条实验都做，可先 URAT1 |
| 双节点计算候选 | Vecabrutinib | 对接+NLRP3 模型同向；BTK，II 期因疗效停。两条实验都做，可先 NLRP3 |
| 对照 | lesinurad @ 9DKB；MCC950 @ 7ALV（类似物对照，非自对接） | 校准，不是新提名 |
| 方法学负例 | EGCG、红霉素类大环 | 可进裸 Pareto，审计降级 |

Preferred 门控：\(S_U\ge 90\) 且 \(S_{N,\mathrm{dock}}\ge 90\)；Veber + Ro5 的 HBD/HBA/logP；MW 200–550；显式降级红霉素/epothilone 骨架。

## 数据三套（禁止混用）

| 集合 | 用途 |
|------|------|
| 临床库 8319 / 对接池 1588 / P2 完整案例 **1580** | 主筛选 |
| TrueDecoy / RandomDecoy | 只选 Π\*，不参与临床库排名 |
| 8973 distill | 仅 URAT1 回顾（可选 SI），不作双靶 Pareto |

## 已退役（不要再打开或引用）

已删除：旧全文稿 / 多份大纲、JMM 双稿策略、TAPE-GATE / MASFL / Teacher 蒸馏设计、历史商业对接（Glide XP）生产分数/短名单/图、三态对接计划。  
`data/repurposing/p2/` 为 gnina P2 生产漏斗归档（完整案例 1,580）。`data/repurposing/pareto/` 不再存放对接分数。

