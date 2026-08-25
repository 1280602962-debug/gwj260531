# 投稿文稿指南（唯一写作入口）

> 2026-08-18 起取代所有旧大纲 / 旧全文稿。  
> 先投 ***Molecular Diversity***；若拒，转投 ***Journal of Computer-Aided Molecular Design***（均为 Springer hybrid，可走非 OA）。

## 一句话贡献

预先锁定的 TrueDecoy/RandomDecoy 比较表明：**现有 URAT1 对接读出（含规则幸存者 P2）不足以在多样化临床库上做活性筛选**（RandomDecoy AUC 0.54，EF@1% 0.22，与随机无法区分）。据此，痛风双节点重定位**不能**从本漏斗的对接百分位鉴定候选。P2 百分位表与化学过滤后的 7 个名字只作启发式湿实验面板 / SI，不是 hit。主张入口：[`PROJECT_CLAIM_REFRAME.md`](PROJECT_CLAIM_REFRAME.md)。

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
| 主张重锁（P2 不能当 VS） | [`PROJECT_CLAIM_REFRAME.md`](PROJECT_CLAIM_REFRAME.md) |
| 对照 Mol Divers 的旧计划（历史备忘，不再指导主产品） | [`MOL_DIVERS_REVISION_PLAN.md`](MOL_DIVERS_REVISION_PLAN.md) |
| MD 应跑哪些体系（协议诊断，不是 lead 验证） | [`MD_RUN_PLAN.md`](MD_RUN_PLAN.md) |

## 正文结构（Mol Divers；转 JCAMD 时把 R1 提前加重）

1. Introduction — 已有中文稿  
2. Methods — [`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)；协议筛选在临床库之前完成；生产对接为 P2（`num_modes=1`）；含 2.12 统计分析小节（配对 bootstrap、超几何检验）；MD 六体系清单（2.11，轨迹未报数值）  
3. Results — [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)  
   - R1 诱饵相似性泄漏审计：RandomDecoy 无近邻泄漏，弱活分子的重叠符合设计意图  
   - R2 双诱饵 P0–P5：**无一达到多样化库可用的 VS 门槛**；P2 仅为“Random 非零”规则下的幸存者，RandomDecoy 上与随机无法区分  
   - R3 把该弱读出迁到临床库：已知 URAT1 药不回收；大环占据 Pareto——作为阴性应用结果，不是 hit 漏斗成功  
   - R4 化学过滤 P2 尾部得到的 7 个名字降为 SI 启发式面板，**正文不称为优选/双节点候选**  
   - R5 姿态 QC 与 lesinurad Arg477 ≈ 14 Å：P2 不是构象金标准  
   - R6 MD 六体系 = 协议诊断（晶体对照 vs P2 姿），不是 lead 验证；轨迹未报  
4. Discussion — [`DISCUSSION_DRAFT_CN.md`](DISCUSSION_DRAFT_CN.md)：P2 不能当活性筛选器；课题主产品改为阴性协议评价；湿实验纳入与百分位脱钩  
5. Conclusions — [`CONCLUSIONS_DRAFT_CN.md`](CONCLUSIONS_DRAFT_CN.md)：对接读出不够用；不从百分位鉴定双节点 hit

## 主张边界（全文禁止）

- 发现 / 鉴定双靶抑制剂、双口袋、协同、1+1>2、临床推荐  
- 对接分或 MM-GBSA = \(K_i\) / 亲和力  
- 富集 = 结合位点证明；P2 Top-1 = 近原生姿  
- EGCG 或 canagliflozin 作为当前主推荐 lead  
- Glide XP / 默认 Vina 作为生产读出  
- 未实现模块：三态 \(S_{\mathrm{trap}}\)、生成式路径、SLC22 迁移主创新、Teacher 蒸馏
- **P2 / 双百分位门控 = 筛选活性或双节点 hit 鉴定**
- 把 7 个化学过滤名字写成优选候选、lead、或待验证的双节点计算候选（主结论）

## 降级名单（SI 启发式面板，不是正文产品）

| 角色 | 分子 | 允许的说法 |
|------|------|------------|
| 启发式 URAT1 摄取面板（因羧酸，非因百分位） | GSK-3008348 | 可进湿实验；不声称 P2 已富集 |
| 启发式 NLRP3/IL-1 面板（因模型高分/炎症文献，非因双门控） | Vecabrutinib | 单端探索；不声称双节点已成立 |
| 结构对照 | lesinurad 晶体 @ 9DKB；NP3-146 @ 7ALV | MD / 姿态诊断 |
| 方法学负例 | 红霉素类 Pareto | 对接偏置 |

## 期刊

主产品改为阴性协议评价后，优先 ***JCAMD*** / ***JCIM***。*Molecular Diversity* 的 hit 模板与本结论冲突。

## 数据三套（禁止混用）

| 集合 | 用途 |
|------|------|
| 临床库 8319 / 对接池 1588 / P2 完整案例 **1580** | 主筛选 |
| TrueDecoy / RandomDecoy | 只选 Π\*，不参与临床库排名 |
| 8973 distill | 仅 URAT1 回顾（可选 SI），不作双靶 Pareto |

## 已退役（不要再打开或引用）

已删除：旧全文稿 / 多份大纲、JMM 双稿策略、TAPE-GATE / MASFL / Teacher 蒸馏设计、历史商业对接（Glide XP）生产分数/短名单/图、三态对接计划。  
`data/repurposing/p2/` 为 gnina P2 生产漏斗归档（完整案例 1,580）。`data/repurposing/pareto/` 不再存放对接分数。

