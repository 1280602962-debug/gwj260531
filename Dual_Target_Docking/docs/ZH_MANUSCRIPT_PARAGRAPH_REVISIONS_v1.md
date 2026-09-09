# 用户中文稿逐段修改说明（nature-skills）

对照用户粘贴的 Introduction–Conclusion。轴：`task=manuscript`，`paper_type=research`，`journal=generic`（JCIM）。技能：`nature-writing`（intro / method / experiments / discussion / conclusion）、`nature-shared`（main-text-discipline、nature-introduction、nature-results-discussion、discussion-argument-language）、`nature-polishing`（failure-modes）。

术语锁定：dual / A-only / B-only / neither；方向性判别；\(\mathrm{summary}_{\min}\)；对接评分；未使用池留出集；对应口袋 / 非对应口袋。文献编号改回仓库 `REFERENCES_JCIM.md`（POLYGON = 20，超大规模对接 = 19，Zhou = 9，Kinase-Bench = 22，Caba = 23）。

SI 编号一律用合并后的 S1–S13（旧 S25/S28/S29/S30/S32/S34/S39–S41/S46/S48–S49/S54 已废）。

---

## 1. 引言

**第 1 段。** 开头三句仍是教科书式“多靶点旨在…复杂疾病网络”，收敛太慢（nature-introduction：第一段末应露出具体瓶颈）。删去通路冗余/耐药铺陈；把 POLYGON 与超大规模对接改写成“前瞻发现能力”，并立刻点明这不是本文的评价问题。用户文献 4,5 应对 20 与 19，不可对调。用户“并不容易。3”在仓库参考文献中是 Kitchen 对接综述，已并入 [1,2]，避免错引。

**第 2 段。** 保留对接定义与 DUD / DUD-E / LIT-PCBA。用户 8–10 在仓库中是 5–7；CASF 未进入本段（不靠堆文献收窄问题）。用户 11–14 收为 [12,13,21]，对应“对照定义可改变方法排序”这一句。

**第 3 段。** 四状态与 Zhou / Kinase-Bench 的任务区分保留。用户 15→[9]，16→[22]。删“进一步表明”式的清单连接，改成 Kinase-Bench 目标与本文任务的对比。

**第 4 段。** 删除“首先 / 随后 / 进一步”四步路线图（intro 只预告问题与研究路线，不预告实验步骤）。末句改为：本文关注评价结论的来源，而不是新打分函数。冻结日期与 DualFourClass-Bench 不进引言。

## 2. 方法

**2.1。** 原段几乎重印引言。压缩为：四状态定义、口袋匹配评分、以及后文将做的三类检验。保留“dual 为正类”。

**2.2。** 数据处理保留。公式改为行内 LaTeX。高置信重标限定为“原冻结已打分配体”，普查后五对未覆盖；结果见 Table S3，不把 352/352 写成八对都做了。

**2.3。** 原稿把“A-only 和 B-only ≥50”写成八对共同准入，不成立。改为：厚面板操作性门槛；EGFR/HER2 未过门槛而保留；八个靶对按供给厚度与结构可行性分别设定候选池、配额和骨架上限（Table 1），全文不写两波查库。Table 1 数字与 9V8H 肽脚注保持冻结表。

**2.4.1–2.4.4。** 技术内容基本保留。用户“Tables S2 和 S3”改为对接盒与共晶 RMSD 同为 Table S2。

**2.4.5。** 独立 GNINA 仅三对，不比较引擎总体性能，指向 Table S9。

**2.5。** 公式排版为 \(S_{\mathrm{Vina}}\)、\(\mathrm{summary}_{\min}\)。负类交换 → Table S4；AND 过滤 → Table S13；文献簇 → Table S10。不把 bootstrap 算法再在结果里复述。

**2.6。** 化学基线 Table S5；对应口袋 Table S6；晶体替换 Table S8；留出集 Table S7；种子 Table S9；BindingDB / 时间切分 Table S11 / S12。原稿 assay-context、MCL1、PM110、J0 49 对等 provenance 不进正文方法。

**2.7。** 只保留 Table S1 与 Table S2。

## 3. 结果

**3.1。** 普查漏斗数字正确，保留。删与 2.3 重复的筛选过程句，只保留八对名单和 EGFR 供给受限。

**3.2。** 原稿以 \(\mathrm{summary}_{\min}\) 范围开场，把旗舰发现埋在后面。改为先写 EGFR 固定评分通道 0.430 对 0.808，再 JAK1/TYK2 重复，再写哪些靶对差值含 0。旧 Table S34 → **Table S4**；S25 / S46 → **Table S13**。\(\mathrm{summary}_{\min}\) 范围和 PPARG/PPARA 区间作为第二段必要支持，不作为段首。

**3.3。** 数字正确。旧 S28 / S19 → Table 2 / **Table S5**。补 16 个方向最大 |Δ| = 0.023 的出处为 Table S5。

**3.4。** 共晶 RMSD 旧 S3 → **Table S2**。写出 5DXT 0.505 与 4JSX 0.639，旧 S30 → **Table S8**。独立 GNINA 旧 S32 → **Table S9**；JAK1/TYK2 dual-versus-neither 按 SI 写作 0.705（该表无区间则不编造）。PPARG 重评分指向 Table S9，留出集 Table S7。删“表明结构特异性”这类跨段综合，只保留当前实验的局部推断。

**3.5。** 阈值网格旧用户 S4 与仓库旧 S4 冲突：阈值在合并后是 **Table S3**。中位数旧 S29 → **Table S3**。种子旧 S54 → **Table S9**。留出集旧用户 S13 与仓库 Table S13（EGFR 操作点）冲突：留出集是 **Table S7**。删除“JAK1/TYK2（0.475）略有上升”，改为由 0.365 升至 0.475。文献簇旧 S39–S40 → **Table S10**；功效旧 S31 → **Table S10**。高置信结果只报原冻结已打分配体。

**3.6。** 删除英文 leftover “pair”。BindingDB 旧 S48–S49 → **Table S11**；时间切分旧 S41 → **Table S12**。

## 4. 讨论

**4.1。** 任务区分保留。删除 EGFR AUROC 重报（Results 已建立）。文献定位改为 Zhou / Kinase-Bench / DUD–LIT-PCBA / Gu，不再复述 Figure 2。

**4.2。** 0.023 可保留一个锚点数字。Caba 等核对为 *J. Cheminform.* 2024, 16, 40，现为参考文献 [23]；用户原稿 [18] 在仓库中是 DOCKSTRING，已改正。

**4.3。** 激酶交叉对接 → [14]。明确独立 GNINA 检验的是评价设定能否换引擎再现，不是对应口袋优势。对应口袋结论保持：仅 EGFR/HER2 与 AChE/BChE 在主集排除 0，留出集全部含 0。

**4.4。** 前瞻 vs 回顾的边界保留。benchmark 改为“回顾性双靶对接评价”，避免把本文说成新方法基准。

**4.5。** 保留样本、异质性、外部不可得、Vina 为主。不写已删除靶对。删除“unused-pool holdout”英文夹杂。

## 5. 结论

保留两个评价任务不是难度差。把“对应口袋、受体替换、替代评分和独立对接均未显示结构对应优势”拆开：对应口袋无跨靶对稳定优势；受体替换改变个别数值；独立对接再现的是设定差异。不引入新数据、不写冻结年表。
