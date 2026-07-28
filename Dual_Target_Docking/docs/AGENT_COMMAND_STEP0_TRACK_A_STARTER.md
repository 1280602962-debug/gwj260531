# Agent command — Step 0 / 第一步（Track A 启动包）

> **本命令是规划 v2 批准后的第一件可执行任务。**  
> 权威规划：[`EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md`](EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md)  
> 对应工作包：**W1 + W3**（零对接）。W2 prep 对照、Track B 大批对接 **不在本命令内**。

把下面整段复制给 agent 即可。

---

```text
【角色】
你是 DualFourClass 项目的文档/分析 agent（云端或本地均可）。
仓库根：含 Dual_Target_Docking/。
权威规划：Dual_Target_Docking/docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md。

【约束 — 强制】
1. 可大批对接，但本步 **禁止任何新对接 / 重对接 / 升 E / 换引擎跑分**。
2. **无湿实验**。
3. **禁止**重启 EGFR 同协议扩样（旧 AGENT_COMMAND_STAGE1 已关闭）。
4. **禁止**把 rtm_min_z 写成已验证通用主臂；禁止调 shortfall/clash/flags 进主分数。
5. **禁止**开启 Track B 第三对靶面板建设或大批对接（除非用户另发命令）。
6. 本步只做：Track A 启动包（W1）+ 预注册臂冻结确认（W3）。
7. 所有新写文字必须标明数据池：Exploration（已有结果）≠ Confirmation。
8. 做完后 git add / commit / push；用 ManagePullRequest 更新 PR（若环境要求）。
9. 结束时用中文给出 ≤10 行总结：完成了什么、文件路径、下一步建议（W2 或等用户批 Track B）。

【任务名称】
STEP0_TRACK_A_STARTER = W1（诊断主文素材包）+ W3（候选臂预注册冻结确认）

【背景（只读，勿改结论）】
- EGFR/HER2 panel40：rtm_min_z 点估计优于 vina_mean（Exploration）。
- EGFR/HER2 panel120 S1：**No-Go**；全面板 Δ≈−0.039；旧40仍正、新70显著负。
- PIK3CA/mTOR panel48：异质失败模式；rtm 提升亦不作确认。
- 方向仍是架构无关双靶四类决策；失败的是「单 pair 调参→同 pair 扩样确认」。

════════════════════════════════════
PART A — W3：预注册臂冻结确认（先做，短）
════════════════════════════════════

【输入】
- Dual_Target_Docking/data/protocols/CANDIDATE_ARMS_V0.yaml（应已存在）
- Dual_Target_Docking/docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md §5

【动作】
1. 打开并核对 YAML 是否含：vina_mean/min、rtm_mean/min/(min_z)、gnina_*、rank_consensus_*、
   forbidden 列表、selection_protocol=leave_one_pair_out、selection_objective=median_delta_auroc_vs_vina_mean。
2. 若缺失关键字段：补全 YAML（仍标 DRAFT / B0），不要新增「看完分数才发明」的臂。
3. 写一页冻结记录：
   Dual_Target_Docking/data/protocols/STEP0_ARM_PREREG_FREEZE.md
   必须包含：
   - 冻结日期与 git commit 意图说明
   - 主基线：vina_mean
   - 选规则：LOTO + median(ΔAUROC)
   - Confirmation：holdout 只评一次
   - 明确：现有 panel40/48/120 分数 **不得**用于改候选臂清单
   - 状态：B0 预注册完成；Track B 未开工

【完成标准 A】
- [ ] CANDIDATE_ARMS_V0.yaml 可读且自洽
- [ ] STEP0_ARM_PREREG_FREEZE.md 已写
- [ ] 未根据旧面板 AUROC 增删主竞赛臂

════════════════════════════════════
PART B — W1：Track A 诊断主文素材包（核心）
════════════════════════════════════

【输入 — 只读引用，路径以仓库实况为准】
- data/egfr_her2_panel120_v0/analysis/STAGE1_VERDICT.md
- data/egfr_her2_panel120_v0/analysis/STAGE1_EXPAND_ANALYSIS.md
- data/egfr_her2_panel120_v0/tables/stage1_bootstrap_delta.csv
- data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv
- data/egfr_her2_panel120_v0/tables/panel_v0_120.csv
- data/egfr_her2_panel40_v0/ （及既有 bootstrap / cognate QC 文档）
- data/pik3ca_mtor_panel48_v0/analysis/failure_typology_v0/FAILURE_TYPOLOGY_V0.md
- data/pik3ca_mtor_panel48_v0/analysis/bootstrap_ci_v0/（若有）
- docs/JCIM_GAP_TWO_SCENARIOS.md（方案甲 ≈ Track A）
- docs/RESEARCH_DIRECTION_REFREEZE.md（方向；执行以 v2 为准）

【产出目录】
创建（若不存在）：
Dual_Target_Docking/data/track_a_starter_v0/
  README.md
  FIGURE_PLAN.md
  CLAIM_CEILING.md
  PAPER_OUTLINE_TRACK_A.md
  tables/   （可从上游复制或 symlink 说明，勿静默改原数）
  scripts/  （仅当需要重画/汇总时；禁止重跑对接）

【B1 — CLAIM_CEILING.md】
写死 Track A 允许/禁止的 claim：

允许：
- C1：朴素 vina_mean 在四类任务上失败/硬负污染（EGFR panel40/120 证据）
- C2：姿态 Top-1 不可靠、重打分可改变姿态排序（cognate + mode 统计，Exploration）
- 同系列扩样失败：panel40 子集 vs 新70 分裂（S1）
- 失败分型 T1/T2/T5 作为诊断框架
- 跨 pair 异质性（EGFR vs PIK3CA/mTOR）为探索性对照

禁止：
- 「rtm_min_z 是经验证的通用双靶决策臂」
- 「协议已在独立数据上确认」
- 任何湿实验/乘客/moiety 主 claim
- 把 shortfall/clash 阴性结果改写成阳性方法

【B2 — FIGURE_PLAN.md】
为 Track A 规划 3–5 张主图/表（每张写：目的、数据文件、横纵轴/分组、一句 caption 草稿）。
至少覆盖：
  Fig1  任务定义：Dual vs A-only/B-only（示意或计数表）
  Fig2  C1：vina_mean AUROC / Top10 硬负（panel40 与/或 panel120）
  Fig3  S1 分裂：旧40 vs 新70 的 ΔAUROC（点+CI）——这是 Track A 核心新图
  Fig4  失败分型：EH40_23（T2）↔ PIK3CA/mTOR T5 对照（引用 typology）
  Fig5（可选） 两对靶森林图（Exploration only，图注写明非确认）

若缺汇总表：用脚本从现有 CSV **只做聚合**，写入
  data/track_a_starter_v0/tables/fig3_old40_new70_delta.csv
（数字必须与 STAGE1_EXPAND_ANALYSIS 一致；不一致则停下来查，禁止「调到好看」。）

【B3 — PAPER_OUTLINE_TRACK_A.md】
按期刊短文/诊断文结构写提纲（中英标题均可，正文提纲用中文亦可）：
  1 Title / Abstract 要点（诚实阴性 + 任务缺口）
  2 Intro：文献缺口（四类硬负决策，非又一个 consensus VS）
  3 Methods：面板、对接协议、指标、bootstrap；**写明探索性**
  4 Results：C1 → C2 线索 → S1 外推失败 → 分型
  5 Discussion：为何单 pair 试点不能当确认；v2 方法论含义（简述）
  6 Limitations：无湿实验、prep 混杂、N、架构 unknown
  7 Data availability：路径清单
每节列「已有文件指针」+「仍缺什么（本步不填坑的标 TODO）」。

【B4 — README.md】
- 本包目的、与 v2 Track A 关系
- 完成清单勾选
- 明确下一步选项：
  (i) 用户批 W2 prep 对照命令
  (ii) 用户批 Track B（另发 AGENT_COMMAND）
  (iii) 直接进入 Track A 写作（人工）

【完成标准 B】
- [ ] 上述 4 个 md 齐全且路径可点
- [ ] Fig3 所需数字有表或明确引用 STAGE1 分析，无新对接
- [ ] CLAIM_CEILING 与 v2 一致
- [ ] 未修改上游 docking 原始分数含义

════════════════════════════════════
PART C — 文档挂钩（轻量）
════════════════════════════════════

1. 在 docs/PROJECT_MASTER_PLAN.md 顶部「现行命令」处增加本文件链接
   （STEP0），并注明 STAGE1 命令已关闭。
2. 在 docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md §11「立即下一步」
   将 W1/W3 标为进行中或完成（按实况）。
3. 不要大改 v2 科学内容；只更新状态勾选。

════════════════════════════════════
PART D — 交付与退出
════════════════════════════════════

【交付物清单】
1. data/protocols/STEP0_ARM_PREREG_FREEZE.md
2. data/protocols/CANDIDATE_ARMS_V0.yaml（核对/必要时小补）
3. data/track_a_starter_v0/{README,FIGURE_PLAN,CLAIM_CEILING,PAPER_OUTLINE_TRACK_A}.md
4. （可选）data/track_a_starter_v0/tables/fig3_*.csv
5. 相关 docs 状态挂钩更新
6. git commit message 示例：
   "Add Track A starter pack and freeze candidate-arm pre-registration (Step 0)."

【退出检查 — 全部为是才能宣布 STEP0 完成】
- [ ] 零新对接
- [ ] 未宣称 C3/C4 成立
- [ ] 未启动第三对靶对接
- [ ] PR/分支已更新（若适用）
- [ ] 中文总结已给出下一步二选一：W2 vs 等 Track B 批准

【明确不在本步】
- LigPrep vs RDKit 20 分子对照（那是 W2，另开命令）
- GNINA / 新面板 / LOTO 实跑（Track B）
- 论文 LaTeX/投稿账号操作
- 重算以「让 RTM 变好」为目的的任何分数
```

---

## 给用户的使用说明

| 项 | 内容 |
|----|------|
| **谁跑** | 云端文档 agent 或本地分析 agent 均可 |
| **算力** | 无对接；几分钟到一小时级文档/表聚合 |
| **成功后** | 你手里有 Track A 可写大纲 + 算臂已预注册 |
| **下一命令** | 二选一：`AGENT_COMMAND_W2_PREP_CONFOUND.md`（待写）或 `AGENT_COMMAND_TRACK_B_B0B1.md`（待你批准后写） |

## 本步成功长什么样

打开 `data/track_a_starter_v0/PAPER_OUTLINE_TRACK_A.md` 能按图写作；打开 `data/protocols/STEP0_ARM_PREREG_FREEZE.md` 能看到「选臂规则已冻、旧分数不能改臂」。
