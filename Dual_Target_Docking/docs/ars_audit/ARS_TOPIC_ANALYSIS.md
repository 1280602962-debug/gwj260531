# ARS（academic-research-skills）课题分析报告

> 依据 [Imbad0202/academic-research-skills](https://github.com/Imbad0202/academic-research-skills) v3.19.0  
> 技能：`deep-research`（RQ + Architect + DA + 3W）· `academic-paper-reviewer`（methodology / devil's advocate）· `academic-pipeline`（mode advisor · AI failure modes · claim verification 预检）  
> 输入：本仓库 `Dual_Target_Docking/docs/` 既有方案（非完稿论文）  
> 日期：2026-07-23 · PR 分支 `cursor/dual-target-docking-standalone-0b1a`

**与 ResearchStudio 审计的关系：** ResearchStudio 偏「idea 质量 / scoop / 可证伪」；本报告偏「RQ 精炼 / 方法蓝图 / 审稿人攻击面 / AI 科研失败模式 / 投稿路由」。二者互补，不互相替代。

---

## 0. Mode Advisor 路由（`mode_advisor.md`）

| 你现在的状态 | ARS 建议入口 | 本报告实际执行 |
|--------------|--------------|----------------|
| 有较清晰方向、要诊断缺口与可行性 | `deep-research` quick / socratic 片段 + architect | ✅ RQ Brief + Methodology Blueprint |
| 要对标核心论文 | `deep-research` three-way-scan | ✅ 3W 短单 |
| 方案当「准稿」做方法审查 | `academic-paper-reviewer` methodology-focus + DA | ✅ 方法卡 + 压力测试 |
| 尚无实验结果 → **不要**进 full paper / full pipeline Stage 2 | pipeline Stage 1 only | ✅ 停在 Scoping + 预诚信闸门 |
| 有实验后再开 | `academic-pipeline` Stage 2→2.5→3 | ⏳ 门控：打假表产出后 |

**常见误路由警告：** 「直接帮我写 NMI 全文」= 写作无证据 → ARS 明确反对。应先实验，再 `academic-paper plan/full`。

---

## 1. Research Question Brief（`research_question_agent`）

### Topic Area

Structure-based dual-target virtual screening：同一配体对靶点 A/B 的兼容性排序与评测协议。

### Candidate RQs（生成后择优）

| # | Candidate | FINER Avg | 未选原因 |
|---|-----------|-----------|----------|
| C1 | 朴素 mean/sum 融合是否系统性抬升 A-only/B-only？ | 4.2 | — **主 RQ** |
| C2 | softmin 校准融合是否优于 rank fusion？ | 3.6 | 易退化成「又一个加权」，缺诊断骨架 |
| C3 | fused vs linked 谁更难打分？ | 3.2 | 重要但宜作子问题，不宜主 RQ |
| C4 | 细胞活性能否验证双靶结合？ | 2.4 | **Feasible/Ethical-claim 弱**：标签语义不对，且你缺亲和力 |

### Primary Research Question（锁定建议）

> Under leakage-controlled paired activity labels, do standard dual-target docking score aggregations (mean / min / rank fusion) systematically elevate A-only and B-only molecules relative to true dual actives, and does a per-target calibrated shortfall fusion (e.g. softmin of \(\hat p_t-\theta_t\)) recover dual-vs-single ranking performance across multiple target pairs?

### FINER Assessment

| Criterion | Score | Justification |
|-----------|-------|---------------|
| **F**easible | **3** | 方法可及（GNINA/RTMScore/公开库），但 **真双靶配对样本稀疏** 是硬约束；私有细胞可做 holdout，不能当主证据。 |
| **I**nteresting | **5** | 与实务默认操作直接冲突；若打假成立，改变双靶 VS 默认协议。 |
| **N**ovel | **4** | 相对 Pérez-Castillo 的 rank fusion 与 VSDS-VD 单靶基准有清晰增量；非「全新算子」级 novelty。 |
| **E**thical | **5** | 计算+公开数据为主；私有分子去标识；不涉及高危双用途。 |
| **R**elevant | **5** | 直接服务双靶设计筛选；对接 NMI/JCIM 社区问题（VS 评测尺子）。 |

**Average: 4.4 / 5** · 无单项 <2 · **通过 FINER 门槛**（≥3.0）

### Scope Boundaries

**IN SCOPE**
- 同一配体、两端独立对接的普通双靶小分子（merged/fused/linked）
- 公开配对活性（ChEMBL/BindingDB）+ 文献 curated design_type
- 校准融合与 Dual-VSDS 协议；TrueNegative / RandomDecoy 两套 decoy
- 私有 NLRP3/JNK1 **时间盲测排序**（细胞 L2/L3）

**OUT OF SCOPE**
- 新对接采样器 / 通用亲和力 SOTA 竞赛
- PROTAC 三元复合物桥接
- 药物–药物协同（DeepSynergy 线）
- PK/成药性预测模型；动物药效主结论
- 把「未测」标为 inactive

**ASSUMPTIONS**
- 两端独立二元复合物近似成立（无强制共定位）
- 对接分含可校准的活性信息（否则所有融合头应一起失败 → 可证伪）
- design_type 可由文献可靠标注

### Sub-questions

1. **SQ1（诊断）：** 在 scaffold / leave-pair-out 下，mean/sum/rank 的 Top-k 中 A-only 污染率是否显著高于校准短板融合？  
2. **SQ2（协议）：** Dual-VSDS 在 TrueNegative vs RandomDecoy 上结论是否翻转（VSDS-VD 同构）？  
3. **SQ3（外推）：** 公开多靶点对上学到的协议，能否在 NLRP3/JNK1 细胞 holdout 上保持排序一致性（非结合金标准）？

---

## 2. Methodology Blueprint（`research_architect_agent`）

| 组件 | 选择 | 与 RQ 的逻辑链 |
|------|------|----------------|
| **Paradigm** | Pragmatist / positivist-leaning computational | 可测决策误差，不解释主观体验 |
| **RQ type** | Comparative + evaluative（尺子审计） | 「默认融合 vs 校准短板」对照 |
| **Method** | Controlled retrospective VS benchmarking + ablation | 对齐 `controlled_diagnostic_design` |
| **Data strategy** | D2-public（规模）+ D2-curated（真双靶）+ D5 holdout | 公开可复现 + 私有外推 |
| **Analytical framework** | dual-vs-single accuracy, PR-AUC, EF@k, ECE；两套 decoy | 指标直接服务 RQ，非只报 RMSD |
| **Validity** | scaffold / leave-series / leave-pair；校准集≠测试集；PoseBusters 门控 | 防泄漏与姿态垃圾分 |

**方法–问题一致性：** ✅ 高（问「聚合是否有害」，方法是「同姿态换聚合头」）。  
**一致性风险：** 若主文变成「我们的 GNN 更强」，则方法与 RQ 脱钩（ARS Mode 7 frame-lock）。

### Validity checklist（设计内建）

- [ ] 同一姿态特征，只换融合头（隔离机制）
- [ ] trivial / null 基线：random、单靶 max、raw mean
- [ ] 校准拟合不得触碰 test scaffolds
- [ ] 断言测试：未测不得进 inactive
- [ ] ≥2 公开靶点对 + leave-pair-out
- [ ] 负对照：\(\tau\to\infty\) / 退回 raw mean → 指标回落

---

## 3. Three-Way Scan（`three-way-scan` / ars-3w）

对标三篇「必比」论文（WHY / HOW / WHAT）：

| Paper | WHY（为何做） | HOW（怎么做） | WHAT（得到什么） |
|-------|---------------|---------------|------------------|
| **VSDS-VD** (NMI 2025) | RMSD 评测偏离真实 VS | True/Random/MassiveDecoy + 多引擎 | 精度≠物理合理≠富集；层级 VS |
| **Pérez-Castillo 2017** | 双靶富集需要融合多打分 | Rank 算术/几何平均 | 双靶对上融合优于单打分 |
| **EquiScore** (NMI 2024) | 打分泛化与泄漏 | 物理先验 + 增强 + 严划分 | 外测稳健、姿态来源不敏感 |
| **本课题（规划）** | 双靶默认聚合有系统偏置 | 校准 + softmin + Dual-VSDS + dual-vs-single | （待证）朴素融合抬升单靶偏倚分子 |

**Cross-paper synthesis**

- **Common WHY：** 现有尺子/协议误导决策。  
- **Divergent HOW：** VSDS-VD=单靶 decoy 分层；Pérez=双靶 rank 融合案例；你=配对标签 + 跨靶校准 + 硬负。  
- **Strongest WHAT（已发表）：** VSDS-VD 的「设定可翻转」。你必须产出同等级别的反直觉表。  
- **Unresolved gap：** 尚无开放、泄漏控制、多靶点对的 **dual-vs-single** 基准。

---

## 4. Devil's Advocate Stress Test（`devils_advocate_*`）

只挑战、不打分。按攻击强度排序：

### A1 — 「你们只是把 Pérez-Castillo 换了个壳」（High）

**反驳负担：** 必须在同一数据上显示 rank-mean **仍**抬升 A-only，而校准短板显著降低污染；并开放 Dual-VSDS。若效应量小 → 建议改投 JCIM，勿称 NMI。

### A2 — 「细胞好 = 双靶结合」是偷换概念（High）

当前方案文字上已降级，但 Abstract/Title 稍有不慎即被杀。  
**要求：** Abstract 主句不得出现 cell validation as binding proof。

### A3 — 「校准用了标签，当然打赢未校准」（High / Mode 4 shortcut）

若校准与评测共享信息，结果是 tautology。  
**要求：** scaffold-held-out 校准；报告校准 ECE；展示未校准 softmin 消融。

### A4 — 「Dual-VSDS 只是 VSDS-VD 换皮」（Med-High）

**要求：** 新任务对象（配对分子 × 双标签 × design_type），不是「再加一个靶跑一遍 VSDS」。

### A5 — 「配对数据太少，统计不稳定」（Med-High）

**要求：** 预先登记最小 dual 正例数；不足则主贡献改为协议/负结果论文。

### A6 — 「fused/linked 专家是后hoc 调参」（Med）

**要求：** 类型条件开关消融；预先定义类型标注协议。

### A7 — 「对接本身无信息，融合头都是噪声整形」（Med）

**可证伪友好：** 若所有头 dual-vs-single≈随机，应诚实报告并讨论对接失败模式（仍可能有诊断价值）。

---

## 5. Methodology Review Card（`methodology_reviewer_agent`，方案级）

> 注意：无实验结果，本卡评估的是 **设计是否足以回答 RQ**，不是结果质量。

| 维度 | 判断 | 证据 / 缺口 |
|------|------|-------------|
| RQ–方法匹配 | **Pass** | 同姿态换融合头直接服务诊断 RQ |
| 数据收集适当性 | **Warn** | 来源清晰，但稀疏性与 curated 规模未量化 |
| 分析指标正确性 | **Pass** | dual-vs-single / PR-AUC / 双 decoy 合适；避免只报 ROC |
| 可重复性 | **Warn→Fail if frozen late** | θ、引擎版本、box、种子、校准划分仍有 open decisions |
| 外部效度 | **Warn** | 依赖 ≥2 靶点对；单对 NLRP3–JNK1 不够 |
| 统计计划 | **Fail（现状）** | 缺预注册式：效应量、置信区间、多重比较、最小 n |

**方法审稿人一句话：** 设计方向正确，但目前仍是「协议草稿」而非「可重复实验方案」；补齐冻结参数与样本量门槛前，不建议进入写作 Stage。

---

## 6. AI Research Failure Modes 预检（`ai_research_failure_modes.md`）

现阶段 **尚无实验数字** → Mode 1/3/5/6 多为 **INSUFFICIENT EVIDENCE（预期）**。按 ARS：这些模式在真写论文进 Stage 2.5 时 **必须有 run log，否则阻断**。

| Mode | 现状 | 对本课题的具体含义 |
|------|------|-------------------|
| 1 Implementation bug | INSUFFICIENT | 对接/融合脚本需保存完整日志与种子；警惕「完美整数 EF」 |
| 2 Hallucinated citation | SUSPECTED 风险 | 文献表已较扎实；定稿前跑 S2/DOI 校验（ARS citation-check） |
| 3 Hallucinated results | INSUFFICIENT | **禁止**在无 CSV 前写「提升 X%」进草稿 |
| 4 Shortcut reliance | **SUSPECTED（设计层）** | 最大风险=校准泄漏 / 系列泄漏冒充机制；必须 scaffold 外校准 |
| 5 Bug as insight | Watch | 「RandomDecoy 上反转」若来自标签错误会很像洞见 |
| 6 Methodology fabrication | Watch | Methods 易写成「理想流程」；必须以 YAML/实际命令为准 |
| 7 Frame-lock | **SUSPECTED if NMI-only** | 过早锁死 NMI 叙事；若打假表失败应允许退回 JCIM / 负结果协议文 |

**Stage 2.5 预阻断条件（写给未来的自己）：**  
无打假表原始表、无校准划分证明、无引擎版本锁 → **BLOCK**。

---

## 7. Claim Registry 预登记（`claim_verification_protocol` Phase E 预备）

投稿前每条主张必须能挂证据。当前规划主张分级：

| ID | Claim（规划） | 允许强度 | 所需证据 | 风险 |
|----|---------------|----------|----------|------|
| K1 | 朴素 mean/sum 系统性抬升 A-only/B-only | **主 claim** | Dual-VSDS 多 split 表 | 效应量不足则降级 |
| K2 | 校准短板融合改善 dual-vs-single | 主 claim | 同姿态消融 + 负对照 | 校准泄漏 |
| K3 | Dual-VSDS 是首个…开放配对基准 | 仅 `SUPPORTED_WITHIN_SEARCH` | 写清检索库与日期 | ADV-E5 绝对「first」 |
| K4 | 细胞 holdout 验证了双靶结合 | **禁止** | — | 标签错误 |
| K5 | 方法预测 PK / 成药性 | **禁止** | — | 超范围 |
| K6 | 新 SOTA docking sampler | **禁止** | — | 与 RQ 冲突 |

**Novelty 写法（ARS E5）：**  
「Based on searches of PubMed/ChEMBL literature and docking VS benchmarks through YYYY-MM, we did not identify a leakage-controlled multi-pair dual-vs-single benchmark with architecture labels…」——禁止裸写 “first”。

---

## 8. 与 ResearchStudio 审计对照

| 维度 | ResearchStudio | 本 ARS 报告 | 合流建议 |
|------|-----------------|-------------|---------|
| Idea 强度 | 67/100 strong | FINER 4.4/5 | 想法够格推进 |
| Scoop | Level 3 | 3W + DA A1/A4 | Delta 必须可测 |
| 最大漏洞 | 数据稀疏、校准循环、缺 kill-switch | 同左 + 统计计划缺失 + NMI frame-lock | 先打假表门控 |
| 可完成性 | JCIM 高 / NMI 条件 | 同；并给出 Stage 路由 | 实验→再写作 |

---

## 9. ARS 推荐执行路径（人机协作，非全自动）

```text
[现在] deep-research scoping ✅（本报告）
    → 冻结 θ / 引擎 / ≥2 靶点对 / 校准划分
    → experiment-agent 或自管实验（ARS 不跑实验）
    → 产出：打假表 + 消融 CSV + YAML + logs
[门控] 若 K1 不成立 → 改 RQ 为「协议/负结果」或降期刊
    → academic-paper plan → full（IMRaD / NMI 结构）
    → academic-pipeline Stage 2.5 integrity（7 modes + citation）
    → academic-paper-reviewer full（含 DA）
    → revision loops → Stage 4.5 → finalize
```

**不要做：** 跳过实验直接 `academic-paper full` 生成「有数字」的 NMI 草稿（触发 Mode 3/6）。

---

## 10. 编辑部综合（模拟 EIC，方案阶段）

**Verdict: Major Revision to the research plan（非拒稿 idea）**

必须修改后才能进入写作：
1. 主 RQ 锁定为诊断式（本报告 §1），子问题服务主 RQ。  
2. 冻结可重复实验包（θ、版本、划分、最小 n）。  
3. 登记 claim K1–K3；删除 K4–K6。  
4. 预写 falsification / 负对照。  
5. 期刊策略：NMI = 条件目标；JCIM = 保底可完成出口。

**Strengths：** 问题真实、与高分评测文同构、自觉切割 sampler/协同/PK。  
**Weaknesses：** 证据未发生；参数未锁；统计计划空白；存在校准捷径与 frame-lock 风险。

---

## 11. 三行行动令（ARS 版）

1. **把主 RQ 改成 §1 的诊断句**，所有实验表格只为回答它。  
2. **先做 K1 打假表**；这是 Stage 2 写作许可证。  
3. **写作启动后** 强制跑 ARS Stage 2.5 七模式 + citation-check；细胞/PK 主张保持禁止列表。
