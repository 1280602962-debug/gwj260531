# 高分文重定思路：skills 复审报告（2026-07-23）

> **调用的 skills（本机复跑）**  
> - ResearchStudio：`idea-quality` · `scoop_check` · `idea_spark` 模式  
>   `assumption_audit_and_pivot` · `controlled_diagnostic_design`  
> - ARS：`research_question_agent`（FINER）· `mode_advisor`（仍处 Scoping，禁止直接写全文）  
> Skills 仓库：`/tmp/skills/ResearchStudio` · `/tmp/skills/academic-research-skills`  
> **输入：** 用户否决「拼分主叙事」后，仍希望冲击高分方法文；需**明确双靶要解决的真问题**。

---

## 0. 一句话判决

**旧主线（校准拼分）应降级或放弃作为冲高分创新。**  
若仍要冲高水平方法文，应把要解决的问题重定为：

> **独立两端对接把「一个分子」当成「两个无关单靶配体」——这一假设在连接型/杂合型双靶上系统性失效；需要的是按设计类型隔离该失效、并给出可复现的结构一致性格验协议（而不是再发明一个分数加权）。**

拼分最多当附录对照；**主问题是对接评价对象错了，不是加权公式差了。**

---

## 1. ResearchStudio：先审计旧思路的「承重假设」

模式：`assumption_audit_and_pivot`

| 项 | 内容 |
|----|------|
| **旧思路承重假设** | 两端独立对接分数里**已经有足够信息**；失败主要来自跨靶分数不好比 / 聚合方式。 |
| **该假设在真实场景是否成立** | **高度可疑。** 用户指出且成立：对接错了，拼分救不了。连接型分子还有构象耦合与连接约束，独立对接根本未建模。 |
| **旧方法失败是否可归因于该假设** | 是。Pérez-Castillo / Jaiteh 已做「两端对接 + 融合/取高排」；再优化聚合 = 在同一假设上打补丁 → idea-quality 的 **深度弱、与问题错配**。 |
| **枢轴（pivot）** | 把承重假设从「聚合」改到「**评价对象**」：双靶对接应检验的是「同一化学图在两个口袋中的**可兼容构象对**」，不是两个独立分数的算术。 |

模式：`controlled_diagnostic_design`（对高分诊断文）

| 要隔离的混淆 | 做法 |
|--------------|------|
| 化学型混淆 | **同一协议下**分 fused / merged / linked 报告，禁止混报一个总分 |
| 尺子混淆 | 同时报：姿态（有共晶时）、物理有效性、筛选富集、**只强一端硬负**；禁止只报平均分 |
| 无效基线 | 独立两端对接 + mean/min/rank 必须作为 null；随机分必须垫底 |
| 外部锚 | 双端共晶（少但硬）+ 文献 curated 有类型分子 +（可选）实验室系列 |

---

## 2. 到底要解决双靶分子的什么问题（问题清单 → 选题）

结合近年综述（Proschak 2024 JMC；Drug Discovery Today 2024 in silico MTDL；Frontiers 2025 STaMPs）与计算前沿（POLYGON / LaMGen / FuseDiff），双靶问题可拆成：

| ID | 问题 | 高分潜力 | Scoop / 拥挤度 | 与你们资源匹配 |
|----|------|----------|----------------|----------------|
| **Q1** | 靶组合是否值得做 | 中（偏系统药理） | 中 | 弱 |
| **Q2** | 融合/连接化学与成药性 | 中高（药化主刊） | 高（大量 JMC） | **强**（合成侧） |
| **Q3** | 独立对接对「真双靶」假阳性高 | 中 | 已有 Zhou 2013 警示 | 中 |
| **Q4** | **设计类型如何使独立双对接失效** | **高（诊断+协议，贴 NMI 尺子文）** | **相对空**（未见按 Morphy 类型的双对接失效审计） | **强**（有文献池+少量双端共晶+实验室 fused/linked） |
| **Q5** | 联合生成双姿态/双靶分子 | 高 | **极挤**（POLYGON、CLM、LaMGen、FuseDiff、MTD2025） | 弱（缺大规模生成+湿实验流水线） |
| **Q6** | 分数怎么拼 | 低–中 | 已有融合案例 | 基建级 |
| **Q7** | 实验验证贵、假双靶浪费合成 | 中（转化叙事） | — | **强** |

**冲高分且可辩护的缺口：Q4（可带 Q3/Q7 作后果）。**  
**不要选 Q5 当主创新（易被 scoop）；不要选 Q6 当主创新（不新、且救不了对接）。**

---

## 3. ARS：FINER 研究问题 Brief（重定后）

### Topic Area

结构基双靶评价：同一配体、两种设计架构（紧凑融合 vs 连接杂合）下，独立两端对接协议的有效性边界。

### Candidate RQs

| # | Candidate | FINER 均分 | 去留 |
|---|-----------|------------|------|
| C1 | 朴素拼分是否抬高只强一端？ | ~3.2 | **降级**：承重假设不稳；新颖性弱 |
| C2 | 独立两端对接在 linked 上是否显著差于 fused？失效能否用「构象/连接一致性」刻画？ | **4.2** | **主 RQ** |
| C3 | 再训一个联合双姿态生成器能否超越 FuseDiff/LaMGen？ | ~2.8 | 否：F 低、N 低（拥挤） |
| C4 | 细胞活性能否证明双结合？ | ~2.0 | 否：标签语义错 |

### Primary Research Question（锁定）

> Under a frozen independent dual-docking protocol, does performance (pose recovery where available, enrichment, and dual-vs-single ranking) **collapse specifically on linked/hybrid dual ligands relative to fused/merged chemotypes**, and can a **connectivity-/ensemble-aware consistency check** (not a new sampler) restore dual-vs-single discrimination on architecture-stratified benchmarks?

### FINER Assessment（主 RQ）

| Criterion | Score | Justification |
|-----------|-------|---------------|
| **F**easible | **3** | 双端共晶少是硬约束；可用「单端共晶+类型分层富集+实验室系列」补；不能承诺大规模双姿态 RMSD 榜。 |
| **I**nteresting | **5** | 直接打中「对接不准 / 类型怎么办」；与拼分叙事切割清晰。 |
| **N**ovel | **4** | 相对 Pérez 融合与生成式双靶，增量是**架构分层的失效诊断 + 一致性协议**；非新扩散模型。 |
| **E**thical | **5** | 公开数据+有限私有盲测。 |
| **R**elevant | **5** | 直接服务双靶药化筛选与合成优先级。 |

**平均 ≈ 4.2 / 5 · 通过门槛。**

### Scope

**IN**  
- 普通小分子双靶（fused / merged / linked）；独立对接协议的**分层审计**  
- 双端共晶个案作姿态锚；文献 curated 类型集；只强一端硬负  
- 可选：连接一致性 / 构象系综门控（协议层，非新采样器）

**OUT**  
- 以「更好 mean/softmin」为唯一创新  
- 与 FuseDiff/LaMGen 比生成 SOTA  
- 细胞 = 双结合金标准  
- PROTAC 三元主线  

---

## 4. Idea-Quality（对新 idea 打分）

### Idea Card（建议主 idea）

```markdown
# Architecture-Stratified Dual Docking Audit

## Motivation
Dual-target docking almost always docks one ligand twice independently, then aggregates scores.
This treats a single connected molecule as two unrelated ligands. That assumption may hold for
compact fused/merged chemotypes but fails for linked hybrids whose viable poses are coupled by
the linker. Score fusion cannot repair missing coupled pose information. The field lacks an
architecture-stratified diagnostic that isolates this failure and a protocol-level consistency
gate before dual ranking.

## Method
1. Curate architecture labels (fused/merged/linked) for paired dual-activity molecules + Tier-A
   both-end co-crystals.
2. Freeze independent dual docking (GNINA/Vina + PoseBusters).
3. Report all metrics stratified by architecture; include A-only/B-only hard negatives.
4. Define connectivity/ensemble consistency checks (linker strain, pocket clash, shared-graph
   pose plausibility) as gates — not as a new generative sampler.
5. Show: independent docking + fusion looks OK on fused, collapses on linked; gating restores
   dual-vs-single ranking; kill-switch = gate off → return to failure.
6. Optional holdout: lab fused vs linked series ranking concordance with dual biochemical labels.
```

### Axis scores（idea-quality）

| Axis | Score | Reason |
|------|-------|--------|
| **A Problem position** | **4** | 缺口非显然：把失效归因于「评价对象/架构」而非「加权」；命名本身有信息量。 |
| **B Method quality** | **3** | depth 3（诊断+协议，非新算子）；soundness 4；feasibility 3（共晶稀缺）。 |
| **C Problem-fit** | **5** | 方法直接对准「对接不准 × 类型」；不再错打到拼分。 |

**Overall：`round(100*(4+3+5-3)/12) = 75` · Verdict: strong**  
（相对旧拼分 idea 的 ~67：提升来自 **C**，不是更大算子。）

### 对照：若仍坚持拼分主 idea

A3 / B2 / C2 → 明显更弱；且承重假设已被否。**不建议冲高分。**

---

## 5. Scoop-Check（对新 novelty）

**Research problem：** 独立双对接在不同双靶化学架构上是否系统性失效，以及如何用协议层一致性门控纠偏。  

**Novelty：** 架构分层的双对接失效诊断 + 连接/系综一致性门控协议 + 开放分层基准（非新生成器、非新加权公式）。

### 四轴分解

- **Framing：** 架构分层双对接审计（fused vs linked）  
- **Mechanism：** 一致性门控 + 分层评测（非 FuseDiff 式联合生成）  
- **Insight：** 独立对接成功可能是 fused 化学型的假象；linked 才暴露协议破产  
- **Domain：** 结构基双靶 / MTDL  

### 高威胁先验

| 先验 | 重叠 | 威胁 |
|------|------|------|
| Pérez-Castillo 2017 融合打分 | framing 部分 · mechanism 不同 | 低–中（你们已主动离开融合主叙事） |
| Zhou 2013 双靶对接假阳性 | insight 近 · 无架构分层 | 中（必须 cite；你们用类型分层深化） |
| TwistDock 2019 | 仅 bivalent 特例 | 低（通解声明需克制） |
| POLYGON / LaMGen / FuseDiff | 都是**生成** | **高若你们改去做生成**；作诊断协议则 **低–中**（Related Work 切清） |
| VSDS-VD 2025 | 评测美学近 · 单靶 | 中（叙事可对标，任务不同） |
| Proschak 2024 DTDL 数据画像 | 化学侧 | 低（互补：他们挖数据，你们挖对接协议失效） |

**Scoop 等级（诊断协议主线）：Level 3 — Medium Overlap（可辩护）**  
**Scoop 等级（若改回联合生成 SOTA）：Level 2 — 高危，不建议。**

**Delta 句（投稿用）：**  
Unlike score-fusion dual docking and dual-target generative models, we isolate **ligand architecture** as the confound that breaks independent dual docking, release an architecture-stratified audit, and show a **connectivity-consistency gate** (not a new sampler) is necessary for dual-vs-single decisions on linked hybrids.

---

## 6. 高分文应长什么样（流程重写）

```text
1. Intro：双靶对接的隐藏假设 =「一分子 = 两次无关单靶对接」
2. 诊断：按 fused / linked 分层后，独立对接+拼分的表现分裂
3. 机制：linked 失败来自构象耦合/连接约束未被建模（用例与共晶锚）
4. 协议：一致性门控 + 硬负（只强一端）+ 双负例设定
5. 杀伤开关：关掉门控 → linked 上 dual-vs-single 回落
6. （可选）实验室 fused vs linked 系列外推
7. 开源：分层数据、YAML、门控脚本
```

**明确不写：** 「我们找到了最好的 softmin。」  
**明确解决的问题：** **双靶对接评价协议在设计类型上的系统性盲区。**

### 期刊预期（诚实）

| 出口 | 条件 |
|------|------|
| **Nat. Mach. Intell. 类** | 诊断反直觉表扎实 + 多靶点对 + 协议可复现 + 主张克制；**仍是条件可达，非默认成功** |
| **Chem. Sci. / JCIM** | 更稳妥的主出口 |
| **J. Med. Chem.** | 若主证据变成实验室 fused/linked 系列 + 计算协议，药化叙事可走，但创新定位变「工作流+对照」 |

---

## 7. 与旧材料怎么处理

| 旧产出 | 新定位 |
|--------|--------|
| 拼分 / Dual-VSDS 融合打假 | **附录或负对照**，证明「只拼分不够」 |
| 三组公开靶点对 | **仍用**，但必须加架构分层标签 |
| 文献 100+/高分刊池 | **优先补 fused/linked 标签**（服务 Q4） |
| 双端共晶 catalog | **姿态锚**，不再假装有大规模双 RMSD 榜 |
| NLRP3/JNK1 细胞 | **系列外推**，不作结合金标准 |

---

## 8. Implementability / 证伪（ResearchStudio）

**可完成：** 是（JCIM/Chem Sci 默认；NMI 条件）。  

**Kill-switch：**  
若 fused 与 linked 在独立双对接上**没有**可重复的分层差异，则「架构混淆」假说破产 → 改发数据报告或退回路线收缩，**禁止**再包装成拼分创新。

**最大风险：** 双端共晶太少导致姿态证据薄 → 应用「富集 + 硬负 + 实验室系列」三角支撑，并在摘要写明姿态锚的 n。

---

## 9. 最终拍板建议

1. **要解决的双靶问题（唯一主问题）：**  
   **独立两端对接协议忽视分子设计类型（尤其连接型）导致的系统性评价失效。**  

2. **不要再当主问题的：**  
   找一个更好的对接分加权公式；追赶联合生成 SOTA。  

3. **下一步执行（Scoping → 实验，仍不要写全文）：**  
   - 给冻结三对 + 文献金标准分子打上 fused/merged/linked；  
   - 先出一张「按类型分层的独立双对接表现」表（有无门控都先跑基线）；  
   - 若分层差异出现 → 再设计一致性门控与杀伤开关；  
   - 若分层差异不出现 → 停止该高分叙事。

**Skills 共识句：**  
旧拼分 idea 因承重假设错误而 **problem-fit 失败**；新 idea 用 `assumption_audit_and_pivot` 转到架构分层诊断，idea-quality 升至 **strong (~75)**，scoop **Level 3 可辩护**——这是在「仍想发高分」约束下，目前最不靠谱程度最低的一条路。
