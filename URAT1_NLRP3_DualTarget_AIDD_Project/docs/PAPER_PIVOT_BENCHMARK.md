# 论文方向调整说明（务实版）

> **状态**：2026-06 定稿方向  
> **替代**：原 TAPE-GATE / 双靶 AIDD 算法论文主线（见 `MANUSCRIPT_OUTLINE.md` 旧版，**不再作为投稿叙事**）  
> **执行计划**：`URAT1_THREE_STATE_BENCHMARK_PLAN.md`  
> **双轨策略（快速发表 + JCIM）**：`TWO_PAPER_STRATEGY.md`

---

## 一、修订后的文章方向

### 1.1 一句话定位

**不是**「我们提出了 URAT1/NLRP3 双靶 AI 发现新算法」，**而是**：

> **系统评估 URAT1 抑制剂在三态 cryo-EM 结构上的刚性对接何时有效、何时失效，并给出可复现的替代打分协议。**

这是 JCIM / J. Cheminformatics 近年更接受的 **benchmark + critical evaluation** 叙事，而非未经验证的 discovery pipeline。

### 1.2 推荐标题（择一）

1. *Benchmarking URAT1 inhibitor docking across inward, occluded, and outward cryo-EM conformations: when rigid ensemble docking fails and what to do instead*
2. *A critical evaluation of three-state Glide protocols for the human urate transporter URAT1*
3. *URAT1 three-state docking benchmark: structure mapping, pose viability, and rescoring limits for SLC transporters*

### 1.3 核心科学问题

| 问题 | 本文回答方式 |
|------|----------------|
| 三态 PDB 应如何正确映射？ | 9DKB / 9B1K / 9B1L（纠正 9JDZ 误用） |
| 刚性 Glide 三态对接是否普适可行？ | **否** — 实测 4 药在 9B1K/9B1L 零 pose |
| 失败时工程上可接受什么替代？ | inward dock + 叠合 pose 转移 + Prime/MM-GBSA / 惩罚项 |
| 替代协议能否区分已知活性药与 decoy？ | Gate + enrichment 定量检验 |
| 单态 inward 对接是否足够？ | 与 Protocol C 对比，诚实讨论信息损失 |

### 1.4 论文贡献（投稿可写 4 条）

1. **结构映射校正**：明确 URAT1 inward/occluded/outward 的权威 PDB 对应关系及文献依据（Dai 2024; Suo 2025; Wu 2025 局限）。
2. **协议对比基准**：Protocol A（单态 9DKB）/ B（刚性三态 Glide）/ C（inward + pose 转移 + rescoring）/ D（IFD 敏感性，仅 4 药）的系统比较。
3. **开放可复现资源**：Gate 四药 + benchmark 化合物 + property-matched decoys 的 SMILES、对接参数、预计算分数与 Maestro 工作流文档。
4. **方法学警示**：结合 Sindt 2025 等对 rescoring 局限的讨论，报告何时 $S_\pi$ 方向性有效、何时仅为启发式。

### 1.5 降级或移出主线的内容

| 原叙事 | 新定位 |
|--------|--------|
| TAPE-GATE 双靶融合 | 不写入主文；最多在 Discussion「未来工作」一句 |
| Teacher M-CPDL / 8973 蒸馏 | **暂停**；Gate 未过前不跑全量 |
| Path B 生成式优化 | 补充材料或独立后续项目 |
| NLRP3 assay-conditioned ML | 若有 rigor 统计可写 **短文/第二篇**，非本篇主线 |
| OAT 迁移学习（+0.004 Spearman） | 技术笔记或 SI，**不作创新点** |
| 「双靶发现漏斗」 | **禁止**作为主 claim |

### 1.6 目标期刊

| 期刊 | 匹配度 | 说明 |
|------|--------|------|
| **JCIM** | 高（若 benchmark 完整 + 开源） | Application / Benchmark 类，非 Methods 新算法 |
| **J. Cheminformatics** | 高 | 资源型、协议型文章更宽松 |
| **J. Chem. Inf. Model. 姊妹刊** | 中 | 需更强 decoy 集与统计 |
| **SLAS Discovery / Mol. Informatics** | 备选 | 若 JCIM 拒稿 |

---

## 二、为什么「新算法」路线在你当前条件下不可行

以下不是泼冷水，而是对照 **JCIM 2025–2026 近期发文标准** 与 **你已有的实证结果** 得出的结论。

### 2.1 数据规模不足以支撑「新深度学习架构」

| 数据集 | 规模 | JCIM 新算法论文常见要求 |
|--------|------|-------------------------|
| URAT1 活性 | 822 SMILES | 多任务 benchmark + 外部验证集 |
| NLRP3 | 513 SMILES | 同上；且 7.2% 跨 assay 冲突 |
| 两靶重叠 | **0** | 无法训练标准双靶 MTL 并宣称协同 |

近期 JCIM 方法论文（如 P-gp GNN、IOMemP）通常包含：**新架构 + 多靶/多构象公开基准 + 与 5+ 基线严格对比**。你的数据量更适合 **评估现有工具**，而非训练新模型。

### 2.2 URAT1 ML 基准已削弱「AI 创新」叙事

回顾性 benchmark：**4 个必回收阳性药中 2 个 NO_GO**（`URAT1_NO_GO`）。在审稿人眼中，这意味着：

- 不能声称「我们的 ML 优于传统对接」
- 任何 fusion / Teacher 标签都建立在 **未验证的弱证据** 上
- 继续包装成算法论文会被要求 **湿实验验证 hit**，而你目前没有

### 2.3 结构瓶颈：三态刚性对接实测失败

你已验证的事实：

- **9DKB**：4 药 SP→XP 均有 pose ✓
- **9B1K / 9B1L 刚性 Glide**：4 药 **零 pose**（`GRID-ENERGY MIN FAILED`, GlideScore=10000）
- **urate 对照**在 B1K/B1L 有 pose → grid 本身可用，**药物 pose 不可行**
- 叠合后配体共定位但 **大量 clash**

因此：

- 「三态 Glide 系综 + Boltzmann $\pi$」**不是**已验证算法，而是 **待 benchmark 的假设**
- Pose 转移 + MM-GBSA（Scheme 2）是 **工程折中**，文献有 ensemble rescoring 先例，但 **无 URAT1 专用先例** — 只能作为 benchmark 的一种协议 C，不能包装成 novel method

### 2.4 算力与许可约束

| 任务 | 规模 | 现实 |
|------|------|------|
| 全量三态 Glide | 8973 × 3 ≈ 2.7 万 job | B1K/B1L 大部分会失败，**浪费许可** |
| MM-GBSA × 三态 × 全库 | 极高 | 无集群自动化 pipeline |
| IFD | 每分子数小时 | 仅适合 4–12 个敏感性案例 |

JCIM 级「新算法」往往需要 **可扩展、可开源的端到端 pipeline**；你目前是 **手工 Maestro + 局部脚本**，适合 benchmark 论文，不适合 claim scalable AI method。

### 2.5 时间与专长结构不匹配

完整新方法论文通常需要：

1. 架构设计与实现（1–3 月）
2. 多靶外部验证（需合作或公开数据）
3. 与 jcim.5c01609 一致的统计 rigor（效应量、scaffold CV、practical significance）
4. 可选湿实验（JCIM application 类几乎必需）

你已有的强项是：**结构生物学事实核对、对接实操、诚实的失败记录** — 这恰好是 **benchmark 论文** 的素材，不是 **算法论文** 的素材。

### 2.6 OAT 迁移与双路径不构成「practically significant」创新

- OAT→URAT1 迁移：Spearman **~0.726 vs ~0.722**（Δ ≈ 0.004）→ 统计上可能显著，**方法学上无意义**
- Path A/B 双路径：**脚本骨架存在，无回收率结果**
- MASFL Teacher：**未实现** `02_teacher_mcpdl.py`

继续强行包装为「新算法」会在审稿第一轮被 **lack of novelty / insufficient validation** 拒稿。

### 2.7 无湿实验 → 不能写 discovery 类 application

JCIM 2025–2026 的 transporter / inflammasome **发现类**文章普遍含：

- 酶学或细胞 IC50 验证至少 1–3 个新 hit，或
- 对已知药的机制新见解 + 多种 biophysical 证据

你没有实验合作时，**唯一诚实且可发表的角度是计算方法比较**，不是「我们发现了双靶候选」。

---

## 三、接下来你需要做什么（分阶段清单）

### Phase 0 — 冻结叙事（立即，1–2 周）

- [ ] 通读本文 + `URAT1_THREE_STATE_BENCHMARK_PLAN.md`，确认团队一致
- [ ] **停止** 8973 全量 B1K/B1L 刚性对接
- [ ] **停止** 把 Teacher M-CPDL / 双靶 funnel 当主故事写摘要
- [ ] 用 `MANUSCRIPT_OUTLINE_BENCHMARK.md` 替换写作大纲
- [ ] README 首段改为 benchmark 项目描述（可选）

### Phase 1 — 完成 Gate 与协议对比（计算核心）

**对象**：4 药 + lesinurad redock + 8–12 benchmark 化合物 + 50–200 property-matched decoys

| 协议 | 内容 | 优先级 |
|------|------|--------|
| **A** | 仅 9DKB Glide SP→XP | P0 — 已有部分结果 |
| **B** | 刚性三态 Glide | P0 — **记录失败率**（这是重要结果） |
| **C** | 9DKB dock → 叠合至 B1K/B1L → Prime Minimize → MM-GBSA → $S_\pi$ | P0 — Scheme 2 |
| **D** | 9B1K IFD（仅 4 药） | P1 — 敏感性上界 |

**Gate 通过标准**（详见 `TEACHER_GATE_QC_DATASETS.md`）：

- Gate 1：lesinurad redock RMSD ≤ 2.0 Å（9DKB）
- Gate 2：四药 Protocol C 下 $S_\pi > 0$（**4/4**）
- Gate 3：活性集 A vs decoy D 的 median($\pi_{in}+\pi_{oc}$) 分离

### Phase 2 — 评估与统计（JCIM 风格）

- [ ] 指标：pose 可行率、redock RMSD、$S_\pi$ 方向、enrichment@k、AUC（对 decoys）
- [ ] Protocol A vs C 在已知活性药上的 rank / enrichment 对比
- [ ] 引用 Sindt 2025（jcim.5c00730）讨论 rescoring 局限 — **诚实写**
- [ ] 按 jcim.5c01609：报告效应量与置信区间，不只 p 值
- [ ] 记录叠合 RMSD（你已有 B1K ~9.3 Å, B1L ~18.9 Å）作为局限

### Phase 3 — 开源与可复现

- [ ] GitHub release 子集：`URAT1-3State-Docking-Benchmark`
- [ ] 包含：SMILES、yaml 配置、预计算分数表、Maestro 操作 SOP
- [ ] 实现 `scripts/utils_three_state_scoring.py`（Boltzmann $\pi$、惩罚规则）
- [ ] 可选：`scripts/run_three_state_benchmark.py` 汇总各协议输出 CSV

### Phase 4 — 撰稿

- [ ] 按 `MANUSCRIPT_OUTLINE_BENCHMARK.md` 写初稿
- [ ] 图表：三态示意图、协议对比表、enrichment 曲线、4 药 clash 案例图
- [ ] SI：完整对接参数、decoy 生成规则、失败 log 摘录
- [ ] **不写** 双靶发现结论；Discussion 明确 transporter vs kinase docking 范式差异

### Phase 5 — 可选升级（有资源再做）

| 升级 | 解锁能力 |
|------|----------|
| 1–3 个化合物湿实验 IC50 | 可投 application 更强叙事 |
| NLRP3 assay-conditioned 严格统计 | 第二篇短文 |
| 真正 transporter-aware GNN + 外部数据 | 方法论文（长期） |

---

## 四、现在明确不要做的事

1. ❌ 跑 822+8000 条 × 3 grid 刚性 Glide（已知 B1K/B1L 对药物失败）
2. ❌ 在 Gate 2 未 4/4 前启动 Teacher 标签与 PC-Student
3. ❌ 投稿摘要写「novel dual-target AIDD algorithm」
4. ❌ 把 OAT +0.004 Spearman 写成主要创新
5. ❌ 用 9JDZ 代表 occluded/outward（已纠正，勿回退）

---

## 五、与仓库其他文档的关系

| 文档 | 角色 |
|------|------|
| `URAT1_THREE_STATE_DOCKING.md` | PDB 映射与 Glide 流程 — **仍有效** |
| `TEACHER_GATE_QC_DATASETS.md` | Gate 数据集 — **Phase 1 直接用** |
| `INNOVATION_POINTS.md` | 旧 TAPE-GATE 创新点 — **参考用，不作投稿主线** |
| `MANUSCRIPT_OUTLINE.md` | 旧大纲 — **已废弃主叙事** |
| `MANUSCRIPT_OUTLINE_BENCHMARK.md` | **新写作大纲** |
| `MASFL_V3_WORKFLOW.md` | 全量蒸馏 — **Gate 通过后 Phase 5** |

---

## 六、心理预期（务实）

- **可发表**：一篇诚实的 URAT1 三态对接 benchmark，对领域有价值（纠正 PDB 误用 + 报告刚性对接失败 + 提供替代协议比较）
- **不可发表（当前）**：双靶 AI 发现、Teacher 蒸馏、生成式 hit 列表
- **JCIM 不是不可能**，但文章类型要从 **Algorithm** 改为 **Benchmark / Application (computational only)**

把「对接失败」写进论文，在 JCIM 是 **加分项**（critical evaluation），在算法包装叙事里是 **致命项** — 这就是必须改方向的根本原因。
