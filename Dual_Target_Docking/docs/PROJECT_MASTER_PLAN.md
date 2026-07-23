# 课题总规划：双靶兼容性打分与 Dual-VSDS

> **现行投稿主线已更新：** 双药效团 / passenger 污染 + moiety-resolved 评测。  
> → **NMI 执行规划：** [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md)  
> 下文保留早期「融合打假」总览；分数融合降为基线，fused/linked 为分层轴。  
> **立项依据专文：** [`TOPIC_ESTABLISHMENT_REPORT.md`](TOPIC_ESTABLISHMENT_REPORT.md)。

| 来源 | 文档 | 核心结论 |
|------|------|----------|
| **现行 NMI 规划** | [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md) | passenger 诊断 + moiety-resolved + Dual-VSDS-Moiety |
| ResearchStudio | [`researchstudio_audit/RESEARCHSTUDIO_AUDIT.md`](researchstudio_audit/RESEARCHSTUDIO_AUDIT.md) | 早期 fusion 路线；需按 moiety 重估 |
| ARS | [`ars_audit/ARS_TOPIC_ANALYSIS.md`](ars_audit/ARS_TOPIC_ANALYSIS.md) | 有打假表前不写作（仍适用） |
| 旧主方案 | [`NMI_DUAL_COMPATIBILITY_PLAN.md`](NMI_DUAL_COMPATIBILITY_PLAN.md) | 融合主线（历史） |

---

## 1. 这个课题的内容是什么（What）

### 1.1 一句话

**现行：** 单分子双靶含两药效团；整分子单口袋打分被 passenger 污染 → **moiety-resolved 打分 + 分靶校准 + 短板双靶决策 + 开放基准**。  
**历史表述：** 把「同一分子两端是否同时够好」做成可校准的双靶兼容性排序 + Dual-VSDS。

### 1.2 研究对象与非对象

| 做 | 不做 |
|----|------|
| 同一配体 → 靶 A、靶 B **各做一次独立对接** → 融合两侧信号 | 新 DiffDock/采样器 |
| 标签：dual / A-only / B-only / inactive（未测≠inactive） | 药物联用协同预测 |
| 校准到每靶 \(p(\mathrm{active})\) 后短板融合（如 softmin） | PROTAC 三元桥接主线 |
| 开放 Dual-VSDS + 泄漏控制划分 | 把细胞活性当双靶结合金标准 |
| 私有 NLRP3/JNK1 细胞数据作 **时间盲测排序** | PK/成药性预测模型 |

### 1.3 主研究问题（ARS 锁定）

> 在泄漏控制的配对活性标签下，**朴素 mean / min / rank 融合是否系统性抬升 A-only/B-only？**  
> 每靶校准后的短板融合，能否在多靶点对上改善 **dual-vs-single** 排序？

子问题：TrueNegative vs RandomDecoy 结论是否翻转；fused/linked 是否需要分治；公开协议能否外推到 NLRP3/JNK1 细胞 holdout（排序相关，非结合证明）。

### 1.4 交付物长什么样

1. **方法**：可插拔对接后端（GNINA/Vina ± RTMScore + PoseBusters）+ 校准短板融合头  
2. **基准**：Dual-VSDS（多靶点对、四类标签、两套 decoy、scaffold/leave-pair 划分）  
3. **证据**：一张「打假表」证明朴素融合有害 + 消融证明短板机制必要  
4. **开源**：数据 schema、YAML、脚本、Zenodo  
5. **可选附录**：私有细胞 holdout；细胞–PK 脱钩案例（不宣称预测 PK）

---

## 2. 为什么要进行这个课题（Why）

### 2.1 现实痛点（你们实验室语境）

- 已有 **NLRP3/JNK1** fused & linked 分子，部分细胞有效，但动物 PK 差。  
- 计算侧若继续「两端打分相加选候选」，会把 **单靶很强、另一端很弱** 的分子推到前列 → 浪费合成与筛选。  
- 短期很难补齐全面亲和力/动物药效 → 适合做 **方法+基准** 文，而不是成药故事。

### 2.2 领域缺口（skills 共同认定）

| 已有 | 仍缺 |
|------|------|
| Pérez-Castillo：双靶 rank 融合案例 | 跨靶**概率校准** + dual-vs-single **硬负** |
| VSDS-VD：单靶 VS 分层 decoy | **配对分子**的双靶任务与基准 |
| EquiScore/CleanSplit：泄漏与打分纪律 | 把纪律迁到**双靶决策** |
| POLYGON 等：双靶**生成** | **判别/排序**尺子 |

### 2.3 为什么值得做（收益）

- **科学：** 纠正双靶 VS 默认协议的系统偏置（CleanSplit/VSDS-VD 同构的「打假」贡献）。  
- **工程：** 给实验室一条「先过兼容性门槛再谈细胞/PK」的筛选 SOP。  
- **发表：** 资源型基准 + 诊断发现，比再卷一个 docking SOTA 更贴 NMI；保底可发 JCIM。

### 2.4 为什么「现在」做

- 方案与文献对标已齐；skills 判定 **idea 够格、未被 scoop 死**。  
- 卡点在执行（数据与打假表），继续空转文档边际收益低。

---

## 3. 怎么开展这个课题（How）

### 3.1 总原则

1. **问题驱动方法**：先证明尺子坏了，再给修正协议。  
2. **零件可已知**：GNINA、RTMScore、PoseBusters、Z-score、softmin 都可以旧。  
3. **门控实验**：没有「mean 抬升 A-only」表 → 不写 NMI 主文、不宣称创新成立。  
4. **期刊双轨**：NMI = 条件目标；JCIM/Chem. Sci. = 保底出口。

### 3.2 六阶段路线（执行版）

| 阶段 | 目标 | 完成定义（Definition of Done） | 建议工具 |
|------|------|--------------------------------|----------|
| **P0 冻结** | 锁死可复现参数 | θ（如 pAct≥6）、主指标、GNINA/Vina 版本、≥2 公开靶点对、校准/测试划分规则书面化 | 文档 |
| **P1 数据** | D2-public + D2-curated + D4 decoy + D5 holdout 入库 | 断言测试：未测∉inactive；InChIKey 去重；scaffold 标签可复现 | K-Dense: database-lookup, bioservices, rdkit/datamol；statistical-power 定最小 n |
| **P2 对接** | 两端独立对接特征表 | YAML+种子；PoseBusters 门控；导出 vina/cnn/rtm/pb | 自管 GNINA；DiffDock 仅消融；modal 可选 |
| **P3 打假★** | 证明朴素融合有害 | Top-1%/5% 中 A-only 污染率：mean/sum ≫ 校准短板；严格 split | statistical-analysis |
| **P4 融合** | 校准 softmin + 硬负 +（可选）类型专家 | 负对照：\(\tau\to\infty\) 或退回 raw mean → 指标回落；fused/linked 开关消融 | 自写融合脚本 |
| **P5 外推+开源** | leave-pair-out + D5 排序盲测 + 发布 | Zenodo/GitHub；细胞主张不超过「排序 holdout」 | schematics / writing（有表后） |

### 3.3 最小充分实验集（不够就不要扩 scope）

1. ≥**2** 个公开靶点对（建议再加 1 对做 leave-pair）  
2. 同一姿态特征上：raw mean / min / z-mean / rank / calibrated softmin  
3. 两套 decoy：TrueNegative + RandomDecoy  
4. 一条 kill-switch 负对照（短板关掉 → 回到基线）  
5. （可选）NLRP3/JNK1 细胞 holdout 排序相关图  

### 3.4 明确不做的并行线

- 再训通用 pose scorer  
- 动物药效 / PK 预测主文  
- PROTAC 三元主线  
- 五线课题并行（JNK1 选择性等可作案例，不作第二主创新）

### 3.5 投稿叙事骨架（对齐 VSDS-VD）

1. Intro：默认 mean/sum 有害（打假预告）  
2. Dual-VSDS 构建与泄漏协议  
3. 对接与校准方法  
4. Results：打假 → 校准短板 → 消融 → 双 decoy → holdout  
5. Discussion：细胞≠结合；PK 脱钩；局限  

---

## 4. 课题的创新点是什么（Novelty）

### 4.1 主创新（只打 1–2 条）

| # | 创新点 | 一句话证据形态 |
|---|--------|----------------|
| **N1** | **任务重定义**：双靶评估 = 校准后的短板约束 + dual-vs-single 决策 | 打假表：mean 抬 A-only；softmin 降污染 |
| **N2** | **Dual-VSDS 开放基准**：配对标签 × 设计类型 × 双 decoy × 泄漏控制 | 别人能复现你们的表 |

### 4.2 次要 / 支撑（不单独撑 NMI）

- fused/linked 条件专家（需消融，防「堆启发式」）  
- 细胞 holdout + PK 脱钩案例（叙事完整，非方法核心）

### 4.3 明确「不是创新」

- softmin / Z-score / GNINA / RTMScore / PoseBusters 本身  
- 「我们对接精度更高」类 RMSD leaderboard  

### 4.4 相对先验的 Delta（cover letter 用）

> Unlike Pérez-Castillo et al. (rank fusion on one dual pair) and VSDS-VD (single-target VS audit), we treat dual assessment as a **calibrated shortfall / dual-vs-single ranking problem**, release **Dual-VSDS**, and show standard aggregations systematically elevate single-target–biased molecules.

Scoop 等级：**Level 3（中等重叠）**——可辩护，但 Delta 后半句必须有表。

### 4.5 Skills 对创新强度的共识

- ResearchStudio：方法深度中等（零件旧），**问题位置强** → 靠诊断+基准上分  
- ARS：FINER Novel=4；禁止写「细胞验证了双靶结合」「全球首个」绝对句  
- 模式对标：`controlled_diagnostic_design`（尺子审计）+ `reframe_as_solvable_object`

---

## 5. 离真正落地还有多远（Gap）

### 5.1 成熟度总览

```text
想法/文献/叙事     ████████████░░░░  ~75%  ← 已完成（本仓库文档）
可复现参数冻结     ████░░░░░░░░░░░░  ~25%  ← θ/引擎/靶点对未锁
公开配对数据       ██░░░░░░░░░░░░░░  ~15%  ← schema 有，表未建
对接流水线         █░░░░░░░░░░░░░░░  ~10%  ← 无 YAML/批跑
打假表/融合代码    ░░░░░░░░░░░░░░░░   ~0%  ← 创新门控未过
开源发布/投稿包    ░░░░░░░░░░░░░░░░   ~0%
────────────────────────────────
整体落地（可投 JCIM 包）约 20–30%
整体落地（可冲 NMI 包）约 15–20%（还差效应量与多靶点对）
```

### 5.2 距离分层

| 层级 | 含义 | 现状 | 还差什么 |
|------|------|------|----------|
| **L1 方案可讲清** | 组会/开题说得通 | ✅ 已达成 | — |
| **L2 可复现实验可跑** | 外人按文档能出表 | ❌ | P0 冻结 + P1–P2 代码与数据 |
| **L3 创新主张可成立** | K1 打假显著 | ❌ | P3 门控表 |
| **L4 可投方法刊** | JCIM/Chem.Sci. | ❌ | L3 + 开源 + 消融 |
| **L5 可冲 NMI** | 反直觉+多对+协议质量≈VSDS-VD | ❌ | L4 + ≥2–3 对泛化 + 叙事克制 |
| **L6 实验室筛选 SOP** | 真用于 NLRP3/JNK1 选分子 | 部分概念可用 | L3 后把协议写进日常筛选；**不依赖**发 NMI |

**「落地」若指发高分文：** 关键路径在 **L3→L4**；NMI 是 L5 附加条件。  
**「落地」若指帮实验室少走弯路：** L3 打假成立即可开始改筛选规则，不必等 L5。

### 5.3 最大风险（按致命度）

1. **真 dual 配对太少** → 打假表方差大 / 创新落空 → 改发「协议+负结果」或降期刊  
2. **校准泄漏** → 假阳性创新（ARS Mode 4）  
3. **过早锁死 NMI** → frame-lock；打假失败仍硬写  
4. **细胞/PK 写进 Abstract** → 审稿直接杀  

### 5.4 近期行动清单（补完规划）

**本周（P0）**
- [ ] 冻结 θ、主指标列表、引擎版本字符串  
- [x] 选定 ≥2 公开靶点对（写出 UniProt/ChEMBL ID）→ **已锁 3 对**：PIK3CA/mTOR、EGFR/HER2、Mcl-1/Bcl-xL（见 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) / [`../data/public_pair_selection/FROZEN_PUBLIC_PAIRS.yaml`](../data/public_pair_selection/FROZEN_PUBLIC_PAIRS.yaml)）  
- [ ] 写死：校准只用 train scaffolds  
- [ ] 用 statistical-power 估 dual 最小 n；达不到则预设「降级期刊」开关  

**下一阶段（P1–P2）**
- [ ] ChEMBL 配对抽取脚本 + 断言测试  
- [ ] 文献 curated 30–80 条 design_type  
- [ ] GNINA YAML 跑通 1 个靶点对小样  

**创新门控（P3）**
- [ ] 产出打假表；**显著才进入 P4/写作**  
- [ ] 不显著 → 启动保底叙事（协议文 / JCIM）  

---

## 6. 三套 skills 结论一页纸

| 问题 | 共识答案 |
|------|----------|
| 课题能不能做？ | **能**；idea 合格，非 abandon |
| 会不会已被做完？ | **不会被 scoop 死**（Level 3）；须写清相对 Pérez 与 VSDS-VD 的差 |
| 最大漏洞？ | 数据稀疏、校准循环、缺 kill-switch、参数未冻 |
| 创新靠什么？ | **打假 + 任务重定义 + 开放基准**，不靠新对接器 |
| 先发哪？ | 条件冲 NMI；**保底 JCIM** |
| 现在最该干什么？ | **冻参数 → 建 D2 → 对接 → 打假表**；不要写全文 |
| K-Dense 角色？ | 执行加速（库/分子/统计）；不替代科学门控 |

---

## 7. 给负责人的五句定调

1. **内容：** 双靶兼容性排序协议 + Dual-VSDS 基准。  
2. **动机：** 默认分数融合会系统性错选单靶偏倚分子；实验室需要可计算的筛选尺子。  
3. **路径：** 冻参数 → 配对数据 → GNINA 特征 → 打假 → 校准短板 → 开源 → 双轨投稿。  
4. **创新：** 证明朴素融合有害，并给出校准短板 + dual-vs-single 评测与开放基准。  
5. **距离：** 叙事已 75%，可投稿包约 20–30%；**离「创新落地」只差一张过关的打假表和前面的数据/对接工程。**
