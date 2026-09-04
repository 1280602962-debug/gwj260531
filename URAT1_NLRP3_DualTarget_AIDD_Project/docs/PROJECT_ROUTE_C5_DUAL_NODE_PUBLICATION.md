# C5 路线：URAT1–NLRP3 双节点候选论文（目标 *Molecular Diversity*，无湿实验）

> 取代 C2 / C3 / C4 作为**唯一现役主线**。C2（转运周期 Δ 打分）、C3（结构时间窗口叙事）、C4（选择性改题）各自的可用零件被吸收进本路线，不再作为独立主线。
> 冻结资产不覆盖：`data/repurposing/p2/`、`data/campaigns/c1/07_clinical_dock/acid_dual_a1_frozen/`、`data/si/`。
> 门控在看到任何新分子名之前锁定于 `config/campaign_c5.yaml`。

---

## 0 一句话产品

> 以两个实验结构臂各自**独立回顾校准**的几何/相互作用门控，从 8,319 个临床阶段分子的酸性等效子集中，提出并分层一组可实验验证的 URAT1–NLRP3 **双节点候选假说**；对接分数全程不作亲和力解释。

标题口径（草案）：
*Dual-anchor geometric screening of clinical-stage acid pharmacophores identifies putative URAT1–NLRP3 dual-node leads*

关键定语必须保留：**putative / hypotheses / prioritized**。不写 dual inhibitors，不写 discovered。

### 为什么这个形态能投这本刊

对照近半年 *Mol. Divers.* 无湿虚筛文（Klebsiella PBP3+β-内酰胺酶萜类双靶、µ 阿片偏倚激动剂、CERS2、PCOS 3β-HSD1、M3 抗毒蕈碱），过审靠五件事：具名候选、疾病/靶点钩子、漏斗厚度、**化学成分**、可复现外观。本项目全部具备，且多出同行普遍没有的两项：

1. **两臂门控各自有回顾性富集统计**（同行几乎只做"对接→MD 稳定→点名"的自指闭环）；
2. **多种子稳定性 + 明确的对照失败停规则**。

期刊 scope 明确排除"没有显著化学成分的计算生物学"。因此化学内容必须显性化：酸性等效基团分布（羧酸 272 / 四唑 9 / 酰基磺酰胺 27）、Murcko 骨架多样性、成药性与 ADMET，不能只有几何数字。

---

## 1 现在就能用的结果（不需重算）

### 1.1 主文级（可直接进正文图表）

| 资产 | 数值 | 文件 | 正文角色 |
|---|---|---|---|
| NLRP3 自对接对照 | NP3-146@7ALV CNNscore 选姿 RMSD **0.82 / 0.67 / 0.68 Å**（seeds 42/43/44） | `05_metrics/`、`RESULTS_DRAFT_CN.md` §3.7.1 | 结构方法有效性的正对照 |
| URAT1 酸门回顾校准（A1 规则） | 228 羧酸 active vs 64 羧酸 true decoy；OR = **3.18**，95% CI **1.77–6.81**，Fisher *p* ≈ **4.5×10⁻⁴**；sens 0.447 / spec 0.797 / PPV 0.887 | `05_metrics/acid_gate_retrospective_benchmark/acid_gate_benchmark_summary.json` | **全文最强可信度锚点**，同行罕有 |
| A2 规则的回顾表现 | OR ≈ **0.97**，*p* = 1.0；sens 0.952 / spec 0.047 | 同上 | 主文如实报告：A2 是**召回型可及性门**，不是判别器 |
| 酸池化学构建 | 1,588 → 酸等效 **303**（羧酸 272 / 四唑 9 / 酰基磺酰胺 27）→ 化学软过滤 **156** | `07_clinical_dock/acid_pool/acid_pool_summary.json` | 化学成分主体，满足 scope |
| A2 多种子双臂通过 | URAT1 121/120/121；NLRP3 宽松 78/77/83；双靶宽松 59/59/61 | `acid_dual_a2/acid_dual_summary_a2_seed4*.json` | 稳定性表 |
| A2b NLRP3 结构门 | 结构 74/74/75；双靶结构 **56/57/53**；≥2/3 ≈ **54**；3/3 = **38** | `acid_dual_a2/nlrp3_structural_summary_seed4*.json` | 漏斗末段 |
| 三种子交集 | dual 交集 **42**；≥2/3 dual **59**；化学审计后 eligible **40** | `08_nomination/acid_a2_eligible_audited.csv` | 漏斗表 |
| 候选逐分子结构指标 | overlap / IFP Jaccard / 关键接触 / Arg477 距离 | `08_nomination/acid_shortlist_a2_competition.csv` | 候选表 T1 |
| NLRP3 缩库模型 | AUROC **0.893**、AUPRC **0.914**、EF@10% 1.57（n = 513 分子 / 609 记录 / 25 assay） | `docs/MODEL_TRAINING_SUMMARY.json` | 仅"生物学缩库"，不作结合证明 |
| lesinurad 姿态失败 | 自由对接 CNNscore 选姿 RMSD **4.30 / 4.31 / 4.88 Å**；生产姿 Arg477 **14.2 Å**；晶体最小 6.7027 Å | `05_metrics/pass_fail.json`、§3.5 | **改写为"为何必须用几何优先门"的动机**，不作头条失败 |

### 1.2 SI 级（进补充材料，不进正文主线）

- 冻结 P2 协议比较 P0–P5 全表（True EF@1% 2.59 / Random 0.22；P5 2.80 / 0；P0 1.94 / 1.94）——`data/si/protocol_enrichment_ci/protocol_ef_ci.csv`
- 配对 bootstrap（P2 vs P1/P3/P4 显著，vs P0/P5 不显著）
- 诱饵相似性泄漏审计（RandomDecoy 骨架重叠 0，无 TC > 0.5）
- legacy P2 audit set（7 人）、裸 Pareto 四大环、τ 阈值敏感性
- URAT1 回归模型评估（RMSE 0.663 / R² 0.508 / Spearman 0.726；lesinurad、dotinurad 未回收 → 不用于主排序）
- 逐种子完整表

### 1.3 必须降级或废弃

- **Rank 轨全部内容不得作为产品**。P2 51/7 名单只以"协议迁移局限"出现在 SI。
- `MANUSCRIPT.md`、`MD_RUN_PLAN.md`、`CONCLUSIONS_DRAFT_CN.md`、`DISCUSSION_DRAFT_CN.md`、`METHODS_DRAFT_CN.md` 相对 HEAD 已过时，投稿前按 C5 重写。
- Vecabrutinib、Zelenirstat、Deucrictibant、Praliciguat、MLN-0415、BI 653048 不得出现在候选表。

---

## 2 投稿前必须修掉的三处（不修会被审出来）

### M1 已知配体面板有重复且缺阴性对照
`nlrp3_structural_panel` 阳性 10 个中 `CHEMBL3183703` 与 `MCC950` SMILES 相同；背景是 20 个临床酸，不是性质匹配阴性。当前 spec 仅 9/20。
**修法**：去重 → 阳性 n = 9；背景改为「临床酸 20 + 性质匹配诱饵 ≥ 40（相对 NLRP3 已知活性物 max Morgan TC ≤ 0.5）」。重算 Fisher / OR / CI。

### M2 A2b 结构门在本面板上没有增加特异性
seed42 宽松门与结构门通过集合完全相同，等于该门只对临床池起作用、对面板不区分。
**修法**：在 M1 扩充面板上**重新标定**四个阈值（overlap、IFP Jaccard、关键接触、clash），要求结构门在面板上相对宽松门**至少提升特异性**；若无法提升，则如实降级为"姿态质控"，不叫"结构兼容门"。阈值必须在看临床名单前锁定。

### M3 primary 候选内部矛盾
PF-04620110 三种子 NLRP3 结构门**全部不过**（overlap 0.48 / IFP 0.45 / 接触 4/7），却被列为 primary。在"候选发现"型论文里这是审稿人第一刀。
**修法**：primary 只保留 3/3 dual-structural 且化学干净者（现有：PSI-697、PF-03882845、Lanifibranor；Admilparant 为 2/3，可列 primary 但须标注 seed44 宽松失败）。PF-04620110 移至"通路证据支持、口袋可及但未恢复共晶模式"的单独类别，或直接移出候选表。GSK-3008348 保持结构对照，不作候选。

---

## 3 需要往下推进的计算（W1–W5，本地对接 + MD 可用）

设计原则：**新算力全部投在"让两臂对称、各自可校准"上**，不投在提高排序名次上。

### W1 URAT1 四联体交叉对接 —— 把自对接失败变成方法学结果
Suo/Fedor/Lee *Nat. Commun.* 2025, **16**:5178 同构建体四结构：

| PDB | 状态 | 分辨率 | 配体 |
|---|---|---|---|
| 9DK9 | apo | 2.68 Å | — |
| 9DKA | holo | 3.00 Å | benzbromarone (R75) |
| 9DKB | holo | 2.74 Å | lesinurad (LES) |
| 9DKC | holo | 2.55 Å | TD-3 (A1A45) |

**做法**：3 个晶体配体 × 4 个受体（含 apo）× 3 种子，`run_gnina_batch.py` + `config/docking_c1.yaml` 同参数；报告 Top-1 / Top-3 / best-of-9 RMSD 与酸根–Arg477 距离矩阵。

**为什么值得做**：lesinurad 的硫醚–乙酸臂是四者中最柔性的配体，4.3 Å 极可能是**配体柔性驱动**而非搜索盒或方法失败。刚性配体（benzbromarone、TD-3）若在同参数下自对接良好，就得到一个可发表的判断：*该口袋的姿态保真度依赖配体柔性，因此柔性酸必须用锚定几何门而非自由 Top-1*——这正好把 §3.5 的失败转成 W2 的方法依据。

**判据**：至少 2/3 刚性配体自对接 Top-1 ≤ 2.0 Å。若刚性配体也失败 → 受体准备/盒子有系统问题，回头查准备流程，不得直接进 W2。

### W2 URAT1 酸锚 IFP 门 —— 让 URAT1 臂与 NLRP3 臂对称
现状不对称：NLRP3 有 IFP 门（overlap + Jaccard + 关键接触），URAT1 只有**单一 Arg477 距离**。单距离门就是 A2 那个 OR≈0.97 的东西。

**做法**：
1. 从 9DKB(lesinurad)、9DKC(TD-3)、9DKA(benzbromarone) 三个**实验尿酸排泄剂姿态**提取共识关键残基集（Arg477、Phe365、Phe449、Ser35、Tyr…，按 `extract_c1_crystal_refs.py` 实测确定，不预设）；
2. 定义 URAT1 门 = 酸根–Arg477 距离 + 口袋重原子重叠 + 关键残基 IFP Jaccard + 无冲突，四项与 NLRP3 侧同构；
3. **在同一回顾集（228 羧酸 active vs 64 羧酸 true decoy）上标定**，复用 `run_acid_gate_benchmark.py`，不重新对接。

**判据（预登记）**：新 URAT1 IFP 门 OR 的 95% CI 下界 > 1。
- 满足 → 成为主文 URAT1 臂门控，与 NLRP3 臂并列，两臂各有一条 OR 统计。这是全文相对同行的核心增量。
- 不满足 → **回落到 A1 规则（OR 3.18）作判别门 + A2 作可及性门**，二者取交集定义 tier-1。不得因为想要新指标而放弃已成立的 A1。

### W3 tier 分层 —— 不需要新对接，只需重打分
现有 SDF 已含 9 个模式，A1 与 A2 只是姿态选择规则不同。

- **tier-1**：A1 判别门通过（CNNscore 首选姿即几何兼容）∩ A2 可及性通过 ∩ ≥2/3 种子 dual-structural ∩ 化学审计通过；
- **tier-2**：仅 A2 + dual-structural + 化学通过（即现有 40 eligible 的余部）；
- **tier-3 / 不入表**：仅宽松门。

A1 临床双靶 keep 为 24（seed42 冻结），A2 ≥2/3 dual 为 59，因此 tier-1 必然是 24 的子集，规模小、门槛硬——正好当 primary。**具体数目由脚本算出，不预设。**

### W4 NLRP3 面板重建（含 M1/M2 修复）
阳性去重后 9 个（NP3-146、MCC950 及 ChEMBL 磺酰脲活性物），新增 ≥40 个性质匹配诱饵；三种子；重标定结构门阈值；报告 OR / CI / Fisher。
可选加强：若能取到第二个 NACHT 抑制剂共晶结构，做 2×2 交叉对接，检验 IFP 门是否跨结构稳健。

### W5 MD（分层，两个对照永远在）
同行标配 100–300 ns + MM-GBSA；本文要做得比"RMSD 平了就说稳定"更硬。

**体系**：
- 对照 A：lesinurad 晶体羧酸根姿 @ 9DKB（POPC 膜，转运体必须建膜）
- 对照 B：NP3-146 共晶姿 @ 7ALV（可溶蛋白）
- Discovery：tier-1 中 2–3 个 × 两靶

**参数**：3 个独立复制 × 200 ns（膜体系需额外平衡），力场按现有 pipeline。

**读出必须与门控同构**（不只是 RMSD）：
- URAT1：酸根–Arg477 盐桥占据率（%帧）、关键残基接触持续率、配体 RMSD
- NLRP3：关键残基 IFP 持续率、口袋重叠时序、配体 RMSD
- MM-GBSA 仅作**同池相对比较**，明确写不作亲和力

**停规则（不变）**：对照 A 失去 Arg477/Phe 笼 → URAT1 侧一律不解释；对照 B 漂出 → NLRP3 侧一律不解释。MD 合格 ≠ 双靶成立。
`md_authorized` 在 W1–W4 完成并冻结短名单后方可置 true。

### W6（可选，推荐但不改题）OAT1 反筛作为**排除性**注释
只做一层：把 tier-1/tier-2 与参照药对接进人 OAT1 结构（Jeon *Structure* 2025 apo/olmesartan/probenecid；Wu & Luo *Sci. Adv.* 2025 cidofovir），报告"预测的分泌型 OAT 反靶风险"。

理由：URAT1 抑制若同时封堵分泌型 OAT，尿酸排泄收益会被抵消——这是**必要滤器**，不是新主题。校准参考 Taniguchi *JPET* 2019 同测定矩阵（dotinurad URAT1/OAT1 ≈ 110×；benzbromarone 16.5×；lesinurad 0.23×、probenecid 0.066× 方向反转）。

**边界**：NLRP3 仍是主线第二节点；OAT1 只出现在一张表 + 讨论一段，标题摘要不出现 selectivity。**不做**选择性排序或 SURI 分类主张。

---

## 4 明确不需要补的（逐条否掉）

| 不做 | 理由 |
|---|---|
| 重开 Rank 轨 / 9DKB 全诱饵对接（L3 full decoy） | 排序不是产品；`allow_L3_full_decoy: false` 已锁；lesinurad 自由自对接三种子失败 |
| 新建 RandomDecoy 或扩充诱饵以"提高 EF" | 用于活性排序的富集不是本文主张；冻结数字不得美化 |
| URAT1 回归模型重训 / OAT 迁移学习 | 已评估：命名药回收失败，主排序不依赖它；再训不改变主张 |
| FEP / RBFE / 绝对结合自由能 | 无同系列、无实验锚点，n 太小，不可辩护 |
| C2 转运周期 Δ 打分（outward/occluded 差值） | Costanzi & Vilar *JCC* 2011、HOLIgraph *J. Cheminform.* 2025 已做过双态对接，作为"方法创新"卖不成立；且外向态模型不可靠 |
| 生成式设计（REINVENT/GFlowNet 等） | 会产出无化学来源的新分子，与"临床阶段重定位"故事冲突，且无合成/验证路径 |
| 网络药理学 + 草药提取物 | 期刊 scope 明确排除草药/多草药制剂与无化学成分的计算生物学 |
| Conformal selection / 保形筛选 | n ≈ 8–40，统计功效不足，已降级 |
| DFT / 量化优化 | 本文无具体电子结构问题，纯装饰 |
| 追 Pareto 前沿四个大环（Idremcinal 等） | 已证为高分子量偏倚产物，仅作 SI 偏倚证据 |
| 用诱饵标签训 ML 再选协议 | 方法学错误，已在项目内禁止 |
| 任何湿实验替代性声明 | 无实验即写 putative；不得写 validated / confirmed |

---

## 5 文章骨架

**图**
- F1 双臂分层工作流（含两个自对接对照位置）
- F2 URAT1 四联体交叉对接保真度矩阵（W1）
- F3 两臂门控回顾校准（W2 URAT1 OR + W4 NLRP3 OR，并列 forest）
- F4 漏斗 + 多种子稳定性（8,319 → 1,588 → 303 → 156 → 59/54 → tier-1）
- F5 候选姿态图：每个候选 URAT1 酸锚 + NLRP3 IFP 双面板
- F6 MD：盐桥占据率 + IFP 持续率（含两对照）

**表**
- T1 候选表（tier、Arg477 距离、overlap、IFP、接触、q_N、MW/QED、临床阶段、原适应症）
- T2 两臂门控校准统计（n、OR、CI、Fisher p、sens/spec）
- T3 化学组成（酸等效基团分布、Murcko 骨架、ADMET 摘要）

**讨论必须自己先说的三条**（同行普遍不说，说了反而加分）
1. 对接分数不解释为亲和力；A2 门 OR≈0.97 已明确其为可及性门；
2. 姿态保真度依赖配体柔性（W1 直接证据），故柔性酸用锚定几何门；
3. 无湿实验，产物是可测假说；URAT1 摄取抑制的 IC50 与单一姿态占据不是同一个量。

**局限**必须写：单一 NLRP3 结构、无膜环境下的对接、NLRP3 缩库模型来自细胞层数据而非直接结合、tier-1 规模小。

---

## 6 执行顺序与止损

```
W1 四联体交叉对接 ──┐
W4 面板重建(M1/M2) ─┤→ 门控冻结 → W2 标定 → W3 分层 → 冻结短名单 → W5 MD → W6 OAT1 注释 → 重写正文
M3 短名单修正 ──────┘
```

止损点：
- W1 刚性配体也失败 → 停，查受体准备，不进 W2；
- W2 新门 OR CI 含 1 → 回落 A1∩A2，不新造指标；
- W4 结构门无法提升特异性 → 如实降级为姿态质控；
- W5 任一对照失败 → 该靶侧不解释，MD 不进正文主张。

任何一步失败都**不改变产品形态**：候选表仍是主结果，只是支撑层数减少。
