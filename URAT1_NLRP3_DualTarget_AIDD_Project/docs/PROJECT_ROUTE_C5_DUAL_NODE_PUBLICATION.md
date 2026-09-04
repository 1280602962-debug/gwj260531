# C5 路线：URAT1–NLRP3 双节点候选论文（目标 *Molecular Diversity*，无湿实验）

> 取代 C2 / C3 / C4 作为**唯一现役主线**。C2（转运周期 Δ 打分）、C3（结构时间窗口叙事）、C4（选择性改题）各自的可用零件被吸收进本路线，不再作为独立主线。
> 冻结资产不覆盖：`data/repurposing/p2/`、`data/campaigns/c1/07_clinical_dock/acid_dual_a1_frozen/`、`data/si/`。
> 门控在看到任何新分子名之前锁定于 `config/campaign_c5.yaml`。

> **2026-09-04 核实修订。** 本次修订用实际检索/计算核对了上一版的关键假设，发现并修正了六处问题（详见 §7）。核实产物：
> `data/campaigns/c5/00_verification/w1_reference_ligand_verification.json`（RCSB 结构/配体/可旋转键核实）、
> `data/campaigns/c5/00_verification/m1_m2_nlrp3_panel_verification.json`（面板去重后 Fisher p 重算）、
> `data/campaigns/c5/03_tiering/{tier1_candidates.csv,tier2_candidates.csv,tier_summary.json}`（W3 分层已用现有脚本 `scripts/build_c5_tier_assignment.py` 真实算出，不再是占位符）。

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
| 候选逐分子结构指标 | overlap / IFP Jaccard / 关键接触 / Arg477 距离 | `03_tiering/tier1_candidates.csv`（**取代**旧 `08_nomination/acid_shortlist_a2_competition.csv`，见 §7） | 候选表 T1 |
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

## 2 投稿前必须修掉的三处（不修会被审出来）——**均已核实并部分修复，见下**

### M1 已知配体面板有重复且缺阴性对照 —— **已修复并重算**
核实：`nlrp3_structural_panel/panel_ligands.csv` 中 `CHEMBL3183703` 与 `MCC950` 的 SMILES **逐字符相同**，确系重复。去重后阳性 n=10→9，用与冻结脚本相同的 `scipy.stats.fisher_exact` 重算：

| | 阳性 n | 背景 n | 通过/阳性 | 通过/背景 | Fisher p（两侧） |
|---|---|---|---|---|---|
| 冻结（含重复） | 10 | 20 | 10/10 | 11/20 | 0.01340 |
| **M1 修复后** | **9** | 20 | 9/9 | 11/20 | **0.02703** |

去重后仍 <0.05，但显著性**弱了约一倍**，不是"无关紧要的小修"。产物：`data/campaigns/c5/00_verification/m1_m2_nlrp3_panel_verification.json`。
背景仍是 20 个随机临床酸，不是性质匹配诱饵——**这部分尚未修复**，仍需补 ≥40 个诱饵（W4）。

### M2 A2b 结构门在本面板上没有增加特异性 —— **已逐行核实为真**
逐行核对 `nlrp3_panel_metrics_seed42.csv` 全部 30 行：`keep_nlrp3_pose` 与 `keep_nlrp3_structural` **完全相同**，无一例外。即结构门（overlap + IFP + 关键接触）在当前背景上贡献的判别力恰好为零——不是约等于零，是逐行相等。这说明当前 20 个背景分子只分成"完全不进口袋"和"完全进口袋且姿态合理"两类，没有"进了口袋但姿态错"的中间情形来考验结构门。
**修法不变**：W4 扩充诱饵后重新标定；若仍无法提升特异性，如实降级为"姿态质控"。

### M3 primary 候选内部矛盾 —— **根因比原判断更严重，已用机械脚本修复**
原判断以为只是"手写表里错标了一个分子"。核实后发现：`08_nomination/acid_shortlist_a2_competition.csv` 的 primary/backup 分层**根本不是门控算出来的**——生成脚本 `build_c1_acid_shortlist_a2.py` 里硬编码了 `PRIMARY_TIER1 = ["PF-04620110"]`、`PRIMARY_TIER2 = ["ADMILPARANT", "RUNCACIGUAT", "LANIFIBRANOR"]`、`BACKUP_TIER = ["PSI-697", "PF-03882845"]` 这样的**人工名单**，门控数据只是拿来给这些预选名字打分，PF-04620110 从未真正被结构门选中过。

**修法**：弃用该硬编码脚本的分层逻辑，改用纯机械交集（见 §3 W3、已跑通并产出真实数字）。

### M4（本次核实新发现）化学软排除是名单而非结构筛
`build_c1_acid_shortlist_a2.py` 的 `SOFT_EXCLUDE_SUBSTR` 只硬编码了三个头孢药名，不是 β-lactam 环的 SMARTS。用 SMARTS 一查，tier-2 的 24 个候选里立刻多出 3 个头孢类（Cefetrizole、Cefazedone、Cefoxazole）。
**修法**：把名单排除换成/追加 β-lactam 环 SMARTS（`[#6]1[#6][#7][#6]1=O`）等结构告警，`scripts/build_c5_tier_assignment.py` 已内置该检查并标记 `beta_lactam_flag` 列，供候选表最终过滤用。

---

## 3 需要往下推进的计算（W1–W5，本地对接 + MD 可用）

设计原则：**新算力全部投在"让两臂对称、各自可校准"上**，不投在提高排序名次上。

### W1 URAT1 四联体交叉对接 —— 把自对接失败变成方法学结果
Suo/Fedor/Lee *Nat. Commun.* 2025, **16**:5178 同构建体四结构（**已用 RCSB API 核实**：DOI 10.1038/s41467-025-60480-3、PubMed 40467597、四个分辨率与文档一致）：

| PDB | 状态 | 分辨率 | 配体 CCD | 配体 |
|---|---|---|---|---|
| 9DK9 | apo | 2.68 Å | — | — |
| 9DKA | holo | 3.00 Å | R75 | benzbromarone |
| 9DKB | holo | 2.74 Å | **A1AIL**（原文档误写 LES，已订正） | lesinurad |
| 9DKC | holo | 2.55 Å | A1A45 | TD-3 |

> **环境澄清（已核实，纠正上一轮判断）**：这台云沙箱其实**有网络**，四个结构都能直接 `curl https://files.rcsb.org/download/*.cif` 抓到（9DKC 的旧版 `.pdb` 格式返回 404，必须用 `.cif`）；`rdkit`/`pandas`/`scipy` 也能 `pip install --user` 装上。**结构下载、配体 SMILES 提取、可旋转键统计、受体准备前处理，这一层现在就能在这里做**，不必等本机。唯一真正卡在"你本机"的，是 gnina 对接本身（这台机器没有二进制、没有 GPU 假设）。核实产物：`data/campaigns/c5/00_verification/w1_reference_ligand_verification.json`。

**做法**：3 个晶体配体 × 4 个受体（含 apo）× 3 种子，`run_gnina_batch.py` + `config/docking_c1.yaml` 同参数；报告 Top-1 / Top-3 / best-of-9 RMSD 与酸根–Arg477 距离矩阵。

**判据需要修正（原判断的"刚性配体"假设部分错误，已用 RDKit 核实）**：

| 配体 | 可旋转键数 | 备注 |
|---|---:|---|
| lesinurad | 5 | 已知失败案例：硫醚–CH₂–COOH 摆臂 |
| benzbromarone | **3** | 酮连接，无摆臂——**唯一真正更刚性的参照** |
| TD-3 | **5** | 硫醚–C(CH₃)₂–COOH 摆臂，**与 lesinurad 同一类柔性负担，可旋转键数完全相同** |

原方案把 TD-3 当"刚性对照"是错的：TD-3 和 lesinurad 共享同一个硫醚–羧酸摆臂化学型，只是端基从 –CH₂– 换成 –C(CH₃)₂–。**修正后的判据**：
- **受体准备是否正常**只看 benzbromarone：自对接 Top-1 RMSD ≤ 2.0 Å 通过即可，这是本三角里唯一独立于柔性假说的检验。
- **TD-3 的角色改为"柔性假说的同类复现"，不是对照**：若 TD-3 自对接也失败（Top-1 RMSD 明显偏高），这**支持**而不是推翻"柔性硫醚–羧酸摆臂导致姿态歧义"的解释，反而加强用几何优先门（W2）的必要性；若 TD-3 意外通过，说明柔性不是唯一因素，需要重新审视 lesinurad 特有的什么（例如溴苯环的取向）。
- 若 benzbromarone 也失败 → 受体准备/盒子有系统问题，回头查准备流程，不得直接进 W2。

**残基编号提醒**：`arg477_coords.json` 已记录"Prepared 9DKB PDBQT renumbers this residue as ARG A 476; same guanidinium"——即准备后的受体残基编号相对文献编号偏移 1，W2 建关键残基图时**每个残基都要单独核对编号**，不能整批套用文献序号。

### W2 URAT1 酸锚 IFP 门 —— 让 URAT1 臂与 NLRP3 臂对称
现状不对称：NLRP3 有 IFP 门（overlap + Jaccard + 关键接触），URAT1 只有**单一 Arg477 距离**。单距离门就是 A2 那个 OR≈0.97 的东西。

**做法（关键残基来源已更正为文献锚定，不是"实测确定"）**：
1. **关键残基不是靠几何反推，而是已有定量突变体文献**（`docs/PROJECT_ROUTE_C3_CRYOEM_DUAL_NODE.md` 已整理，直接复用）：
   Tan 等 *Sci. Rep.* 2017, 7:665（verinurad/苯溴马隆/磺吡酮/丙磺舒的 fold-change）+ Guo/Chen *Nat. Commun.* 2025, 16:1512（S35Q、R477N、D389A、F360T 等）+ Dai & Lee *Cell Res.* 2024（五 Phe 笼 F241/F360/F364/F365/F449）给出残基集：
   `S35, M214, F241, F360, F364, F365, D389, K393, Q437, F449, R477`（+ Q473）。
   这与 NLRP3 侧关键残基集（Dekker 2021 文献锚定 Ala227/228、Arg351、Met408、Tyr443、Phe575、Arg578）是**同一建法**，不是两套逻辑。
2. **编号必须逐个核对**（见 W1 残基编号提醒），不能把文献序号直接套进准备后的受体 PDBQT。
3. 定义 URAT1 门 = 酸根–Arg477 距离 + 口袋重原子重叠 + 关键残基 IFP Jaccard + 无冲突，四项与 NLRP3 侧同构。

**阈值标定方法需要修正（原方案有数据窥视风险）**：原计划"网格搜索让 OR 的 CI 下界最大"，这等于在同一 228 vs 64 回顾集上先调参再报告显著性，属于用同一份数据既选阈值又检验阈值——这正是项目自己在别处明令禁止的"用诱饵标签训 ML 再选协议"的同构错误。

**修正后的标定方法（锚定法，不调参）**：仿照 NLRP3 门"锚定在自对接表现之下"的建法（overlap≈1.0、IFP≈0.84–1.0、关键接触 6–7/7，是从 NP3-146 自对接直接读出来的，不是在背景集上调出来的）：URAT1 侧四个阈值同样从**晶体自对接的实测值**（benzbromarone、lesinurad、TD-3 三者的 overlap/IFP/接触数）取一个保守下界（如最低值或略低于最低值），锚定后**只评估一次** 228 vs 64 回顾集，不回头调整。

**判据（预登记，不变）**：新 URAT1 IFP 门 OR 的 95% CI 下界 > 1。
- 满足 → 成为主文 URAT1 臂门控，与 NLRP3 臂并列，两臂各有一条 OR 统计。这是全文相对同行的核心增量。
- 不满足 → **回落到 A1 规则（OR 3.18）作判别门 + A2 作可及性门**，二者取交集定义 tier-1。不得因为想要新指标而放弃已成立的 A1。

### W3 tier 分层 —— **已用现有冻结文件机械算出，不再是占位符**
现有 SDF 已含 9 个模式，A1 与 A2 只是姿态选择规则不同。零新对接，纯重打分/重连接，已用 `scripts/build_c5_tier_assignment.py` 跑通：

- **tier-1** = A1(seed42) 判别门通过（24）∩ A2 ≥2/3 种子 dual-structural（54）∩ 化学审计通过（40）→ **n = 13**；
- **tier-2** = A2 ≥2/3 dual-structural ∩ 化学审计通过，去掉 tier-1 → **n = 24**；
- 产物：`data/campaigns/c5/03_tiering/{tier1_candidates.csv, tier2_candidates.csv, tier_summary.json}`。

**tier-1 的 13 个中**：
- GSK-3008348 FREE BASE 按 M3 结论标记为"结构对照，非候选"，真正候选 **12** 个；
- 只有 5 个（Lanifibranor、Admilparant、PF-03882845、PSI-697、GSK-3008348）与旧的手写 `acid_shortlist_a2_competition.csv` 重叠；PF-04620110 不在其中（与 M3 结论一致，机械重算自动排除了它，不需要手动"降级"）；
- 另有 8 个新出现的候选未曾在旧表中出现：Tonapofylline、Caficrestat、Lintitript、Spiroglumide、Cavosonstat、Runcaciguat、Fulimetibant、CR-3465 Free Acid（其中 Runcaciguat 此前只出现在硬编码脚本的 `PRIMARY_TIER2` 名单里，从未真正被门控选中，现在是**门控本身**选出来的，性质不同）。
- Tier-1 中头孢类结构告警数 = 0（用 β-lactam SMARTS 核对）。

**tier-2 的 24 个里发现新问题（M4，见 §7）**：3 个头孢类抗生素（Cefetrizole、Cefazedone、Cefoxazole）混在里面，因为项目现有的化学软排除名单（`SOFT_EXCLUDE_SUBSTR`）只硬编码了 3 个具体头孢药名（`CEFCANEL`、`CEFAZAFLUR`、`CEFOVECIN`），不是结构筛（β-lactam SMARTS），漏掉了名字不在列表里的同类结构。这三个已在 tier-2 输出里被 `beta_lactam_flag=True` 标出，若把 tier-2 拿去当 W5 的 backup MD 槽位来源，投稿前必须先排除它们。

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

---

## 7 本次核实修订记录（2026-09-04）

上一版方案的方向判断（产品形态、两臂对称、tier 分层逻辑、止损点）**站得住**，但五处具体细节此前是猜测或未经计算的占位符，本次逐一用实际检索/计算核实：

| 编号 | 上一版怎么说 | 核实后发现 | 处置 |
|---|---|---|---|
| E1 环境能力 | "9DK9/9DKA/9DKC 本仓库没有，你本机下载" | 这台云沙箱**有网络**，`curl files.rcsb.org` 直接可抓；`pip install rdkit/pandas/scipy` 可装。结构/配体/化学层面的准备工作现在就能做，只有 gnina 对接本身仍卡在本机 | 已下载三个结构核实分辨率/配体码，见 W1 |
| E2 9DKB 配体代码 | 文档/config 全部写 `LES` | 本地 `9DKB.cif` 与刚下载的新副本一致，实际 CCD 是 `A1AIL`；`LES` 是错的（几何数据本身没问题，只是引用代码错） | 本文档与 `config/campaign_c5.yaml` 均已订正 |
| E3 W1"刚性配体"假设 | benzbromarone、TD-3 都算刚性对照 | RDKit 核实可旋转键：benzbromarone=3，TD-3=**5（与 lesinurad 相同）**。TD-3 和 lesinurad 同属硫醚–羧酸摆臂化学型 | W1 判据改为：仅 benzbromarone 测受体准备；TD-3 失败反而支持柔性假说 |
| M1/M2 面板 | "去重、重标定"写成一句话待办 | 实际重算：去重后 Fisher p 从 0.0134 变 **0.0270**（弱了一倍，仍显著）；逐行核对确认结构门与宽松门 30/30 行完全相同 | 数字已写入正文；背景诱饵扩充仍是待办 |
| M3 短名单 | "PF-04620110 该降级" | 根因更严重：整张短名单来自脚本里硬编码的人工名单，不是门控产物 | 弃用硬编码脚本，改用机械交集重算，见 W3 |
| M4（新发现） | 未提及 | 化学软排除是三个头孢药名的名单，不是结构筛，tier-2 混进 3 个头孢类 | 加 β-lactam SMARTS 检查，已内置进新脚本 |

**核实产物清单**：
- `data/campaigns/c5/00_verification/w1_reference_ligand_verification.json`
- `data/campaigns/c5/00_verification/m1_m2_nlrp3_panel_verification.json`
- `data/campaigns/c5/03_tiering/tier1_candidates.csv`（n=13，含 1 个结构对照）
- `data/campaigns/c5/03_tiering/tier2_candidates.csv`（n=24，含 3 个待排除头孢类）
- `data/campaigns/c5/03_tiering/tier_summary.json`
- `scripts/build_c5_tier_assignment.py`（可重跑，无硬编码分子名）

**结论**：核实没有推翻方案的整体逻辑，但把三处"计划要做"的事提前变成了"已经做完的验证"（W1 结构可获得性、W3 tier 分层、M1 面板重算），同时发现并修了两处真实缺陷（LES 引用错误、旧短名单靠硬编码而非门控）。唯一没有改变的判断是：**gnina 对接本身仍必须在你本机完成**，这一点上一轮的结论是对的。
