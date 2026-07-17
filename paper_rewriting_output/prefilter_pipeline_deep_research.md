# 对接前初筛漏斗深度调研（Step0–7）

**日期：** 2026-07-17  
**问题：** 527k 胺 → 丙烯酰胺 → ECFP 四锚点 → Murcko/hinge 双轨 → 5–20k 对接就绪，是否有创新？如何评判、量化不确定性、用何算法升级？  
**检索方式：** PaperSpine Research 规范（并行 SOTA/方法调研 + WebSearch；无 Scholar/PubMed MCP 时用 Web 回退）  
**与已确认贡献对齐：** `confirmed_contribution.md` 的 Framework（必达）+ Chemistry（升级）

---

## 0. 总判决（先读）

| 主张 | 判决 |
|------|------|
| Step0–7 **单独**作为方法学创新 | **弱 / 基本不成立** |
| Step0–7 作为 **JNK2/8ELC 可复现工程 SOP** | **合理，可写 SI/方法** |
| Step0–7 + **校准指标 + 共价专属过滤 + 湿实验闭环** | **中–强，可支撑贡献 B/C** |
| 把 ECFP 双轨写成「近三年非 OA 标准流程」 | **过度声称**（标准的是模块，不是这套命名与阈值） |

**一句话：** 你们现在的指令是**合格 triage 工程**；创新不在「又写了一套 RDKit 流水线」，而在是否把 **反应性 / 选择性 / scaffold-hop 表征 / 不确定性** 接到 **8ELC–Cys116 校准框架** 上，并拿出可证伪证据。

---

## 1. 文献地图（代表性，非穷尽）

### 1.1 共价 VS / 库 / AF3（方法底座）

| ID | 文献 | 出处 | 与本漏斗关系 |
|----|------|------|--------------|
| L1 | London et al., DOCKovalent large-library covalent docking | Nat Chem Biol / PMC4232467 | 大库共价对接范式；强调 warhead 反应性不进打分 |
| L2 | Zhu et al., CovDock-VS | JCIM 2014, 10.1021/ci500175r | 共价对接 VS + decoy 富集；对接后仍需相互作用过滤 |
| L3 | Shamir/London, Discovery of Covalent Ligands with AF3 | JACS 2025/26, 10.1021/jacs.5c22222 | mPAE 主排；胺/丙烯酰胺大库前瞻；**对接前 triage 须服务 AF3 配额** |
| L4 | CovalentInDB 2.0 | NAR 2024, 10.1093/nar/gkae946 | 2M 级共价可筛库 + 与已知共价药 Tc；外部对照 |
| L5 | Kin-Cov integrated covalent kinase workflow | JMC 2023, 10.1021/acs.jmedchem.3c00088 | 残基可及性/反应性/结构验证一体；比纯 ligand Tc 更「共价」 |
| L6 | Reactive Docking | JCIM 2023, 10.1021/acs.jcim.3c00832 | 反应几何进入 ranking |
| L7 | HASTEN ML-boosted giga-scale docking | JCIM 2023, 10.1021/acs.jcim.3c01239 | 大库用 ML 代理减对接量——类比你们用 ligand 预筛减 AF3 |
| L8 | Enamine acrylamide / REAL acrylamides | 供应商库 | 商业丙烯酰胺 vs ChEMBL 胺枚举的双源对照 |

### 1.2 JNK 共价化学（锚点与目标边界）

| ID | 文献 | 关系 |
|----|------|------|
| J1 | Lu et al., YL5084 / 8ELC | JMC 2023, 10.1021/acs.jmedchem.2c01834 — 选择性锚点 + 共晶模板 |
| J2 | Wydra et al., 56d ligand-first | JMC 2025, 10.1021/acs.jmedchem.5c00884 — 正交骨架 |
| J3 | Zhang et al., JNK-IN-8 | Chem Biol 2012 — pan 校准，非选择性终点 |
| J4 | Nat Commun 2024 reversible/cyclic warheads for JNK Cys116 | 10.1038/s41467-024-52573-2 — **说明「只枚举标准丙烯酰胺」已不是前沿唯一路径** |
| J5 | Axcelead MCS2024 covalent JNK design talk | 从可逆骨架装 warhead + docking 优先 — MedChem 主流叙事 |

### 1.3 Scaffold hop / 表征 / 不确定性

| ID | 文献 | 关系 |
|----|------|------|
| H1 | ErG | JCIM 2006, 10.1021/ci050457y — scaffold hop 经典 2D 药效团图 |
| H2 | SwissSimilarity 2021 (ECFP/ErG/pharm2D) | IJMS 2022, 10.3390/ijms23020811 |
| H3 | SHAFTS | JCIM 2011 + JMC prospective RSK2 |
| H4 | HybridSim-VS | Bioinformatics 2017, 10.1093/bioinformatics/btx418 |
| H5 | GATNN-VS neural scaffold hop | JCIM 2020, 10.1021/acs.jcim.0c00622 |
| H6 | JCIM 2024 scaffold-hopped ID by LBVS | 10.1021/acs.jcim.4c00342 — ECFP 局限 |
| H7 | Enrichment confidence bands | J Cheminform 2022, 10.1186/s13321-022-00629-0 |
| H8 | UMAP clustering split for VS eval | J Cheminform 2025 — Murcko split 仍可能高估 |
| H9 | MolPAL active learning | Chem Sci 2021, 10.1039/D0SC06805E |
| H10 | KDBNet calibrated geometric DL for kinase | Nat Mach Intell 2023 |
| H11 | BIreactive / GSH reactivity ML | 反应性预测用于共价库 |
| H12 | PAINS utility & limits | ACS Chem Biol 2018 — 勿黑箱硬杀 |

---

## 2. Step0–7 逐模块：标准性 / 薄弱点 / 创新潜力

| Step | 内容 | 在非OA/领域中的地位 | 薄弱点 | 能否单点创新 |
|------|------|---------------------|--------|--------------|
| 0 | 体检、流式 | 工程标配 | 无 | 否 |
| 1 | 胺位数/MW 预筛 | 标配 | 芳香胺规则粗 | 否 |
| 2 | 胺→丙烯酰胺 | **非常常见**（REAL/共价 VS） | 忽略反应性窗、第二 warhead 仅作变体 | 否（组合可有） |
| 3 | MW/clogP/PAINS/SA | 标配 | PAINS 硬过滤有争议；缺 GSH/过反应 | 否 |
| 4 | ECFP4→四锚点 | **默认基线** | warhead 刷分；4 锚过窄；无 decoy 校准 | 否 |
| 5 | Murcko+hinge 双轨 | 思路常见、正式命名少 | Murcko≠真正 hop；hinge SMARTS 脆弱；阈值任意 | **弱**（需校准才升格） |
| 6 | Top 5–20k | 资源预算 | 无 Recall@k / EF 支撑 | 否 |
| 7 | 验收漏斗表 | 好工程 | 不等于科学验证 | 否 |

**与近三年非 OA 真实叙事的差距：**  
- JMC（Lu/Wydra）：**MedChem SAR + kinact/KI**，计算是辅助。  
- JACS COValid：**AF3 主排 + 前瞻实验**，ligand ECFP 双轨不是其卖点。  
- 工业大库：**Enamine 已有丙烯酰胺 / REAL**，或 ML-boosted docking（HASTEN），不是「先喊活性初筛再 Tc」。

---

## 3. 合理性与可解释性：评审级检查表

### 3.1 合理（保留）

1. 527k 不能全对接 → 必须 triage。  
2. 单丙烯酰胺 + 物化窗对齐已知共价 lead。  
3. Sim / Novel **分文件分排序**（避免一个榜淹没新骨架）——与你们 Phase0「方法不可互换」一致。  
4. pan（JNK-IN-8）单独成仓 —— 防止选择性目标被 pan 化学劫持。  
5. 全程可追溯列（tc、anchor、scaffold、discard_reason）——可解释性优于黑盒。

### 3.2 需改口 / 需补证据

1. 「活性初筛」应改称 **ligand-based enrichment / chemotype triage**。  
2. Tc 0.25–0.55 / <0.25 在未网格校准前是 **hyperparameter**，不是科学常数。  
3. 装弹头后再算 ECFP，丙烯酰胺子结构会抬高与锚点 Tc → 建议并行 **warhead-stripped / core Tc**。  
4. hinge SMARTS ≠ hinge binding；Novel 仓假阳性风险高。  
5. 缺 **JNK1 负筛 / 选择性维度** —— 与已确认贡献「JNK2>JNK1」不完全对齐。  
6. 缺 **反应性窗** —— Nat Commun 2024 已把 warhead 精度推到可逆/环状；只枚举标准丙烯酰胺偏保守（可接受为 Phase1，不可吹成前沿）。

---

## 4. 不确定性：可量化实验清单（对接前就能做）

### 4.1 必做（低成本，WSL+RDKit）

| 实验 | 指标 | 用途 |
|------|------|------|
| Tc 阈值网格 | 仓规模、Murcko 数、两仓 Jaccard | 证明阈值落在稳定平台 |
| Leave-one-anchor-out | 回收率、Top-k overlap | 检验是否被单一锚点绑架 |
| Warhead-stripped vs full Tc | Spearman、分仓重分配比例 | 量化 warhead 刷分 |
| Random / 无关激酶共价锚点对照 | 假「Sim」比例 | 基线 enrichment |
| Hinge 消融 | Novel 规模与噪声样例 | 证明 SMARTS 有无贡献 |

### 4.2 应用回顾性富集（中成本）

构建：JNK 共价阳性（文献+ChEMBL）+ property-matched decoy（或 CovalentInDB 非 JNK 丙烯酰胺）。  

指标：**EF@0.5/1/5%、BEDROC、PR-AUC、Recall@10k**；Sim/Novel **分开报**。  
不确定性：enrichment 曲线置信带（J Cheminform 2022）。  

### 4.3 校准与适用域（若上 ML 分数）

Reliability curve、ECE、Brier；conformal coverage；每个分子输出 AD 旗标（max Tc to actives、cluster 距离）。

### 4.4 与下游一致性（对接/AF3 后）

Sim 高分 ∩ mPAE 差、Novel 低 Tc ∩ mPAE 好 —— 量化初筛与 Phase0 主门控的分歧（延续你们 Glide↔AF3 弱相关教训）。

---

## 5. 现有算法：在 Step0–7 上「加一层」而非推倒

按 **WSL 可实现性 × 与贡献对齐** 排序：

### P0 — 最小可发表补丁（强烈建议）

1. **Core/warhead-stripped Tc** + product Tc 双列。  
2. **阈值网格 + leave-one-out + 漏斗敏感性表** 写入方法。  
3. **反应性/过反应 SMARTS 分桶**（双亲电、过活 Michael、已知 alert）；PAINS 改为 flag 非硬删。  
4. **JNK1 轻量负筛计划**（哪怕先用同源口袋描述符/粗对接 Δscore，后置也可）。

### P1 — Scaffold-hop 补强（Track-Novel 真要新）

1. **ErG + Gobbi 2D pharmacophore**（RDKit）作 Novel 召回。  
2. **8ELC/YL2056 三维药效团 + shape**（USRCAT / ROSHAMBO / 商业 ROCS·SHAFTS）只对 Novel 子集。  
3. Novel 入选改为：**低 ECFP + 高 ErG/pharm + Murcko 新 + Cys 矢量几何过线**。

### P2 — 算力友好的大库策略

1. **MolPAL / active learning**：用 Glide 或 AF3 作昂贵 oracle，surrogate 用 ECFP+ErG+shape。  
2. 对照轨道：**商业 Enamine acrylamide** vs **ChEMBL 胺枚举**（双源，已在漏斗文档中）。

### P3 — 谨慎当作主创新

生成模型造共价分子、端到端 DL 活性回归（JNK 标签少、易泄漏）、宣称「新相似度算法」。

---

## 6. 与 `confirmed_contribution` 的映射（避免跑偏）

| 贡献层 | Step0–7 能支撑什么 | 不能支撑什么 |
|--------|-------------------|--------------|
| **Framework 必达** | 可复现建库→分仓→缩库 SOP；与 8ELC/Cys116 规则衔接；校准附录 | 「发现了活性分子」 |
| **Chemistry 升级** | 交付 Sim/Novel 可采购或可合成名单 | 无 kinact/KI 的效价宣称 |
| **Analysis/benchmark** | JNK2 位点 prefilter enrichment + 不确定性 | 通用共价 VS 新理论 |

**审稿人最可能的两刀：**  
1. “这不就是 ECFP 近邻搜索吗？” → 用 **校准表 + ErG/3Dpharm 消融 + Phase0 衔接** 回答。  
2. “没有新分子/没有选择性数据” → 明确 Framework vs Chemistry 分期；疾病不作主 claim。

---

## 7. 可实现创新包（Minimal Viable Innovation）

在不推翻现有 agent 指令前提下，把 Step4–5 升级为：

```text
特征：
  tc_full_* , tc_core_* , erg_* , pharm2d_*
  murcko, hinge_hits, reactivity_bucket, sa, purchasable
分仓：
  Sim  : tc_core 网格标定窗 + (YL/56d 优先; JNK-IN-8→pan)
  Novel: tc_core 低 + Murcko 新 + (ErG/pharm 高 或 8ELC-3Dpharm 过线)
  硬排除: 过反应桶、近重复 tc_full>0.7
交付：
  阈值敏感性 + LOO + EF/Recall@10k 附录
  06_dock_ready_{sim,novel,pan}.csv
下游：
  Glide 松阈值 → AF3 mPAE 门控 → JNK1 counterscreen 设计
```

**可写进引言的诚实卖点：**  
不是发明 ECFP，而是在 **JNK2 Cys116 / 8ELC** 上给出 **经校准、分轨、可量化不确定、并与 AF3 门控衔接的丙烯酰胺枚举 triage**，用于推进区别于 YL5084/56d 的可检验候选。

---

## 8. 打开思路的提问清单（给自己用）

1. 若去掉全部 Tc，只留 8ELC 药效团 + 反应性，还能不能富集已知阳性？  
2. 56d 在「仅 YL 锚点」时是否被误扔？（正交骨架压力测试）  
3. Novel Top100 目视：有多少根本不像激酶配体？  
4. 与 Enamine 在库丙烯酰胺重叠多少？重叠高则「枚举创新」弱、应强调选择性/几何校准。  
5. 首轮 96 化合物分层抽样湿实验：高 Sim / 高 Novel / 边界 / 随机 —— 这才是 enrichment 的金标准。

---

## 9. 检索局限（诚实披露）

- 本环境 **无** CNKI/PubMed/Semantic Scholar MCP；以 WebSearch + 全文抓取 + 并行研究子代理为主。  
- 未声称读完所有共价 VS 论文；上表为 **与本漏斗决策相关的高权重样本**。  
- 「ARS Codex」：工作区仅安装 **PaperSpine**（含 Codex 宿主约束）；本次按 PaperSpine Research（SOTA map + 并行专题）执行，而非独立名为 ARS 的 skill。
