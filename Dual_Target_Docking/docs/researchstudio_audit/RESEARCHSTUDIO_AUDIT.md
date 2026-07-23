# ResearchStudio 审计报告：双靶兼容性打分 / Dual-VSDS

> 依据 [microsoft/ResearchStudio](https://github.com/microsoft/ResearchStudio) 的 **Idea** 侧技能执行：  
> `idea-quality` · `scoop-check`（+ `paper-search`）· `idea_spark` 的 critique / falsification / implementability / anti-patterns。  
> **未**跑完整 IdeaSpark 五阶段「生成新 idea」（你已有定稿方向）；本报告是对**现有课题**的质量与可完成性审计。  
> 想法原文：[`idea_card.md`](idea_card.md) · 主方案：[`../NMI_DUAL_COMPATIBILITY_PLAN.md`](../NMI_DUAL_COMPATIBILITY_PLAN.md)

审计日期：2026-07-23 · 分支 `cursor/dual-target-docking-standalone-0b1a`

---

## 总判决（先看这里）

| 维度 | 结论 |
|------|------|
| **Idea-Quality** | **67 / 100 · Verdict: strong（贴边）** |
| **Scoop-Check** | **Level 3 — Medium Overlap**（可辩护，但 delta 必须写死） |
| **Critique（软判断）** | **revise**（非 abandon）：可证伪结构缺失 + 若干实现空洞 |
| **可完成性** | **可完成 JCIM / Chem. Sci. 级方法+基准文；冲击 NMI 条件可达成但非默认成功** |
| **最大致命漏洞** | **配对真双靶标签稀疏 + 校准/评测循环风险 + 细胞 holdout 被误读成 binding 金标准** |

一句话：思路方向正确、与 VSDS-VD/CleanSplit 同构，**不是被 scoop 死的 idea**；真正卡脖子的是**数据密度与评测因果设计**，不是再发明一个对接器。

---

## 1. Idea-Quality（`idea-quality` skill）

### Decomposition

- **Problem / gap：** 双靶筛选实务 = 两端独立对接 + raw mean/min/rank；跨靶分数不可比 + 药效不对称 → 系统性抬升 A-only/B-only。
- **Method：** Dual-VSDS 配对基准 + 每靶校准 \(p(\mathrm{active})\) + softmin/短板融合 + dual-vs-single 硬负 + fused/linked 条件 + 私有细胞时间盲测。
- **Why it should work：** 若失败主因是「不可比分数的朴素聚合」，则校准 + 短板敏感目标应在 dual-vs-single 上压过 mean/sum；若失败主因是对接本身无信息，则所有融合头应一起垮（可证伪）。

### Axis scores

| Axis | Score | Evidence（引自 idea） | Reason |
|------|-------|----------------------|--------|
| **A — Problem position** | **4** | “systematically promotes A-only or B-only”；“do not treat cross-target calibration… as a first-class ML problem” | 缺口真实且非显然；命名本身有信息量。略扣：领域内已有 Perez-Castillo 式融合案例，需证明「任务级」而非「又一次加权」。 |
| **B — Method quality** | **3** | 校准 + softmin + hard negatives + Dual-VSDS | **depth 3**：零件已知，深度靠任务重定义与基准，非新算子。**soundness 4**：短板逻辑自洽。**feasibility 3**：依赖稀缺配对标签与大量对接；若干阈值未锁。 |
| **C — Problem-fit** | **4** | Method steps 1–5 直接对准「朴素融合偏置」与「开放评测」 | 方法打在声称缺口上；私有细胞/PK 若喧宾夺主会伤 C（当前写法已自觉降级）。 |

**Overall: 67 / 100 · Verdict: strong**  
（公式：`round(100*(4+3+4-3)/12)=67`；A/C 均 >2，未触发 gate。）

**Strongest point：** 把双靶评估重定义为 **dual-vs-single 决策问题**，与「再训一个 docking scorer」刻意切割。  
**Most fixable weakness：** 补一条 **ResearchStudio 式 falsification kill-switch**（见 §4），并把 θ / 靶点对 / 校准划分写成作者决策而非散文。

---

## 2. Scoop-Check（`scoop-check` + `paper-search`）

### 2.1 Decomposed claim

- **Problem framing：** 同一配体对 A/B 的双靶虚拟筛选排序与评测（非三元 PROTAC，非联用协同）。
- **Core mechanism：** 每靶概率校准 + 短板敏感融合（softmin/阈值边距）+ dual-vs-single 硬负目标 + 泄漏控制配对基准。
- **Key insight：** 朴素 mean/sum/rank 在跨靶不可比与药效不对称下系统性偏袒单靶强分子；兼容性是校准后的短板约束，不是分数算术平均。
- **Application domain：** 小分子双靶 / polypharmacology structure-based VS（含 fused/linked）。

### 2.2 检索

Queries：

1. `dual-target docking score fusion virtual screening`
2. `polypharmacology docking consensus scoring benchmark`
3. `calibrated softmin dual inhibitor ranking hard negatives`
4. `dual target ligand docking virtual screening enrichment`（二次）

来源：arXiv / OpenAlex / Crossref（ResearchStudio `paper-search`）。原始输出见 `steps/paper_search_raw.json`。

另用模型知识补入检索未顶置但直接相关的锚点：VSDS-VD、EquiScore、CleanSplit、PoseBusters、POLYGON。

### 2.3 高威胁候选（深度对照，非全量 67 篇）

| 先验 | 轴重叠 | 判定 |
|------|--------|------|
| **Pérez-Castillo et al. 2017** — *Fusing Docking Scoring Functions… Dual Target Ligands* | framing ✓ · mechanism **partial**（rank mean，无校准/硬负/开放配对基准）· insight partial · domain ✓ | **最接近祖先**；3 轴级威胁若你写不清增量 → Level 2 |
| **VSDS-VD (Gu et al. NMI 2025)** | framing 评测美学 ✓ · mechanism **differ**（单靶）· insight 部分同构（设定可翻转）· domain docking VS ✓ | **叙事孪生、任务不同**；必须 cite 并声明 Dual-VSDS = 双靶配对推广 |
| **ML / Exponential consensus scoring (Ericksen 2017; Palacio-Rodríguez 2019)** | mechanism 融合 ✓ · framing 多为**单靶多打分/多构象** · domain VS ✓ | 零件先验；不是双靶兼容性任务 |
| **Ensemble docking fusion (Bajusz 2019)** | GEOM/HARM 等融合规则 | 单靶多结构；Related Work 必备 |
| **EquiScore / CleanSplit (NMI)** | 泄漏与打分纪律 | 方法学盟友，非 scoop |
| **POLYGON / CLM multi-target (Nat. Commun.)** | domain 双靶 ✓ · mechanism **生成** ≠ 判别排序 | 体裁不同 |

### 2.4 Comparison（最差先验封顶）

对 **Pérez-Castillo 2017**：约 **2 轴匹配**（framing + domain）→ **Level 3**。  
对 **VSDS-VD**：约 **2 轴匹配**（评测哲学 + domain）→ **Level 3**。  
**总体 Verdict：Level 3 — Medium Overlap**。

### 2.5 Delta（必须写进 cover letter / Intro）

> Unlike Pérez-Castillo et al. (2017), which fuse docking **ranks** to enrich dual ligands for one disease pair without cross-target probability calibration or A-only/B-only hard negatives, and unlike VSDS-VD, which audits **single-target** docking under True/RandomDecoy regimes, this work reframes dual assessment as a **calibrated shortfall / dual-vs-single ranking problem**, releases a **leakage-controlled paired Dual-VSDS** benchmark spanning multiple target pairs and design types, and shows that standard mean/sum fusion systematically elevates single-target–biased molecules.

若写不出这句的可测后半句（「shows that…」必须有表），则 novelty 塌缩为 Level 2。

---

## 3. Critique 五检（`critique.txt` 逻辑，软审计）

| Check | Verdict | 说明 |
|-------|---------|------|
| **gap_closure / Reject lessons**（对照 `controlled_diagnostic_design`） | **borderline** | 有「打假朴素融合」意图，但若只有单靶点对、无 trivial baseline、或用同一校准集自证，会撞上 *circular validation* / *single-task overclaim* / *missing trivial baselines*。 |
| **recipe_application** | **borderline→applied 若补齐** | 诊断型成功条件要求：**只动一根轴的对照**（mean vs softmin，同对接姿态）、**trivial 基线**、**外部锚**（binding 标签 ≠ 对接分）。目前方案有方向，实验设计表未钉死。 |
| **anti_pattern** | **watch：`heterogeneous_decomposition` + prior** | fused/linked 分治 = 异质分解；再叠物理先验/类型专家时，必须做「关分解 / 关先验」独立消融（见 ResearchStudio `anti-patterns.md`）。否则 reviewer 问「到底谁在干活」。 |
| **paper_pointed_threat** | Pérez-Castillo + VSDS-VD | **addressable**（见 Delta）；非 exact-mechanism → **不触发 hard-floor abandon**。 |
| **falsification_structure** | **deficient（现状）** | 主方案有阶段与指标，但缺 ResearchStudio 要求的完整 kill-switch 四件套（见下）。 |

**软判决：`revise`（非 abandon）** — 方向可进；投稿前必须修 falsification + 校准划分 + 多靶点对。

**parametric_family_concern（软信号）：**  
`consensus docking` · `rank fusion for polypharmacology` · `multi-objective docking score aggregation` — 投稿前再扫一轮 2024–2026 ChemRxiv/JCIM。

---

## 4. 可证伪条款（建议锁定的 kill-switch）

按 `falsification_structure_check` 补全：

**最小实验：** 同一 Dual-VSDS split、同一 GNINA 姿态；仅替换融合头：`raw mean` vs `z-mean` vs `calibrated softmin`。  
**结果指标（方向）：** dual-vs-single pairwise accuracy **上升**；Top-1% 中 A-only/B-only 占比 **下降**。  
**负载变量：** 短板敏感度 / softmin temperature \(\tau\)（或阈值边距 \(m_t=\hat p_t-\theta_t\)）。  
**负对照（打在结果指标上，非同义反复）：** 将 \(\tau\to\infty\)（softmin→均值）或强制 \(\theta_A=\theta_B=0\) 且改回未校准 raw mean——预期 **dual-vs-single 精度回到 mean 基线**，且 A-only 污染回升。若负对照后指标仍高，则增益来自对接/校准数据泄漏而非短板机制 → **主张失败**。

当前文档 **缺此条** → 按 ResearchStudio 标准属 **deficient**；补上后可 re-audit 为 sound。

---

## 5. Implementability 空洞（`implementability_audit`）

工程师冷启动时仍无法直接编码的 **open 作者决策**：

| ID | 空洞 | 风险 |
|----|------|------|
| O1 | **θ / pAct 阈值**未冻结 | 标签定义漂移，整表不可复现 |
| O2 | **公开靶点对清单**未冻结（2–3 对？） | 做完一对会被打 single-task overclaim |
| O3 | **校准集 vs 评测集划分**未写死 | **循环验证**：用评测标签校准再评同一批 |
| O4 | softmin \(\tau\) / 是否可学习 | 方法不可复现 |
| O5 | A-only/B-only **最小样本数**门槛 | 硬负太少则 dual-vs-single 指标不稳 |
| O6 | design_type 标注协议（谁标、双人？） | fused/linked 专家不可审计 |
| O7 | 细胞 L2/L3 与 binding 混用边界 | NMI reviewer 最爱打的「标签语义」 |
| O8 | 对接失败 / PB-fail 分子如何进入排序 | 选择偏倚 |

**可编码部分（已够清楚）：** 独立两端对接 → 特征表 → 融合脚本 → 两套 decoy 报表。  
**不可假装已清楚：** O1–O3、O5。不锁死则「可完成」只停留在 PPT。

---

## 6. 流程漏洞清单（按严重度）

### P0 — 不修则论文主张不成立

1. **配对活性稀疏：** ChEMBL 两端同时有可靠实测的分子远少于「看起来很多」；真 dual 更少。若 D2-public 的 dual 正例 < 临界规模，打假表方差爆炸。  
2. **校准泄漏：** 必须 **按 scaffold / 分子簇** 做校准拟合，禁止用 test 分子拟合 isotonic。  
3. **未测 ≠ inactive：** 方案已写，但工程上极易在负例脚本里违约——需断言测试。  
4. **「打假表」必须先出：** 没有 mean/sum 系统性抬升 A-only 的 Figure 1，NMI 叙事空心（CleanSplit 同构要求）。

### P1 — 不修则难上 NMI，仍可能上 JCIM

5. **靶点对过少 / 仅 NLRP3–JNK1：** 私有盲测不能当主证据；公开至少 2–3 对 + leave-pair-out。  
6. **细胞 holdout 过宣：** 细胞活性 ≠ 双靶结合；只能报排序相关/富集，不能报「验证了双靶结合」。  
7. **与 VSDS-VD 区分不足：** reviewer 会说「你们是双靶版附录」。必须强调 **配对标签 + dual-vs-single** 是新任务，不是多跑一个靶。  
8. **fused/linked 专家无消融：** 撞 anti-pattern 风险。

### P2 — 工程质量

9. 尚无对接 YAML / 融合代码（文档超前于实现）。  
10. 引擎版本、随机种子、box 定义未冻结。  
11. PK 案例若写成预测模型 = 超范围（方案已警告，执行时易滑）。

---

## 7. 可完成性判断

### 7.1 什么叫「可完成」

| 目标 | 可行性 | 条件 |
|------|--------|------|
| **完成一套可投稿计算包（JCIM / J. Cheminform. / Chem. Sci.）** | **高** | 锁 O1–O3；做出 D2≥2 公开靶点对；跑通打假表+校准 softmin 消融；开源 |
| **完成 Dual-VSDS 开放基准（资源型）** | **中高** | 劳动密集（文献 curated + ChEMBL 配对清洗）；技术无奇迹，要人月 |
| **NMI 接收** | **中 / 条件可达成** | 必须同时具备：① 显著且可复现的「朴素融合失败」；② 多靶点对泛化；③ 协议质量对齐 VSDS-VD；④ 清晰 delta 对抗 Pérez-Castillo；⑤ 私有细胞仅作 holdout |
| **NMI + 湿实验成药故事** | **低（按你当前约束）** | 缺亲和力/动物时不要并行这条线 |

### 7.2 ResearchStudio 模式匹配

| 模式 | 匹配度 | 含义 |
|------|--------|------|
| `controlled_diagnostic_design` | **高** | 主叙事应是「诊断尺子坏了 + 给出修正协议」 |
| `characterize_limit_then_surpass` | **中** | 你有「limit」直觉（mean 不可比），但缺**紧刻画**；实证打假可代替形式化，但别装定理 |
| `reframe_as_solvable_object` | **高** | 把「双靶好不好」重定义为 dual-compatibility score |
| 纯 benchmark construction | **风险** | IdeaSpark 对「只建库」偏苛；NMI 的 VSDS-VD 证明 **benchmark + 可执行协议 + 反直觉发现** 可过——你必须有反直觉发现，不能只有表头 |

### 7.3 完成路径（审计建议的最小充分集）

```text
锁 θ / 引擎 / 2–3 公开靶点对
    → 建 D2-public（断言：无未测当负）
    → 固定对接 → 导出姿态分
    → scaffold 外校准 → 打假表（mean 抬 A-only）【门控：无此表则降级 JCIM】
    → softmin + 硬负消融
    → leave-pair-out
    → D5 细胞排序盲测（附录级）
    → 开源 → 投 NMI；若打假效应量弱 → 改投 JCIM/Chem.Sci. 不硬冲
```

---

## 8. 与「放弃 / 改写」的边界

**不建议放弃** 整个双靶兼容性方向（scoop-check 未到 Level 1/2 硬封顶）。  

**建议改写/收缩的宣称：**

- 删掉任何「新对接算法 / 新通用 scorer」暗示。  
- 细胞数据不得写进 Abstract 主贡献句。  
- softmin 本身不宣称理论创新；宣称 **任务 + 诊断发现 + 基准**。  
- 若公开 dual 正例实在太少：主贡献改为 **「双靶 VS 评测协议 + 负结果（朴素融合有害）」**，方法头降为协议组件——仍可发，但 NMI 更难。

---

## 9. 审计产物索引

| 文件 | 内容 |
|------|------|
| `idea_card.md` | ResearchStudio idea-quality 输入卡 |
| `RESEARCHSTUDIO_AUDIT.md` | 本报告 |
| `steps/paper_search_raw.json` | scoop-check 检索原始输出 |

---

## 10. 给作者的三行行动令

1. **本周：** 冻结 θ、主引擎版本、≥2 个公开靶点对；写死校准/测试划分（O1–O3）。  
2. **下阶段唯一门控实验：** 打出「mean/sum 抬升 A-only」表；失败则立即降级期刊预期。  
3. **补 falsification 段落后再谈 NMI**；同时对 fused/linked 做开关消融以避开 anti-pattern 质疑。
