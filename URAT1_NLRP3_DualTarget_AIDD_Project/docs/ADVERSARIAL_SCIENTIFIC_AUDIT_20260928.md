# URAT1–NLRP3 Dual-Target AIDD：对抗性科学审计

**审计对象：** `origin/cursor/urat1-nlrp3-dualtarget-aidd-e43d`（tip `bc6d5559`，证据锁基点 `cf5fd758`）  
**审计日期：** 2026-09-28  
**性质：** 作者授权的投稿前对抗审计（working draft）。不是期刊编辑决定。  
**原则：** 先读真实数据、脚本和配置，再独立判断；形式完整不能替代科学有效。项目已有的 `CONSOLIDATION_AUDIT.md` 只作交叉核对，不作免责金牌。

---

## 0. 调用的 skills 与方法边界

从 GitHub 读取并按当前版本执行：

| 包 | 仓库 | 本次实际使用的 skills |
|---|---|---|
| scientific-agent-skills | [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills) | `literature-review`，`hypothesis-generation`，`experimental-design`，`peer-review`（claim–evidence / reproducibility），`scientific-critical-thinking`，`scientific-writing`，`scholar-evaluation` |
| academic-research-skills | [1Dimension/academic-research-skills](https://github.com/1Dimension/academic-research-skills)（v1.8） | `academic-paper-reviewer`（EIC / methodology / domain / Devil’s Advocate / editorial synthesizer），`integrity_verification_agent`，`argument_builder` |

方法约束（来自上述 skills，不是装饰）：

- 文献缺口只报告“在本次检索边界内未定位”，不写“没有任何前人工作”。
- 假设、机制、预测、证据必须分开；HARKing 必须标成探索性。
- 报告完整性 ≠ 设计质量；数值可复现 ≠ 科学有效。
- Devil’s Advocate 的 CRITICAL 一旦成立，模拟决定不能是 Accept。
- 检索日期 2026-09-28；查询清单见第 1.4 节。

Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065

---

## 独立总判（先给结论）

**课题“URAT1+NLRP3 双节点”作为靶组合已经不是文献缺口。** 2025–2026 已有湿实验双靶/多靶分子（Zhang *Nat. Commun.* 2025 compound 32；Sun *Eur. J. Med. Chem.* 2025 HNW005）和临床双作用路径（arhalofenate；dotinurad / benzbromarone / probenecid 的炎症表型）。把临床库对接短名单写成“双靶发现”会超过证据。

**当前仓库的真实产物**是：一条可审计的临床酸性分子几何漏斗，冻结 12 个 primary / 21 个 reserve 的**口袋相容性假设**。主稿 `MANUSCRIPT_MASTER.md` 对“无双靶活性、无 MD、无新靶点确认”写得很克制，数字也与表一致。这只能说明**文案没有明显撒谎**，不能说明**课题现在值得作为研究论文发表**。

| 发表路径 | 判断 |
|---|---|
| 现在投 JCIM / JMC / EJMC，当双靶发现或新方法 | **不值得。** JCIM 明文不收无实验验证的常规对接筛选；EJMC/JMC 已有湿实验双靶文章。 |
| 现在投 JCAMD / *Molecular Diversity*，当 12 个 lead | **不值得。** A1 漏掉已知 URAT1 药、A2 无区分力、W4 阳性为同一化学家族、无 MD。 |
| 缩成 2–3 个具名假设 + 协议失败审计的短通讯，投计算化学中低档期刊 | **勉强可以考虑**，前提是执行下面 Critical/Major 的最低修复。 |
| 等 URAT1 摄取 + NLRP3 结合/ATPase（PF 须加 MR 对照）或至少 MD 对照通过 | **这才是该课题重新变得可发表的最低科学门槛。** |

**模拟编辑决定（academic-paper-reviewer / editorial standards，仅供作者规划）：Reject — Premature / Insufficient Contribution**，针对任何把 12/21 当成双靶候选发现的稿件。若改写成“几何门失效审计 + 2–3 个待验假设”，可再评估为 Major Revision 级别的短文，而不是 Article。

---

## 1. 课题创新性

**Skills：** literature-review，hypothesis-generation，peer-review / Devil’s Advocate，EIC originality，scholar-evaluation（只评作品，不评人）。

### 1.1 仓库自己声称的问题与假设

| 对象 | 仓库文本 | 审计意见 |
|---|---|---|
| 研究问题 | “临床阶段及已上市小分子中，哪些酸性分子具有值得实验验证的 URAT1 与 NLRP3 口袋结构相容性？”（`README.md`） | 这是**筛选问题**，不是可证伪的机制假设。 |
| 隐含假设 | NLRP3 相关分类器缩库 + GNINA 几何门 → 可检验的双口袋结构假设（`campaign_final.yaml`，`MANUSCRIPT_MASTER.md`） | 假设把“相关活性分类”偷换成“NACHT 口袋相容”，把“COOH–Arg 邻近”偷换成“URAT1 抑制”。 |
| 贡献 | 12 个冻结结构 primary；PF-03882845 为后续优先对象 | 贡献是名单与协议记录，不是新靶、新方法或新活性。 |
| 已否认 | 首次双靶、新靶点确认、临床用途、MD 轨迹 | 否认正确；**标题仍卖“双节点结构候选”**，读者会按发现文来读。 |

`hypothesis-generation` 要求：观察、问题、假设、机制、预测、替代解释必须分开。当前设计**没有能区分“口袋几何碰巧通过”与“双靶活性”的判别预测**。对 PF-03882845，已知 MR 拮抗 + 醛固酮肾病模型炎症基因下降（Meyers 2010；Orena 2013）是更简约的替代机制；主稿承认了这一点，但提名优先级仍把它写成 follow-up 首选，缺少预先登记的否证规则（例如：MR 对照下 IL-1β 仍降才算 NLRP3 臂阳性）。

### 1.2 文献缺口（2026-09-28 检索边界）

**已定位、足以占用“双节点”新颖性的前人工作：**

1. **Zhang et al., *Nat. Commun.* 2025, 16:7430.** DOI [10.1038/s41467-025-62645-6](https://doi.org/10.1038/s41467-025-62645-6)。按 URAT1+NLRP3 药效团设计 64 个类似物；compound 32（NLRP3/URAT1-IN-1）URAT1 IC50 3.81 µM，IL-1β IC50 2.61 µM，SPR KD 27.8 µM，体内降尿酸 + 痛风关节炎。**湿实验双靶已经做成。** 仓库 SI `published_dual_structures.csv` 已收录。
2. **Sun et al., *Eur. J. Med. Chem.* 2025.** DOI [10.1016/j.ejmech.2025.117644](https://doi.org/10.1016/j.ejmech.2025.117644)。HNW005：NLRP3 KD 204.6 nM / IC50 1.7 µM，URAT1 IC50 6.4 µM。仓库已收录，但 HNW005 SMILES 仍标厂商页，SI 未取到。
3. **Shi et al., *J. Med. Chem.* 2026.** DOI [10.1021/acs.jmedchem.6c00058](https://doi.org/10.1021/acs.jmedchem.6c00058)。URAT1 IC50 0.19 µM + IL-1β IC50 3.39 µM 的 dual-acting 候选。**仓库 SI 未收录。** IL-1β 抑制 ≠ 已证明的 NLRP3 NACHT 结合，但足够说明“URAT1 + 抗炎”已是 JMC 主线。
4. **Arhalofenate（临床双节点）：** URAT1/OAT4/OAT10 促尿酸排泄 + 抑制 MSU 驱动的 IL-1β。仓库主文未定位为最近临床对照。
5. **已上市 URAT1 药的炎症/NLRP3 表型：** dotinurad、benzbromarone、probenecid 均有巨噬细胞 IL-1β / caspase-1 文献。Lesinurad 专利曾主张 “dual inhibitor of URAT1 and inflammasome”。**“URAT1 药兼抗炎”不是意外发现。**
6. **Agarwal et al., *J. Biomol. Struct. Dyn.* 2024.** FDA 库 + 7ALV + ML/MD，无湿实验。说明“临床库对接 NLRP3 7ALV”也不是方法缺口。

**本次检索边界内未定位的、也是唯一还说得过去的增量：** 针对**实验 URAT1 冷冻电镜 9DKB** 的临床库双靶对接短名单，以及 PF-03882845 作为 URAT1/NLRP3 配体的湿实验证实。那是**候选级**增量，不是靶组合级或方法级增量。

FINER 的 Novel 项：在有记录的检索完成前必须标为 unresolved。本次检索之后，**靶组合 Novel = 不成立**。

### 1.3 原创性来源判定（EIC / ScholarEval 作品向）

| 声称的原创性来源 | 判定 |
|---|---|
| 新靶组合 | **否。** 2025–2026 已被湿实验和临床路径占据。 |
| 新方法 | **否。** GNINA CNN rescore + XGBoost 骨架 CV + 几何门是常规管线；P2 排序轨自己锁定为失败迁移。 |
| 新数据 | **弱。** 临床库公开；NLRP3 标签来自 ChEMBL 异质 assay；对接姿势是自产，但 C1 生产 SDF 无 GNINA 版本日志。 |
| 新候选 | **最多是未验证假设。** 12 个名字来自几何 ∩ 历史药名排除；去掉药名规则会另入 MK-5108 / Piromidic acid / Pradofloxacin。 |
| 新视角（失败审计） | **这是仓库里真正有价值的东西**：lesinurad Top-1 失败、A2 LR+≈1、W1 选择失败、1588 不可比特重放。目前被埋在“双节点候选”叙事下面。 |

**So what 测试：** 若结论成立，领域不会多一个双靶策略，只会多一份“临床酸 + 9DKB/7ALV 几何通过名单”。对痛风药物化学，这份名单在没有 IC50 之前几乎不改变任何人的下一步合成或实验优先级——除非缩到 1–2 个名字并立刻做转运/结合实验。

### 1.4 检索查询（文献 skill：可复现边界）

日期 2026-09-28。数据库/界面：WebSearch + 出版商页 WebFetch（Nature、ACS JCIM about、PMC）。不是系统评价，不是 PRISMA。

代表性查询：`dual-target URAT1 NLRP3 inhibitor gout 2020..2026`；`HNW005 tranilast URAT1 NLRP3`；`Zhang Eurycoma longifolia Nature Communications 2025`；`Pyridoimidazolyl Sulfonamide Compound 2 URAT1 Journal of Medicinal Chemistry 2026`；`arhalofenate URAT1 NLRP3`；`FDA approved drug library virtual screening NLRP3 7ALV`；`9DKB URAT1 docking virtual screening`；`JCIM virtual screening without experimental validation`；`PF-03882845 URAT1 NLRP3`。

---

## 2. 对接审计

**Skills：** experimental-design，reproducibility（peer-review statistical_reproducibility），ablation（作为设计完整性，不是单独 skill 名），methodology reviewer，Devil’s Advocate。

### 2.1 生产协议（从配置与表读取，不从旧 Markdown 宣传读取）

| 项 | 锁定值 | 证据 |
|---|---|---|
| 引擎 | GNINA，`cnn_scoring=rescore` | `config/docking_final.yaml` |
| exhaustiveness | 32 | 同上 |
| seeds | 42/43/44；A1 提名只用 42 | `urat1_gate_definition.yaml` |
| URAT1 生产受体 | **仅 9DKB** | `MANUSCRIPT_FACTS.yaml` |
| NLRP3 生产受体 | 7ALV | 同上 |
| 盒 | URAT1 22³ @ `[99.966, 102.967, 105.699]`；NLRP3 20³ | `docking_final.yaml` |
| 姿势选择（A1，决定 12 人名单） | **先 CNNscore Top-1，再几何** | `urat1_gate_definition.yaml` A1 |
| 质子化 | 受体 pH 7.4；配体 Dimorphite→Meeko，单选微观状态 | `docking_final.yaml` |

### 2.2 致命设计问题：采样成功、选择失败，生产仍用失败的选择器

Lesinurad@9DKB，exhaustiveness 32（`data/redock_smoke/redock_results_lesinurad_9DKB.csv`）：

| 读出 | Top-1 RMSD (Å) | best-of-9 (Å) | Top-1 ≤2 Å |
|---|---|---|---|
| Vina P1 | 0.856 | 0.856 | 是 |
| **CNNscore P0（A1 同选择器）** | **4.163** | 0.994 | **否** |
| CNNaffinity P2 | 4.163 | 0.994 | 否 |

Benzbromarone@9DKA W1 门：三 seed Top-1 RMSD 3.585 / 3.603 / 3.595 Å，best-of-9 ≈1.11 Å，全部 `pass: false`。主稿 Results 也写了同一事实。

**实验设计结论：** 这不是“对接没跑起来”，而是 **Fisher 意义上的选择器失效**：近晶体姿势在 ensemble 里，CNNscore Top-1 不取它。A1 却把同一选择器当成 primary 入选规则。已知临床 URAT1 药（lesinurad / verinurad / puliginurad）在同一逻辑下 A1 失败。因此 12 个 primary 的可辩护解释只有一句：**“CNNscore Top-1 酸锚几何相容”**，不能写成 URAT1 活性富集或优于已知降尿酸药。

W1 失败后仍授权 W4（`run_c5_w4_nlrp3_panel.py`，`campaign_c5.yaml` `search_ok_selection_fail`）。科学上 W4 是另一个靶，不依赖 W1；叙事上这是**失败后改规则继续往下走**，不能叫预注册漏斗。

### 2.3 双臂“结构交集”几乎是单臂

A1（seed 42，228 vs 64）：sens 0.447，spec 0.797，LR+ **2.20**。  
A2：sens 0.952，spec **0.047**，LR+ **≈1.00**，Fisher p=1.0。  
漏斗第 5 步的 54 个“dual structural ≥2/3 seeds”用的是 **A2 ∧ NLRP3**，不是 A1。A2 对 URAT1 活性没有判别力，则 54 这个数主要是 NLRP3 结构门 + 几乎不过滤的 URAT1 可容纳性。

12 个 primary 另加 **A1 seed42**。同一漏斗里两套 URAT1 规则：宽的用来做“双结构”分母，严的用来做名单。这是 post hoc 分层，不是单一预指定分析。

### 2.4 NLRP3 门：家族内自检，失败对接当阴性

W4 seed42 结构门 `intent_to_score`：TP=9，FN=0，FP=11，TN=**29**，spec 0.725。  
`complete_case`（去掉 5 个全 seed 失败诱饵）：TN=24，spec **0.686**。  
临床酸背景切片上，structural 与 loose **完全相同**（spec 0.40）。

阳性 9/9 来自 MCC950/NP3-146 磺酰脲家族（配置已写 `validation_family_limitation`）。诱饵来自 **URAT1 true-decoy 性质匹配**，不是实验 NLRP3 阴性。`failed_docking_interpretation` 写了“计算不通过 ≠ 实验无活性”，但主数字仍把失败算进 TN=29。

这是 `missing-as-fail` 直接进入特异度。按 experimental-design：阴性对照必须不能走目标机制、同时共享偏倚通路。URAT1 decoy 当 NLRP3 阴性不满足这一条。

### 2.5 消融完整性

**做了的：** OAT 迁移；NLRP3 阈值 5.5/6.0/6.5；P0–P5 协议；decoy 泄漏审计；药名排除反事实（MK-5108 等）；Arg ±0.5 Å 敏感性。

**必须做、没做、因此不能做因果声称的：**

| 缺失消融 | 为什么是最低必要而不是锦上添花 |
|---|---|
| 漏斗顺序（先 ML 再酸再对接 vs 先酸再对接） | 否则“双节点”与“NLRP3 富集酸”无法分开 |
| 结构-only 完整重建 12/21 | 现有反事实只列出 3 个名字，不重建短名单 |
| A1 距离门 vs W2 IFP **作为成员规则** | 4/12 primary W2 ≥2/3 失败（Lintitript 0/3）但仍是 primary |
| ≥2/3 seed vs 3/3 vs 单 seed | A1 提名单 seed 42 |
| 9DKB vs 9DKA/9DKC 对同一 156 人重对接 | 只有刚体转移敏感性 |
| num_modes 1 vs 9 | 遗产 P2 导出 1 mode，C1 门 9 modes |
| 临床背景 vs 磺酰脲阳性作为 NLRP3 验证主表 | 临床背景上结构门零增量 |

### 2.6 可重复性

| 问题 | 证据 | 对结论的影响 |
|---|---|---|
| 1588 池不能比特重放 | 2026-07-01 同交 joblib+CSV，无当时 pip 清单；同阈值重训得 1377 | ML 缩库是漏斗第一刀，成员不可重算 |
| 模型哈希 ≠ 历史生成器证明 | `MANUSCRIPT_FACTS.yaml` | 只能以冻结表为成员依据 |
| C1 生产酸对接 0 条 GNINA 版本日志 | `gnina_execution_log_summary.csv` | 12/21 所依据的姿势版本 **unknown** |
| 文档曾写 GNINA 1.3.1，沉积日志为 1.3.2 | 同上 | 中等，但生产姿势本身无日志 |
| `n_jobs=-1` | pickle / SI 重训 | 未证实是 1377 vs 1588 的原因，禁止写成已验证因果 |

**伪重复：** 三 seed 不是三个配体。验证表按 seed 分列是对的；任何把 3×N 当样本量的句子都是错的。

---

## 3. 投稿前检查

**Skills：** peer-review claim–evidence，scientific-writing evidence binding，academic-paper argument_builder，EIC journal-fit，integrity 的数据层（本次不是 100% 文献逐条 DOI 核验）。

### 3.1 中心论题能否成立

主稿可辩护论题只有这一句：

> 在已存临床酸性分子上，用 NLRP3 相关分类器缩库和 GNINA 几何门，可以得到一份**待实验检验**的口袋相容名单；该名单**不是**已验证双靶抑制剂。

标题“临床酸性化学空间中的 URAT1–NLRP3 双节点结构候选”比论题强一档。Abstract 工作稿前半仍按发现文节奏写漏斗，最后一句才收回。按 scientific-writing：**title keyword 必须各有一项度量支持**。“双节点候选”没有活性度量。

### 3.2 Claim–evidence（摘要级；全表见 CSV）

| ID | 声称 | 证据 | 强度 | 超证？ |
|---|---|---|---|---|
| C1 | 8319→1588→303→156→12/21 | `screening_funnel.csv`，`checks.json` | 强（计数） | 否 |
| C2 | NLRP3 CV AUROC 0.893 | `model_validation.json`，513 分子骨架 CV | 中 | **是，若暗示对 12 人的 NACHT 结合** |
| C3 | A1 LR+ 2.20 | 228 vs 64 羧酸回顾 | 中 | **是，若写成 URAT1 检索器**（sens 0.45，已知药失败） |
| C4 | 54 个双结构分子 | A2∧NLRP3 ≥2/3 | 弱 | **是**：A2 spec 0.047 |
| C5 | W4 spec 0.725 | intent-to-score，含 5 失败诱饵 | 弱 | **是，若写成 NLRP3 非活性识别** |
| C6 | 12 个结构候选值得验证 | 冻结表 + 几何 | 弱 | **是，若“值得”隐含优先级高于已知双靶分子** |
| C7 | PF-03882845 优先后续 | 多 seed 几何 + 药化注释 | 弱 | **是，若暗示 URAT1/NLRP3 新靶**（主稿已否认，提名仍在） |
| C8 | 非首次双靶 | HNW005 + compound 32 | 强 | 否；但漏 Shi 2026 / arhalofenate |
| C9 | 无 MD / 无湿实验 | `md_authorized: false` | 强 | 否 |
| C10 | 结构过滤自动得到 12/21 | 药名排除反事实为否 | — | **数据与这句话相反**；主稿已改口，冻结名单未改 |

逻辑链断裂点（Devil’s Advocate CRITICAL 候选）：**从“几何通过”推到“双节点候选”缺少活性端点。** 这不是措辞问题，是论题问题。

### 3.3 稿件逻辑

1. 归档 `paper_spine_ars_analysis/confirmed_contribution.md` 仍说主贡献是 **P2 对接不能当 URAT1 检索器、双百分位名单不是双节点候选**。现行主稿相反。仓库里两套贡献声明并存。
2. 54∩40=37 是交集不是顺序筛选；40 含药名规则。漏斗图若画成 54→40→12 就是错的。现行表已写对，任何 PPT/图必须跟表。
3. Lintitript Arg 7.676 Å，门限 7.7027 Å（晶体 6.70+1.0）。+0.5 Å 会踢掉它。阈值不是回顾面板上 FP 最优。这是典型 post hoc margin。
4. Runcaciguat cLogP 6.88。软门没有 logP，主稿已承认；仍须防止“类药过滤后的 12 个”这种句子。
5. 英文稿、目标期刊、CONSORT/TRIPOD 类报告清单均未锁定。状态“尚未投稿定稿”与 README“并非投稿就绪”一致。

### 3.4 期刊匹配（EIC）

| 期刊 | 匹配 | 原因 |
|---|---|---|
| **JCIM** | 不匹配 | Aims：不收无充分实验验证的常规对接应用；2015 编辑信把“virtual screening without experimental validation”列为不太可能过审。本工作无新算法、无湿实验、无大基准方法贡献。 |
| **EJMC / JMC** | 不匹配 | 只要体外活性；同领域已有 HNW005、compound 32、Shi 2022/2026 dual-acting。 |
| **JCAMD** | 差 | 优先实验验证；可作为失败协议短文，不能当 12 lead 发现。 |
| **Molecular Diversity** | 差到边缘 | 社区常预期具名 2–5 个分子 + MD；12 个方法漏斗过重、产品过弱。 |
| **JBSD / BMC Chemistry / Molecules 一类** | 才是计算-only 短名单的现实出口 | 即便如此，也必须先把论题降到假设名单，并修 Critical 项。 |

---

## 4. 代码 / 人为因素审计

**重点：** hard-coded selection，fallback，missing-as-fail，post hoc threshold，candidate-name exclusion，frozen-output provenance。加上模型泄漏。

### 4.1 仍决定 12/21 成员的人为规则

| 类别 | 代码位置 | 是否仍进入 12/21 | 披露？ |
|---|---|---|---|
| 药名子串排除 FLOXACIN/PIROMIDIC/MK-5108/… | `config/chemistry_final.yaml` L20–29；历史脚本 `build_c1_acid_shortlist_a2.py` | **是**（40 eligible 的 provenance） | 是；决定是不改 12/21 |
| β-lactam SMARTS 末级排除 | `chemistry_final.yaml`；`freeze_c5_shortlist.py` | 是（21 reserve） | 是 |
| GSK 结构对照排除 `REP_07907` | 同配置 | 是 | 是 |
| A1 距离 7.7027 Å = 晶体 min+1.0 | `c1_acid_pose_selection.py`；`urat1_gate_definition.yaml` | **是（primary）** | 部分 |
| A1 只用 seed 42 | 同上 | 是 | 是 |
| NLRP3 结构地板 overlap/IFP 0.50、5/7 接触 | `c1_nlrp3_pose_metrics.py` L35–39，“upgrade 2026-08-27” | 是（A2 臂） | 阈值锚定在 NP3-146 自对接之下，有前瞻声明也有看过阳性后加严的痕迹 |
| 分类器 p≥0.5 生产池 | `screen_repurposing_library.py` | 是（第一刀） | 是 |
| 酸几何匹配器只认 COOH 氧 | `chemistry_final.yaml` L36 | 影响解释，156 中 21 个无该药效团 | 是 |

**现行机械路径** `build_c5_tier_assignment.py` 不再硬编码药物商品名，这是改进。但 **40 eligible 仍是带药名排除的历史集合**，12/21 是在这个集合上切出来的。结构-only 反事实已经证明名单会变。`pre_submission_decision: accept_as_limitation_do_not_replace_frozen_12_or_21` 把人为规则冻成论文事实。

### 4.2 已退役但仍在仓库里的硬编码挑选

`scripts/build_c1_acid_shortlist_a2.py`：

```python
PRIMARY_TIER1 = ["PF-04620110"]
PRIMARY_TIER2 = ["ADMILPARANT", "RUNCACIGUAT", "LANIFIBRANOR"]
BACKUP_TIER = ["PSI-697", "PF-03882845"]
```

`__main__` 会 SystemExit。**不进入现行 12/21。** 仍构成“曾经按名字点将”的 provenance。现行 12 人里仍有 LANIFIBRANOR、ADMILPARANT、RUNCACIGUAT、PSI-697、PF-03882845——与旧名单高度重叠。机械规则是否只是把旧偏好翻译成阈值，需要结构-only 重建来否证，而作者选择不重建。

W4 阳性 SMILES 硬编码在 `run_nlrp3_structural_panel.py` 的 `HARDCODED_POSITIVES`（NP3-146、MCC950）。影响验证面板，不影响 12/21 成员。

W4 失败诱饵 ID 硬编码：`FAILED_DECOYS = {"C5W4D_002", ...}`（`score_c5_w4_nlrp3_panel.py`）。这是事后登记的失败集合，再用于 intent-to-score。

### 4.3 Fallback / missing-as-fail

- 无 pose → `keep_*=False`（NLRP3 / W2）：缺失当不通过。对验证特异度，这把计算失败变成 TN。
- `prepare_ligands_c1.py` 中 Meeko 路径 `except Exception: pass`：配体可能被静默丢掉，进不了 12/21。这是 **missing-as-exclude**，会改变分母。
- `freeze_c5_shortlist.py`：SDF 缺失仍调用 `evaluate_sdf`，W2 标失败，**不踢出名单**（注释轨）。
- `merge_docking_pareto.py` inner join：缺对接则丢行。C5 主路径不用 Pareto 选 12/21。

### 4.4 Post hoc / HARKing 时间线（`campaign_c5.yaml` 已记录，仍算 HARKing）

1. A1 绝对 4.0 Å 因严于晶体 6.70 被退役，改晶体+1 Å。  
2. A2 修正（2026-08-27）：先几何再 CNNscore。  
3. W1 Top-1 失败 → 解释改成 `search_ok_selection_fail` → 授权 W4。  
4. W2 锚点从 CNN 重对接受体改为沉积晶体坐标。  
5. M3：名字分层改为机械分层。  
6. 整理阶段明确 **不是新的前瞻预注册**（`campaign_final.yaml`）。

这些记录比瞒报好。它们不能把探索性漏斗升级成确证性发现。

### 4.5 冻结产物溯源

- `freeze_a1_exploratory.py`：`shutil.copy2`，**不从 SDF 重算**。A1 腿是复制冻结。
- `freeze_c5_shortlist.py`：0 次新对接，重打 W2 注释，不洗牌名字。
- `build_evidence.py`：SHA-256 + 成员重算检查。这能证明**从表到表一致**，不能证明**从原始对接作业到表一致**（C1 无版本日志）。
- 模型：`nlrp3_model.joblib` 与 1588 表同日提交；**没有同时生成的环境哈希**。冻结输出的 provenance 在 ML 第一刀处断开。

### 4.6 泄漏

| 泄漏 | 状态 |
|---|---|
| NLRP3 训练分子骨架 CV | 做了；不能证明临床库 ID 从未出现在训练折 |
| 8319 ∩ 训练 = 3 个临床 NLRP3 药 | 归档 presubmission 有；**现行主稿未写** |
| 12/21 与训练 ChEMBL ID 精确交集 | 证据层检查为空（好） |
| 同一分类器既训练又给 8319 打分再切 p≥0.5 | **筛选泄漏/循环使用**，不是标签泄漏；12 人全部来自该池 |
| W4 诱饵 = URAT1 decoy | 标签错误，不是训练集泄漏 |
| W4 阳性从 `nlrp3_records.csv` 挖磺酰脲 | 与训练集同语料，家族自检 |
| ECFP4 相对已发表双靶 Tc ≤0.33 | 骨架不撞车，不恢复靶组合新颖性 |

---

## 5. 分级问题与最低必要修复

只列**改变能否发表或改变 12/21 含义**的修复。润色、补引用、补图不算最低必要。

### Critical

| ID | 问题 | 为何是 Critical | 最低必要修复 |
|---|---|---|---|
| K1 | 靶组合当创新卖点 | 2025–2026 湿实验双靶已存在；论题“双节点发现”地基塌掉 | **删掉发现式标题/贡献。** 最多写“临床库在 9DKB/7ALV 上的几何短名单”；或改成对接选择器失败审计。 |
| K2 | A1 使用已证明失败的 CNNscore Top-1 | 共晶配体 Top-1 RMSD>4 Å，已知 URAT1 药 A1 失败；12 人由该规则定义 | **停止用 A1 Top-1 定义 primary。** 最低替换：best-of-9 几何 或 Vina Top-1（lesinurad 自对接通过的那条），并重算名单；若不愿重算，**把 12 人降为探索性附录，正文只报告选择器失败。** |
| K3 | 无活性证据仍维持“候选”产品 | 逻辑链在几何处断开 | **正文候选数降到 0 或 2–3 个具名假设**，并写清否证实验；12/21 只许出现在 SI 作为漏斗审计。 |
| K4 | 1588 ML 池不可重放，却是漏斗第一刀 | 成员生成器与沉积物无法闭合 | **要么给出可比特重放的训练/打分环境，要么把 1588 改为“历史冻结输入”并补结构-only 对照漏斗作为主分析。** |

### Major

| ID | 问题 | 最低必要修复 |
|---|---|---|
| M1 | 药名排除仍塑造 40→12/21，却声称结构规则 | 主分析改用结构-only 集合，或把 12/21 标为“含历史药名规则的遗产名单”，不以之投稿 |
| M2 | A2 LR+≈1 仍用于“双结构 54” | 删除 54 作为双靶证据；双臂必须两边都有判别力 |
| M3 | W4 失败对接计入 TN；诱饵非 NLRP3 阴性；阳性同一家族 | 主表只用 complete-case；特异度不得解释为 NLRP3 识别；验证改用非磺酰脲阳性或承认无跨家族证据 |
| M4 | W2 不是成员规则，4/12 失败 | 要么 W2 进入成员，要么主表不要把 IFP 画成质量保证 |
| M5 | C1 生产 SDF 无 GNINA 版本日志 | 重跑 156×2×3 并写版本，或声明 12/21 姿势不可版本化复现、不得当主结果 |
| M6 | PF 优先与 MR 混杂 | 预注册：IL-1β 实验必须带 MR 拮抗对照；做不到就把 PF 从“优先对象”降为“需排除 MR 解释的案例” |
| M7 | 文献缺口陈述不完整 | 补 Shi 2026 JMC、arhalofenate、Agarwal 2024 7ALV FDA VS、dotinurad/probenecid 炎症表型 |
| M8 | 漏斗顺序 / ML vs 结构 未消融 | 至少做一次“跳过 p≥0.5、仅酸+几何”的完整短名单，报告与 12 人的差集 |

### Moderate

| ID | 问题 | 最低必要修复 |
|---|---|---|
| D1 | Arg +1.0 Å 保住 Lintitript | 主文报告 +0.5/+1.0 两个名单，或论证 1.0 的事前理由（晶体误差，不是保住某人） |
| D2 | 归档 spine 与现行贡献相反 | 删除或醒目标记内部 spine，避免投稿包带两套故事 |
| D3 | 盒子 20 Å vs 22³、P2 `num_modes=1` vs C1 9 modes | 方法只锁定生产配置，遗产导出标 closed |
| D4 | 8319∩训练=3 未进主稿 | 一句话写入 Methods |
| D5 | `except Exception: pass` 配体准备 | 失败必须写进 manifest，禁止静默丢分子 |
| D6 | 无英文定稿、无目标期刊 | 若仍要投，先锁期刊再锁论题，不要先写中文发现稿 |

### Minor

| ID | 问题 | 最低必要修复 |
|---|---|---|
| N1 | HNW005 SMILES 来自厂商页 | 拿到论文 SI 或标 “structure not used for any activity claim” |
| N2 | 文档曾写 GNINA 1.3.1 | 已部分纠正；清掉残余 |
| N3 | 退役硬编码名单仍可被误跑 | 已 SystemExit；再加测试禁止 import 侧效应即可 |
| N4 | Runcaciguat cLogP>5 | 已承认无 logP 硬门；图表不要画成 Ro5 通过 |

---

## 6. 模拟评审摘要（academic-paper-reviewer）

**Field：** 计算药物发现 / 痛风双靶，应用对接，非新方法。  
**成熟度：** 探索性筛选记录，未到投稿 Article。

| 角色 | 分数意向（相对 JCIM/JMC 标准） | 要点 |
|---|---|---|
| EIC | Originality 弱；journal fit 差 | 靶组合不新；无实验；像内部竞选记录 |
| Methodology | Rigor 弱到中 | 选择器失效仍用于入选；验证面板偏倚；1588 不可重放 |
| Domain | 文献覆盖中 | 已引 HNW005/32，缺临床双节点与 2026 JMC |
| Perspective | 实际影响低 | 不改变合成或临床决策 |
| Devil’s Advocate | **CRITICAL：几何 ≠ 双靶活性** | 更简约解释：CNNscore Top-1 酸几何过滤器 + NLRP3 磺酰脲自相似 + 药名规则 |

加权远低于 Reject 线（rubric <50）。**DA CRITICAL 存在 ⇒ 不能 Accept。**

---

## 7. 什么情况下这个课题重新“值得做/值得发”

按 hypothesis-generation：先冻结判别预测，再收结果。

1. **论题改成二选一，不要混装：**  
   - (A) 开源对接选择器在 URAT1 共晶上失败的可重复审计；或  
   - (B) 1–2 个临床分子的双靶湿实验（URAT1 摄取 IC50 + NLRP3 SPR/ATPase 或带 MR 对照的 IL-1β）。
2. **若走 (B)：** 用 best-of-9 或实验上说得通的姿势规则重算短名单；去掉药名排除；不要把失败对接当阴性。
3. **MD 不是装饰，但是第二位。** 没有功能实验时，MD 最多支持“姿势不立刻爆炸”，支撑不了双靶声称。有功能实验时，再按已锁的 `md_protocol.yaml` 跑对照。
4. **不要再扩 12 个 lead。** 名单越长，越像没做选择。

当前仓库作为**内部证据层**是有价值的：漏斗可重建、限制写得比多数对接稿诚实。作为**科学论文产品**，在修完 K1–K4 之前不应投稿。

---

## 8. 证据索引

- 稿件与锁：`docs/MANUSCRIPT_MASTER.md`，`MANUSCRIPT_FACTS.yaml`，`README.md`，`config/campaign_final.yaml`
- 漏斗与门：`data/manuscript/tables/screening_funnel.csv`，`urat1_gate_validation.csv`，`nlrp3_gate_validation.csv`，`urat1_redocking_summary.csv`
- 冻结名单：`data/frozen/primary_candidates_frozen.csv`
- 反事实：`data/manuscript/si/structure_only_chemistry_counterfactual.csv`，`published_dual_structures.csv`，`chemistry_limitation_acceptance.md`
- 对接失败：`data/redock_smoke/redock_results_lesinurad_9DKB.csv`，`data/campaigns/c5/01_crossdock/gate_benzbromarone_9dka.json`
- 代码：`scripts/build_c5_tier_assignment.py`，`c1_acid_pose_selection.py`，`freeze_c5_shortlist.py`，`freeze_a1_exploratory.py`，`score_c5_w4_nlrp3_panel.py`，`manuscript_pipeline/build_evidence.py`，`config/chemistry_final.yaml`，`config/urat1_gate_definition.yaml`，`config/nlrp3_gate_definition.yaml`
- 对照矩阵：`data/manuscript/audit/claim_evidence_matrix_20260928.csv`

本次**没有**重跑 GNINA、没有重训 XGBoost、没有做 100% 参考文献 DOI 逐条核验（integrity Phase A）。文献判断基于 2026-09-28 的检索边界。湿实验有效性无法从本仓库证实，也未声称证实。
