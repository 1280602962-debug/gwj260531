# PaperSpine 文章思路规划：分子从哪来 + 全文叙事

**Scene:** journal（偏 *JCIM* 方法叙事）  
**Contribution:** Option A（已锁定）  
**Purchase:** 690（家族/pose 锚）+ 2231（MD 偏好假说，grade C 风险明示）  
**软件现实：** 历史 Glide 筛选不可“假装开源可复现”；正文用 **Selection / Confirmation** 双层叙述。

---

## 0. 核心贡献句（全文只守这一句）

> 我们给出一条可追溯的商业库→酶活板短名单管线，并证明在 JNK1/2/3 上，常用计算 isoform 选择性过滤器（Δsel / Gly87 / ML 选择性标签）**不足以作为采购依据**；采购分子用于检验**家族结合富集**，对 JNK1 偏好仅作预注册次级观察。

审稿人应记住的 payoff：  
**“别再用便宜 Δsel 买 JNK 选择性；但可以用诚实管线富集家族 binder。”**

---

## 1. 分子“从哪来”——必须这样讲（最重要）

### 1.1 禁止的讲法（一写就被动）

| 错误叙述 | 为何危险 |
|----------|----------|
| “我们用 Gnina/Vina 筛选发现了 690/2231” | 事实是 Glide 历史漏斗选出；开源是后验验证 |
| “MD 证明 2231 具有 JNK1 选择性故采购” | MD≠酶选择性；且 2231 overall MD fail |
| “通过选择性过滤器筛选得到高选择性候选” | Δsel/Gly87/ML 已失败，且采购已解耦 |
| “开源重现了全部 Glide 排序” | 做不到也不需要；引擎本就会不一致 |

### 1.2 推荐总叙事：两层漏斗

```text
Layer A — SELECTION（发现/缩库，历史）
  商业库 → ML 家族活性召回 → Glide XP 对接短名单
  → ADMET → MD pose-QC → 人工采购决策
  ※ 选择性过滤器只做回顾失败分析，不进采购硬门

Layer B — CONFIRMATION（验证，开源、可复现）
  对已购 690/2231（+对照）→ Vina/Gnina 多 seed
  → 无约束 MD replicas →（到货后）JNK1/2/3 IC50
```

**一句话写进 Abstract/Intro：**  
*Candidates were prioritized by a documented ML–docking–ADMET–MD triage; isoform-selectivity filters were evaluated separately and not used as purchase gates. Open-source docking/MD then assessed pose consensus for the purchased set prior to enzymatic testing.*

### 1.3 两个分子各自的“出身证明”（Results 里分开写）

| 分子 | 在文中的角色 | 如何描述“怎么来的” | 必须同时承认的风险 |
|------|--------------|-------------------|-------------------|
| **690** | RQ-A 锚：家族活性 / pose 可信 | 通过 ML 召回 → Glide 短名单 → ADMET → **MD overall pass、grade A**；铰链三亚型都高 → **pan-leaning 活性锚**，不是选择性代表 | 可能无 isoform 偏好；活性可能弱于 E1/CC-90001 |
| **2231** | RQ-B 假说：最强 MD JNK1 bias | 同漏斗内 **MD bias 排序第 1**、Δsel_dock 偏向 JNK1、score_JNK1 最强；为检验“MD 不对称能否翻译成酶学偏好”而采购 | **pose grade C、pass_md_overall 否**；Vina 显示 JNK2 pose 不稳；旧 200 ns 有配体约束 → 正文降级为 historical |

**配对逻辑（写进 Purchase rationale 小节）：**  
- 690 = “管线在 pose-QC 及格区能否给出可测活性？”  
- 2231 = “即便 MD overall 不及格，最强 bias 信号是否仍值得酶学证伪？”  
- 二者合起来服务 Option A，而不是“两个都是优质 selective hit”。

### 1.4 Methods 里固定模板段落（可直接改英文）

**2.x Candidate selection (historical triage).**  
Commercial-library compounds were processed through a cascade of (i) ML family-activity recall (p_family≥6.0), (ii) Glide XP docking against JNK1/2/3 structures (3ELJ/3E7O/3TTI), (iii) ADMET/physchem filters, and (iv) MD-based pose quality control (ligand RMSD and hinge H-bond occupancy). Docking-derived Δsel, Gly87 occupancy, and ML selective-class labels were computed for analysis but **were not applied as hard purchase filters** after failing a literature benchmark (see Results). Two compounds were purchased for prospective assays: **690** (MD overall pass, grade A; pan-leaning hinge profile) and **2231** (strongest MD JNK1-bias score among the shortlist, purchased as a directional hypothesis despite grade-C / overall-MD-fail flags).

**2.y Open-source pose confirmation.**  
Independent of the historical Glide campaign, purchased ligands were re-evaluated with AutoDock Vina and Gnina (multi-seed) and unrestrained MD replicas to assess pose consensus. These calculations validate geometric stability for reporting; they do **not** redefine the purchase set.

**2.z Enzymatic panel.**  
JNK1/2/3 IC50 for 690, 2231, E1, and CC-90001 under pre-registered endpoints (any-isoform activity primary; SI≥3 vs both off-isoforms for preference).

---

## 2. Introduction 论证阶梯（PaperSpine checklist）

| 阶梯 | 写什么 | 落到贡献 |
|------|--------|----------|
| Problem | JNK1 bias 有转化动机（CC-90001/E1）；近同源 ATP 口袋 | 领域问题 |
| Progress | 对接差/IFP/FEP/ML 在其他激酶上有成功先例 | 公平承认先验 |
| Gap | JNK1/2/3 上缺乏对“便宜选择性过滤器”的硬负校准；采购仍常误用 | **具体缺口** |
| RQ | RQ-C 过滤器可信吗？RQ-A 管线能否富集家族活性？RQ-B 次级 | 问题精确化 |
| Contribution promise | 负校准 + 诚实管线 + 预注册酶活板 | Core contribution |
| Evidence preview | 文献 benchmark 失败表；漏斗；690/2231 设计；开源验证；IC50 | 证据预告 |
| Reader payoff | 以后别拿 Δsel 买 JNK 选择性；可用管线喂酶活板 | reviewer payoff |

Intro **不要**以“我们找到了两个抑制剂”起笔，而以“选择性计算在近同源位点会失效”起笔。

---

## 3. Results 单元 ↔ 贡献验证（results_validation）

| Results 单元 | 验证哪条贡献承诺 | 允许解读 | 禁止解读 |
|--------------|------------------|----------|----------|
| R1 漏斗图 | 管线可追溯、可操作 | 缩库路径清楚 | 宣称高 hit rate |
| R2 选择性尸检（C5） | **RQ-C 主贡献** | 过滤器不适合采购 | “所有对接都无用” |
| R3 采购设计（690 vs 2231） | 采购与选择性过滤解耦；次级假说可检验 | 设计合理 | MD=选择性 |
| R4 新颖性/PAINS | 非琐碎已知物/assay 假象风险可控 | 指纹远离已知 | kinome 干净 |
| R5 Vina/Gnina 共识 | 开源确认层 | pose 稳/不稳如实 | 开源重筛发现分子 |
| R6 无约束 MD replicas | pose QC 可重复 | 稳定性均值±方差 | 酶选择性证明 |
| R7 IC50（到货后） | RQ-A；条件满足才 RQ-B | 按 C4 规则 | 事后改 SI；MD 证实选择性 |

---

## 4. 推荐全文骨架（更新版，替换旧 2157 叙事）

### Title（工作标题）
Computational Isoform-Selectivity Filters Fail for JNK1/2/3: A Documented Family-Binder Enrichment Pipeline with Prospective Purchase of Pose-Credible and Bias-Hypothesis Candidates

### Abstract 六句
1. 近同源 JNK 口袋 → 选择性难  
2. 便宜过滤器常被当采购依据，缺负校准  
3. 我们做 ML→Glide→ADMET→MD 管线，并在文献集上测过滤器  
4. Δsel/Gly87/ML 失败；采购解耦；购 690（QC 锚）与 2231（bias 假说）  
5. 开源对接/MD 做 pose 确认；同批 JNK1/2/3 IC50（预注册）  
6. 结论：过滤器不可买选择性；管线可服务家族活性检验  

### 章节顺序（建议）
1. Introduction（阶梯）  
2. Methods  
   - 2.1 Data/benchmark  
   - 2.2 ML（召回 vs 选择性分类失败）  
   - 2.3 Historical Glide triage（Selection）  
   - 2.4 Selectivity metrics evaluated (not purchase gates)  
   - 2.5 ADMET + MD QC definitions  
   - 2.6 Purchase rule for 690 & 2231  
   - 2.7 Open-source confirmation (Vina/Gnina/MD)  
   - 2.8 Assay + pre-registered analysis  
3. Results（上表 R1–R7）  
4. Discussion：负结果意义；2231 风险设计；开源 vs Glide 不一致正常；无 kinome 局限  
5. Conclusions  

---

## 5. 三种到货结局，Discussion 怎么接（预写好）

| IC50 结局 | 分子来源叙述不变？ | Discussion 重点 |
|-----------|-------------------|-----------------|
| 仅 690 有活性，2231 无 | **不变** | 管线可富集；MD bias #1 为假阳性 → 强化“MD≠选择性” |
| 2231 有偏好（SI 达标） | **不变** | 次级发现可报；仍保留过滤器失败主线；勿升格成“选择性发现文” |
| 两者皆无，对照正常 | **不变** | enrichment failure 也可发（方法+负结果）；检查 assay/chemotype | 

**关键点：** 分子“出身故事”在到货前就写死；IC50 只改变 Results 末节与 Discussion 语气，不改 Selection 叙事。

---

## 6. 与开源/Gnina 叙事的衔接句（防审稿人抓）

> Although initial triage used Glide XP within an institutional docking environment, all pose-consensus and MD stability analyses reported for the purchased compounds were repeated with openly documented tools (AutoDock Vina, Gnina, and unrestrained MD). Score disagreements across engines are expected and are interpreted only as sensitivity of ranking, not as a mandate to revise the locked purchase set.

---

## 7. 你现在写作时的执行顺序（PaperSpine）

1. 守住 `confirmed_contribution.md` 一句（已锁；购买集已是 690+2231）  
2. 用本文 **Selection/Confirmation** 改 Methods 草稿  
3. Results 先写 R2（尸检）+ R3（采购出身）——这两节不依赖到货  
4. R5/R6 等 Gnina/MD 补完再填  
5. R7 到货后按 C4 填，禁止改终点  

---

## 8. 一句话文章观

这篇文章不是“讲述两个神奇分子如何被 AI 发现”，而是：  
**讲述一条诚实缩库管线如何选出可测分子，并证明那些诱人的选择性分数不该主导采购；690/2231 是管线上两种互补角色的代表，不是两个已证实的高选择性先导物。**
