# NMI 及高分文章：流程最可参考清单

面向本课题（**双靶兼容性打分 + Dual-VSDS 基准**，冲击 *Nature Machine Intelligence*）的**叙事与实验设计对标**，不是泛泛文献列表。  
判定标准：是否「纠偏主流做法 → 重定义任务/评测 → 泄漏控制基准 → 消融证明朴素基线失败 → 开放资源 ± 有限湿实验」。

相关主方案见 [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md)。

---

## 0. 一句话结论

| 用途 | 首选文章 | 期刊 |
|------|----------|------|
| **整篇实验骨架最像你要写的** | Gu et al. **VSDS-VD** | **Nat. Mach. Intell. 2025** |
| **打分方法叙事最像（校准/泄漏/外测）** | Cao et al. **EquiScore** | **Nat. Mach. Intell. 2024** |
| **「打假+泄漏」cover letter 语气** | Graber et al. **PDBbind CleanSplit** | **Nat. Mach. Intell. 2025** |
| **姿态合理性评测语言** | Buttenschoen et al. **PoseBusters** | Chem. Sci. 2024（NMI 常引） |
| **双靶化学叙事（非对接打分主线）** | Isigkeit / POLYGON 等 | Nat. Commun. |

**不要把 POLYGON / CLM multi-target 生成文当主骨架照抄**——它们是「生成+合成验证」，你是「判别/排序+基准」。

---

## 1. Tier S：流程几乎可直接映射（优先精读）

### 1.1 Gu et al. — VSDS-VD（**最像**）

- **题名**：*Benchmarking AI-powered docking methods from the perspective of virtual screening*  
- **期刊**：*Nat. Mach. Intell.* **2025**, 7, 509–520  
- **DOI**：[10.1038/s42256-025-00993-0](https://doi.org/10.1038/s42256-025-00993-0)  
- **资源**：[GitHub VSDS-VD](https://github.com/shukai1997/VSDS-VD) · Zenodo 数据/代码

**叙事弧（建议你按同构写）：**

1. 指出主流评测偏「再对接 RMSD」，**忽略物理合理性与真实 VS**  
2. 构建分层基准：**TrueDecoy / RandomDecoy / MassiveDecoy**  
3. 系统比较：AI 对接 vs 物理对接 vs AI 重打分  
4. 关键**结论可翻转**的现象（TrueDecoy 上 Glide 强，RandomDecoy 上 AI 更强）  
5. 给出可执行协议：**层级 VS**（快筛→重打分→精筛）  
6. 开放数据与代码

**映射到 Dual-VSDS：**

| VSDS-VD | Dual-VSDS（你） |
|---------|-----------------|
| TrueDecoy / RandomDecoy | dual 正例 + A-only/B-only 硬负 + RandomDecoy |
| RMSD vs PB-valid vs EF | 单靶姿态/PB + **dual-vs-single** + PR-AUC |
| 「精度≠物理合理≠富集」 | 「单端高分≠双靶兼容」 |
| 层级 VS | 两端阈值门控 → softmin/校准融合 → 可选精筛 |
| 不发明新 sampler | **同样：复用 GNINA/Vina/RTMScore** |

**应照抄的实验设计习惯：** 同一分子多协议并报表；至少两套 decoy；消融「原生分 vs 重打分」；结论写清适用边界。

---

### 1.2 Cao et al. — EquiScore

- **题名**：*Generic protein–ligand interaction scoring by integrating physical prior knowledge and data augmentation modelling*  
- **期刊**：*Nat. Mach. Intell.* **2024**, 6, 688–700  
- **DOI**：[10.1038/s42256-024-00849-z](https://doi.org/10.1038/s42256-024-00849-z)

**可参考的结构：**

1. 问题：打分泛化差、对对接姿态来源敏感  
2. 方法：物理先验 + 数据增强（零件已知，组合新）  
3. **泄漏分析 + 更严冗余去除**  
4. 外测：DEKOIS2.0 / DUD-E 等，对比十余方法  
5. **对姿态生成器不敏感** 作为卖点  
6. 类似物排序 / 可解释性作次要支撑  

**映射：** 你的「校准到 \(p(\mathrm{active})\) + softmin」≈ 他们的「让打分在真实筛选设定下可泛化」；你的 Dual-VSDS ≈ 他们的 PDBscreen/外测纪律。不必学他们再训一个大 GNN 作主创新。

---

### 1.3 Graber et al. — PDBbind CleanSplit（**打假型**）

- **题名**：*Resolving data bias improves generalization in binding affinity prediction*  
- **期刊**：*Nat. Mach. Intell.* **2025**  
- **DOI**：[10.1038/s42256-025-01124-5](https://doi.org/10.1038/s42256-025-01124-5)

**可参考的结构：**

1. 揭露 CASF/PDBbind **train–test 泄漏**导致 DL 亲和预测虚高  
2. 发布 CleanSplit 过滤协议  
3. **重训 SOTA → 性能断崖**（证明问题真实）  
4. 再展示自己的模型在严划分下仍稳  

**映射：** 你的对应句是——「在 random split / 朴素 mean-fusion 上看起来很好；换 **scaffold + dual-vs-single 硬负** 后 mean/sum **系统性推高 A-only/B-only**」。这类「先打假再给协议」是 NMI 极吃的 cover letter。

---

## 2. Tier A：评测语言与图表套路（强烈建议对齐）

| 文章 | 期刊 | 学什么 | 不学什么 |
|------|------|--------|----------|
| **PoseBusters** (Buttenschoen 2024) | Chem. Sci. | PB-valid 门控；「RMSD 合格≠物理合格」话术 | 再做一个姿态 RMSD leaderboard |
| **PoseBench** (Morehead et al.) | *Nat. Mach. Intell.* 2025 · [10.1038/s42256-025-01160-1](https://doi.org/10.1038/s42256-025-01160-1) | 统一 toolkit；多设定（apo、多配体、盲对接）；开放可复现 | 主贡献写成 co-folding 竞赛 |
| **Decoding limits of DL docking** (Li/Cao 2025) | Chem. Sci. · [10.1039/D5SC05395A](https://doi.org/10.1039/D5SC05395A) | 五维诊断表：精度/物理/相互作用/VS/泛化 | 综述口吻盖过方法贡献 |
| **PoseX** | arXiv 2505.01700 | 时间切割 + cross-docking；松弛后处理消融 | 与 Dual-VSDS 任务不同，仅作协议补充 |
| **JCIM 2025** DiffDock-L + Vina/GNINA/RTMScore | JCIM | **采样与打分解耦** 的消融表设计 | 停留在单靶 |

---

## 3. Tier B：高分但体裁不同（参考叙事边界，勿当主模板）

### 3.1 双靶/多靶 **生成 + 湿实验**（Nat. Commun.）

| 文章 | DOI | 可参考 | 为何不是主模板 |
|------|-----|--------|----------------|
| **POLYGON** (Nat. Commun. 2024) | [10.1038/s41467-024-47120-y](https://doi.org/10.1038/s41467-024-47120-y) | 公开数据识别双活性 → 多靶对生成 → docking  sanity → **合成 32 个 + 生化/细胞** | 主贡献是生成；你不做 de novo 主线；湿实验量远超你当前能力 |
| **Automated multi-target CLM** (Nat. Commun. 2024) | [10.1038/s41467-024-52060-8](https://doi.org/10.1038/s41467-024-52060-8) | 「pooled selective ligands → 生成双靶 → 对接预筛 → 合成验证」闭环 | 同上；对接只是后验检查 |

**你可借鉴的 20%：** 应用章用「公开靶点对案例 + 私有 NLRP3/JNK1 盲测排序」代替大规模合成；细胞结果作 **L2 holdout**，不宣称 binding gold standard。

### 3.2 超快 VS / 资源库（Science / NeurIPS）

| 文章 | 期刊 | 可参考 |
|------|------|--------|
| **DrugCLIP** genome-wide VS | *Science* 2025 · [10.1126/science.ads9530](https://doi.org/10.1126/science.ads9530) | 基准 + 湿实验 hit rate + **开放 GenomeScreenDB**；资源型闭环 |
| DrugCLIP (NeurIPS 2023) | NeurIPS | zero-shot VS、EF 早期富集报表 |

**边界：** 他们卖速度与覆盖；你卖 **双靶决策正确性**。可学「开放数据库/排行榜」包装 Dual-VSDS，勿转对比学习主创新。

### 3.3 方法引擎型（Nat. Comput. Sci. 等）

| 文章 | 用途 |
|------|------|
| **KarmaDock** (*Nat. Comput. Sci.* 2023) | 消融对照引擎；真实 VS 案例写法 |
| **DeepDock** (*Nat. Mach. Intell.* 2021) | 早期几何深度学习打分；历史引用 |
| **MISATO** (*Nat. Commun.* 等) | 数据资源叙事 |

---

## 4. Tier C：双靶分数融合「祖先文」（必须 cite，证明缺口）

这些不是 NMI，但是 reviewer 会问的 prior art：

| 文章 | 要点 | 你的增量怎么写 |
|------|------|----------------|
| Perez-Castillo et al. *Molecules* 2017 · [PMC5725543](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5725543/) | MAO-B / A2A **rank fusion** 提高双靶富集 | 无跨靶校准、无 A-only 硬负、无 fused/linked、无开放配对基准 |
| Ensemble docking fusion (*Molecules* 2019) | GEOM/HARM 优于 MIN/MAX | 单靶多构象融合 ≠ 跨靶兼容性 |
| Consensus docking surveys | Z-score / rank-by-rank | 零件已知；任务未重定义 |

**Cover letter 一句：** Prior dual-target docking papers fuse ranks/scores; they do not treat cross-target incommensurability, dual-vs-single hard negatives, or architecture-conditioned evaluation as a first-class ML problem with a leakage-controlled public benchmark.

---

## 5. 建议你「按周精读」的顺序（执行用）

| 顺序 | 文章 | 产出（读完应能写出） |
|------|------|----------------------|
| 1 | **VSDS-VD** | Methods 里 decoy 分层 + 表格骨架 + 层级协议图 |
| 2 | **EquiScore** | Introduction「泛化/泄漏」段 + 外测对比表 |
| 3 | **CleanSplit** | Results 第一张「打假」图：naive fusion 崩溃 |
| 4 | **PoseBusters** + VSDS-VD SI | PB-valid 门控写法 |
| 5 | **Perez-Castillo 2017** | Related Work 差距表 |
| 6 | **PoseBench**（选读） | 开放 toolkit / 多设定报告方式 |
| 7 | POLYGON（选读 Discussion） | 私有细胞数据如何当「有限前瞻验证」而不喧宾夺主 |

---

## 6. 可直接套用的文章章节骨架（对齐 VSDS-VD × CleanSplit）

```text
1. Introduction
   - Dual inhibitors matter; practice = two independent dockings + naive mean/sum
   - We show this is systematically biased toward A-only/B-only (Figure 1 punchline)
   - Contribution: dual-compatibility objective + Dual-VSDS + holdout

2. Dual-VSDS construction  (≈ VSDS-VD §dataset)
   - Target-pair selection; activity pairing rules
   - Labels: dual / A-only / B-only / inactive (never untested→inactive)
   - Decoys: TrueNegative + RandomDecoy
   - Splits: scaffold / target-pair / time; leakage audit

3. Docking & calibration protocol  (≈ EquiScore 可插拔)
   - Engine: GNINA or Vina; optional RTMScore; PoseBusters gate
   - Per-target calibration to p(active); softmin / threshold-margin fusion
   - Baselines: raw mean, min, z-mean, rank fusion

4. Results
   - 4.1 Naive fusion fails on dual-vs-single (CleanSplit-style cliff)
   - 4.2 Calibrated softmin recovers dual enrichment
   - 4.3 Ablations: calibration, hard negatives, design_type experts
   - 4.4 TrueDecoy vs RandomDecoy (VSDS-VD-style flip possible)
   - 4.5 External: public pairs + private NLRP3/JNK1 cell holdout (ranking only)

5. Discussion / Limitations
   - Cell ≠ dual binding; PK decoupled; no new sampler claim

6. Data & Code availability
```

---

## 7. 与本仓库主张的对齐检查

| NMI 成功模式 | 本课题是否具备 |
|--------------|----------------|
| A 纠偏 | 有：mean/sum → 假双靶 |
| B 任务重定义 | 有：dual-compatibility |
| C 组合方法 | 有：校准+softmin+硬负+类型条件 |
| D 开放基准 | 规划中：Dual-VSDS |
| E 湿实验闭环 | 弱但可用：细胞 holdout（勿夸大） |

**最值得模仿的两篇合在一起 = 你的 ideal paper：**  
**VSDS-VD（评测与协议美学） × CleanSplit（打假力度） × EquiScore（打分/泄漏纪律）**，主题换成双靶兼容性。

---

## 8. 链接速查

| ID | 文章 | 链接 |
|----|------|------|
| S1 | VSDS-VD | https://doi.org/10.1038/s42256-025-00993-0 |
| S2 | EquiScore | https://doi.org/10.1038/s42256-024-00849-z |
| S3 | CleanSplit | https://doi.org/10.1038/s42256-025-01124-5 |
| A1 | PoseBench NMI | https://doi.org/10.1038/s42256-025-01160-1 |
| A2 | PoseBusters | https://doi.org/10.1039/D3SC04185A |
| A3 | Decoding limits | https://doi.org/10.1039/D5SC05395A |
| B1 | POLYGON | https://doi.org/10.1038/s41467-024-47120-y |
| B2 | CLM multi-target | https://doi.org/10.1038/s41467-024-52060-8 |
| B3 | DrugCLIP *Science* | https://doi.org/10.1126/science.ads9530 |
| C1 | Perez-Castillo dual fusion | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5725543/ |
