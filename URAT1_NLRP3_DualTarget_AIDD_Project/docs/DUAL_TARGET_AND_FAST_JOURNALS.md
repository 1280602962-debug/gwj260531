# 双靶叙事 + 快速期刊选择（JMGM / 非 OA）

> 补充：`TWO_PAPER_STRATEGY.md`（双轨总策略）  
> 问题来源：双靶是否完全放弃？JMGM 近半年标准？非 OA 且审稿快的期刊？

---

## 一、双靶不是「完全不能提」，但不能按原 TAPE-GATE 方式提

### 1.1 什么算「创新」、什么不算（对照你现有数据）

| 叙事 | 有无数据支撑 | 能否作为创新点 |
|------|--------------|----------------|
| 痛风「代谢 URAT1 + 炎症 NLRP3」疾病网络 | ✅ 文献 + 临床逻辑 | ✅ **Introduction / Discussion 动机** |
| 对两靶 **代表抑制剂** 做平行结构计算 | ⚠️ URAT1 有；NLRP3 **需补 2 药对接+MD** | ✅ **可写成快速线主角度**（见下） |
| TAPE-GATE / Pareto 融合 / 双路径漏斗 | ❌ 无回收结果 | ❌ 不能写 |
| 8973 蒸馏 / Teacher M-CPDL | ❌ 未跑通 | ❌ 不能写 |
| 「发现双靶新 hit」 | ❌ 无湿实验、无双靶化合物列表 | ❌ 不能写 |
| 「新 AI 双靶算法」 | ❌ 0 重叠 SMILES；URAT1 ML 2/4 | ❌ 不能写 |
| ChEMBL 0 重叠导致双靶 ML 困难 | ✅ 可量化 | ✅ **Discussion 方法学观察**（不是算法贡献） |

**结论**：双靶可以保留为 **疾病网络驱动的 paired computational pharmacology**（双节点平行表征），**不是**「双靶 AI 发现算法」。

### 1.2 可诚实的「双靶创新点」（快速线最多写 2 条）

1. **疾病网络框架**：在痛风代谢–炎症双节点背景下，**并行**表征临床 URAT1 抑制剂（9DKB）与 NLRP3 抑制剂（8ETR/7ALV）的结合模式与动力学（**各做各的**，不假装有统一漏斗）。
2. **靶标类别对比**：转运体 inward 结构（9DKB）与 NACHT 变构口袋（8ETR）在 **对接验证、MD 稳定性、可外推性** 上的差异——这是计算药理学观察，不是新软件。

**不能写**：「首次提出双靶虚拟筛选平台」「dual-target AIDD」「融合排序优于单靶」。

### 1.3 快速线升级方案：论文 A′（双靶平行版）

在原 **论文 A**（仅 URAT1 @ 9DKB）基础上，**增加约 1–2 周计算** 即可合为一篇：

| 部分 | 内容 | 数据状态 |
|------|------|----------|
| **URAT1** | 四药 @ 9DKB + redock + MD + MM-GBSA | ✅ 对接已有 |
| **NLRP3** | MCC950 @ 7ALV；GDC-2394 @ 8ETR redock + MD | ❌ **需补** |
| **ML** | URAT1 2/4 → **不写**；NLRP3 AUROC 0.89 → 仅 **SI 一句**（训练集内 benchmark） | 弱化 |
| **双靶融合** | **不写** | — |

**推荐标题方向**：

*Computational pharmacology of URAT1 and NLRP3 inhibitors in a metabolic-inflammatory dual-target framework for gout: structure-based docking and molecular dynamics studies*

**与论文 B（JCIM）分工**：
- A′：双靶 **疾病动机 + 两节平行结构计算**（单态 9DKB + NLRP3 单结构）
- B：URAT1 **三态 benchmark**（刚性失败、decoy）——**不再重复 NLRP3**

---

## 二、JMGM 近半年发文画像与门槛（2025-12 — 2026-06）

**期刊**：*Journal of Molecular Graphics and Modelling*（Elsevier，ISSN 1093-3263）  
**指标**：IF ≈ 3.0；CiteScore ≈ 5.9；**官方录用率约 13%**（较严）

### 2.1 官方时间线（Elsevier Insights，订阅路线无 APC）

| 阶段 | 官方中位/典型 |
|------|----------------|
| 投稿 → 首次决定 | **5 天** |
| 投稿 → 审稿后决定 | **51 天** |
| 投稿 → 录用 | **101 天**（约 **3.4 个月**） |
| 录用 → 在线发表 | **3 天** |

**结论**：JMGM **首决快，但录用并不快**；比 JCAMD / J. Mol. Modeling 整体周期偏长。

### 2.2 出版模式

- **Subscription（非 OA）**：**不向作者收 APC**，文章订阅可读——符合你「非 OA」要求。
- OA 可选，APC 约 **USD 3,790**——投稿时选 **Subscription** 即可。

### 2.3 编辑部明确拒稿标准（Guide for Authors）

> *Routine applications of standard modelling approaches, providing only very limited new scientific insight, will not meet our criteria for publication.*

即：**常规 Glide + MD 四药表征**，若无明确新见解，JMGM **可能直接拒**。

### 2.4 近半年典型内容（Vol 137–138, 2025；Vol 145–146, 2026 在刊）

| 类型 | 代表 | 特点 |
|------|------|------|
| 靶点抑制剂设计 | VEGFR2 oxindole（Vol 138）；EGFR L858R/T790M/C797S | 对接 + MD + **DFT** |
| 药效团/QSAR | PAK1 3D-QSAR + pharmacophore（2026 in press） | 有新化学系列 |
| 方法/ML | sQC ChemBERTa+ECFP VS（Vol 146）；FocusLG 激酶亲和力 | **新模型 + 基准** |
| 实验+计算 | Transthyretin mangiferin（Vol 137） | 对接 + **体外验证** |
| 物理化学 | DFT、MXene、超分子 | 非药化主战场 |

**JMGM 上少见**：纯「网络药理学双靶」无实验、无新分子系列；双靶 MTDL 近期更多见于 *Pharmaceuticals* 等（QSAR 设计 16 类似物 + 双靶对接 + MD）。

### 2.5 你的工作投 JMGM 的匹配度（诚实）

| 优势 | 劣势 |
|------|------|
| 9DKB 新 cryo-EM；URAT1 转运体 | 仅 4+2 已知药物，易被视作 routine |
| 可加双靶疾病框架 | 无 DFT、无新类似物、无湿实验 |
| redock + MD 规范 | 官方 13% 录用率 |

**建议**：JMGM **可作冲刺目标**，但 **不要作为「最快录用」首选**；若投，需在 Discussion 强化 **转运体 vs 炎症小分子的结构药理学差异**，并避免「我们又做了一遍 Glide」的语气。

---

## 三、非 OA + 审稿较快期刊排序（针对你的快速线）

以下均为 **可选 Subscription、不收 APC**（投稿时勿选 Open Access）。时间来自官网/期刊页，实际因审稿人而异。

### 第一梯队：优先投（快 + 接受案例研究）

| 优先级 | 期刊 | 首决（官方） | 录用周期（参考） | IF≈ | 为何适合你 |
|--------|------|--------------|------------------|-----|------------|
| **1** | **Chemical Biology & Drug Design** (Wiley) | 宣传 **5–30 天** | 审稿文化偏快，少催补实验 | 3.3 | 对接+MD+疾病故事；**明确可选非 OA** |
| **2** | **Journal of Molecular Modeling** (Springer) | **4 天** | LetPub 约 **2.3 月** | 2.5 | **默认非 OA**；MD 案例友好；录用相对容易 |
| **3** | **J. Comput.-Aided Mol. Des.** (Springer) | **10 天** | 约 **2–3 月** | 3.0 | 案例+redock；2026 刊文重验证 |

### 第二梯队：可投但周期或门槛更高

| 期刊 | 首决 | 录用（官方） | 备注 |
|------|------|--------------|------|
| **Molecular Diversity** | — | 常 **3–5 月** | 非 OA 可选；爱 MD+ADMET；无实验时要收敛 claim |
| **JMGM** | 5 天 | **101 天** | 非 OA 可选；**拒 routine**；录用率低 |
| **Computational Biology and Chemistry** | 42 天 | 84 天+ | 明示拒「不成熟建模」；需 MD 支撑结论 |

### 第三梯队：不建议为「快」而投

- 纯 OA 期刊、高额 APC 期刊  
- 新办期刊（如部分 2025 年新刊）：审稿可能快，但认可度风险大  

### 推荐投稿顺序（务实）

```
方案 1（最稳最快）：
  J. Mol. Modeling → 若拒 → Chemical Biology & Drug Design → 若拒 → JCAMD

方案 2（你想冲略高 IF）：
  Chemical Biology & Drug Design → JCAMD → J. Mol. Modeling

方案 3（坚持 JMGM）：
  仅在 MD+双靶 Discussion 写扎实后投；接受 3–4 月+ 周期和较高拒稿风险
```

**JMGM 不适合作为「尽快发表」的第一选择**，只适合作为 IF 与品牌略高时的 **第二轮**。

---

## 四、三条线最终分工（含双靶）

| 线 | 期刊档 | 双靶怎么出现 | 预计周期 |
|----|--------|--------------|----------|
| **A′ 快速** | J. Mol. Modeling / CBDD / JCAMD | Intro 疾病网络 + Results **两节平行**（URAT1 + NLRP3） | **2–4 月** |
| **B JCIM** | JCIM / J. Cheminformatics | **不写双靶**；仅 URAT1 三态 benchmark | **6–12 月** |
| **搁置** | — | TAPE-GATE 融合、8973 蒸馏 | Gate 后 |

---

## 五、快速线 A′ 还需补什么（双靶版）

### URAT1（已有基础）
- [x] 9DKB 四药 SP→XP pose  
- [ ] lesinurad redock RMSD  
- [ ] MD 100 ns × 4 + MM-GBSA  

### NLRP3（补双靶平行节）
- [ ] GDC-2394 @ **8ETR** redock（RMSD ≤ 2 Å）  
- [ ] MCC950 @ **7ALV**（共晶为 MCC950 类类似物，Methods 中说明）对接 + 与文献药效团比较  
- [ ] MD 50–100 ns × 2 + MM-GBSA  
- [ ] **不要**写 NLRP3 ML 为主结果（MCC950/GDC 在训练集）

### 写作
- [ ] Abstract 用 *computational pharmacology* / *paired characterization*  
- [ ] 不出现 fusion、Teacher、virtual screening of millions  
- [ ] Discussion 一段：ChEMBL **0 SMILES 重叠** → 双靶 ML 需分靶证据后再整合（引向未来工作）

---

## 六、措辞对照（双靶专用）

| 说法 | 可否 |
|------|------|
| 「痛风代谢–炎症双靶治疗策略的计算药理学评估」 | ✅ |
| 「并行分析 URAT1 与 NLRP3 临床抑制剂结合模式」 | ✅（补 NLRP3 计算后） |
| 「提出双靶 AI 发现算法」 | ❌ |
| 「虚拟筛选鉴定双靶先导化合物」 | ❌ |
| 「ML 模型成功支持双靶筛选」 | ❌ |
| 「为双靶先导发现奠定基础」 | ⚠️ 仅 Discussion 未来时态 |

---

## 七、相关文档

| 文件 | 内容 |
|------|------|
| `MANUSCRIPT_OUTLINE_FAST_9DKB.md` | 原单靶 URAT1 快速大纲（可扩展为 A′） |
| `TWO_PAPER_STRATEGY.md` | 双轨总表 |
| `config/docking_ensemble.yaml` | NLRP3：7ALV、8ETR |
