# JCIM 可行性复核：先前工作量估算哪里错了，以及现在怎么走

> 复核对象：会话中给出的「冲 JCIM 需补 2–3 对新靶 ×100 ×2 端 + 统一 prep + GNINA + 打包，约 400–600 次对接」这一估算。  
> 数据池：**Exploration**（全部由已冻结分数表与已抓取 ChEMBL 字典重算，**零新对接**）  
> 复现脚本：[`../data/jcim_feasibility_v0/scripts/`](../data/jcim_feasibility_v0/scripts/)  
> 日期：2026-07-29

---

## 0. 一句话结论

先前估算**方向对、数量级也大致对，但把瓶颈找错了**：真正的限制不是对接算力，而是 **公开数据里 strict 硬负例的供给**。同时它漏了一件已经免费完成的事（EGFR 统一 prep），也漏了一个新增负担（PIK3CA/mTOR 同样是 LigPrep）。

---

## 1. 两个新事实（本次核算得到）

### 1.1 EGFR 统一 prep 其实不用再对接，而且结论没变

M4-min 已把旧 40 用 RDKit+meeko 重跑，新 70 本来就是 RDKit → **全 RDKit 的 EH110 可直接拼出来，零新对接**。重算 z 后：

| arm | family | D/A | D/B | `min(D/A,D/B)` |
|-----|--------|-----|-----|----------------|
| `vina_mean` | docking | 0.680 | 0.282 | **0.282** |
| `vina_min` | docking | 0.648 | 0.253 | 0.253 |
| `rtm_min` | docking | 0.530 | 0.256 | 0.256 |
| `rtm_min_z` | docking | 0.533 | 0.253 | 0.253 |
| `heavy_atoms` | baseline | 0.700 | 0.369 | 0.369 |
| `MW` | baseline | 0.648 | 0.416 | 0.416 |
| **`cLogP`** | baseline | 0.526 | 0.482 | **0.482** |
| `TPSA` | baseline | 0.703 | 0.427 | 0.427 |

**结论：** 在完全统一 prep 下，EGFR/HER2 上**所有**对接臂仍全部低于最好的非对接基线（cLogP 0.482）。  
→ M3 的 EGFR No-Go **不是 prep 造成的**；这条 Track A 头条现在是 prep-clean 的，可以直接写进论文，不需要 M4-full 再花算力。

脚本/表：`data/jcim_feasibility_v0/scripts/unified_prep_eh110.py`、`tables/eh110_unified_prep_*.csv`

### 1.2 真正的瓶颈：strict 硬负例在公开数据里极度稀缺

用已抓取的 ChEMBL pChEMBL 字典，对 12 个候选对统计 strict 标签（dual 双端 ≥6.5；A_only A≥6.5 且 B≤5.5；B_only 反之）：

| pair | 两端都测 | θ=6 dual/A/B | **strict** dual/A/B | gray | 可支撑 strict 面板 |
|------|----------|--------------|---------------------|------|--------------------|
| PIK3CA/mTOR | 2713 | 2002/266/236 | 1552/**80**/**81** | 0.34 | **是** |
| AChE/BChE | 2537 | 986/483/225 | 687/**189**/**78** | 0.42 | **是** |
| PIK3CA/PIK3CB | 1990 | 988/299/213 | 602/**56**/**67** | 0.54 | **是**（同工酶，叙事偏近） |
| Mcl-1/Bcl-xL | 305 | 82/77/24 | 39/41/**12** | 0.45 | 否 |
| **EGFR/HER2** | 1751 | 1182/207/46 | 951/39/**7** | 0.31 | **否** |
| Mcl-1/Bcl-2 | 371 | 160/87/34 | 65/24/**2** | 0.64 | 否 |
| AKT1/p70S6K | 601 | 453/11/100 | 373/**1**/35 | 0.30 | 否 |
| JAK2/HDAC1、BRD4/HDAC1、BRD4/HDAC6、PARP1/MET、CDK6/BRD4 | 6–82 | — | 硬负 ≤3 | — | 否 |

**门槛设为「两侧 strict 硬负各 ≥50」时，12 对里只有 3 对通过。**

关键含义：

1. **EGFR/HER2 全库只有 7 个 strict B_only** —— 这正是 M2 里 `B_only_strict n=7` 的来源。也就是说，EGFR 面板**无论加多少对接**都不可能支撑 strict 四类分析；M2=Weak 有一部分是**数据供给上限**，不是采样不够。  
2. 「再加 2–3 对厚硬负的新靶」这个前提，在已审计集合里**最多只剩 2 对**（AChE/BChE、PIK3CA/PIK3CB），且其中一对是同工酶。  
3. 想要更多对靶，必须**先扩大 ChEMBL 层面的候选对审计**（零对接），而不是先买算力。

脚本/表：`scripts/audit_strict_label_supply.py`、`tables/strict_label_supply.csv`

---

## 2. 先前估算逐条复核

| 先前说法 | 判决 | 更正 |
|----------|------|------|
| 需要补 2–3 对新靶 | **部分成立** | 已审计集合中只有 2 对合格（含 1 对同工酶）；需先扩审计再定 K |
| 「EH110 统一 prep」是待做工作 | **错** | 已可零成本拼出，且结论不变（见 §1.1） |
| GNINA 属 P1（加分项） | **错，应为 P0** | 「你只是打分函数差」是对核心 claim 的第一质疑；且 GNINA 以 **rescore 模式复用现有姿态**，几乎不增加采样成本 |
| 约 400–600 次对接是主要负担 | **数量级对，定位错** | 对接便宜；真正贵的是 **每对靶的受体选择 + cognate QC + 面板策展**，以及**数据供给根本不足** |
| 未提 PIK3CA/mTOR 的 prep | **漏项** | panel48 也是 **LigPrep**；跨对靶要统一 prep，需 48×2≈96 次重跑 |
| 每对 N≈80–120 | **需修正** | Track A 的估计量是「现象能否跨对靶复现」，不是单对显著性；应**多对靶 × 中等 N**，并在建面板时就按 strict 定额抽样（而非先按 θ=6 建再重贴标签，否则 gray 必然 30–60%） |
| Track A 两对靶不够 JCIM Article | **成立** | 保留该判断 |

---

## 3. 现有结果已经足够支撑的 claim（无需新对接）

1. **池化指标会自我抵消**：EGFR D/A 0.680 与 D/B 0.282 平均成 ~0.52 的假随机（M1；统一 prep 后仍然）。  
2. **EGFR 上对接不敌平凡基线**，且 **prep 已排除**（§1.1）。  
3. **PIK3CA/mTOR 上对接有真实信号**（D/A 0.698、D/B 0.597，体积基线 ~0.46）→ 结论按对靶分化，不可平均。  
4. **RTM 绝对表现强依赖配体准备**（M4：D/A 0.80→0.66，D/B 0.61→0.47）→ 协议层要求。  
5. **标签任务原则上可辨**（oracle @σ=0.5 两端 >0.90），失败是端不对称/分数问题，不是「标签全噪声」。  
6. **四类任务存在公开数据供给上限**（§1.2）→ 这是**新的、可写进论文的结构性发现**。

---

## 4. 一条现实的 JCIM 路线

### 4.1 论文形态

**评测 + 基准型 Article（非方法型，不需要赢）**

> How the four-class dual-target docking decision task should be evaluated — and why current practice hides its failures.

| 贡献 | 内容 | 成本 |
|------|------|------|
| **C1** | 形式化四类双靶决策任务；证明池化 AUROC 是两个方向 AUROC 的加权平均、可相互抵消（代数 + 实证） | 0 |
| **C2** | **公开数据供给审计**：扩到 ≥50 对候选靶，量化 strict 硬负供给与 gray 比例；结论「该任务在公开数据上是供给受限的，仅少数对靶可支撑严格评测」 | 0（ChEMBL API） |
| **C3** | 在 **K=4 对已对接靶**（统一 prep、**≥3 打分通道**：Vina / GNINA-CNN / RTMScore）上报告方向分解 + **必报平凡基线** | 中 |
| **C4** | prep 敏感性作为协议要求（M4 证据） | 已有 |
| **C5** | 发布 **DualFourClass-Bench**：面板、标签规则、分数、脚本、方向分解评测代码 | 低 |

这个形态**不要求方法胜出**，因此不受 Stage M 的 Weak 阻断；且正面回应 JCIM「不接受单靶对接应用且无实验验证」的条款——它不是应用文，是任务定义 + 基准 + 评测方法学。

### 4.2 最小充分规模

| 项 | 目标 |
|----|------|
| 审计对靶（零对接） | **≥50** |
| 已对接对靶 | **4**：PIK3CA/mTOR、EGFR/HER2、AChE/BChE、+1（PIK3CA/PIK3CB 或 Mcl-1/Bcl-xL 作「薄硬负 / 异质折叠」案例） |
| 每对 N | ~100–120，**按 strict 定额抽样** |
| 打分通道 | Vina + GNINA(CNN rescore) + RTMScore |
| 平凡基线 | heavy_atoms / MW / cLogP / TPSA / Morgan-to-dual |
| 统计 | 每对方向分解 AUROC + 配对 bootstrap CI；跨对靶森林图 |
| prep | 全库 RDKit ETKDG + meeko |

### 4.3 新增对接预算（诚实计数）

| 任务 | 新 Vina 采样 |
|------|--------------|
| EGFR/HER2 统一 prep | **0**（已完成） |
| PIK3CA/mTOR 转 RDKit（48） | ~96 |
| PIK3CA/mTOR 扩至 strict 定额 ~110 | ~220（若扩） |
| AChE/BChE 新面板 ~110 | ~220 |
| 第 4 对 ~110 | ~220 |
| GNINA CNN **rescore**（复用姿态） | 0 新采样（~880 次重打分） |
| **合计** | **约 540–760 次 Vina 作业** |

对接本身不重；**主要工程量在 2–3 套新受体的结构冻结 + cognate QC + 面板策展**。

### 4.4 执行顺序（先零成本，再买算力）

```
J0  扩大 ChEMBL 供给审计（≥50 对）           零对接  ← **完成** `data/jcim_j0j1_v0/`（49 对已审计；API 宕机时 fetch queue=22）
J1  按 strict 定额 + holo 结构可得性选定 K=4  零对接  ← **完成草案**（待用户批准才对接）
J2  统一 prep 重跑 PM48（+ 可选扩面板）        小     ← **需用户批准**
J3  新受体结构冻结 + cognate QC（2–3 套）      中（判断密集）
J4  新对靶面板对接 + RTM                       中
J5  GNINA rescore 全部对靶                     小（复用姿态）
J6  方向分解 + 基线 + CI + 森林图 + 打包发布    零对接
```

**J0+J1+Track A（零对接）已完成。**  
**冲 JCIM 白话下一步：** [`JCIM_NEXT_ROUTE_PLAIN.md`](JCIM_NEXT_ROUTE_PLAIN.md)  
**对接阶段命令：** [`AGENT_COMMAND_JCIM_DOCKING_PHASE.md`](AGENT_COMMAND_JCIM_DOCKING_PHASE.md) — **本轮已执行**。  
下一本地对接命令（**勿自行执行**）：用户批准后的 J2 PM48→RDKit + J3 新受体冻结。

任一步 No-Go 都能退到 **Mol. Inf. / JCAMD 版**（现有 2 对 + Stage M 已够），不至于全盘落空。

---

## 5. 如果不做 J0–J5

**不要硬投 JCIM Article。** 现有 2 对靶 + Stage M 的正确归宿是 **Molecular Informatics / JCAMD**（hybrid，非 OA），叙事仍用测量洞见 + 供给审计；将来补足 K=4 与多通道后再写 JCIM 基准文，并引用第一篇。

---

## 6. 与其他文档的关系

| 文档 | 关系 |
|------|------|
| `PLAN_V2_REDTEAM_AND_REDESIGN.md` | Stage M 起源；本文补「供给受限」这一层 |
| `data/stage_m_v0/analysis/STAGE_M_VERDICT.md` | Track B=Weak 仍有效；本文说明 Track A 路线不受其阻断 |
| `JCIM_GAP_TWO_SCENARIOS.md` | 规模锚点仍可用；本文取代其「先扩 EGFR」式建议 |
| `PUBLIC_TARGET_PAIR_SELECTION_REPORT.md` | 其 θ=6 审计需按 strict 标签重做（J0） |
