# 方向重定：PaperSpine × ResearchStudio × Deep Research 综合报告

> **Skills 调用记录（2026-07-24）**  
> - **PaperSpine**：Problem → SOTA → Gap → Claim spine → Kill criteria  
> - **ResearchStudio**：`paper_search` + `scoop_check` + `idea_quality`（内嵌）  
> - **academic-research-skills**：`deep-research`（FINER RQ + 证据综合）/ Nature 政策锚点（不冲 NMI 默认）  
> - **scientific-agent-skills**：`medchem` 仅作面板 triaging 提醒；本轮不跑生成式设计  
> - 检索产物：`research_runs/direction_refreeze_20260724/`  
> - 对接证据：`evidence/docking_v1_metrics.json`

**一句话结论：**  
冻结方向为 **架构无关的双靶对接决策尺子（Dual-VSDS-Decision）**；乘客/moiety **出局**；默认目标刊 **JCIM / J. Cheminform.**（非 NMI）。

---

## 0. 第一次对接结果回顾（证据层）

### 0.1 Cognate QC（TAK-285 / EH40_01，3POZ / 3RCD）

| 端 | Vina Top-1 | RTM 选中 | 含义 |
|----|------------|----------|------|
| EGFR 3POZ | RMSD **9.51 Å**（翻面错构） | **1.02 Å**（原 mode 2 → #1） | **采样够、Top-1 打分不可靠** |
| HER2 3RCD | ~1.87 Å 过关 | ~1.97 Å 仍过关 | 重打分不破坏已合格端 |

共晶配体在 40 面板按 RTM mean/min 均排 **#1**。

### 0.2 四类面板排序（Dual vs A∪B）

| 规则 | AUROC | Top10 dual | Top10 硬负例 |
|------|-------|------------|--------------|
| Vina mean | **0.55** | 4 | **6** |
| RTM mean | **0.69** | 6 | 4 |
| RTM min | **0.69** | **7** | **3** |
| RTM min(z) | **0.71** | 7 | 3 |

姿态换位：最佳 RTM mode ≠ Vina mode1 — EGFR **36/40**，HER2 **20/40**。

### 0.3 硬负例分型（决定下一步协议）

| 分子 | 类 | Vina→RTM | 解读 |
|------|----|----------|------|
| **EH40_18** | A_only | #2 → **#27** | 一端虚高；**重打分即可压下** |
| **EH40_33** | B_only | mean#6 → min#19 | **弱端短板规则有效** |
| **EH40_23** | A_only | 仍 **#2**（min/z 仍 #2） | 两端 RTM 都高；**短板不够，需姿态/相互作用门控** |
| EH40_28 | B_only | 仍 Top5 | 同源激酶 ATP 骨架污染 |

**已能断言：** 朴素双端 Vina 融合近随机；ensemble/RTM 必要且有效；协议尚未闭环（顽固假阳性仍在）。  
**不能断言：** 乘客机制；活性回归（Spearman 仍弱）；可外推到任意靶对。

---

## 1. PaperSpine：问题脊骨

### Problem（可操作）

> 对靶 A、B 独立对接后，如何用**可复现决策规则**把 Dual 从 A-only/B-only 硬负例中排出来——且**不依赖**分子是 merged / linked / 其它？

失败机制（已观测）：

1. Top-1 姿态错（TAK-285 EGFR）  
2. 两靶分数尺度不可比  
3. mean 允许「一端极强掩盖另一端」  
4. 真正难负例是 **A-only/B-only**，不是随机 decoy  

### SOTA 地图（scoop 关键威胁）

| 层级 | 代表 | 已做 | 未做 |
|------|------|------|------|
| 双靶对接可行性 | Zhou 2013 JCIM | 激酶对对接；高 FP、低富集 | 无四类硬负例决策协议；无现代重打分/QC |
| 分数融合 | Pérez-Castillo 2017 | 双端多函数融合提升已知 dual 富集 | 评 decoy/已知 dual，非 Dual vs A/B-only |
| 前瞻双靶 VS | Jaiteh 2018 JMC | 异质靶对共识排名 + 实验命中 | 个案发现，非可复现尺子 |
| 多目标 VS | Fromer 2024 Pareto | Pareto 优于标量化；EGFR/IGF1R 例 | 优化采办，非姿态审计+校准决策 |
| 生成式双靶 | DualDiff / AIxFuse / CombiMOTS | 生成同时高双端分 | **不解决**筛选决策；且常标量化 Vina |
| 单靶评测生态 | VSDS-VD / RTMScore / PoseBusters / CASF | 重打分、物理合理性 | **不是**双端四类决策问题 |
| 化学分类 | Morphy；DTDL 2024 | linked/fused/merged | 设计分类 ≠ 对接决策 |

**Scoop 等级：3.5 / 5（部分拥挤，未全 scoop）**

### Gap（可写进 Intro）

文献回答的是：「对接/融合能不能富集已知 dual / 发现命中？」  
缺口是：「给定两端独立对接证据，怎样做**可复现决策**，专门惩罚单端硬负例，并跨化学架构评测？」

### DELTA（对外主张，必须窄）

> 不宣称新对接器或乘客机制。宣称：把双靶对接形式化为 **四类硬负例决策问题**，并交付 **姿态 QC → 可选重打分 → 分靶校准 → 弱端短板/门控** 的架构无关协议 + 泄漏控制基准。

### Claim spine

| ID | Claim | 现状 |
|----|-------|------|
| C1 | 朴素 Top-1 mean 近随机 / 硬负例污染重 | ✅ EGFR/HER2 40 面板 |
| C2 | ensemble/重打分改变姿态并改善诊断 | ✅ TAK-285 + mode 统计 |
| C3 | 校准 + min/shortfall/门控优于裸 mean | 🟡 min/z 略升；**EH40_23 未解** |
| C4 | 第二对公开靶外推 | ❌ 未做 |
| C5 | 开放四类基准 + 协议 YAML | ❌ 未打包 |

### Kill criteria（事先写死）

1. 公平复现 Pérez 式融合后，你的协议无增益  
2. min/shortfall/门控在 **≥3 靶对**上不优于 mean/几何/谐波/Pareto  
3. 只在 EGFR/HER2 有效  
4. 标签噪声大于方法信号  
5. 改善只来自「换更好的单靶打分」，决策规则无额外贡献（消融失败）

---

## 2. ResearchStudio idea_quality（方向选择）

比较三案（绝对分 0–100；A/C 门控生效）：

| 方向 | A 问题位 | B 方法 | C 契合 | 总分 | 判决 |
|------|---------|--------|--------|------|------|
| **A 决策尺子** | 4（四类硬负例真实且开放） | 3（组件已知，组合+基准是增量） | **5**（直接打你要的「覆盖所有构型」） | **~67 strong 边缘** | **采纳** |
| B 乘客/moiety | 2（覆盖面窄，与目标冲突） | 3 | 1 | 门控 → weak/borderline | **拒绝** |
| C 纯基准无协议 | 4 | 2（资源文） | 3 | ~50 borderline | **备用**（协议增益弱时降级） |

**Nature / NMI 适配（academic nature policy）：**  
NMI 买的是任务纠偏 + 大规模开放基准 + 强方法洞察。当前证据与增量深度 **不够默认冲 NMI**；Nature 政策亦要求 AI 使用披露、人类对内容负责——本课题应把复现脚本与协议 YAML 当一等公民，而不是「AI 生成叙事」。  
**默认目标刊：JCIM / J. Cheminform. / Digital Discovery；NMI 仅当 C3+C4+C5 全绿。**

---

## 3. 重新冻结的课题方向

### 是什么

**Dual-VSDS-Decision：架构无关的双靶对接决策尺子**

输入：任意小分子 + 两端口袋结构  
过程：对接（固定协议）→ top-K →（可选）RTM/同类重打分 → PoseBusters/相互作用门控 → 分靶校准 → 弱端决策（min / shortfall）  
输出：双靶优先级分数 + 可审计中间量  
评测：四类标签 Dual / A_only / B_only / neither；主指标 AUROC、Top-k enrichment、硬负例 top%  
报告：按 merged / linked / other **分层**，不假设单一机理

### 不是什么

- 乘客 / moiety 封面  
- 新对接采样器  
- 生成式双靶设计主线（DualDiff 等是邻域，不是本题）  
- PROTAC 三元；表型=结合证明  

### 为什么现在定这个

1. 你的约束：必须覆盖所有双靶构型 → 只能卖协议/评测，不能卖一种分子机理  
2. 第一次对接已击中真实缺口：朴素融合失败、重打分部分修复、顽固假阳性暴露协议缺口  
3. Scoop 显示融合本身不新，但 **四类硬负例决策 + 可复现尺子** 仍有空间  

---

## 4. 执行路线（按证据驱动，不再换故事）

```text
Phase 0  [DONE] 方向冻结；EGFR/HER2 cognate QC；Vina/RTM 基线
    │
Phase 1  [NOW]  协议闭环（只在 EGFR/HER2 40±扩面板）
    │           - 冻结 YAML：引擎/盒子/top-K/RTM/门控/融合
    │           - 消融：Vina mean/min | RTM mean/min | z-校准 | shortfall
    │                   | PoseBusters/PLIF 门控
    │           - 专治 EH40_23：姿态与相互作用一页纸
    │           - Go：硬负例 Top10 比例相对 Vina 下降 ≥30% 且 dual 回收不崩
    │
Phase 2  [GO后] 第二对靶外推（优先 PIK3CA/mTOR）
    │           - 同一 YAML；重复 C1–C3
    │
Phase 3         第三对可选（Mcl-1/Bcl-xL）+ 架构分层标注
    │
Phase 4         打包 Dual-VSDS-Decision（数据+脚本+协议）→ 写 JCIM 稿
```

### 本阶段唯一允许的实验问题

> 在 RTM 之上，**校准 + 弱端规则 + 姿态/相互作用门控**能否把 EH40_23 类压出前列，且不打崩 TAK-285 / lapatinib？

若否 → 降级为「诊断+基准」资源文（方向 C），仍不回乘客。

---

## 5. 与旧文档关系

| 文档 | 状态 |
|------|------|
| `NMI_SUBMISSION_PLAN_DECISION_RULER.md` | 现行投稿骨架（与本报告对齐） |
| `PROJECT_MASTER_PLAN.md` | 现行总览 |
| `CRITIQUE_AND_NEXT_STEPS.md` | 红队（乘客前置条件已废） |
| `NMI_SUBMISSION_PLAN_MOIETY.md` | 废弃 |
| 本文件 | **方向重定权威记录** |

---

*本报告由对接实证 + scoop_check + PaperSpine 缺口图共同约束；后续实验只服务 C3→C4，禁止再开新封面故事。*
