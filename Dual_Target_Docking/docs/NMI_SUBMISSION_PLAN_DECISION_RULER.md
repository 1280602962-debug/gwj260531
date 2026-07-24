# 投稿规划：架构无关的双靶对接决策尺子

> **主线定调（2026-07-24 冻结）**  
> 双靶配体化学架构多样（merged / linked / 其它）。课题不假设「乘客药效团」机制，而提出一套**覆盖所有架构**的双靶虚拟筛选**决策协议**：姿态可审计 → 可选重打分 → 分靶校准 → 短板/门控融合，并在四类标签上证明优于朴素双端融合。  
> 乘客 / moiety 不进主 claim；linked 子集最多作分层附录。

总览见 [`PROJECT_MASTER_PLAN.md`](PROJECT_MASTER_PLAN.md)。  
旧稿 [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md) 已废弃。

---

## 1. 要解决的问题（Problem statement）

**实践痛点：** 做双靶 VS 时，人们常对靶 A、靶 B 各跑一遍对接，再把分数 mean/min/乘积拼起来。  
**已观测失败：** 在 EGFR/HER2 40 分子四类面板上，Vina Top-1 mean 的 Dual vs A/B-only AUROC ≈ 0.55（近随机）；RTM 重打分后升至 ≈ 0.69，但仍有顽固硬负例（如 A-only 冲到 #2）。

**科学问题：**

> 在不依赖特定双靶分子架构假设的前提下，怎样定义并评测一套**可复现的双端决策尺子**，使真正的 dual 相对 A-only/B-only 可分？

**卖点结构：** 任务纠偏 + 协议 + 开放基准；**不是**新采样器，**不是**单一机理故事。

| 项 | 内容 |
|----|------|
| 默认目标刊 | **JCIM / Briefings in Bioinformatics** |
| NMI 仅当 | 三对靶 + 强诊断 + 开放基准同时成立 |
| 一句话 | Naive dual-end docking fusion fails to rank dual binders over single-target hard negatives; an architecture-agnostic decision protocol (ensemble/rescoring + calibration + shortfall) recovers discrimination on leakage-controlled public pairs |

---

## 2. Claim ladder（只承诺能支撑的）

| 层级 | Claim | 证据 |
|------|-------|------|
| C1 诊断 | 朴素双端 Top-1 融合在四类分子上失败 / 近随机 | EGFR/HER2 面板 + Vina 指标（已有） |
| C2 协议必要性 | 姿态 ensemble / 重打分改变排序且改善 cognate QC | TAK-285 RMSD；mode≠1 比例；RTM AUROC（已有雏形） |
| C3 决策规则 | 校准 + shortfall/min/门控优于裸 mean；降低硬负例 top% | 待补：EH40_23 类消融 |
| C4 外推 | 同一协议在第二对公开靶上仍优于朴素融合 | 待做：PIK3CA/mTOR |
| C5 资源 | 公开基准 + 协议 YAML + 脚本 | 待打包 |
| C6（附录，非必须） | linked 子集上 moiety 是否额外有用 | 可选；失败也不伤主文 |

**不做的 claim：** 乘客是根因；moiety 适用于所有双靶；新对接引擎；协同；PROTAC；表型=结合。

---

## 3. 为什么「能涵盖所有双靶」

| 做法 | 作用 |
|------|------|
| 输入只要求「分子 + 两端口袋」 | 不依赖是否可切成两个药效团 |
| 主指标按 **四类活性标签** | Dual / A_only / B_only / neither，与化学架构正交 |
| **架构分层报告**（merged / linked / unknown） | 覆盖全体，同时诚实展示差异，不把一种机理说成全体 |
| 决策函数定义在分数上 | 对任何对接引擎可插拔 |

「涵盖」= **同一协议评所有人**；不是「一种生物学解释解释所有人」。

---

## 4. 实验骨架

1. **姿态 QC**：共晶配体自对接；Top-1 与 top-K；必要时 RTM/GNINA 重打分  
2. **四类面板**：两端都有实测；θ=6；强制含硬负例  
3. **基线**：Vina Top-1 mean / min  
4. **协议臂**：top-K 重打分 → 分靶校准 → shortfall / 门控  
5. **主读数（冻结）**：AUROC Dual vs A∪B；Top10 enrichment；硬负例进入 Top10 的比例  
6. **第二对靶重复** 后才升期刊预期

---

## 5. 与现有结果的对齐

| 结果 | 对主线的含义 |
|------|----------------|
| TAK-285：EGFR Vina Top-1 失败，RTM 救回 | C2：协议必须含姿态审计/重打分 |
| Vina mean AUROC 0.55 | C1：朴素融合不够 |
| RTM mean/min AUROC ~0.69 | 协议方向对，尚未完成 C3 |
| EH40_23 仍 #2 | 必须做短板/门控，否则尺子不可交付 |
| 与 pChEMBL Spearman 弱 | 主文卖排序决策，不卖活性回归 |

---

## 6. Go / No-Go

| 条件 | 动作 |
|------|------|
| C3：硬负例 top% 相对 Vina 明显下降，且 dual 回收不崩 | 扩第二对靶 |
| C4：第二对仍优于朴素融合 | 打包基准；默认 JCIM/Brief Bioinform |
| 仅 EGFR/HER2 有效 | 缩小为「同源激酶双靶协议笔记」，不冲高刊 |
| 任何时候又想回乘客封面 | **拒绝**；与「涵盖所有双靶」目标冲突 |

---

*本文替代 moiety 投稿规划，作为现行主张。*
