# 课题总览：架构无关的双靶对接决策尺子

> **先读：** [`RESEARCH_DIRECTION_REFREEZE.md`](RESEARCH_DIRECTION_REFREEZE.md)（PaperSpine × ResearchStudio × 第一次对接证据）  
> **后续实验总规划（现行权威）：** [`EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md`](EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md)  
> **历史规划（v1，已废止权威）：** [`EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md`](EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md)  
> **阶段1命令：** [`AGENT_COMMAND_STAGE1_EGFR_EXPAND.md`](AGENT_COMMAND_STAGE1_EGFR_EXPAND.md)  
> **红队：** [`CRITIQUE_AND_NEXT_STEPS.md`](CRITIQUE_AND_NEXT_STEPS.md)  
> **投稿骨架：** [`NMI_SUBMISSION_PLAN_DECISION_RULER.md`](NMI_SUBMISSION_PLAN_DECISION_RULER.md)  
> **诊断操作：** [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md)

**方向已冻结（2026-07-24）：架构无关双靶对接决策尺子；不做乘客 / moiety 主线。**  
默认目标刊 **JCIM / J. Cheminform.**（非 NMI）。双靶分子结构种类很多（merged / linked / 其它），课题必须用**同一套决策协议**覆盖。

---

## 1. 到底要解决什么问题

### 一句话

> 给定一对靶点，如何用**可复现的对接 + 重打分 + 双端决策规则**，把真正的双靶活性分子从「只对一端强、或两端都虚高」的分子里排出来？

### 问题拆开

| 层 | 问题 | 不是什么 |
|----|------|----------|
| **任务** | Dual vs A-only / B-only / neither 的**排序与筛选决策** | 不是新对接采样器 |
| **失败模式** | 两端独立对接后做 mean/min 朴素融合，在同源激酶上接近随机（本面板 AUROC≈0.55） | 不是「发现协同」或「预测细胞表型」 |
| **覆盖面** | 同一协议适用于多种双靶化学架构；报告时**按架构分层**，不假设单一机制 | 不是只服务 linked / 双药效团分子 |
| **交付** | 公开「双靶 VS 决策尺子」：协议 YAML + 四类标签基准 + 复现脚本 | 不是私有靶点专属黑盒打分 |

### 为什么这个题成立

主流双靶（EGFR/HER2、PI3K/mTOR 等）大量是 **merged / 紧凑 ATP 双抑制剂**，不是清一色 linker 双头。  
因此「乘客药效团污染」**不能**当封面故事——它最多是 linked 子集的附录假说。  
真正跨构型共通、且你已用数据碰到的，是：

**整分子双端对接分数 ≠ 双靶活性排序；需要可审计的决策协议（姿态 QC → 可选重打分 → 校准 → 短板/门控融合）。**

---

## 2. 是什么 / 不是什么

### 是什么（交付）

1. **诊断**：证明朴素双端融合在四类分子上失败（已有 EGFR/HER2 40 面板 + Vina/RTM 证据）
2. **协议**：固定对接 → top-K →（可选）RTM/同类重打分 → 分靶校准 → shortfall / min 决策
3. **基准**：多公开靶点对、四类标签、架构分层报告（Dual-VSDS-Decision 或同等命名）
4. **外推**：第二对靶（优先 PIK3CA/mTOR）重复同一协议

### 不是什么

| 不做 | 原因 |
|------|------|
| 乘客 / moiety 主线 | 无法覆盖所有双靶构型；你已明确放弃 |
| 新 DiffDock / 自研采样器 | 创新在决策与评测，不在采样 |
| 分数融合当唯一卖点 | 融合是组件；卖的是「失败诊断 + 完整尺子」 |
| 药物联用协同、DTI GNN | 任务不同 |
| PROTAC 三元主线 | 另表；协议层可注明 out of scope |
| 细胞表型 = 双靶结合证明 | 标签语义错误 |

---

## 3. 公开靶点对（已冻结）

见 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) / `FROZEN_PUBLIC_PAIRS.yaml`：

1. **PIK3CA / mTOR** — 主规模、多为 merged 双 ATP；测尺子外推
2. **EGFR / HER2** — 姿态金标准（TAK-285）+ 第一张失败/纠偏诊断表
3. **Mcl-1 / Bcl-xL** — 异质口袋；可选 linked 子集作**附录分层**，非主线

标签：pChEMBL ≥ 6；未测 ≠ inactive；主评测必须含 A-only / B-only。

---

## 4. 执行顺序（与乘客脱钩）

1. **冻结读数**：并列报告 `vina_mean` **与** `rtm_min_z`（禁止只报 RTM）  
2. **对靶特异 E**：EGFR/HER2 `exhaustiveness_v0_1=8`；PIK3CA/mTOR `=16`  
3. **失败分型 + 化学型警告层**（flags 不改分）：T2/T5 写入 Limitations  
4. **冻结阈值决策消融**：shortfall / consensus 若无法同时压硬负并保护误伤 dual → 诚实写边界（已见 panel48）  
5. **第二对靶证据已在**；扩面板 / 第三对排在分型与协议回写之后  
6. linked / moiety 仅当附录假说，**不进主 claim**

**禁止：** 为打掉 T2 假阳性拧 clash；宣称 C4 已成功外推；乘客主线回归。

---

## 5. 相关文档

| 文档 | 角色 |
|------|------|
| [`CRITIQUE_AND_NEXT_STEPS.md`](CRITIQUE_AND_NEXT_STEPS.md) | 红队 + 下一步（已改决策尺子） |
| [`NMI_SUBMISSION_PLAN_DECISION_RULER.md`](NMI_SUBMISSION_PLAN_DECISION_RULER.md) | 现行投稿主张 |
| [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md) | 诊断操作（whole-mol 协议；moiety 降为可选附录） |
| [`DUAL_TARGET_SCORING_IMPLEMENTATION.md`](DUAL_TARGET_SCORING_IMPLEMENTATION.md) | 校准/短板组件 |
| `NMI_SUBMISSION_PLAN_MOIETY.md` | **已废弃主线**，仅作历史 |

---

*主线：架构无关的双靶对接决策尺子。乘客路线正式退出封面。*
