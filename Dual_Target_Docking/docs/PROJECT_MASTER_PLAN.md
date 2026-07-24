# 课题总览：Moiety-resolved 双靶对接评测

> **现行唯一主线：** [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md)  
> **第一张诊断表：** [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md)

---

## 1. 是什么

单分子双靶含两个药效团。独立对接 + **整分子**打分会把「乘客」药效团算进单口袋分数，双靶排序（Dual vs A-only/B-only）会系统性出错。

本课题交付：

1. **诊断**：证明整分子融合会抬高假双靶  
2. **改法**：moiety 计分 + 分靶校准 + 短板融合（同一姿态上改决策，不换采样器）  
3. **基准**：Dual-VSDS-Moiety（公开多靶点对、四类标签、可复现协议）  
4. **案例**：NLRP3/JNK1 外部锚点后 holdout（非循环）

## 2. 不做什么

| 不做 | 原因 |
|------|------|
| 新 DiffDock / 自研采样器 | 创新不在采样 |
| 分数融合当封面创新 | 降为基线对照 |
| 药物联用协同、DTI GNN | 任务不同 |
| PROTAC 三元 | 另表 |
| 细胞表型 = 双靶结合证明 | 标签语义错误 |

## 3. 公开靶点对（已冻结）

见 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) / `data/public_pair_selection/FROZEN_PUBLIC_PAIRS.yaml`：

1. PIK3CA / mTOR — 主规模  
2. EGFR / HER2 — 诊断 demo + 姿态金标准（TAK-285）  
3. Mcl-1 / Bcl-xL — 异质口袋外推  

标签：pChEMBL ≥ 6；未测 ≠ inactive；主评测必须含 A-only/B-only。

## 4. 执行顺序

1. EGFR/HER2 小面板：TAK-285 QC → whole-mol vs moiety 诊断表  
2. 扩到三对公开靶 + 泄漏控制划分  
3. 开放基准与脚本  
4. NLRP3/JNK1 锚点案例（可选附录）  
5. 有打假表前不写全文  

## 5. 相关文档

| 文档 | 角色 |
|------|------|
| [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md) | 投稿主张与 WP |
| [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md) | 诊断操作细则 |
| [`NMI_REFERENCE_PAPER_PLAYBOOK.md`](NMI_REFERENCE_PAPER_PLAYBOOK.md) | 高分文流程对标 |
| [`DUAL_TARGET_SCORING_IMPLEMENTATION.md`](DUAL_TARGET_SCORING_IMPLEMENTATION.md) | 校准/短板等组件笔记 |
| [`REFERENCES_AND_MOLECULES.md`](REFERENCES_AND_MOLECULES.md) | 文献与分子 |
| [`DUAL_MULTI_TARGET_DOCKING_SURVEY.md`](DUAL_MULTI_TARGET_DOCKING_SURVEY.md) | 方法综述 |

---

*旧融合总方案、skills 审计、路线重审等历史文档已删除；以本文 + moiety 规划为准。*
