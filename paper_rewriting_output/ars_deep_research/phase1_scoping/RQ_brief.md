# Phase 1 — Research Question Brief（ARS deep-research）

**Mode:** full / lit-review → plan  
**Date:** 2026-07-17  
**Skill:** `academic-research-suite` (ARS-Codex v0.1.18) → `ars/deep-research`

## Primary RQ

在约 **52.8 万 ChEMBL 唯一胺** 的条件下，对 **JNK2 Cys116 / PDB 8ELC 丙烯酰胺共价筛选**，对接前 triage（胺预筛→装弹头→物化 QC→锚点相似→双轨分仓→缩库）应如何设计，才能在证据上区别于「任意工程过滤」，并与已确认贡献（可复现 Framework + 可选新化学）对齐？

## Sub-questions

1. 文献中对接前 ligand-based / 反应枚举 / 双轨（类似物 vs 新骨架）哪些是 **标准模块**，哪些有 **位点特异性校准** 先例？  
2. ECFP Tanimoto + Murcko + hinge SMARTS 的主要失败模式是什么？可用哪些可实现算法补救？  
3. 不确定性和“有效性”应如何 **可证伪地量化**（在湿实验之前）？  
4. 最小可行初筛计划（MVS）应包含哪些硬步骤与硬门控？

## FINER

| Criterion | Score | Note |
|-----------|-------|------|
| Feasible | 高 | WSL+RDKit；对接/AF3 后置 |
| Interesting | 高 | 直接决定 527k→对接名单质量 |
| Novel | 中 | 模块不新；**位点校准+选择性/反应性** 有空间 |
| Ethical | 高 | 计算 triage，无人体数据 |
| Relevant | 高 | 绑定 confirmed Framework contribution |

## Scope

**In:** 对接前 triage；丙烯酰胺主 warhead；四锚点化学；与 Glide/AF3 分层衔接；可量化校准。  
**Out:** 疾病模型主 claim；声称新相似度算法发明；无校准的“活性预测 IC50”。
