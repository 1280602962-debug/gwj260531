# 中英文润色修改摘要（2026-09-12，审稿风险消除轮）

本轮按投稿前审稿风险清单改中英文，不改冻结 AUROC，不加新对接/MD/ML，不改官方论文 PR #32。

## 已改

1. **摘要。** ECFP4 改为 receptor-free class information；结论句改为 experimental measurements at both targets。
2. **引言首段末句。** 改为 docking-score interpretation when ligands have experimental measurements at both targets。三问未改。
3. **2.3。** 删除 hard negatives / 硬负样本；改为 A-only / B-only 或 strict selective class。
4. **2.3。** 删除 new docking / 新对接；G5 改为 common rigid-receptor Vina protocol，并指向 Table S14。不加 prespecified。
5. **Table S14。** 新增最后 17 对 + G3 两对 + EGFR/HER2 的对级纳入/排除审计。删除 “Other pairs that passed earlier gates were not docked.”
6. **EGFR/HER2。** 改为 supply-limited case；删 thick-supply / 厚供给。4.1 补 JAK1/TYK2 同量级现象。
7. **Table S15。** 未新增。现有冻结表不能按 8 对 × 类别给出 IC50/Ki/Kd/EC50/Potency 与 source-document 完整交叉，避免虚构。
8. **3.6。** 压成四句；删三次重复的 “not a claim that no external data exist”。
9. **F2/F10。** 从 4.5 移到 4.1，作为 pair-dependent 方向反转。
10. **4.3。** GNINA 结论收窄到 EGFR/HER2 与 JAK1/TYK2。
11. **4.5。** 写明 θ=6.0 对 6.5/5.5 与骨架上限差异；绝对 AUROC 不用于跨对排序。
12. **分层 bootstrap。** Table S10 未重算 class-stratified `summary_min`，2.5.3 如实写未重算，不虚构。
13. **4.2。** 与摘要对齐；保留 “does not establish that ECFP4 is a formally better predictor than docking.”
14. **结论首段。** experimental-state definition；measurements at both targets。
15. **结论次段。** matched-pocket advantages were not consistently reproduced；高 AUROC 不能证明 pocket-specific 3D information。
16. **Data Availability。** GitHub + 可由分数/元数据表复现。无 Zenodo DOI，不写 pending。SI 改为 S1–S14。
17. **Figure 1 图注。** 末句改为 shown separately from the supply census。
18. **3.5。** 未再展开。
19. **4.5 样本量模拟。** 未把 0.109/0.119/0.128 放回主文。
20. **术语。** 终稿去掉 hard negative、new docking、thick-supply、both-end experimental labels、original/added pair。

## 数字未改

Table 2 / Table 3 冻结 AUROC 与区间未重算。
