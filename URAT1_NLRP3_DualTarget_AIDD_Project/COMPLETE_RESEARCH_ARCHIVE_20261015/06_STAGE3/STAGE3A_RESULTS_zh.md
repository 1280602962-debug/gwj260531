# 第三阶段 A

脚本 `scripts/nlrp3_v2_stage3a_audit.py` 没有运行 GNINA。它只划定谁可以进入下一次对接。

v2 池并且通过化学软过滤的有 157 个。其中 118 个已经在历史 A2 对接集里，39 个没有。118 个里 103 个有羧酸、15 个没有。39 个里 32 个有羧酸、7 个没有。

对接决定只看羧酸是否存在。32 个有羧酸，记为 dock。7 个没有，记为 exclude，原因是酰磺胺或四唑，羧酸–Arg477 规则不适用。这 7 个是 TASISULAM SODIUM、TASISULAM、REPINOTAN、QUINOTOLAST、GALOSEMIDE、SULOFENUR、DU 125530。β-内酰胺、PAINS 和 Brenk 在这一步只记录，不从对接名单里删掉。

157 个里 PAINS 5 个、Brenk 38 个、β-内酰胺 11 个、PAINS 与 Brenk 都干净的 118 个。这些标记沿用 RDKit FilterCatalog，没有在第三阶段重调。

这一步没有修改冻结的 12 和 21。
