# 中英文润色修改摘要（2026-09-12，nature-polishing 文字轮）

按仓库先前使用的 `Yuan1z0825/nature-skills`（`nature-polishing` + `nature-shared`）审英文稿。轴：`paper_type=research`，`section=whole manuscript`，`language=en`，`journal=generic`（JCIM）。未改冻结 AUROC，未改文章框架。

## 本轮文字修改

1. **摘要破折号。** `pairs—dual...neither—` 改为括号，去掉 em dash。
2. **术语账本。** 活性数据统一为 experimental measurements at both targets；both-end 仅保留给对接分数和 neither 类。
3. **3.1 标题。** Both-end experimental supply → Paired experimental supply。
4. **3.3。** 删除 matched or exceeded docking AUROC，与摘要/4.2 的 competing-explanation 对齐。
5. **2.3。** Table S14 由三处压成一处；EGFR/HER2 长句拆开。
6. **3.6。** 两句 “no pair” 合并为一句准入标准 + 一句未形成外部集。
7. **4.1。** 一段拆成三段（机制、异质性、与 Zhou 的定位）；DUD-E- 连字符去掉。
8. **4.3。** 删除 “particularly for homologous targets” 同义重复。
9. **4.4/4.1。** the present study/analysis 统一为 this study。
10. **4.5。** 去掉与 3.6 重复的 BindingDB/2018 复述；双 “main limitation” 开场合并。

中文稿做了对应修改。

## PR 37 language/data-alignment revision (2026-09-13)

1. Table 2 now reports both directional AUROCs with pointwise bootstrap intervals; the descriptive summary_min interval remains the locked pooled estimate.
2. Added the class-stratified summary_min sensitivity and its provenance table without replacing the primary bootstrap.
3. Replaced vague “drug-like filter” wording with the executable molecular-identity, molecular-weight, heavy-atom, and metal-element rules, and separated manual scope decisions from executable filters.
4. Defined cognate redocking top-1, top-3, and all-saved-pose RMSD; clarified that GNINA was a targeted sensitivity analysis.
5. Added GroupKFold and logistic-regression settings, conditional best-descriptor wording, and operational external-gate wording.
6. Reassembled both manuscripts and regenerated the checksum manifest and submission pack. Numeric audit: 114 PASS, 0 FAIL, 3 NOTE.
