# Agent command — Stage 1 EGFR/HER2 expand (S1 gate)

```text
【约束】可大批对接；无湿实验。权威规划见 docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md。
本阶段只做阶段1（S1），不要开第三对靶，不要上复杂ML，不要调clash/shortfall刷分。

【任务】EGFR/HER2 panel40 → N≈100–120 扩展并对接重评。

【难点条款（强制）】
- 主标签=四类活性；架构不作为选取硬过滤（可后标 compact_ATP/clear_linker/unknown）
- 同Murcko骨架配额上限≤5
- 重点补充 A_only/B_only 硬负；两端皆有人类靶定量活性；规则与panel40一致（θ等写进manifest）
- 产出端水平分数表；抽样更新T2/T5个案（L2）；禁止flags进gated score

【对接协议】
- 结构：3POZ / 3RCD；E=8；seed=20260727；n_modes=9；RTM best-of-9
- 与panel40 as-run一致，不升E化妆

【分析】
- 臂：vina_mean、rtm_min_z（及已有可比臂）
- 配对ligand-bootstrap 95% CI：ΔAUROC(rtm_min_z − vina_mean)
- 预冻结判定：CI排除0 → Go(S1)；否则 No-Go（写明降级建议）

【交付】
- 扩展名单CSV、四类计数、对接/RTM表、bootstrap CI、STAGE1_VERDICT.md
- commit/push；中文5行总结Go/No-Go与下一步
```
