# Hard-negative case notes — EH40_23

- **panel_id**: EH40_23
- **chembl**: CHEMBL3633938
- **class**: A_only (EGFR active / HER2 measured weak) — hard negative for dual ranking
- **pref_name**: (none)

## Rank before/after gating

| arm | rank |
|-----|------|
| vina_mean | 4 |
| vina_min | 4 |
| rtm_mean | 2 |
| rtm_min | 2 |
| rtm_min_z | 2 |
| gated_rtm_min | 2 |

## Clash gate (conservative)

- clash_distance_A = 2.2; fail if clash_atom_pairs >= 3
- 3POZ RTM-best mode 6: clash_count = 0; fail=0
- 3RCD RTM-best mode 2: clash_count = 0; fail=0
- Vina mode1 clashes: 3POZ=0, 3RCD=0

## Conclusion

**门控未能打掉**：EH40_23 在 RTM 最优 pose 上未触发 clash 门控（gated_rtm_min 名次仍为 2，与 rtm_min=2 相同）。
观察：该分子为 anilinoquinazoline 类骨架，与共晶配体/Lapatinib 化学空间接近，Vina+RTM 易给出高分伪双重命中；简单几何 clash 不足以识别 HER2 端假阳性。未改动标签或分数刷榜。下一步建议：PoseBusters / PLIP 关键氢键与 hinge 几何、或 HER2 特异性相互作用指纹 / MM-GBSA 再打分。
