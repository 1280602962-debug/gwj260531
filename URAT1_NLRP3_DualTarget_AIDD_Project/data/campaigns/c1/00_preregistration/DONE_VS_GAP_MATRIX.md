# C1 vs 用户审阅意见：完成度对照（更新）

| # | 审阅要求 | 状态 | 证据 |
|---|----------|------|------|
| 1 | Rank 关闭，只走 Acid | **已完成** | `campaign_c1.yaml` |
| 2 | 旧 7 不预承诺 / 不选 MD | **已完成** | `do_not_precommit`；未开 discovery MD |
| 3 | A1 Arg ≤7.7027 Å | **已完成** | Amendment A1；脚本 assert |
| 4–5 | L2 lesinurad 败 / NP3 过 | **已完成** | `02_selfdock/` |
| 6 | 禁 L3 Rank | **已完成** | yaml + 无 GPU |
| 7 | Acid 子集 303→156 | **已完成** | ligand prep |
| 8 | Acid 双靶重算 | **进行中** | 烟雾 4 酸→watcher 启全量 156×2 |
| 9 | 冻结 C1 shortlist | **未做** | 等批对接 |
| 10 | shortlist 后才 MD | **延后** | Results §3.6 |
| 11 | Methods/Results 双阶段 | **正文主文已同步** | §2.1/2.9/2.14；Results §3.0–3.7 |
| 12 | Discussion/Conclusions | **顶栏+开头已改**；正文中段仍旧 |
| 13 | 口袋占位硬门 | **已补进脚本** | COM≤6 Å vs 晶体配体 |
| 14 | 不为留候选改门 | **遵守** | |

## 烟雾即时核对（自由对接，seed42；sanitize=False 重算后）

| ID | 分子 | Arg (Å) | A1 | COM (Å) | keep_URAT1 | keep_NLRP3 | dual |
|----|------|--------:|:--:|--------:|:----------:|:----------:|:----:|
| REP_00207 | lesinurad | 14.11 | fail | 1.63 | no | yes | no |
| REP_05846 | verinurad | 12.30 | fail | 0.75 | no | yes | no |
| REP_07907 | GSK-3008348 | ~2.98 | **pass** | (recalc) | **yes** | yes | **yes** |
| REP_07015 | puliginurad | 14.55 | fail | 1.62 | no | yes | no |

- GSK 曾因 RDKit 默认 sanitize 丢掉全部 9 姿 → 假 `no_poses`；已改 `load_poses(..., sanitize=False)` 并重启全量批。
- lesinurad/verinurad/puliginurad 自由对接口袋可占但 Arg 远：与 Rank 关闭一致；Acid 阳性对照仍是晶体 `local_only`。
- **不**因烟雾里 GSK 通过就把它预锁定为 MD；须等全量 shortlist + 审计。

## 自动化

- 全量批：`data/campaigns/c1/07_clinical_dock/acid_dual/` + `acid_dual_batch.log`（可 SKIP 恢复）
- 预计 CPU：156×2×~2–4 min ≈ 10–20 h
