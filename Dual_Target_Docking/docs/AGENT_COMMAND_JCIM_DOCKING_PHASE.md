# Agent command — 冲 JCIM 对接阶段（已批准路线）

> 白话路线：[`JCIM_NEXT_ROUTE_PLAIN.md`](JCIM_NEXT_ROUTE_PLAIN.md)  
> 依据：`data/jcim_j0j1_v0/`（供给审计 + K=4 草案）+ Stage M（Track B=Weak，本文不做方法竞赛）  
> **用户意图：冲 JCIM。** 对接仅服务评测/基准，不服务「调出通吃决策臂」。

复制下面整段给 agent。若本轮只做文档授权、不做本地对接，开头加：  
`本轮只做第一步；第二步及以后跳过。`  
若本地有对接算力，用：  
`第一步到第三步全做；第四步 GNINA 若环境无则 Skip 并写明。`

---

```text
【角色】
你是 DualFourClass 冲 JCIM 的执行 agent。
权威白话路线：Dual_Target_Docking/docs/JCIM_NEXT_ROUTE_PLAIN.md
K=4 草案：data/jcim_j0j1_v0/analysis/J1_PAIR_SELECTION_DRAFT.md
协议角色：data/protocols/PAIR_ROLES_DRAFT_J1.yaml
Stage M：Track B=Weak — 禁止把本阶段做成「选通用主臂」竞赛。

【约束 — 强制】
1. 无湿实验。
2. 论文形态=评测/基准；不宣称通用决策尺子已验证。
3. 禁止再扩 EGFR 对接；EGFR 只用现有统一 prep EH110 作供给受限案例。
4. 全库配体准备写死：RDKit ETKDG + meeko（不得与 LigPrep 混用进主表）。
5. 主指标=方向分解 AUROC(dual vs A_only) 与 (dual vs B_only)；必报平凡基线；池化只作附录。
6. 新面板按 strict 标签（dual 两端≥6.5；A_only：A≥6.5且B≤5.5；B_only 对称）定额抽样；灰区分子不进主面板或单独表。
7. 新受体必须先结构冻结 + cognate QC，通过后才批量对接。
8. 不做 shortfall/clash 调参进主分数；不把 flags/架构灌进分数。
9. 可复现：每步 protocol.yaml + 分数表 + README；git commit/push；更新 PR。
10. 结束中文总结：完成到第几步、路径、是否需用户补 GNINA/扩 PM。

【K=4 冻结名单】
1. PIK3CA/mTOR — 开发对；需 RDKit 重跑（现有 panel48=LigPrep）
2. EGFR/HER2 — 案例；零新对接
3. AChE/BChE — 新开发对
4. PIK3CA/PIK3CB — 同工酶对照（叙事写过近）

════════════════════════════════════
第一步 — 文档授权冻结（云端可做）
════════════════════════════════════

1. 新建 data/protocols/PAIR_ROLES_APPROVED_JCIM.yaml
   - 从 PAIR_ROLES_DRAFT_J1.yaml 复制并设 docking_authorized: true
   - 写明批准日期、白话路线文档路径、禁止项（不扩 EGFR 等）
2. 更新 docs/JCIM_NEXT_ROUTE_PLAIN.md 或 PROJECT_MASTER_PLAN：状态=对接阶段已批准
3. 产出目录约定：
   - data/pik3ca_mtor_panel48_rdkit_v0/   （第二步）
   - data/ache_bche_panel_v0/             （第二/三步）
   - data/pik3ca_pik3cb_panel_v0/         （第二/三步）
   - data/jcim_bench_v0/                  （第五步汇总，可先建 README 占位）

【第一步完成标准】授权 YAML 存在且 authorized=true；未开始对接也算本步完成。

════════════════════════════════════
第二步 — PM 统一 prep + 新受体冻结（本地）
════════════════════════════════════

### 2A. PIK3CA/mTOR panel48 → RDKit
- 输入：现有 panel48 SMILES/名单；受体冻结 4L23 / 4JT6
- 动作：RDKit ETKDG+meeko → Vina（与冻结 E/seed/n_modes 对齐，写进 protocol）→ RTM best-of-K
- 输出：data/pik3ca_mtor_panel48_rdkit_v0/tables/… + MANIFEST + 与 LigPrep 旧分对照表（prep 敏感性）
- 可选：按 strict 定额扩到 N≈110（若做，名单单独 CSV，勿 silently 混进 48）

### 2B. 新受体冻结（AChE、BChE、PIK3CB）
- AChE 候选：4EY4 / 6O5V 等；BChE：6ZWI / 1P0I 等；PIK3CB：选 holo，PIK3CA 复用 4L23
- 每端：准备蛋白、盒子、cognate 配体重对接 QC（RMSD/关键相互作用写一页）
- QC 不过禁止进入第三步该端

【第二步完成标准】
- [ ] PM48 RDKit 全成功分数表
- [ ] 三个新蛋白（或明确复用）冻结文件 + cognate QC md
- [ ] EGFR 无新作业

════════════════════════════════════
第三步 — 新两对面板对接（本地）
════════════════════════════════════

对 AChE/BChE 与 PIK3CA/PIK3CB：

1. 从 ChEMBL（已有字典优先）按 **strict 定额** 抽 ~100–120：
   建议下限：dual≥25，A_only≥20，B_only≥20，neither≥10（不足则如实降 N 并在 MANIFEST 说明供给上限）
2. 骨架配额、排除共价/PROTAC 特例（规则写入 MANIFEST）
3. RDKit prep → 双端 Vina → RTM
4. 产出每对：panel CSV、scores、job_status、protocol.yaml、四类计数

【第三步完成标准】两对均有完整分数表；失败配体有清单；未调阈值刷分。

════════════════════════════════════
第四步 — GNINA 重打分（本地；可 Skip）
════════════════════════════════════

- 对已有 Vina 姿态做 CNN rescore（不强制重新全局对接）
- 覆盖：EH110（若姿态在）、PM RDKit、AChE/BChE、PIK3CA/PIK3CB
- 若无 GNINA：analysis 写 STATUS: SKIPPED + 原因；不阻断第五步，但 JCIM 文需在 Limitations 写「单引擎采样」

════════════════════════════════════
第五步 — 汇总分析与基准打包（云端可做）
════════════════════════════════════

1. 四对靶统一脚本：方向分解 + 基线 + bootstrap CI + 森林图数据表
2. 写入 data/jcim_bench_v0/：README、CLAIM_CEILING、主表、复现命令
3. Track A 素材升级为 JCIM 提纲（Evaluation Article）：C1–C5 对号入座
4. 更新 docs/JCIM_NEXT_ROUTE_PLAIN.md 状态为「对接完成 / 待写作」

【第五步完成标准】别人能按 README 复现主表；claim 不超过评测/基准。

════════════════════════════════════
退出
════════════════════════════════════

- 未批准的靶对零对接
- 未把 rtm_min_z 封为通用主臂
- commit message 示例：
  "Start JCIM docking phase: authorize K=4 and [PM RDKit / receptor freeze / panels]."
- 中文总结告诉用户：下一步是写作投稿，还是补 GNINA/扩 PM
```

---

## 你怎么用

| 你的环境 | 复制时加一句 |
|----------|----------------|
| 只有云端 | `本轮只做第一步；第二步及以后跳过。` |
| 本地能对接 | `第一步到第三步全做；第四步有 GNINA 就做否则 Skip。` |
| 对接已齐只分析 | `只做第五步；前步已完成则核对路径。` |

## 成功后你拥有什么

K=4 里除 EGFR 案例外的对接分数 +（可选）GNINA + 可投稿 JCIM 的基准包骨架。
