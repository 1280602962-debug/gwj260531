# Agent command — Stage M（M1–M5 测量审计）

> **这是当前真正的第一步。**  
> 权威：[`PLAN_V2_REDTEAM_AND_REDESIGN.md`](PLAN_V2_REDTEAM_AND_REDESIGN.md) · [`EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md`](EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md)  
> 已有 v0 产物：`data/plan_v2_redteam_v0/`（M1/M3 初版 + M2 边界普查）  
> **M1–M3 + M5 = 零对接（云端可跑）；M4 = 本地对接（可分开发）。**

把下面整段复制给 agent。若只跑零对接部分，在开头加一句：`本轮只做 PART A–C 与 PART E；跳过 PART D（M4）。`

---

```text
【角色】
你是 DualFourClass 的测量审计 agent。
仓库含 Dual_Target_Docking/。
权威红队文：Dual_Target_Docking/docs/PLAN_V2_REDTEAM_AND_REDESIGN.md。
规划 v2：Dual_Target_Docking/docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md（Track B 被 Stage M 门控）。

【约束 — 强制】
1. 无湿实验。
2. PART A/B/C/E：**禁止新对接**；只读已冻结分数表与面板 CSV，可写新分析脚本与表。
3. PART D（M4）：仅当用户明确批准（或本命令未加「跳过 M4」）才做；对接协议必须与 panel120 as-run 一致（3POZ/3RCD, E=8, seed=20260727, n_modes=9），**统一 RDKit ETKDG+meeko**（或全程 LigPrep，二选一写死，不得混用）。
4. 禁止把 rtm_min_z / rank_consensus 写成已验证主臂；禁止调 shortfall/clash/flags 进主分数。
5. 禁止重启同协议 EGFR「扩样求显著」；禁止未经 Stage M 门控开 Track B。
6. 主读数必须是**方向分解**：AUROC(dual vs A_only) 与 AUROC(dual vs B_only) 分报；池化 Dual vs A∪B 只作附录。
7. 平凡基线（至少 heavy_atoms, MW, cLogP；建议加 TPSA）为必报表。
8. 所有数字必须可复现：脚本入库、随机种子写死、输出 CSV + STAGE_M_VERDICT.md。
9. 做完 git add/commit/push；更新 PR（若环境要求）。
10. 结束用中文 ≤15 行：M1–M5 各 Go/No-Go/Skip、关键数字、是否允许开 Track B。

【任务名称】
STAGE_M_MEASUREMENT_AUDIT = M1 + M2 + M3 + (可选 M4) + M5 → STAGE_M_VERDICT

【产出根目录】
Dual_Target_Docking/data/stage_m_v0/
  README.md
  analysis/STAGE_M_VERDICT.md          # 总门控（必须）
  analysis/M1_DIRECTIONAL.md
  analysis/M2_LABELS.md
  analysis/M3_BASELINES.md
  analysis/M4_UNIFIED_PREP.md          # 若跳过 M4：写 SKIPPED + 原因
  analysis/M5_ARM_LIST.md
  tables/                              # 所有数值表
  scripts/                             # 可复现脚本
  protocol/                            # M4 时写 prep 冻结说明

【只读输入】
- data/egfr_her2_panel120_v0/tables/panel_v0_120.csv
- data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv
- data/egfr_her2_panel120_v0/analysis/STAGE1_*.md
- data/pik3ca_mtor_panel48_v0/tables/panel_v0_48.csv
- data/pik3ca_mtor_panel48_v0/tables/ablation_ligand_scores.csv
- data/plan_v2_redteam_v0/             # 可复用/扩展，勿静默改上游分数
- data/protocols/CANDIDATE_ARMS_V0.yaml
- docs/PLAN_V2_REDTEAM_AND_REDESIGN.md

【已完成可核对、勿重复发明结论】
plan_v2_redteam_v0 已显示（须在 M1/M3 中复核写入 stage_m_v0）：
- EGFR/HER2 vina_mean：D/A≈0.689，D/B≈0.311，池化≈0.516
- 重原子数池化≈0.549 > vina_mean > rtm_min_z
- PIK3CA/mTOR 无反转；对接臂 0.65–0.69 > 体积基线 ~0.46
- 硬负贴边：EGFR A_only 60%、B_only 78% 在 pChEMBL 6±0.5

════════════════════════════════════
PART A — M1 方向分解（零对接）
════════════════════════════════════

【目标】把方向分解固化为正式主指标定义与主表。

【动作】
1. 扩展/迁移 plan_v2_redteam 脚本到 stage_m_v0/scripts/run_m1_directional.py
2. 对两对靶、各决策臂 + 平凡基线，输出：
   tables/m1_directional_auroc.csv
   列至少：pair, subset, arm, n_dual, n_A_only, n_B_only,
           auroc_D_vs_A, auroc_D_vs_B, auroc_pooled,
           top10_A_only, top10_B_only, top10_dual
3. EGFR 子集：all / old40 / new70（标明 prep 混杂，结论措辞谨慎）
4. 写 analysis/M1_DIRECTIONAL.md：定义主指标、展示表、明确废弃「仅报池化」

【M1 完成标准】
- [ ] 主表存在且与红队 v0 关键数字一致（允许 ±0.005 舍入）
- [ ] Top10 按 A_only / B_only 分列
- [ ] 文档写明：池化 AUROC 不作主结论

【M1 门控】完成即 Go（定义性工作，无统计 Go/No-Go）

════════════════════════════════════
PART B — M2 标签效度（零对接）
════════════════════════════════════

【目标】量化标签噪声与阈值敏感性；决定四类任务是否可支撑 Track B。

【B1 — margin 面板】
规则（写死进脚本与文档）：
- dual_strict：两端 pChEMBL ≥ 6.5
- A_only_strict：A≥6.5 且 B≤5.5（已测）
- B_only_strict：B≥6.5 且 A≤5.5（已测）
- neither_strict：两端均 ≤5.5（已测）
- gray：其余两端都测过的分子 → 主分析排除，单独计数

输出：
- tables/m2_margin_panel_counts.csv（每对靶各类 n、gray n）
- tables/m2_directional_on_margin.csv（在 strict 子集上重算 M1 方向 AUROC；若某类 n<8 标记 underpowered）

【B2 — 连续量】
- 目标：min_pchembl（两端最小值；缺失规则写明）
- 指标：Spearman(score, min_pchembl) 全体 + 按类分层描述
- 输出：tables/m2_continuous_spearman.csv

【B3 — 阈值敏感性】
- cutoff ∈ {5.5, 6.0, 6.5} 重新二值化四类（与 panel 构建规则一致：未测≠阴）
- 每个 cutoff 重算方向 AUROC（至少 vina_mean, rtm_min_z, heavy_atoms）
- 输出：tables/m2_threshold_sensitivity.csv
- 判定：结论符号（D/B 是否反转；臂是否 > 体积基线）是否随 cutoff 翻转

【B4 — 噪声天花板】
- 对每端 pChEMBL 加 N(0, σ²)，σ∈{0.3, 0.5, 0.7}；种子=20260728；B=500 次
- 每次扰动后按 cutoff=6.0 重贴四类标签，用**真值分数**（vina_mean）算方向 AUROC，取分布均值/分位作为「标签噪声下可达性」参照
- 更硬的上界（可选）：用 min_pchembl 本身当 oracle 分数，在扰动标签上算 AUROC → 标签可辨性天花板
- 输出：tables/m2_noise_ceiling.csv
- 写清：若 oracle 在 σ=0.5 下 D/A 或 D/B 中位数 <0.65，则四类任务信息量不足，倾向 Track A

【文档】analysis/M2_LABELS.md（含 gray 比例、敏感性是否翻转、天花板数字）

【M2 门控】
- Go：margin 子集上方向结论与主面板同号；oracle/噪声天花板显示任务仍可辨（至少一端方向中位 AUROC≥0.65 @σ=0.5）
- Weak：天花板边缘或 margin n 过小 → Track B 高风险，须在 VERDICT 写明
- No-Go：阈值一变结论就翻，或天花板显示任务不可辨 → **禁止 Track B**，收 Track A

════════════════════════════════════
PART C — M3 平凡基线（零对接）
════════════════════════════════════

【目标】证明/证伪对接臂相对非对接基线的增量。

【动作】
1. 基线组（全部必报）：heavy_atoms, MW, cLogP, TPSA
2. 可选：Morgan fingerprint Tanimoto 到「已知 dual 原型」中位数相似度
   （原型 = 该 pair 面板内 role_note 含 gold/classic 的 dual，或全部 dual 的 leave-one-out 中位相似 — 写死一种，防泄漏说明）
3. 在方向分解指标下，比较每个对接臂 vs 最强体积基线（通常 heavy_atoms 或 MW）
4. 输出：tables/m3_baselines_vs_arms.csv
5. 规则写入文档：臂若不能在 **D/A 与 D/B 的汇总**（建议用 min(D/A,D/B) 或平均，写死）上超过体积基线，标记 fail_baseline

【文档】analysis/M3_BASELINES.md

【M3 门控】
- Go：≥1 对接臂在 ≥1 对靶上稳定超过体积基线（方向汇总），且另一对不严重更差
- No-Go：两对靶上所有对接臂均 ≤ 体积基线 → 与红队一致，方法轨无增量可卖

（注：EGFR 上预期 No-Go 倾向、PIK3CA/mTOR 上预期有信号 — 分别报告，勿平均掉。）

════════════════════════════════════
PART D — M4 统一 prep（本地对接；可 Skip）
════════════════════════════════════

【若本轮跳过】
在 analysis/M4_UNIFIED_PREP.md 写：
  STATUS: SKIPPED
  原因：用户未批准 / 本轮零对接
  影响：旧40 vs 新70 的 RTM 分裂不得写成方法结论；VERDICT 中 M4=Skip，Track B 最高只能 Weak 且须重跑 prep 后才能升 Go

【若执行 — 最小充分方案（推荐先做）】
方案 M4-min：仅重跑 **旧 40**，统一为与新 70 相同的 RDKit ETKDG+meeko
- 配体：panel120 中 from_panel40=yes 的 40 个 SMILES
- 受体/盒子：与 panel120 protocol 相同（3POZ/3RCD）
- E=8；seed=20260727；n_modes=9；随后 RTM best-of-9（与现脚本一致）
- 产出目录：data/egfr_her2_panel40_reprep_rdkit_v0/（或 stage_m_v0/m4_reprep/）
- 比较表：同一 40 分子 LigPrep 旧分 vs RDKit 新分
  tables/m4_old40_prep_delta.csv
  指标：方向 AUROC、per-ligand |Δscore|、姿态 mode 一致性（若可）

方案 M4-full（可选加强）：全部 110 统一 RDKit 重跑（2×110）
- 仅当 M4-min 显示 prep 改变主结论时再升级

【对接完成后分析】
- 若统一 prep 后「新70上 RTM 显著差于 vina」仍成立 → 可升级为方法学结论
- 若消失 → S1 分裂主因是 prep，写入 VERDICT，纠正 Track A 措辞

【M4 门控】
- Go：统一 prep 对照完成并写入结论（无论 prep 是否为主因）
- Skip：允许，但 Track B 不能 Full-Go
- Fail：对接未按协议 / 混用 prep → 重做

════════════════════════════════════
PART E — M5 清洗预注册臂 + 总门控
════════════════════════════════════

【动作】
1. 更新 data/protocols/CANDIDATE_ARMS_V0.yaml 或新建 CANDIDATE_ARMS_V1_STAGE_M.yaml：
   - 主竞赛臂 ≤4：例如
     (i) vina_mean
     (ii) heavy_atoms（对照，非「方法」）
     (iii) 至多一个 rescoring 臂（rtm_min 或 rtm_min_z；校准仅 LOTO 训练折）
     (iv) 至多一个机制臂占位：size_deconfounded_vina 或 endpoint_comparability（可先定义公式，分数可后算）
   - rank_consensus_*：保持 EXPLORATION_DERIVED，移出主竞赛
   - 删除/降级「看完 panel40 才推荐」的措辞
2. analysis/M5_ARM_LIST.md：列出最终主竞赛臂、排除理由、与 Stage M 数字的关系
3. **写总门控** analysis/STAGE_M_VERDICT.md，表格：

| 项 | 状态 | 关键数字 | 含义 |
|----|------|----------|------|
| M1 | Go/Fail | … | |
| M2 | Go/Weak/No-Go | … | |
| M3 | Go/No-Go（按 pair） | … | |
| M4 | Go/Skip/Fail | … | |
| M5 | Go/Fail | 臂数= | |

【总门控规则 — 写死】
- **Track B Full-Go**：M1 Go ∧ M2 Go ∧ M3 至少一对靶 Go ∧ M4 Go ∧ M5 Go
- **Track B Weak（仅文档规划，不启动大批对接）**：M4=Skip 但 M1–M3 非 No-Go；或 M2=Weak
- **Track B No-Go → 收 Track A**：M2 No-Go 或（两对靶 M3 均为 No-Go）或 M1 Fail
- Track A 写作头条必须采用红队 §5 的测量发现（方向抵消 + 基线 + 标签边界），不得以「同系列外推失败」为唯一头条（除非 M4 已排除 prep 混杂）

【README.md】
说明如何复现、哪部分跳过、与 plan_v2_redteam_v0 的关系。

【文档挂钩】
- 更新 docs/EXPERIMENTAL_PLAN_DUALFOURCLASS_V2.md §11 状态：Stage M 进行中/完成
- 更新 docs/PROJECT_MASTER_PLAN.md 现行命令指向本文
- 勿大改科学正文，只更新状态与指向

【退出检查】
- [ ] STAGE_M_VERDICT.md 存在且给出 Track B Full-Go / Weak / No-Go
- [ ] 所有 tables/ 有对应 scripts/
- [ ] 未宣称 C3/C4 成立
- [ ] 未在 M4=Skip 时把旧40/新70 RTM 分裂写成已确认方法结论
- [ ] commit/push 完成
- [ ] 中文总结含：是否允许开 Track B；若 No-Go，下一步是 Track A（STEP0）还是补 M4

【明确不在本命令】
- 第三/四对靶大批对接（Track B）
- STEP0 论文写作包（Stage M 门控后再跑 AGENT_COMMAND_STEP0）
- 调参 shortfall/clash、发明新融合公式刷分
- 湿实验
```

---

## 你怎么用（给用户）

| 场景 | 复制时加的一句话 |
|------|------------------|
| **现在就跑（推荐）** | `本轮只做 PART A–C 与 PART E；跳过 PART D（M4）。` |
| **本地有对接算力** | 不加跳过句，或写 `执行 M4-min（旧40统一 RDKit 重跑）。` |
| **只要门控结论** | 同上零对接版；看 `STAGE_M_VERDICT.md` 的 Track B Full-Go/Weak/No-Go |

**命令文件路径：** `Dual_Target_Docking/docs/AGENT_COMMAND_STAGE_M.md`  
**成功标志：** 打开 `data/stage_m_v0/analysis/STAGE_M_VERDICT.md` 能看到五门状态与是否允许 Track B。
