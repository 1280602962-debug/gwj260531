# Agent command — 下一步立即可跑（J0+J1 + Track A 素材，零对接）

> **现在就可以给 agent 的命令。** Stage M 已完成（Track B=Weak）；本命令只做**零对接**工作。  
> 权威：[`JCIM_ROUTE_ASSESSMENT_V1.md`](JCIM_ROUTE_ASSESSMENT_V1.md) · Stage M 结论：`data/stage_m_v0/analysis/STAGE_M_VERDICT.md`  
> 已有可行性包：`data/jcim_feasibility_v0/`（12 对供给审计 + 统一 prep EH110）

把下面整段复制给 agent。

---

```text
【角色】
你是 DualFourClass 的文档/数据分析 agent（云端即可）。
仓库含 Dual_Target_Docking/。
权威路线：docs/JCIM_ROUTE_ASSESSMENT_V1.md
Stage M 总门控：data/stage_m_v0/analysis/STAGE_M_VERDICT.md → Track B = Weak（禁止开大批对接）
可行性包：data/jcim_feasibility_v0/

【约束 — 强制】
1. 无湿实验。
2. **禁止任何新对接 / 重对接 / 升 E / 换引擎跑分**（本命令零对接）。
3. 禁止开 Track B；禁止同协议 EGFR 扩样赌显著；禁止把 rtm_min_z 写成已验证主臂。
4. 主读数=方向分解 AUROC；平凡基线必报；池化 Dual vs A∪B 只作附录。
5. 引用统一 prep EH110 时用 jcim_feasibility_v0 的数字（或本命令重算），不得把混 prep panel120 的 RTM 分裂写成方法结论。
6. 产出可复现：脚本 + CSV + md；git add/commit/push；更新 PR（若环境要求）。
7. 结束中文 ≤12 行：完成了什么、文件路径、下一步是否需要用户批准本地对接（J2+）。

【任务名称】
NEXT_IMMEDIATE = J0（扩供给审计）+ J1（选定对接对靶草案）+ TA（Track A 素材包升级）

【产出根目录】
Dual_Target_Docking/data/jcim_j0j1_v0/
  README.md
  analysis/J0_SUPPLY_AUDIT.md
  analysis/J1_PAIR_SELECTION_DRAFT.md
  analysis/TRACK_A_FIGURE_CLAIM_PACK.md
  tables/
  scripts/
另：升级或新建
  Dual_Target_Docking/data/track_a_starter_v0/   （若已有则更新，勿与旧 STEP0 冲突）
  Dual_Target_Docking/data/protocols/  （若需补 PAIR_ROLES_DRAFT.yaml）

════════════════════════════════════
PART A — J0：扩大 ChEMBL 供给审计（零对接，核心）
════════════════════════════════════

【目标】把「strict 硬负供给上限」从 12 对扩到 **≥40–50 候选对**，形成论文 C2 主干。

【输入】
- data/public_pair_selection/（现有 mols_*.json、chembl_pair_fourclass.csv、FROZEN_PUBLIC_PAIRS.yaml）
- data/jcim_feasibility_v0/scripts/audit_strict_label_supply.py（可扩展）
- docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md、DUAL_TARGET 文献目录（挖候选对）

【动作】
1. 拟定 ≥40–50 个候选靶对名单（覆盖：激酶同源、异质酶、PPI、转运体/GPCR 若数据够、文献高频双靶）。
   - 必须包含已审计的通过/失败对作对照。
   - 明确排除：NLRP3/JNK1（私有 holdout）；金属酶-only 作唯一主对需标注风险。
2. 若本地已有 pChEMBL 字典则复用；若缺靶：
   - **优先**用 ChEMBL API/已有脚本抓取（若环境禁网或失败：把「待抓取靶列表」写入 tables/j0_fetch_queue.csv，并对已有靶先出完整表，勿假装完成 50 对）。
3. 对每对计算（与 feasibility 脚本一致）：
   - n_both_measured
   - θ=6.0：dual / A_only / B_only
   - strict 6.5/5.5：dual / A_only / B_only / neither / gray / gray_frac
   - min_strict_hardneg；supports_strict_panel（建议阈值：两侧硬负各 ≥50；另报 ≥20 的「薄面板」档）
4. 输出：
   - tables/j0_strict_label_supply.csv
   - tables/j0_candidate_pairs.csv（动机、文献锚点、是否金属酶/PPI）
   - analysis/J0_SUPPLY_AUDIT.md（森林/排序表 + 一句话：有多少对可支撑严格四类面板）

【完成标准 A】
- [ ] 至少对「字典已齐全」的全部对给出 strict 表；目标总数 ≥40（不足则写明缺字典队列）
- [ ] 明确列出 supports_strict_panel=Y 的对与「薄面板」对
- [ ] 复现脚本入库

════════════════════════════════════
PART B — J1：对接对靶草案（零对接，选不采）
════════════════════════════════════

【目标】在不对接的前提下，冻结「若用户批准算力，下一轮对接哪 K=4」。

【动作】
1. 结合 J0 + holo PDB 可得性（可用已有 pdb_holo_counts.csv / RCSB 粗计数；无网则用已有表）选出：
   - **Tier S（建议对接）**：supports_strict_panel=Y，且常规可对接口袋
   - **Tier T（薄/案例）**：硬负不足但有共晶/异质折叠价值（如 Mcl-1/Bcl-xL、EGFR/HER2）
2. 默认草案（可被 J0 结果修改，改了要写理由）：
   - 必留：PIK3CA/mTOR（已有；需统一 prep）
   - 必留：EGFR/HER2（已有；统一 prep 已齐；作「供给受限 + 方向反转」案例，不要求 strict 厚硬负）
   - 优先新对：AChE/BChE（若 J0 仍通过）
   - 第 4 席：PIK3CA/PIK3CB（同工酶对照）或审计中新出现的合格对；Mcl-1/Bcl-xL 作 Tier T
3. 写 analysis/J1_PAIR_SELECTION_DRAFT.md：
   - K=4 名单与角色（development / case / holdout 预留）
   - 每对：建议 N、strict 定额、结构候选 PDB、风险
   - **对接预算表**（仅估算，不执行）
4. 写 data/protocols/PAIR_ROLES_DRAFT_J1.yaml（机器可读）

【完成标准 B】
- [ ] 有明确 K=4 草案 + 预算，且标明「待用户批准后才对接」
- [ ] EGFR 标为 case/supply-limited，不假装能做厚 strict 面板

════════════════════════════════════
PART C — Track A 素材包（马上可成文的部分）
════════════════════════════════════

【目标】把现有 Stage M + 统一 prep EH110 + 供给上限 收成投稿素材（Mol. Inf. / JCAMD 可直接用；JCIM 作 C1/C2/C4）。

【动作】
1. 固化统一 prep 表到本包或 track_a_starter（可复制 jcim_feasibility_v0 结果，勿改上游原始分数含义）：
   - tables/eh110_unified_prep_directional.csv（引用或复制）
2. 写 analysis/TRACK_A_FIGURE_CLAIM_PACK.md，固定：
   **Claim 天花板**
   - 允许：方向抵消；EGFR 对接≤平凡基线（prep-clean）；PM 有对接信号；prep 敏感性；公开数据 strict 硬负供给受限
   - 禁止：rtm 通吃；混 prep 方法结论；通用决策尺子已验证；湿实验/乘客
   **主图计划（每张：数据文件、轴、一句 caption）**
   - Fig1 任务定义（四类）
   - Fig2 方向分解（EH vs PM；含统一 prep EH110）
   - Fig3 平凡基线 vs 对接臂
   - Fig4 标签灰区 / strict 供给（含 J0 若已有；至少含 12 对 feasibility）
   - Fig5（可选）M4 prep 敏感性
   - Fig6（可选）失败分型 T1/T2/T5 指针
3. 写/更新 track_a_starter_v0/PAPER_OUTLINE_TRACK_A.md（章节 + 文件指针 + TODO）
4. 写 track_a_starter_v0/CLAIM_CEILING.md 与 FIGURE_PLAN.md（若与 STEP0 重复则合并更新，以本命令头条为准）

【完成标准 C】
- [ ] 一张图计划能支撑 Mol. Inf./JCAMD 诊断文
- [ ] 头条不再是「同系列外推失败」，而是测量 + 供给上限

════════════════════════════════════
PART D — 文档挂钩与退出
════════════════════════════════════

1. 更新 docs/PROJECT_MASTER_PLAN.md：现行命令指向本文；Stage M 标完成；下一步 J2+ 需批准。
2. 在 docs/JCIM_ROUTE_ASSESSMENT_V1.md §4.4 将 J0/J1 标为进行中或完成。
3. README：如何复现、哪些未做（对接）。

【退出检查】
- [ ] 零新对接
- [ ] Track B 未启动
- [ ] J0 表 + J1 草案 + Track A 素材三件套存在
- [ ] commit/push
- [ ] 中文总结给出：用户若批准本地对接，下一命令应是「J2 PM48 转 RDKit + J3 新受体冻结」（勿自行执行）

【明确不在本命令】
- J2/J3/J4/J5 任何对接或 GNINA
- 第三对靶结构准备以外的湿实验/合成
- 调参 shortfall/clash
- 直接写投稿 PDF/LaTeX（提纲即可）
```

---

## 你怎么用

| 场景 | 做什么 |
|------|--------|
| **现在立刻** | 复制上面整段给云端 agent（本命令） |
| **本命令跑完后** | 你手里有：更大供给审计表、K=4 对接草案、Track A 图/claim/提纲 |
| **需要你再批准才开** | 本地对接命令（另写）：J2 PM48→RDKit、J3 新受体冻结、J4 新面板对接、J5 GNINA rescore |

## 马上能完成 vs 不能

| 马上能完成（零对接） | 不能假装马上完成 |
|----------------------|-------------------|
| J0 供给审计扩展 | 新靶对接与 RTM |
| J1 选对草案 + 预算 | GNINA（需姿态/本地） |
| 统一 prep EH110 固化进素材包 | Track B Full-Go |
| Track A 图计划 / claim / 论文提纲 | JCIM 完整 Article（还差 K=4 对接） |
| Mol. Inf. / JCAMD 级诊断文骨架 | 湿实验 |

## 成功标志

打开 `data/jcim_j0j1_v0/analysis/J0_SUPPLY_AUDIT.md` 能看到「多少对可支撑严格面板」；  
打开 `TRACK_A_FIGURE_CLAIM_PACK.md` 能按图写诊断文；  
且**没有任何新对接产物**。
