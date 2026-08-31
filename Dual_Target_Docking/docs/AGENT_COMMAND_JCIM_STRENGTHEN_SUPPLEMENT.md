# Agent 命令 — 审稿审计后补实验（给本地 agent）

> 权威审计：[`../data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md`](../data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md)  
> 加厚规划：[`JCIM_STRENGTHENING_PLAN_V1.md`](JCIM_STRENGTHENING_PLAN_V1.md)  
> Claim：[`../data/jcim_bench_v0/CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)  
> 已有口袋匹配表：`data/jcim_bench_v0/tables/pocket_matched_*.csv`（T0.1/T0.2 已部分落地）

把下面 **整段** 复制给本地 agent。  
默认跑：**A 组（零对接）+ B 组最小对接核**。若只要分析：开头加 `本轮只做 A 组`。若算力够冲刺：加 `A+B+C 全做`。

---

## 需要补什么（一页汇总）

### 为什么补（三件事）
1. **主指标测错口袋**：dual vs A_only 的标签差只在 B 端，但旧主表用 `vina_mean` 池化分做两端对比 → 必须改成**口袋匹配**，并保留池化为对照。  
2. **信号大多是混淆**：错口袋对照远偏离 0.5；尺寸/效价/TPSA 捷径明显 → 必须做匹配子集 + LE 归一。  
3. **唯一成功对（PM）恰好 E=16，其余 E=8** → 必须做 exhaustiveness 对照 + 单靶 enrichment sanity，否则审稿人会说“不是对接赢了，是采样更足”。

### 任务清单

| 组 | ID | 内容 | 对接？ | 状态 |
|----|-----|------|--------|------|
| A | T0.1–0.2 | 口袋匹配主指标 + 错口袋对照 | 否 | **脚本已有**；需升格主表/图 |
| A | T0.3–0.5 | 效价/尺寸匹配子集；协变量 logistic；LE 主文化 | 否 | **待做** |
| A | T0.6–0.9 | 统一 θ；ChEMBL 聚合敏感；Murcko 分层；scaffold bootstrap | 否 | **待做** |
| A | T0.10 | 版本固定、文档矛盾修复、姿态上传清单 | 否 | **待做** |
| B | T1.1 | exhaustiveness 对照（优先：PM 全员 E=8 **或** AChE 全员 E=16） | 是 | **待做 · 第一优先对接** |
| B | T1.3 | 单靶 enrichment sanity（先 PM 两端） | 是 | **待做** |
| B | T1.4 | PM strict 扩面 ~110 | 是 | **待做** |
| C | T1.2/1.5/T2.* | decoy 受体 / 换结构 / holdout / 补面板 | 可选 | 冲刺 |

**禁止：** 扩 EGFR 对接；Track B 选臂；LigPrep 混进主表；湿实验；把 flags 灌进主分数。

---

```text
【角色】
你是 DualFourClass 本地执行 agent。依据：
  Dual_Target_Docking/data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md
  Dual_Target_Docking/docs/AGENT_COMMAND_JCIM_STRENGTHEN_SUPPLEMENT.md（本文件）
  Dual_Target_Docking/data/jcim_bench_v0/CLAIM_CEILING.md
目标：补齐 JCIM 评测文被审稿人打穿的缺口。不做“通用决策臂胜利”叙事。

【强制约束】
1. 无湿实验；不扩 EGFR 新对接；不混 LigPrep 进主表。
2. 新主指标 = 口袋匹配方向 AUROC：
     D vs A_only → 用口袋 B 分数
     D vs B_only → 用口袋 A 分数
   池化 vina_mean / worst-pocket / LE 归一 / 错口袋对照 = 并列报告，不得删掉旧对照。
3. 必报平凡基线（heavy/MW/cLogP/TPSA）+ bootstrap CI。
4. Prep = RDKit ETKDG + meeko；seed=20260727（面板抽样可用 20260729）；n_modes=9。
5. 对接作业写 protocol.yaml + job_status + 分数表；git commit/push；更新 PR。
6. 结束中文总结：完成了哪些 ID、路径、未完成项、建议用户下一步。

【工作目录】
Dual_Target_Docking/
产出总包：data/jcim_strengthen_t0t1_v0/
  tables/  analysis/  scripts/  protocols/

════════════════════════════════════════════════
A 组 — 零对接（先做完再开 B；云端/本地均可）
════════════════════════════════════════════════

产出目录：data/jcim_strengthen_t0t1_v0/

A1. 升格主表（基于已有 pocket_matched_*）
- 输入：
  data/jcim_bench_v0/tables/pocket_matched_directional_v1.csv
  data/jcim_bench_v0/tables/pocket_specificity_gap_v1.csv
  data/jcim_bench_v0/tables/pocket_matched_size_strata_v1.csv
  data/jcim_bench_v0/tables/forest_summary_min_ci_v1.csv
- 动作：写 analysis/PRIMARY_METRIC_V2.md
  - 明确：主指标改为口袋匹配；旧池化表降为对照
  - 给出 K=4 主结果表（口袋匹配 ± CI，含错口袋对照、LE、尺寸基线）
- 复现命令（已存在）：
  python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py

A2. 效价匹配 + 尺寸匹配子集（T0.3）— 新脚本
- 对每对靶：
  - potency-matched：D vs A_only 时按 pA 最近邻匹配（|ΔpA|≤0.5）；
    D vs B_only 时按 pB 匹配（|ΔpB|≤0.5）
  - size-matched：|Δheavy|≤2
- 输出表：tables/matched_subset_directional_v1.csv
- 在匹配子集上重算口袋匹配 AUROC + CI（B=2000）
- 写清匹配后每臂剩余 n；n<8 标 underpowered

A3. 协变量调整（T0.4）
- logistic：label(dual vs A_only) ~ score_B + heavy + tpsa（对称做 dual vs B_only ~ score_A + …）
- 报告 score 的系数 / OR / p；以及仅 score vs score+covariates 的 ΔAUROC
- 输出：tables/covariate_adjusted_v1.csv

A4. LE / 聚合敏感（T0.5 + 聚合）
- 并列：pocket-matched / wrong-pocket / worst-pocket / pooled-mean / LE-pocket-matched
- 输出：tables/aggregation_sensitivity_v1.csv（可复用 pocket_matched 表扩展）

A5. 统一 θ 重算（T0.6）
- 对四对靶用 θ∈{5.5,6.0,6.5} 及 strict(6.5/5.5) 重贴标签（需双侧 pChEMBL）
- 输出：tables/unified_threshold_sensitivity_v2.csv
- 主文建议：主表统一用一种规则，另一种进 SI

A6. ChEMBL 聚合敏感性（T0.7）— 若本地有 chembl 缓存
- 输入：data/public_pair_selection/mols_*.json
- 对比：现用 max pChEMBL vs median；可选 confidence≥8、Homo sapiens（若字段存在）
- 输出：tables/chembl_aggregation_sensitivity_v1.csv
- 若缺字段：写 SKIP 原因到 analysis/T0_SKIPS.md，不要假装完成

A7. Murcko 分层 + scaffold bootstrap（T0.8/T0.9）
- 真 Murcko（不要 chembl_id[:8] 代理）
- 表：tables/scaffold_inventory_v1.csv
- cluster/scaffold bootstrap：按支架重采样后再算 summary_min CI
- 输出：tables/scaffold_bootstrap_ci_v1.csv

A8. 配体 ML 基线（加分，原 W1.9）
- ECFP4 + Logistic/RF；嵌套 CV；仅 exploration 面板；禁止 peek holdout
- 与口袋匹配对接同表对比
- 输出：tables/ligand_ml_baseline_v1.csv

A9. 协议/文档修复（T0.10）
- 修 data/ache_bche_panel_v0/MANIFEST.md 与 GNINA 状态矛盾（以 jcim_bench_v0/analysis/GNINA_STATUS.md 为准）
- 新建 data/jcim_strengthen_t0t1_v0/ENV_PIN.md：记录 rdkit/meeko/vina/gnina/rtm 版本（用实际 `pip show`/`vina --version`）
- 新建 POSE_UPLOAD_CHECKLIST.md：列出本地 /mnt/d/... 姿态哪些必须进 Zenodo（至少 top1 pdbqt/sdf + 分数）
- 更新 CLAIM_CEILING 一句：主指标=口袋匹配方向 AUROC

【A 组完成标准】
- [ ] PRIMARY_METRIC_V2.md 存在且主表数字与 pocket_matched 一致
- [ ] matched_subset / covariate / unified_threshold / scaffold_bootstrap 四张核心表存在
- [ ] ENV_PIN + 文档矛盾已修
- [ ] analysis/A_GROUP_VERDICT.md：中文一页，说明主结论是否因口袋匹配而改变

════════════════════════════════════════════════
B 组 — 最小对接核（本地；约 380–540 Vina 端次）
════════════════════════════════════════════════

务必先完成 A1（主指标定义冻结），再开对接，避免用错汇总脚本。

────────────────────────────────
B1. T1.1 exhaustiveness 对照（第一优先）
────────────────────────────────
目的：排除「PM 成功只因为 E=16」。

优先方案（二选一，推荐 B1a）：
  B1a. 把 PIK3CA/mTOR panel48_rdkit 全员在 E=8 重跑两端
       - 输入：data/pik3ca_mtor_panel48_rdkit_v0/ 的配体名单 + 4L23/4JT6 盒子
       - 协议：E=8, n_modes=9, seed=20260727, RDKit meeko（与现 RDKit 配体可复用，不必重 prep）
       - 输出：data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_*.csv
                + analysis/EXHAUSTIVENESS_E8_VS_E16.md
       - 分析：同一配体集合上对比 E8 vs E16 的口袋匹配 summary_min 与 cognate RMSD
  B1b. 备选：AChE/BChE 现面板全员 E=16 重跑（~190 端）
       - 若 B1a 算力不足再做

完成标准：有并列表 + 明确结论「E 差异是否解释 PM 优势」。

────────────────────────────────
B2. T1.3 单靶 enrichment sanity（先 PM 两端）
────────────────────────────────
目的：证明对接在单靶上没坏，四类任务失败不是“Vina 废了”。

对受体 4L23、4JT6（必须）；可选再加 4EY7：
1. 从 ChEMBL 抽该靶 actives（建议 pChEMBL≥6.5，n≈50–100）与 property-matched decoys（n≈200–400；匹配 MW/logP/TPSA 或用现成 DUD-E/LIT-PCBA 若靶可得）
2. 同一盒子、E 与主协议一致（PM 用 E=16；对照写清）
3. 报 AUROC / EF1%/EF5%（活性 vs decoy）
4. 输出：data/jcim_strengthen_t0t1_v0/tables/single_target_enrichment_v1.csv
         + analysis/SINGLE_TARGET_SANITY.md

完成标准：PM 两端都有 enrichment 表；若 AUROC≪0.6 须在文中降权对接结论。

────────────────────────────────
B3. T1.4 PM strict 扩面 → N≈110–120
────────────────────────────────
目的：收窄 CI；尺寸分层目前每格仅 2–7 分子。

1. 按 strict 6.5/5.5 定额从 PIK3CA∩mTOR 配对库抽样（勿 silently 混 θ=6）
2. 配额建议：dual/A_only/B_only/neither ≈ 与现比例接近，总 N≈110–120；Murcko 真支架去冗余
3. 新名单：data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv（新目录，勿覆盖 48）
4. RDKit+meeko → Vina(E=16) → RTM best-of-9 → GNINA mode_01 rescore
5. 用口袋匹配指标重算主表 + 与 panel48 对照
6. 输出 analysis/PM110_VS_PM48.md

完成标准：panel110 分数齐；口袋匹配 CI 更新；声称仍受 CLAIM_CEILING 约束。

【B 组完成标准】
- [ ] B1 exhaustiveness 对照结论 md
- [ ] B2 PM 两端 enrichment
- [ ] B3 panel110 或书面说明跳过原因（算力）
- [ ] 汇总 analysis/B_GROUP_VERDICT.md

════════════════════════════════════════════════
C 组 — 可选冲刺（仅当 A+B 完成且算力充足）
════════════════════════════════════════════════

C1. T1.2 decoy 受体：PM48 配体对接到无关口袋（如 4EY7）→ 验证错口袋=分子固有属性
C2. T1.5 PM 换结构子集 40（替代 PDB 或换 4JT6 构象）
C3. T2.2 文献/临床双靶 holdout 名单冻结（可先不定对接）
C4. 禁止：扩 EGFR；无供给硬开第 5 对松标签

════════════════════════════════════════════════
交付与 git
════════════════════════════════════════════════

每完成一大组：
  git add Dual_Target_Docking/data/jcim_strengthen_t0t1_v0 \
          Dual_Target_Docking/data/jcim_bench_v0 \
          Dual_Target_Docking/data/pik3ca_mtor_* \
          Dual_Target_Docking/docs/
  git commit -m "meaningful message"
  git push

结束时中文报告：
1. 完成了哪些 T0/T1 ID
2. 主指标从池化改为口袋匹配后，四对靶结论如何变
3. exhaustiveness / enrichment / PM110 的关键数字
4. 还缺什么才能开写英文稿
```

---

## 本地一键检查（agent 开跑前）

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
ls data/jcim_bench_v0/tables/pocket_matched_directional_v1.csv
ls data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md
```

## 推荐执行口令（给用户选）

| 你对本地 agent 说 | 含义 |
|-------------------|------|
| `本轮只做 A 组` | 零对接，先改主表与混淆对照 |
| `A 组 + B1 + B2` | 分析 + exhaustiveness + 单靶 sanity（性价比最高） |
| `A+B 全做` | 含 PM110 扩面（推荐最小可投） |
| `A+B+C 全做` | 冲刺 |

**默认建议口令：** `A 组 + B1 + B2；B3 若本机夜间空闲再做`
