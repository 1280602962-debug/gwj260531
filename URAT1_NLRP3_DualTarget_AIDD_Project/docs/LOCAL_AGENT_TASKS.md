# 本地 Agent 执行任务书：URAT1–NLRP3 双靶重定位 · 阶段二（漏斗 + 短名单 + MD 文件导出 + 写作）

> 面向本地漏斗执行。生产协议已锁定为 **Π\* = P2**。  
> 写作入口：`docs/MANUSCRIPT.md`。当前 MD 假说分子见该文件（GSK-3008348、Vecabrutinib + 对照）；**不要**用仓库里 Glide 时代的 `pareto_shortlist.csv` / `md_candidate_selection.csv`。

---

## 背景（只读，不要重新决策）

- 课题：痛风 URAT1（代谢，转运体）+ NLRP3（炎症，NACHT）双节点，临床药物库重定位计算筛选。
- 已完成并**锁定**：在 URAT1 TrueDecoy（主判）+ RandomDecoy（否决对照）基准上比较了 P0–P5 六种开源对接/重打分协议，选定生产协议：

  **Π\* = P2（gnina, CNNaffinity, cnn_scoring=rescore, exhaustiveness=32）**

  理由（完整数据见 `docs/PROTOCOL_SELECTION_RESULT.md`）：
  - P2 在 TrueDecoy 上早期富集统计显著（EF@1%=2.54，超几何检验 p≈0.002），且在 RandomDecoy 上非零（EF@1%=0.21）。
  - P5（RTMScore/gnina 构象）虽然 TrueDecoy EF@1% 更高（2.80），但 **RandomDecoy EF@1% = 0**，且该失败经检验是统计显著的真实失败（纯随机排序下出现 0 命中的概率仅约 0.76%），不是噪声。临床库比 RandomDecoy 更像"多样、非匹配"场景，因此 **不选 P5 作生产协议**，只作敏感性分析。
  - P0（gnina CNNscore）两侧最均衡但预注册为负对照，不提为主协议。
  - P1（Vina affinity）、P3（gnina minimizedAffinity）、P4（RTMScore/Vina 构象）与随机无统计显著差异或覆盖不全，排除。
- **不要重新跑 TrueDecoy/RandomDecoy 协议筛选**，也不要重新讨论选哪个协议——这一步已经结束。
- **不要**把结果包装成"发现了双靶抑制剂"；全篇的主张必须是 **computational hypothesis-generating screen**（计算假说生成），需要实验验证。

---

## 任务 0：环境自检（5 分钟）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
git fetch origin cursor/urat1-nlrp3-dualtarget-aidd-e43d
git checkout cursor/urat1-nlrp3-dualtarget-aidd-e43d
git pull origin cursor/urat1-nlrp3-dualtarget-aidd-e43d

which vina || echo "MISSING vina"
which gnina || ls tools/gnina || echo "MISSING gnina"
python3 -c "import rdkit, meeko, pandas, numpy, yaml, sklearn, scipy; print('deps OK')"

ls data/repurposing/screening/docking_pool_p05.csv
ls data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv
```

若任一二进制或依赖缺失，先按 `docs/OPEN_SOURCE_DOCKING.md` 和 `scripts/setup_gnina_wsl_cpu.sh` 安装，再继续。**不要跳过检查直接开跑大批量对接。**

---

## 任务 1：1588 临床池双靶对接（P2，9DKB + 7ALV）

一键脚本已就绪：`scripts/run_funnel_p2.sh`。它会：
1. 准备受体（若未准备）
2. 准备配体（`docking_pool_p05.csv`，约 1588 个）
3. 用 gnina + CNNaffinity 对 9DKB、7ALV 各对接一次
4. 合并打分、算 Pareto 短名单，写入 `data/repurposing/pareto/`

```bash
JOBS=8 bash scripts/run_funnel_p2.sh
```

- `JOBS` 按你机器核数调整（建议核数的一半到全部）。
- 单次运行即可覆盖两个靶点；**不要额外对 True/Random 基准重跑**。
- CPU-only 粗估：单分子 gnina exh=32 约 20–60 秒/靶点，1588×2 建议规划数小时到一晚，视核数。

**产出检查（必须存在且非空）：**
```
results/repurposing/docking_p2/9dkb/docking_9dkb_gnina.csv
results/repurposing/docking_p2/7alv/docking_7alv_gnina.csv
data/repurposing/pareto/pareto_merged_scores.csv
data/repurposing/pareto/pareto_shortlist.csv            # 原始对接 Pareto（仅审计）
data/repurposing/pareto/pareto_shortlist_druglike.csv   # Pareto ∩ MW 200–550
```

若 `docking_status` 列里失败（非 `docked`）比例超过 ~10%，先排查（常见原因：配体准备失败、gnina 超时），修复后重跑，不要带着大面积失败继续下一步。

**重要：** 不要把 `pareto_shortlist.csv` 当成故事分子/MD 候选。对接 Pareto 常被大环内酯等高分子量分子占据；跟进名单以任务 2 的化学优先提名为准。

---

## 任务 2：审计流水线（PAINS / ADMET / 化学空间 / 稳健性 / 提名）

按顺序执行，全部使用默认路径（已与任务 1 的输出对齐）：

```bash
python3 scripts/10_admet_druglikeness.py
python3 scripts/11_chemical_space_novelty.py
python3 scripts/13_pareto_robustness.py
python3 scripts/14_candidate_nomination.py --tau 90 --mw-min 200 --mw-max 550 --top-diverse 12
```

**产出检查：**
```
results/cheminformatics/...                              # ADMET / drug-likeness / novelty 输出
results/pareto_robustness/...                            # bootstrap 敏感性
results/candidates/nominated_candidates.csv              # 双阈值全集 + 化学标记
results/candidates/nominated_shortlist_diverse.csv       # 骨架去冗余后的跟进短名单（主读此表）
results/candidates/candidate_nomination_summary.json
```

打开提名表时确认：
- **跟进分子读 `nominated_shortlist_diverse.csv`**，不是原始 Pareto。
- `preferred_candidate=True`：无 PAINS/Brenk、Lipinski+Veber、MW∈[200,550]、口服吸收替代标志通过。
- 大环/高分子量对接优势分子应落在 `mw_oral_ok=False` 或非 preferred，并在 summary 的 `demoted_high_mw_in_gate` 可见。
- 已知对照药（lesinurad、benzbromarone、verinurad、dotinurad、MCC950、GDC-2394、allopurinol、colchicine）被正确标记为"已知对照"而非"新提名"。
- PAINS/Brenk 命中分子已标记降级原因，不与 clean/preferred candidate 混为一谈。

---

## 任务 3：MD 输入文件（受体/配体导出）

> 本任务只导出起始构象文件。轨迹在有算力的机器上跑。  
> **当前跟进分子（P2 化学提名，非 Glide 裸 Pareto）：** GSK-3008348（URAT1 侧）、Vecabrutinib（NLRP3 侧）；对照 lesinurad @ 9DKB、MCC950 @ 7ALV（若有姿）。  
> **不要 MD：** Zelenirstat、MLN-0415、BI 653048、Deucrictibant、Praliciguat，以及仓库 Glide 短名单中的 EGCG / canagliflozin / 大环内酯。  
> URAT1 必须按 **膜+脂双层** 体系；7ALV 用水盒子。

若要用脚本从 **P2 提名表** 自动挑选（不要喂 Glide 时代 `pareto_shortlist.csv`）：

```bash
python3 scripts/select_md_candidates.py \
  --n-novel 2 \
  --n-controls 2 \
  --output data/md_candidates/md_candidate_selection.csv
```

然后：

```bash
python3 scripts/export_md_ready_candidates.py \
  --selection data/md_candidates/md_candidate_selection.csv \
  --output-dir data/md_candidates
```

交付时写明每个文件夹含义；**不要编造尚未完成的 MD 数值**。已在跑的轨迹填入 Results 时再写数字。

---

## 任务 4：文稿（按 `docs/MANUSCRIPT.md`）

1. 引言：`docs/INTRO_DRAFT_CN.md`；Methods：`docs/METHODS_DRAFT_CN.md`。生产协议 = P2。
2. Results 记录：协议表、漏斗计数、裸 Pareto vs 化学提名、MD 候选及理由。
3. URAT1 MD 按膜蛋白体系写方法；不得把对接分写成亲和力，不得写已发现双靶抑制剂。

---

## 提交规范（每个任务完成后执行）

```bash
git add <本任务涉及的文件>
git commit -m "<清晰描述本任务做了什么>"
git push -u origin cursor/urat1-nlrp3-dualtarget-aidd-e43d
```

- 任务 1 完成后单独提交一次（对接结果 + Pareto 输出）。
- 任务 2 完成后单独提交一次（审计结果）。
- 任务 3 完成后单独提交一次（`data/md_candidates/` 下的选择表、受体/配体文件、manifest、`HANDOFF_NOTES.md`）。
- 任务 4 完成后单独提交一次（文档更新）。
- **不要把任务 0–4 合并成一个大提交**；每个逻辑步骤单独提交，方便回溯。

---

## 禁止事项清单

- 不要重新跑或重新讨论 TrueDecoy/RandomDecoy 协议筛选（P0–P5 已锁定为 P2）。
- 不要把 P5 提升为生产协议。
- 不要对 `true_decoy_benchmark.csv` / `random_decoy_benchmark.csv` 重复对接。
- 不要把仓库 Glide 时代 Pareto / EGCG / canagliflozin 写成当前 lead。
- 不要编造尚未完成的 MD 数值。
- 不要使用确认性发现语言（"we identified dual-target inhibitors"）。
- 不要删除或覆盖 `docs/PROTOCOL_SELECTION_RESULT.md` 中的既有结果表。

---

## 完成标准（Definition of Done）

- [ ] 任务 1：`pareto_shortlist.csv` 存在且非空，对接失败率 <10%
- [ ] 任务 2：`results/candidates/` 下有最终提名表，含 clean candidate 标记
- [ ] 任务 3：MD 输入与当前跟进分子一致（GSK-3008348、Vecabrutinib + 对照），不是 Glide 短名单
- [ ] 任务 4：按 `docs/MANUSCRIPT.md` 更新；无确认性发现语言、无编造 MD 数值
- [ ] 所有任务已分别提交并推送到 `cursor/urat1-nlrp3-dualtarget-aidd-e43d`
