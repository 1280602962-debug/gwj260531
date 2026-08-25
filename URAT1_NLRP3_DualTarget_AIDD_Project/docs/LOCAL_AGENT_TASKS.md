# 本地 Agent 执行任务书：URAT1–NLRP3 双靶重定位 · 阶段二（漏斗 + 短名单 + MD 文件导出 + 写作）

> 面向本地漏斗执行。生产协议已锁定为 **Π\* = P2**。  
> 写作入口：`docs/MANUSCRIPT.md`。当前 MD 假说分子见该文件（GSK-3008348、Vecabrutinib + 对照）。仓库不再提交自动挑选的 MD CSV。

---

## 背景（只读，不要重新决策）

- 课题：痛风 URAT1（代谢，转运体）+ NLRP3（炎症，NACHT）双节点，临床药物库重定位计算筛选。
- 已完成并**锁定**：在 URAT1 TrueDecoy（主判）+ RandomDecoy（否决对照）基准上比较了 P0–P5 六种开源对接/重打分协议，选定生产协议：

  **Π\* = P2（gnina, CNNaffinity, cnn_scoring=rescore, exhaustiveness=32）**

  理由（完整数据见 `docs/PROTOCOL_SELECTION_RESULT.md`）：
  - P2 在 TrueDecoy 上早期富集统计显著（锁定点估计 EF@1%=2.54，12/52；归档 \(\lfloor 0.01N\rfloor\) 重算 2.59，12/51，bootstrap 95% CI 1.31–4.07，超几何 p≈0.0016），且在 RandomDecoy 上非零（EF@1%≈0.22）。
  - P5（RTMScore/gnina 构象）虽然 TrueDecoy EF@1% 更高（2.80），但 **RandomDecoy EF@1% = 0**，且该失败经检验是统计显著的真实失败（纯随机排序下出现 0 命中的概率仅约 0.76%），不是噪声。临床库比 RandomDecoy 更像"多样、非匹配"场景，因此 **不选 P5 作生产协议**，只作敏感性分析。
  - P0（gnina CNNscore）两侧最均衡但预注册为负对照，不提为主协议。
  - P1（Vina affinity）、P3（gnina minimizedAffinity）、P4（RTMScore/Vina 构象）与随机无统计显著差异或覆盖不全，排除。
- **不要重新跑 TrueDecoy/RandomDecoy 协议筛选**，也不要重新讨论选哪个协议——这一步已经结束。
- 投稿形态（τ 敏感性、对照强制对接、MD 对照优先）见 [`MOL_DIVERS_REVISION_PLAN.md`](MOL_DIVERS_REVISION_PLAN.md)。τ=80/85 表在 `data/repurposing/p2/sensitivity_tau80/` 与 `sensitivity_tau85/`，**不要覆盖**生产 τ=90 提名文件。
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

## 任务 1：1588 临床池双靶对接（P2，9DKB + 7ALV）— **已归档**

生产 P2 分数已写入 `data/repurposing/p2/`（及完整姿态包 `docking_export_20260820/`）。不要重对接 1588。

**产出（已存在）：**
```
data/repurposing/p2/docking_9dkb_gnina.csv          # 1582 docked / 1583
data/repurposing/p2/docking_7alv_gnina.csv
data/repurposing/p2/pareto_merged_scores.csv        # 1580 完整案例
data/repurposing/p2/pareto_shortlist.csv            # 4，大环，仅审计
data/repurposing/p2/nominated_shortlist_diverse.csv # 7 优选；跟进 GSK-3008348 + Vecabrutinib
```

若需在新机器重跑（非默认）：`JOBS=8 bash scripts/run_funnel_p2.sh`。**不要**把 `pareto_shortlist.csv` 当成故事分子；跟进读提名表。

---

## 任务 2：审计流水线 — **已归档**

表在 `data/repurposing/p2/`。复现：

```bash
python3 scripts/11_chemical_space_novelty.py --pool data/repurposing/p2/pareto_merged_scores.csv --shortlist data/repurposing/p2/pareto_shortlist.csv --output-dir data/repurposing/p2
python3 scripts/13_pareto_robustness.py --pool data/repurposing/p2/pareto_merged_scores.csv --output-dir data/repurposing/p2/pareto_robustness
python3 scripts/14_candidate_nomination.py --tau 90 --mw-min 200 --mw-max 550 --top-diverse 12
```

已核对：双结构门控 51；优选 7（Veber + Ro5 HBD/HBA/logP + MW 200–550）；裸 Pareto 为 4 个大环；lesinurad/verinurad/colchicine 不在门控内。跟进 GSK-3008348 与 Vecabrutinib。

---

## 任务 3：MD 输入文件（受体/配体导出）

> 本任务只导出起始构象文件。轨迹在有算力的机器上跑。  
> **当前跟进分子（P2 化学提名）：** GSK-3008348（URAT1 侧）、Vecabrutinib（NLRP3 侧）；对照 lesinurad @ 9DKB、MCC950 @ 7ALV（类似物对照姿见 `data/si/mcc950_7alv/`，非自对接）。  
> **不要 MD：** Zelenirstat、MLN-0415、BI 653048、Deucrictibant、Praliciguat，以及大环内酯 / 多酚类对接优势分子。  
> URAT1 必须按 **膜+脂双层** 体系；7ALV 用水盒子。

若要用脚本从 **P2 提名表** 自动挑选：

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

1. 引言：`docs/INTRO_DRAFT_CN.md`；Methods：`docs/METHODS_DRAFT_CN.md`；Results：`docs/RESULTS_DRAFT_CN.md`。生产协议 = P2。
2. Results 已写入：协议表、漏斗 8319→1588→1580、裸 Pareto 4 vs 双结构门控 51 / 优选 7、姿态 QC；MD 数值不报。
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
- 不要把已删除的历史对接短名单或 EGCG / canagliflozin 写成当前 lead。
- 不要编造尚未完成的 MD 数值。
- 不要使用确认性发现语言（"we identified dual-target inhibitors"）。
- 不要删除或覆盖 `docs/PROTOCOL_SELECTION_RESULT.md` 中的既有结果表。

---

## 完成标准（Definition of Done）

- [x] 任务 1：`data/repurposing/p2/pareto_merged_scores.csv` 为 1,580 完整案例
- [x] 任务 2：双结构门控 51 / 优选 7，跟进 GSK-3008348 与 Vecabrutinib
- [x] 任务 3：姿态 QC 已归档；MD 轨迹数值未报
- [x] 任务 4：按 `docs/MANUSCRIPT.md` 把 R1–R4 写入 `docs/RESULTS_DRAFT_CN.md`
- [x] 所有任务已分别提交并推送到 `cursor/urat1-nlrp3-dualtarget-aidd-e43d`
