# 本地 Agent 执行任务书：URAT1–NLRP3 双靶重定位 · 阶段二（漏斗 + 短名单 + MD + 写作）

> 面向：在本地机器上运行的 coding agent（有 shell、文件读写、git 权限）。
> 前提：`vina`、`gnina` 已安装并可在 PATH 或 `tools/gnina` 找到；Python 依赖（rdkit、meeko、pandas、numpy、pyyaml、scikit-learn、xgboost、scipy）已安装；仓库为 `URAT1_NLRP3_DualTarget_AIDD_Project`（Git 仓库根目录的子目录）。
> 你不需要读其他对话记录，本文档自包含。**严格按顺序执行任务 0 → 6**，每个任务完成后按“提交”一节要求提交并推送。

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
data/repurposing/pareto/pareto_shortlist.csv
```

若 `docking_status` 列里失败（非 `docked`）比例超过 ~10%，先排查（常见原因：配体准备失败、gnina 超时），修复后重跑，不要带着大面积失败继续下一步。

---

## 任务 2：审计流水线（PAINS / ADMET / 化学空间 / 稳健性 / 提名）

按顺序执行，全部使用默认路径（已与任务 1 的输出对齐）：

```bash
python3 scripts/10_admet_druglikeness.py
python3 scripts/11_chemical_space_novelty.py
python3 scripts/13_pareto_robustness.py
python3 scripts/14_candidate_nomination.py
```

**产出检查：**
```
results/cheminformatics/...        # ADMET / drug-likeness / novelty 输出
results/pareto_robustness/...      # bootstrap 敏感性
results/candidates/...             # 最终提名表（含 clean candidate 标记）
```

打开 `results/candidates/` 下的提名表，确认：
- 已知对照药（lesinurad、benzbromarone、verinurad、dotinurad、MCC950、GDC-2394、allopurinol、colchicine）被正确标记为"已知对照"而非"新提名"。
- PAINS/Brenk 命中分子已标记降级原因，不与 clean candidate 混为一谈。

---

## 任务 3：挑选 2–4 个代表 lead 做 MD

从 `results/candidates/` 的最终提名表中，按以下优先级挑 **2–4 个**分子：

1. 双轴百分位（S_U 与 S_N）都较高
2. 无 PAINS、无 Brenk
3. 同时通过 Lipinski 与 Veber
4. 非已知对照药（新颖重定位候选优先；若名单里对照药排名很靠前，也可选 1 个对照药作为阳性基线）

对每个 lead：

- **URAT1 侧**：以 9DKB 复合物为体系，配体姿态优先用晶体（若为已知共晶药物）或 gnina/RTMScore 选出的构象（不要用富集协议里已知失败的 Top-1 姿态硬凑）。
- **NLRP3 侧**：以 7ALV 复合物为体系，同样规则选择起始构象。
- MD 设置：**目标 50–100 ns 生产阶段**（按你本地算力可行范围，优先保证质量而非硬凑到 100ns）；报告：
  - RMSD（蛋白骨架 + 配体）随时间曲线
  - 关键残基距离/氢键/盐桥（列出具体残基编号）
  - 可选：MM-GBSA（**仅作同一批复合物之间的相对比较**，不得报告为绝对结合自由能）
- **对照体系**（每个靶点至少 1 条，用于校准解读）：
  - URAT1：lesinurad @ 9DKB
  - NLRP3：MCC950（或 7ALV 共晶类似物）@ 7ALV

MD 工具（GROMACS/AMBER 等）本仓库未内置自动化脚本，按你本地环境常规流程执行；产出轨迹分析图和数值表，保存到：

```
results/md/<pdb>_<compound_name>/rmsd.csv
results/md/<pdb>_<compound_name>/key_contacts.csv
results/md/<pdb>_<compound_name>/summary.md   # 简述体系、力场、时长、关键发现
figures/generated/md/<pdb>_<compound_name>_rmsd.png
```

---

## 任务 4：更新文档与图表（Methods / Results 草稿）

1. 在 `docs/METHODS_DRAFT_CN.md` 中补充：
   - 生产协议 = P2（已在 `docs/PROTOCOL_SELECTION_RESULT.md` 定义，直接引用即可，不要重写协议筛选逻辑）。
   - 1588 双靶对接的具体执行细节（若与已写内容不一致，以本次实际运行参数为准）。
   - MD 体系、力场、时长、生产阶段设置（**必须是真实跑出来的数值，不得编造**）。
2. 在 `docs/RESULTS_DOCKING_9DKB_7ALV.md`（或新建 `docs/RESULTS_FUNNEL_P2.md`）中记录：
   - Pareto 短名单规模、提名后剩余分子数
   - 2–4 个 lead 的 MD 结果摘要
   - 已知对照药在漏斗中的位置（回收情况）
3. 更新 `README.md` 中"实现状态"表，把"重定位库双靶对接"从 `⏳` 改为 `✅`，并链接新文档。

**写作口径（强制）：**
- 通篇使用 "computational dual-node repurposing hypotheses"、"candidate nominations pending experimental validation" 一类措辞。
- **禁止**出现 "identified dual-target inhibitors"、"validated hits"、"potent dual inhibitors" 等确认性表述。
- 局限段必须包含：TrueDecoy/RandomDecoy 富集中等（AUC≈0.58–0.65）；P5 在 RandomDecoy 上失败（写明数值和显著性）；诱饵为库分子而非实验无活；MM-GBSA（如有）仅相对比较；单一起始构象/短时程 MD 的局限。

---

## 提交规范（每个任务完成后执行）

```bash
git add <本任务涉及的文件>
git commit -m "<清晰描述本任务做了什么>"
git push -u origin cursor/urat1-nlrp3-dualtarget-aidd-e43d
```

- 任务 1 完成后单独提交一次（对接结果 + Pareto 输出）。
- 任务 2 完成后单独提交一次（审计结果）。
- 任务 3 完成后单独提交一次（MD 结果，注意大轨迹文件不要直接提交到 git，只提交分析 CSV/图/摘要；若轨迹很大，写清存放路径但加入 `.gitignore`）。
- 任务 4 完成后单独提交一次（文档更新）。
- **不要把任务 0–4 合并成一个大提交**；每个逻辑步骤单独提交，方便回溯。

---

## 禁止事项清单

- 不要重新跑或重新讨论 TrueDecoy/RandomDecoy 协议筛选（P0–P5 已锁定为 P2）。
- 不要把 P5 提升为生产协议。
- 不要对 `true_decoy_benchmark.csv` / `random_decoy_benchmark.csv` 重复对接。
- 不要编造 MD 参数、时长或结果数值；没跑完就如实标注"进行中/待补"，不要在文档里写假数字。
- 不要在任何文档或摘要中使用确认性发现语言（"we identified/discovered dual-target inhibitors"）。
- 不要删除或覆盖 `docs/PROTOCOL_SELECTION_RESULT.md` 中的既有结果表。

---

## 完成标准（Definition of Done）

- [ ] 任务 1：`pareto_shortlist.csv` 存在且非空，对接失败率 <10%
- [ ] 任务 2：`results/candidates/` 下有最终提名表，含 clean candidate 标记
- [ ] 任务 3：至少 2 个 lead + 2 个对照体系完成 MD，产出 RMSD/相互作用摘要
- [ ] 任务 4：Methods/Results 文档更新，README 状态表同步，全篇无确认性发现语言
- [ ] 所有任务已分别提交并推送到 `cursor/urat1-nlrp3-dualtarget-aidd-e43d`
