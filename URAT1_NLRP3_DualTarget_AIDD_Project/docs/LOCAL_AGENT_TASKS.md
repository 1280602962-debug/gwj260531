# 本地 Agent 执行任务书：URAT1–NLRP3 双靶重定位 · 阶段二（漏斗 + 短名单 + MD 文件导出 + 写作）

> 面向：在本地机器上运行的 coding agent（有 shell、文件读写、git 权限）。
> 前提：`vina`、`gnina` 已安装并可在 PATH 或 `tools/gnina` 找到；Python 依赖（rdkit、meeko、gemmi、openbabel、pandas、numpy、pyyaml、scikit-learn、xgboost、scipy）已安装；仓库为 `URAT1_NLRP3_DualTarget_AIDD_Project`（Git 仓库根目录的子目录）。
> **本机算力不足以跑 MD**：任务 3 只产出 MD 输入文件，不在本地运行任何 MD 引擎。
> 你不需要读其他对话记录，本文档自包含。**严格按顺序执行任务 0 → 4**，每个任务完成后按"提交"一节要求提交并推送。

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

## 任务 3：挑选 MD 候选分子，导出可用于 MD 的蛋白/配体文件（**不在本地跑 MD**）

> 本机算力不足以跑 MD。本任务只做两件事：**挑分子** + **产出标准化的受体/配体文件**，
> 交给有算力的机器（云端/工作站/HPC）去跑 MD。**不要在本地尝试运行任何 MD 引擎。**

### 3.1 挑选候选分子

```bash
python3 scripts/select_md_candidates.py \
  --n-novel 4 \
  --n-controls 2 \
  --output data/md_candidates/md_candidate_selection.csv
```

- 默认从 `results/candidates/nominated_candidates.csv`（任务 2 产物）挑选：
  - 最多 4 个 **novel_candidate**：`clean_candidate=True`（无 PAINS/Brenk，通过 Lipinski+Veber），按 `dual_structure_balance` 降序
  - 最多 2 个 **known_control**：已知对照药（优先 lesinurad、MCC950 等），用于校准解读
- 脚本会自动把候选分子关联回 `repurposing_id` / `canonical_smiles` / 对接状态（通过 `data/repurposing/pareto/pareto_merged_scores.csv` 联表），不需要手工处理。
- 若某候选缺 `repurposing_id`（联表失败）会打印 WARNING 并自动剔除，不会中断执行。

**检查**：`data/md_candidates/md_candidate_selection.csv` 存在，行数在 4–6 之间，`has_9dkb_pose` / `has_7alv_pose` 至少一列为 True。

### 3.2 导出 MD-ready 文件（受体 PDB + 配体 SDF/PDB/SMILES）

```bash
python3 scripts/export_md_ready_candidates.py \
  --selection data/md_candidates/md_candidate_selection.csv \
  --output-dir data/md_candidates
```

对每个（化合物, 靶点）组合，只要该化合物在该靶点对接成功，就会生成一个文件夹：

```
data/md_candidates/
  _receptors/
    9DKB_receptor.pdb        # 受体只生成一次，两个靶点共用同一份
    7ALV_receptor.pdb
  9DKB_<repurposing_id>/
    receptor.pdb             # 与该靶点对接时使用的同一受体（蛋白质、去水、去异原子、pH 7.4 加氢）
    ligand.sdf                # gnina P2 产出的对接姿态原始文件
    ligand.pdb                # 同一姿态转成 PDB，便于可视化/组装复合物
    ligand.smi                # 真实标准 SMILES（供后续力场参数化用）
    README.md                 # 来源、对接分数、蛋白质子化 pH、后续 MD 建议步骤
  7ALV_<repurposing_id>/
    ...（同上）
  md_ready_manifest.csv       # 全部导出文件的清单（含 warning 列）
```

**检查**：
- `data/md_candidates/md_ready_manifest.csv` 存在，行数与任务 3.1 选出的（化合物×成功靶点）组合数一致。
- 逐个打开 `warning` 列，若非空需要处理（常见原因：gnina 姿态 SDF 缺失或格式异常，需要回到任务 1 的对接输出核查）。
- 每个文件夹下 `receptor.pdb`、`ligand.sdf`、`ligand.smi` 必须存在；`ligand.pdb` 若因分子解析问题缺失，不阻塞交付，但需在提交说明里注明。

### 3.3 交付说明（写清楚，方便你或云端 agent接手跑 MD）

在 `data/md_candidates/` 下新建 `HANDOFF_NOTES.md`，写明：
- 本批候选分子列表（复制 `md_candidate_selection.csv` 的关键列）
- 每个文件夹里文件的含义（可直接引用上面的目录结构说明）
- **明确声明：这里只有起始构象/文件，MD 尚未运行**，力场选择、复合物组装、溶剂化、平衡、生产阶段都留给下一步执行者
- 建议目标：每个体系 50–100 ns 生产阶段；报告 RMSD/RMSF、关键残基相互作用；MM-GBSA（如做）仅同批相对比较

**不要**在这一步编造任何 MD 数值结果——本任务只产出输入文件，不产出轨迹或能量数字。

---

## 任务 4：更新文档与图表（Methods / Results 草稿）

1. 在 `docs/METHODS_DRAFT_CN.md` 中补充：
   - 生产协议 = P2（已在 `docs/PROTOCOL_SELECTION_RESULT.md` 定义，直接引用即可，不要重写协议筛选逻辑）。
   - 1588 双靶对接的具体执行细节（若与已写内容不一致，以本次实际运行参数为准）。
   - 说明 MD 体系已完成"输入文件准备"（受体/配体导出，见 `data/md_candidates/`），MD 本身**尚未运行**，将在外部算力环境执行；不得编造 MD 参数或结果数值。
2. 在 `docs/RESULTS_DOCKING_9DKB_7ALV.md`（或新建 `docs/RESULTS_FUNNEL_P2.md`）中记录：
   - Pareto 短名单规模、提名后剩余分子数
   - MD 候选分子选择结果（`md_candidate_selection.csv` 摘要：哪几个 novel_candidate、哪几个 known_control，及入选理由）
   - 已知对照药在漏斗中的位置（回收情况）
3. 更新 `README.md` 中"实现状态"表：
   - 把"重定位库双靶对接"从 `⏳` 改为 `✅`
   - 新增一行"MD 候选筛选 + 输入文件导出"标记 `✅`，"MD 模拟本身"标记 `⏳ 待外部算力执行`
   - 链接新文档

**写作口径（强制）：**
- 通篇使用 "computational dual-node repurposing hypotheses"、"candidate nominations pending experimental validation" 一类措辞。
- **禁止**出现 "identified dual-target inhibitors"、"validated hits"、"potent dual inhibitors" 等确认性表述。
- **禁止**声称 MD 已完成或报告任何 MD 数值（RMSD、MM-GBSA 等）——本阶段只完成到"MD 输入文件已备妥"。
- 局限段必须包含：TrueDecoy/RandomDecoy 富集中等（AUC≈0.58–0.65）；P5 在 RandomDecoy 上失败（写明数值和显著性）；诱饵为库分子而非实验无活；MD 尚待外部算力执行，当前无构象稳定性证据。

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
- **不要在本地运行任何 MD 引擎**（GROMACS/AMBER/NAMD 等）；本机算力不支持，任务 3 只导出输入文件。
- 不要编造 MD 参数、时长或结果数值；MD 尚未运行时，文档一律标注"待外部算力执行"，不得写假数字。
- 不要在任何文档或摘要中使用确认性发现语言（"we identified/discovered dual-target inhibitors"）。
- 不要删除或覆盖 `docs/PROTOCOL_SELECTION_RESULT.md` 中的既有结果表。

---

## 完成标准（Definition of Done）

- [ ] 任务 1：`pareto_shortlist.csv` 存在且非空，对接失败率 <10%
- [ ] 任务 2：`results/candidates/` 下有最终提名表，含 clean candidate 标记
- [ ] 任务 3：`data/md_candidates/md_candidate_selection.csv`（4–6 个候选）+ `md_ready_manifest.csv` + 每个文件夹下 `receptor.pdb`/`ligand.sdf`/`ligand.smi` 齐全，`HANDOFF_NOTES.md` 已写明"MD 尚未运行"
- [ ] 任务 4：Methods/Results 文档更新，README 状态表同步，全篇无确认性发现语言、无编造的 MD 数值
- [ ] 所有任务已分别提交并推送到 `cursor/urat1-nlrp3-dualtarget-aidd-e43d`
