# 本地 Agent：8973 XP 对接结果整理 + URAT1 模型/证据提升

**前提**：你已完成 `distill_manifest.csv`（8973 条）在 **9DKB** 上的 Glide **XP** 对接（Maestro 导出 CSV/Excel）。

**核心原则**（项目已定稿）：
- **URAT1 主证据 = 9DKB XP 分数 $S_U$**，不是重训 ML
- ML 当前 **benchmark 2/4（URAT1_NO_GO）**；对接回顾性分析才是主文 Part A
- **不要**把 8973 对接结果与 NLRP3 做几何平均融合当主创新
- **不要**再跑 OAT 迁移（Δρ≈0，已从叙事删除）

---

## 一、人工准备

1. 从 Maestro 导出 **8973 @ 9DKB XP** 分数表（可一个或多个 CSV/Excel）  
2. 放到例如：

```text
results/docking/raw/
  9DKB_xp_scores.csv          # 或你实际的 Maestro 导出名
```

导出列至少包含：**SMILES** + **Glide XP GScore**（或 `r_i_glide_XP`）+ **pose 状态**（可选）

3. `git pull` 拿到脚本：
   - `scripts/merge_8973_docking_results.py`
   - `scripts/analyze_urat1_docking_vs_ml.py`
   - `scripts/00c_augment_urat1_literature_benchmarks.py`

---

## 二、直接复制给本地 Agent 的整段命令

```
你在离线环境工作，不要 git clone / push。

项目根目录：URAT1_NLRP3_DualTarget_AIDD_Project/

背景：8973 条（distill_manifest.csv）9DKB Glide XP 对接已完成。
目标：(1) 清洗合并对接分数与 manifest；(2) 做 URAT1 回顾性验证（活性 vs decoy、四药 benchmark）；
(3) 在诚实前提下尝试有限 ML 提升（文献标签增补 ablation），但主结论仍用对接 S_U。

═══════════════════════════════════════
Phase 1 — 对接结果整理（必做）
═══════════════════════════════════════

【输入】
- data/distill/distill_manifest.csv（8973）
- Maestro 导出：results/docking/raw/9DKB_xp_scores.csv
  （若分多个文件，全部传给 --glide-csv）

【执行】
python3 scripts/merge_8973_docking_results.py \
  --glide-csv results/docking/raw/9DKB_xp_scores.csv \
  --pdb 9DKB \
  --output-dir data/docking

【验收】
- data/docking/8973_9DKB_with_manifest.csv 生成
- data/docking/8973_docking_qc_summary.json 中：
  * coverage_pct ≥ 95%（若低，查 SMILES 规范化/失败 pose）
  * subset A(822) 与 D(8000) 均有覆盖
  * benchmark_panel_9dkb 中四药 lesinurad/benzbromarone/verinurad/dotinurad 尽量 docked=True

【若列名不匹配】
打印 CSV 列名，在 merge_8973_docking_results.py 的 SMILES_ALIASES / SCORE_ALIASES 中补充后重跑。

═══════════════════════════════════════
Phase 2 — 对接 vs ML 回顾性分析（必做，主文 Part A）
═══════════════════════════════════════

【前置】若尚无训练模型，先跑：
python3 scripts/00_prepare_data.py
python3 scripts/02_train_asymmetric_models.py --no-oat-transfer

【执行】
python3 scripts/analyze_urat1_docking_vs_ml.py \
  --merged data/docking/8973_9DKB_with_manifest.csv \
  --hybrid-cv

【输出】
- results/docking/urat1_docking_vs_ml_summary.json
- results/docking/urat1_benchmark_rankings_docking.csv
  （四药：ML 百分位 vs 对接百分位）

【解读写入 SI/主文】
- subset A vs D：ROC-AUC、EF@5%/10%（对接能否分活性/decoy）
- 四药 benchmark：对接是否比 ML 回收更多（预期 docking ≥ ML）
- hybrid_cv：指纹+对接特征是否显著优于指纹 alone（通常小幅，仅 SI）

═══════════════════════════════════════
Phase 3 — ML 有限提升尝试（可选 ablation，非主创新）
═══════════════════════════════════════

3a. 文献活性增补（针对 lesinurad/dotinurad 被 ChEMBL 清洗剔除）

python3 scripts/00c_augment_urat1_literature_benchmarks.py

# 用增补集重训（ablation only）
python3 scripts/02_train_asymmetric_models.py --no-oat-transfer \
  --data-dir data/processed
# 注意：需临时让训练脚本读 urat1_curated_literature_augmented.csv
# 或：cp data/processed/urat1_curated_literature_augmented.csv data/processed/urat1_curated.csv

python3 scripts/07_benchmark_backtest.py

对比重训前后 benchmark 2/4 是否改善；在文稿中写明「数据增补 ablation」，不夸大。

3b. 模型族对比（若未做过）

python3 scripts/08_urat1_model_comparison.py

预期：RF/Chemprop 仍 ~2/4 benchmark，不能替代对接。

【不要做】
- 不要指望 OAT 迁移提升（--oat-transfer 可关掉）
- 不要把 ML 重训写成主贡献
- 不要用 8973 对接分数直接当全库 NLRP3 标签

═══════════════════════════════════════
Phase 4 — 为双靶 Pareto 准备（A′+ 主图，对接整理后立即做）
═══════════════════════════════════════

对同一 8973 manifest 批量 NLRP3 ML（快，不需对接）：

# 若 02_train 尚无 predict 模式，写一小脚本或用 joblib 加载 nlrp3_model.joblib
# 对 data/docking/8973_9DKB_with_manifest.csv 的 canonical_smiles 预测 P_active
# 输出 results/docking/8973_nlrp3_ml_scores.csv

然后合并：
- S_U  = s_u_percentile（来自对接 merged 表）
- S_N  = NLRP3 P_active 百分位
→ Pareto 前沿 + 标注六药（见 docs/PAPER_A_PRIME_PLUS_LOGIC.md）

═══════════════════════════════════════
Phase 5 — 四药/对照补充（若 Phase 1 QC 有漏）
═══════════════════════════════════════

若 verinurad 不在 manifest（已知缺口）：
- 单独 dock @ 9DKB，append 到 raw CSV 后重跑 merge
- 或确认其在 distill_subset_e.csv，合并进分析表

lesinurad redock RMSD：
- 填 results/paper_a_prime/urat1_docking_scores.csv（主文 Table）

═══════════════════════════════════════
完成后汇报
═══════════════════════════════════════

贴出：
1. 8973_docking_qc_summary.json
2. urat1_docking_vs_ml_summary.json
3. urat1_benchmark_rankings_docking.csv
4. （若做了）重训后 benchmark_backtest 四药表

并给一句话结论：对接相对 ML 的 benchmark 回收与 A vs D 富集是否支持「URAT1 结构优先」叙事。
```

---

## 三、输出文件速查

| 文件 | 用途 |
|------|------|
| `data/docking/8973_9DKB_with_manifest.csv` | 8973 + XP 分 + `s_u_percentile` + subset 标签 |
| `data/docking/8973_docking_qc_summary.json` | 覆盖率 QC |
| `results/docking/urat1_benchmark_rankings_docking.csv` | 四药 ML vs 对接排名 |
| `results/docking/urat1_docking_vs_ml_summary.json` | 富集 + ablation 汇总 |

---

## 四、对「提升 URAT1 模型质量」的诚实预期

| 手段 | 能否当主方案 | 预期 |
|------|-------------|------|
| **8973 对接 → $S_U$ 主证据** | ✅ 是 | 主文核心；改善的是**筛选证据**而非 QSAR |
| **A vs D 富集 + 四药 benchmark** | ✅ 是 | 证明对接优于 ML（回顾性） |
| **文献标签增补重训** | △ SI ablation | 可能改善 lesinurad/dotinurad，难达 4/4 |
| **Hybrid（指纹+对接特征）** | △ SI only | 推理仍需先对接，不能替代 $S_U$ |
| **换 Chemprop/RF** | ❌ | 已证 ~2/4，不换 |
| **OAT 迁移** | ❌ | Δρ≈0，不写 |

**论文表述**：*URAT1 evidence was upgraded from fingerprint ML (2/4 benchmark recovery) to 9DKB XP ensemble scoring on the shared 8973 library; ML retained only as auxiliary.*
