# 到货前操作手册：要准备什么文件、做什么、为什么

适用：**无 Schrödinger**；采购 **690 + 2231**；湿实验仅 **JNK1/2/3 IC50**；主叙事 **Option A**（管线 + 选择性预测器失败）。  
配套总览：`WAIT_WINDOW_PLAN_690_2231_zh.md`、`OPENSOURCE_GAPLIST_AND_GNINA_POLICY_zh.md`。

---

# 一、先建立工作目录（一次性）

在项目根 `JNK1_Selectivity_Project/` 下，确保有：

```text
data/structures/pdb/          # 放 3ELJ.pdb 3E7O.pdb 3TTI.pdb（可 gitignore）
results/pose_consensus/       # Vina/Gnina 对接验证
results/md_replicas/          # 无约束 MD（待建）
results/assay/                # 到货后填 IC50
results/assay_analysis/       # C4 预注册输出（已有）
results/chemotype_novelty/    # 已有
results/selectivity_autopsy/  # 已有
results/purchase_risk/        # 已有
paper_rewriting_output/       # 写作与计划
docs/protocols/               # 协议
```

**为什么：** 审稿人要的是“可追溯路径”。文件乱，后面 Results/SI 对不上号。

---

# 二、已经有的（一般不用重做，但要会用）

| 文件/目录 | 你要知道它是什么 | 后面怎么用 |
|-----------|------------------|------------|
| `data/shortlist/md_shortlist_final.csv` | 16 个 MD 短名单 + SMILES + Glide 分 | 查 690/2231 的 SMILES、历史分数 |
| `docs/popular_science/data_tables/27_MD16_选择性排序与报价.csv` | HIT ID、hinge、RMSD、报价 | 写购买理由、核对 HIT100544184=2231、HIT103871685=690 |
| `results/selectivity_autopsy/` | Δsel/Gly87/ML 失败表 | **主文 RQ-C**，与采购无关 |
| `results/chemotype_novelty/` | 690/2231 相对已知 JNK 的 Tc | 主文/SI 新颖性表 |
| `results/pose_consensus/` | Vina 多 seed 结果 | 690 稳；2231-JNK2 不稳 |
| `results/assay/ic50_raw.csv` | 空的 IC50 填表模板 | **到货后唯一湿数据入口** |
| `results/assay_analysis/c4_analysis_lock.json` | 预注册规则 v2 | 到货后禁止改阈值 |
| `results/md_2231_200ns/` | 旧 2231 长 MD（**有配体约束**） | 只能当“历史受限模拟”，不能当无约束证明 |
| `docs/JNK1_PROJECT_REPORT.md` | 漏斗数字与叙事 | Intro/Methods 背景 |

---

# 三、现在要准备的输入文件（开源路径）

## 3.1 受体结构（对接 + MD 共用）

### 要准备的文件
```text
data/structures/pdb/3ELJ.pdb   # JNK1
data/structures/pdb/3E7O.pdb   # JNK2
data/structures/pdb/3TTI.pdb   # JNK3
```

### 怎么准备
```bash
mkdir -p data/structures/pdb
cd data/structures/pdb
for id in 3ELJ 3E7O 3TTI; do
  curl -fsSL -o ${id}.pdb "https://files.rcsb.org/download/${id}.pdb"
done
```
说明见：`data/structures/PDB_DOWNLOAD.md`。

### 为什么
- 这是项目归档的三亚型主 PDB。  
- 无 Schrödinger 时，Vina/Gnina/OpenMM 都要从这些坐标出发。  
- 审稿人会问：对接/MD 用的是哪个结构、哪条链。

### 额外建议准备（强烈）
为每个 PDB 写一份 **1 页笔记**（可放 `docs/protocols/receptor_notes_zh.md`）：
- 用哪条链（通常 A）
- 共晶配体三字母码是什么
- 盒子打算以共晶配体为中心
- 已知质子化/突变注意点（若有）

**为什么：** 避免后面 meeko/`-a` 删残基时自己说不清。

---

## 3.2 配体文件（已购 + 对照）

### 要准备的文件（建议统一建成一张表）
`data/purchase/purchase_panel_smiles.csv`（若没有就新建）：

| compound_id | HIT_ID | role | smiles |
|-------------|--------|------|--------|
| 690 | HIT103871685 | purchased_anchor | `Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1` |
| 2231 | HIT100544184 | purchased_bias_hypothesis | `COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F` |
| E1 | — | positive_control | （文献 SMILES，与 benchmark 表一致） |
| CC-90001 | — | positive_control | （同上） |
| SP600125 | — | optional_tool | （可选，若实验室有货） |

也可直接从已有 CSV 复制 SMILES，但**建议单独一张“湿实验面板表”**，避免和 16 人短名单混淆。

### 为什么
- 到货核对、对接、MD、写文章必须同一套 ID↔SMILES。  
- 2231/690 写错 HIT 是最惨的硬伤。

---

## 3.3 软件与环境（本机/服务器）

### 要准备什么
| 组件 | 用途 | 最低要求 |
|------|------|----------|
| AutoDock Vina 1.2.x | 多 seed 对接（已部分完成） | 可运行二进制 |
| **Gnina** | CNN 精修/重打分 | GPU 更佳，CPU 也可 |
| meeko + gemmi + RDKit | 配体/受体 PDBQT | pip 可装 |
| **OpenMM 或 GROMACS** | 无约束 MD replicas | 有 GPU 最好 |
| OpenFF 或 AmberTools(antechamber) | 配体力场参数 | 二选一 |
| Python + pandas | 分析脚本 | 已有 |

### 要准备的“环境说明文件”
新建：`docs/protocols/OPENSOURCE_ENV.md`，写清：
- 各软件版本号
- 安装方式
- 随机种子约定（如 1,2,3）

**为什么：** Methods/SI 可复现；换机器不丢参数。

---

# 四、按任务拆解：做什么、产出什么文件、为什么

---

## 任务 1 — 共晶 redock（C2c）【建议先做，1–2 天】

### 做什么
对每个 PDB（3ELJ/3E7O/3TTI）：
1. 取出共晶配体  
2. 用 Vina（再加 Gnina）对接回原口袋  
3. 算重原子 RMSD vs 晶体坐标  

### 要准备/产出的文件
```text
results/pose_consensus/redock/
  3ELJ_vina_redock.pdbqt / .sdf
  3E7O_...
  3TTI_...
  redock_rmsd_summary.csv          # PDB, engine, RMSD_A, score
  REDOCK_REPORT.md
```

### 为什么
- 审稿人第一问往往是：“你这套开源协议在已知结构上靠谱吗？”  
- RMSD &lt; 2 Å → 可以说口袋与参数基本可用。  
- 不过关 → 先改盒子/准备，再跑 690/2231，避免后面全废。

---

## 任务 2 — Gnina 验证已购分子（C2b）【2–4 天】

### 做什么
对 **690、2231、E1、CC-90001** × **三亚型**：
1. **固定与 Vina 相同的盒子**（不要重新乱定盒子）  
2. Gnina 对接或对 Vina pose 做 CNN rescore + 局部精修  
3. ≥3 个 seed  
4. 与 Vina top pose 算 RMSD；记录 CNN affinity / CNNscore  

### 要准备的输入
- `data/structures/pdb/*.pdb`  
- 配体 SMILES 表  
- 已有 Vina 盒子参数：`results/pose_consensus/c2_vina_meta.json`（里面有 center/size）

### 产出文件
```text
results/pose_consensus/gnina/
  gnina_scores_by_seed.csv
  gnina_vs_vina_pose_rmsd.csv
  GNINA_CONSENSUS_REPORT.md
```

### 为什么
- 无 Schrödinger 时，需要**第二引擎**交叉验证。  
- 你已有 Vina：690 全过、2231-JNK2 不过；Gnina 用来确认这是不是引擎偶然。  
- **不是**为了重筛库、换化合物。

### 若结果和 Vina/Glide 不一致，你要存什么
在报告里固定三列：
- Vina score / pose  
- Gnina score / pose  
- pairwise pose RMSD  

并写结论句式：
- RMSD≤2 Å → “开源共识 pose”  
- RMSD&gt;2 Å → “该亚型 pose ambiguous，主图不画死，靠 MD/IC50”

**不要**因为 Gnina 更看好别的未购分子就改主线。

---

## 任务 3 — 无约束 MD replicas（C3）【本月主战场，1–3 周】

### 做什么
体系：
```text
{690, 2231} × {JNK1/3ELJ, JNK2/3E7O, JNK3/3TTI} × {seed1, seed2[, seed3]}
= 12～18 条轨迹
```
每条：
- 最小化 → 升温 → NPT 平衡 → **生产 20–50 ns**  
- **生产阶段配体不加位置约束**（与旧 2231 200 ns 带约束不同）  
- 保存轨迹、能量日志、关键参数文件  

### 要准备的输入文件（每条体系一套）
```text
results/md_replicas/systems/690_JNK1_seed1/
  protein.pdb / .gro
  ligand.sdf / .mol2
  ligand.frcmod / openff params   # 力场
  solvated.pdb / system.xml      # 溶剂化体系
  mdin / openmm_script.py        # 运行脚本
  README.md                      # 种子、温度、时长
```

### 产出分析文件
```text
results/md_replicas/
  c3_replica_summary.csv
    # compound, isoform, seed, ligand_RMSD_mean/sd, hinge_HB_occ, pass_qc
  c3_pass_flip_table.csv         # 多种子间 pass/fail 是否翻转
  C3_MD_REPLICA_REPORT.md
  figures/rmsd_*.png  hinge_*.png
```

### 分析指标（与项目旧 QC 对齐）
- 配体 RMSD（蛋白对齐后）  
- 铰链氢键占用率  
- （可选）口袋残基接触频率  

建议沿用归档阈值做“是否通过”的对照（具体数字以 `docs` 里 MD QC 为准，写入报告）：例如 RMSD≤3 Å、hinge≥30%——**先核对旧文档再锁定**，锁了就不要改。

### 为什么
1. 旧筛选大量依赖**单次 MD**；2231 还有**带约束**长轨迹 → 审稿人会说不可靠。  
2. 2231 是 grade C + Vina 显示 JNK2 不稳 → 必须用无约束多 seed 回答“假说还在不在”。  
3. 到货后若 IC50 与 MD 不一致，你需要“多 seed 方差”来写 Discussion，而不是事后编。  
4. 这是无 Schrödinger 条件下，**对得起 JCIM 类期刊**的最关键湿实验前计算。

### 和旧 `results/md_2231_200ns/` 的关系
- 旧文件：**保留**，在报告里写 “ligand-restrained historical run”。  
- 新 C3：**单独目录**，不要覆盖旧结果。  
- 正文引用以 **C3 无约束** 为准。

---

## 任务 4 —（可选但推荐）2231 结构假说表（C8′）【3–5 天，可与 MD 并行】

### 做什么
整理 2231 在三亚型中与以下位点的接触：
- JNK1 特有体积相关位点（项目中 Ile106 叙事）  
- JNK2/3 对应 Leu  
- 铰链氢键  

来源：Gnina/Vina pose + MD 接触频率。

### 产出
```text
results/structural_hypothesis_2231/
  contact_frequency_table.csv
  HYPOTHESIS_2231.md   # 可证伪假说：若偏好存在，预期接触模式是…
```

### 为什么
- RQ-B 需要对“为什么猜 2231 偏 JNK1”有一句结构故事。  
- IC50 若打脸，假说被证伪也是结果（符合 Option A）。

---

## 任务 5 — 回顾富集与脱靶（C9/C10）【各 1–2 天】

### C9 产出
```text
results/enrichment_prior/
  random_vs_shortlist_scores.csv
  ENRICHMENT_NOTE.md
```
**为什么：** n=2 太小，需要告诉审稿人“我们不宣称 hit rate”。

### C10 产出
```text
results/offtarget_prediction/
  690_sea_or_stp.json/csv
  2231_sea_or_stp.json/csv
  OFFTARGET_CAVEAT.md   # 写清：不能替代 kinome
```
**为什么：** 只有三亚型酶活时，审稿人常问脱靶；用预测顶一下 Discussion。

---

## 任务 6 — 写作文件（贯穿整月）

### 要准备/更新的文稿文件
| 文件 | 你要干什么 | 为什么 |
|------|------------|--------|
| `paper_rewriting_output/draft_intro_methods_rqc_en.md` | 扩写 Intro + Methods（开源路径）+ RQ-C | 到货前就能完成一半正文 |
| `paper_rewriting_output/section_blueprints.md` | 按实际采购改 690/2231 | 防写偏成 2157 |
| `paper_rewriting_output/confirmed_*.md` | 已更新购买集，保持一致 | 贡献不漂移 |
| `docs/protocols/SOFTWARE_LICENSE_NOTE.md` | 明确不写 Schrödinger 可复现细节 | 合规 |
| 图注草稿 `paper_rewriting_output/figure_captions_draft_zh.md`（建议新建） | 漏斗、C5、C2/C3、690 vs 2231 | 到货后只插 IC50 图 |

### Methods 必须写清的两层（建议直接贴进草稿）
1. **Selection：** 既往 Glide 商业库漏斗 → 归档短名单 → 采购 690/2231  
2. **Confirmation：** Vina/Gnina 多 seed + 无约束 MD replicas（本文可复现部分）

---

## 任务 7 — 到货当天/当周（湿实验接口）

### 要准备的文件（现在就空着等填）
```text
results/assay/ic50_raw.csv
```
列：`compound_id, isoform, ic50_uM, ic50_nM, pct_inh_10uM, n_replicates, assay_date, notes`

填完后运行：
```bash
python3 scripts/c4_preregistered_ic50_analysis.py
```

### 产出
```text
results/assay_analysis/c4_ic50_si_table.csv
results/assay_analysis/c4_endpoint_verdict.json
```

### 为什么
- 终点在到货前已锁（≤10 µM；SI≥3 才称 preference）。  
- 防止看到 2231 好看/难看就改标准（HARKing）。

### 实验上你还要准备的东西（非计算，但常被忘）
- DMSO 储液浓度、稀释梯度方案  
- 同板：690、2231、E1、CC-90001  
- 至少双复孔  
- 板布局记录（建议另存 `results/assay/plate_map.xlsx`）

---

# 五、推荐时间线（文件视角）

| 周 | 你桌上应出现的新文件 | 核心目的 |
|----|----------------------|----------|
| 第 1 周 | `redock_rmsd_summary.csv`；`gnina/` 初步结果；更新 Methods 草稿 | 证明开源协议可用 + 交叉引擎 |
| 第 2–3 周 | `results/md_replicas/**` 轨迹与 `c3_replica_summary.csv` | 回答 2231 稳不稳 |
| 第 4 周 | 图注、Intro/RQ-C 定稿、C4 冒烟测试 | 到货即测即算 |
| 到货周 | 填 `ic50_raw.csv` → C4 表 | 闭合 RQ-A/B |

---

# 六、最小“可投稿计算包”检查表（到货前勾完）

- [ ] 三 PDB 在本地且笔记写清链/盒子  
- [ ] 购买面板 SMILES/HIT 核对无误  
- [ ] Redock RMSD 表  
- [ ] Vina（有）+ Gnina（补）对 690/2231 的共识报告  
- [ ] 无约束 MD multi-seed 汇总表（即使用 20 ns×2 seed 也比没有强）  
- [ ] C5 尸检表可进主文  
- [ ] C4 lock 仍是 690+2231，未改阈值  
- [ ] Methods 已写成“历史 Glide + 开源验证”  
- [ ] 旧 2231 约束 MD 已降级标注  

---

# 七、用一句话记住

**准备三套东西：**  
1）**结构与配体文件**（PDB + SMILES 面板）；  
2）**开源验证产物**（redock、Gnina、无约束 MD）；  
3）**预注册与写作文件**（C4 空表 + Methods/RQ-C 草稿）。  

**为什么：** 没有 Schrödinger 时，审稿人买账的是“你能复现的验证”，不是“你曾经用 Glide 筛过”。Gnina 不一样也不换货，只把差异写进验证报告。
