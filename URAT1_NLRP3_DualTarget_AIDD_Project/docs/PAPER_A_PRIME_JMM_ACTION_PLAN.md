# 论文 A′ 执行清单 — 双靶平行 · 首投 J. Molecular Modeling

> **叙事**：痛风代谢–炎症双节点下，URAT1 四药 @ 9DKB 与 NLRP3 两药 @ 8ETR/7ALV 的 **paired computational pharmacology**  
> **期刊**：*Journal of Molecular Modeling*（Springer，**Subscription / 非 OA**）  
> **与 8000 条 XP 对接的关系**：见下文 §3 —— **主文不依赖全库结果**

---

## 1. 现有数据总览（截至仓库状态）

### 1.1 论文 A′ **主文必需**（6 个化合物）

| 化合物 | 靶点 | PDB / Grid | 对接 | Redock | MD | MM-GBSA |
|--------|------|------------|------|--------|-----|---------|
| lesinurad | URAT1 | **9DKB** | ✅ 你已完成 SP→XP | ⏳ 需算 RMSD | ❌ | ❌ |
| benzbromarone | URAT1 | 9DKB | ✅ | — | ❌ | ❌ |
| verinurad | URAT1 | 9DKB | ✅ | — | ❌ | ❌ |
| dotinurad | URAT1 | 9DKB | ✅ | — | ❌ | ❌ |
| GDC-2394 | NLRP3 | **8ETR** | ❌ | ⏳ 需做 | ❌ | ❌ |
| MCC950 | NLRP3 | **7ALV** | ❌ | —（共晶为类似物） | ❌ | ❌ |

SMILES 与 ID：`data/distill/teacher_gate_qc_panel_b_direction.csv`（URAT1 四药）+ `data/benchmarks/literature_benchmarks_summary.csv`（NLRP3 两药）

### 1.2 结构 / 配置（已就绪，Methods 直接引用）

| 资源 | 路径 | 用途 |
|------|------|------|
| 三态 PDB 映射（A′ 仅用 9DKB） | `docs/URAT1_THREE_STATE_DOCKING.md` | Introduction 结构选型 |
| NLRP3 7ALV / 8ETR | `config/docking_ensemble.yaml` → `nlrp3_ensemble` | NLRP3 grid 参数 |
| 关键残基 URAT1 | Phe241, Arg477, Tyr527 | 相互作用表 |
| Gate 四药面板 | `data/distill/teacher_gate_qc_panel_b_direction.csv` | 主文化合物列表 |

### 1.3 背景数据（**主文不写或仅 Discussion 一句**）

| 数据 | 规模 | 文件 | A′ 用法 |
|------|------|------|---------|
| URAT1 ChEMBL  curated | 822 | `data/processed/urat1_curated.csv` | Discussion：0 重叠、标签噪声 |
| NLRP3 ChEMBL | 513 化合物 / 39 assays | `data/processed/nlrp3_records.csv` | Discussion：assay 异质性 7.2% |
| URAT1–NLRP3 SMILES 重叠 | **0** | `data/processed/data_summary.json` | Discussion：双靶 ML 困难 |
| URAT1 ML（XGBoost+conformal） | OOF Spearman≈0.73；benchmark **2/4** | `results/benchmark_backtest/` | **不写主文** |
| NLRP3 assay-conditioned ML | AUROC≈0.89；benchmark 2/2（均在训练集） | 同上 | **仅 SI 可选 1 段** |
| OAT 迁移 | Δρ≈0.004 | `results/training/oat_transfer_ablation.json` | **不写** |
| 9B1K/9B1L 刚性失败 | 四药零 pose | 你的 Maestro 记录 | Discussion **1 段局限**；详留给论文 B |

### 1.4 正在运行：8000+ @ 9DKB XP

| 项目 | 说明 |
|------|------|
| 对应文件 | 大概率 `distill_subset_d.csv`（**8000** 条无活性 decoy）或子集 A+D |
| 与 A′ 关系 | **不是主文必需**；属于原 MASFL / 论文 B decoy 管线 |
| 建议 | **后台跑完即可**；不要等它再开 MD / NLRP3 / 写稿 |

---

## 2. 8000 条对接：该怎么用、不该怎么用

### ✅ 跑完后值得做的（1–2 小时，可选 SI）

1. 导出统一 CSV：`compound_id, smiles, glide_score_xp, pose_status`  
2. 查四药在 8000 池中的 **百分位排名**（若四药不在 D 集，单独 dock 一次即可）  
3. 可选：从子集 A（822 活性）抽已 dock 的 top/bottom 与四药对比 —— **仅 SI 图，不进主文结论**

### ❌ 不要做的

- 不要把「筛了 8000」写成主文贡献或虚拟筛选发现  
- 不要等 8000 完成才开始 MD / NLRP3 / 写作  
- 不要把 8000 结果与 NLRP3 做任何「双靶融合排序」

---

## 3. 接下来干什么（按优先级）

### Phase 0 — 本周同步进行（不依赖 8000）

#### P0-1 冻结 Methods 参数（半天）

建 `results/paper_a_prime/methods_snapshot.yaml`（或 Maestro 项目备注），记录：

- Schrödinger 版本号  
- 9DKB：Protein Prep 选项、pH、grid 中心、box 22 Å、SP→XP 参数  
- 将为 NLRP3 使用的 8ETR / 7ALV 同样模板  

#### P0-2 URAT1 主文数据整理（1 天）

- [ ] **lesinurad redock**：晶体 vs XP pose → **RMSD**（Gate ≤ 2.0 Å）  
- [ ] 四药导出：**最佳 XP pose**（.maegz）、GlideScore、2D interaction diagram  
- [ ] 填表 `results/paper_a_prime/urat1_docking_scores.csv`  

| 列建议 |
|--------|
| compound_name, smiles, glide_score_sp, glide_score_xp, key_contacts, redock_rmsd_A |

- [ ] 目检四药 pose：无明显 steric clash  

#### P0-3 启动 NLRP3 对接（1–2 天，与 P0-2 并行）

| 步骤 | 操作 |
|------|------|
| 1 | 下载 **8ETR**（GDC-2394 共晶）、**7ALV**（MCC950 类类似物 NP3-146） |
| 2 | Protein Prep → Grid（Walker B 变构口袋，box ~20 Å，见 yaml） |
| 3 | **GDC-2394** @ 8ETR：SP→XP → **redock RMSD** |
| 4 | **MCC950** @ 7ALV：SP→XP（Methods 写明共晶配体为类似物，非 MCC950 本身） |
| 5 | 导出 pose + 分数 → `results/paper_a_prime/nlrp3_docking_scores.csv` |

#### P0-4 交叉结构对照（可选 SI，半天）

- URAT1：四药与 **9DKA / 9JDY / 9JE1** 共晶模式 **视觉比较**（非 redock）  
- NLRP3：与 8ETR 共晶 GDC-2394 比较  

---

### Phase 1 — MD + MM-GBSA（主文核心，约 1–2 周）

**在对接 pose 目检通过后立即开始**，与 8000 无关。

| 体系 | 数量 | 时长建议 | 软件 |
|------|------|----------|------|
| URAT1–抑制剂 @ 9DKB | 4 | **100 ns** | Desmond 或 GROMACS（全文统一） |
| NLRP3–抑制剂 @ 8ETR/7ALV | 2 | **50–100 ns** | 同上 |

#### 每个复合物流程

1. 复合物构建（蛋白+配体，质子化与对接一致）  
2. 平衡：NVT → NPT  
3. 生产段：记录 RMSD、RMSF、关键距离（如 Arg477–配体）  
4. **MM-GBSA**：末 20 ns 或每 10 ns 一帧  

输出目录：

```
results/paper_a_prime/
├── md/
│   ├── urat1_lesinurad_9DKB/
│   ├── ...
│   └── nlrp3_gdc2394_8ETR/
├── mmgbsa/
│   └── summary.csv
└── figures/   # 后期作图
```

#### MM-GBSA 汇总表（主文 Table）

| compound | target | pdb | dg_bind_mean | dg_bind_std | md_rmsd_plateau |

---

### Phase 2 — 图表与写作（MD 跑起来后即可穿插）

#### 主文图表（8 个）

| Fig | 内容 |
|-----|------|
| 1 | 痛风双节点示意图 + 9DKB / 8ETR 口袋 |
| 2 | URAT1 四药 9DKB 口袋 overlay |
| 3 | lesinurad redock + GDC-2394 redock |
| 4 | URAT1 四药 MD RMSD |
| 5 | NLRP3 两药 MD RMSD |
| 6 | MM-GBSA 六化合物对比条形图 |
| Table 1 | 六药对接分数 + 关键相互作用 |
| Table 2 | MD / MM-GBSA 汇总 + 实验 IC50（趋势对照，**不强求线性相关**） |

#### SI（可选）

- 9DKB vs 9B1H 叠合  
- QikProp / SwissADME 六药  
- NLRP3 ML AUROC 一句 + ROC 图（注明训练集内）  
- 四药在 decoy 池百分位（**仅当 8000 跑完且你愿做**）

#### 写作顺序

1. **Methods**（对接+MD 参数已冻结）  
2. **Results 3.1–3.2** URAT1  
3. **Results 3.3–3.4** NLRP3  
4. **Introduction**（疾病双靶动机）  
5. **Discussion**（转运体 vs NACHT；0 重叠；9B1K/B1L 一句局限）  
6. **Abstract** 最后写  

大纲：将 `MANUSCRIPT_OUTLINE_FAST_9DKB.md` 扩展 NLRP3 两节（或见 `DUAL_TARGET_AND_FAST_JOURNALS.md` §1.3）

---

### Phase 3 — 投稿 J. Molecular Modeling

- [ ] 投稿系统选 **Subscription**（非 OA，无 APC）  
- [ ] Cover letter：强调 **paired computational pharmacology**、9DKB cryo-EM、**非** hit discovery  
- [ ] 上传 SI：对接参数、轨迹分析脚本、可选 decoy 排名  
- [ ] 建议审稿人：URAT1 转运体结构药理学 + NLRP3 变构抑制剂各 1 人  

---

## 4. 四周甘特（与 8000 解耦）

```
周 1  │ P0-2 URAT1 表图 + P0-3 NLRP3 对接 + Methods 冻结
      │ （8000 XP 后台继续）
周 2  │ 提交 MD 6 体系；开始 Methods/Results URAT1 文字
周 3  │ MM-GBSA；作图；Results NLRP3 + Discussion 初稿
周 4  │ 全文润色；8000 如有结果仅补 SI；投稿 JMM
```

**关键路径 = 6 药 MD**，不是 8000 XP。

---

## 5. 每日/每周检查清单

### 本周必须完成

- [ ] `urat1_docking_scores.csv` 四行 + lesinurad RMSD  
- [ ] NLRP3 两药 pose + GDC-2394 redock RMSD  
- [ ] Methods 参数快照文档  
- [ ] 6 个 MD 任务已提交队列  

### 不应花时间的事（A′ 阶段）

- [ ] ~~8973 / 三态 B1K B1L 对接~~  
- [ ] ~~Teacher M-CPDL~~  
- [ ] ~~双靶融合排序~~  
- [ ] ~~等 8000 完成再写稿~~  

---

## 6. 结果文件模板

### `results/paper_a_prime/urat1_docking_scores.csv`

```csv
compound_name,smiles,pdb,glide_score_xp,redock_rmsd_A,key_residues,pose_notes
lesinurad,O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12,9DKB,,,Phe241;Arg477,
...
```

### `results/paper_a_prime/nlrp3_docking_scores.csv`

```csv
compound_name,smiles,pdb,glide_score_xp,redock_rmsd_A,cocrystal_note,key_residues
GDC-2394,...,8ETR,,,native ligand,
MCC950,...,7ALV,,,NP3-146 analog in 7ALV,
```

---

## 7. 相关文档

| 文件 | 用途 |
|------|------|
| `PAPER_A_PRIME_JMM_ACTION_PLAN.md` | 本文件 |
| `DUAL_TARGET_AND_FAST_JOURNALS.md` | 双靶措辞与选刊 |
| `data/benchmarks/literature_benchmarks_summary.csv` | 六药 IC50 / PDB |
| `config/docking_ensemble.yaml` | Grid 与结合位点 |
