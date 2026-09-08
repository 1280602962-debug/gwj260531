# NLRP3 / URAT1 双靶项目 —— 对接数据导出包

导出日期：2026-08-20 ｜ 总计 **27,694 个文件 / 583 MB** ｜ 全部为真实文件拷贝，**无符号链接**

本包汇总了本项目从 2026-07-24 到 2026-08-17 期间所有**对接相关**的产出数据，可直接整体拷贝带走。

---

## 一、做过哪些对接

| 阶段 | 时间 | 靶点 | 配体规模 | 目的 |
|------|------|------|----------|------|
| **阶段一** 协议基准验证 | 07-24 ~ 07-29 | 9DKB（URAT1） | 9,839 | 用 active/decoy 基准比较 6 套打分协议，选出生产协议 |
| **阶段二** 双靶重定位生产对接 | 07-30 ~ 08-02 | 9DKB（URAT1）+ 7ALV（NLRP3） | 1,583 临床药物池 | 用锁定协议筛双靶候选，Pareto 合并 + ADMET + 候选提名 |
| **阶段三** 姿态 QC 与 MD 选型 | 08-17 | 双靶 | 8 个候选 | 对接姿态质控，挑出 3 个化合物 × 2 靶点 = 6 个 MD 体系 |

### 阶段一结论：协议排名（True decoy 基准，按 EF@1% 排序）

| ID | 协议 | 覆盖 | EF@1% | EF@5% | AUC |
|----|------|------|-------|-------|-----|
| **P5** | **RTMScore（打在 gnina 姿态上）** | 9827/9839 | **2.804** | 2.440 | 0.590 |
| **P2** | **gnina CNNaffinity** | 9839/9839 | **2.537** | 1.875 | 0.580 |
| P0 | gnina CNNscore | 9839/9839 | 1.903 | 2.386 | 0.647 |
| P4 | RTMScore（打在 Vina 姿态上） | 8983/9839 | 0.647 | 1.026 | 0.625 |
| P1 | Vina affinity | 9750/9839 | 0.431 | 0.513 | 0.531 |
| P3 | gnina minimizedAffinity | 9839/9839 | 0.423 | 0.511 | 0.503 |

EF=1 / AUC=0.5 为随机基线。P4 覆盖不全（缺 ~765 分子），与 P5 的比较不宜过度解读。
**生产协议最终锁定为 P2**（gnina，`cnn_scoring=rescore`，`exhaustiveness=32`），因为它覆盖完整、无需额外 RTMScore 步骤且 EF 接近最优。

### 阶段二结论

1,583 个临床池配体对双靶各完成 1,582 个（`REP_05842` 配体准备失败，PDBQT 为空文件，两靶的对应 pose 也是空的——这是**原始数据既有情况**，非拷贝错误）。
双靶评分合并后 1,580 个分子进入 Pareto 分析，Pareto 前沿 4 个、shortlist 4 个；经 ADMET/结构过滤后提名候选 **50 个**（默认阈值）或 **85 个**（tau85 宽松阈值），多样性精选各 7 / 8 个。
已知对照药 lesinurad、verinurad、colchicine 均**未**落在 Pareto 前沿。

### 阶段三结论

8 个候选做了双靶姿态 QC（口袋内占位、关键接触、氢键、冲突数）。
最终建 MD 体系的 3 个化合物：**GSK-3008348**、**VECABRUTINIB**、**ZELENIRSTAT**（各含 NLRP3 与 URAT1 两个体系）。
双靶姿态综合分最高的其实是 DEUCRICTIBANT（31.20），但未进入 MD 批次。

---

## 二、目录说明

```
docking_export_20260820/
├── README.md            本文件
├── MANIFEST.txt         完整文件清单（含每个文件大小）与目录统计
├── 00_overview/         四份原始报告：运行汇总、预检审计、项目说明、服务器说明
├── 01_phase1_benchmark_URAT1_9DKB/
├── 02_phase2_dualtarget_9DKB_7ALV/
└── 03_pose_qc_md_selection/
```

### 01 阶段一（541 MB）

| 路径 | 内容 |
|------|------|
| `scores/mol_protocol_scores.csv` | **主表**，9,839 行 × 6 协议得分 |
| `scores/rtm_{vina,gnina}_pose_scores.csv` | RTMScore 逐 pose 原始分 |
| `metrics/protocol_metrics.csv` | 上面那张协议排名表的数据源 |
| `benchmarks/` | true_decoy / random_decoy 基准集 + 去重对接池 |
| `meta/mol_index_map.csv` | mol_id ↔ SMILES ↔ role（active/true_decoy/random_decoy） |
| `receptor/` | 9DKB 受体 PDB/PDBQT、RTM 口袋、lesinurad 晶体参考配体 |
| `config/` | 对接盒子（center 99.980/102.958/105.657，size 20³）与运行配置 |
| `poses/vina_pdbqt/` | **9,750** 个 Vina 姿态（`mol_XXXXX_out.pdbqt`） |
| `poses/gnina_sdf/` | **9,839** 个 gnina 姿态（`mol_XXXXX_out.sdf`，含 CNNscore/CNNaffinity 属性） |
| `logs/` | 汇总日志、超时与跳过清单 |
| `scripts/` | 全套对接脚本 + SLURM 作业模板 |

`mol_protocol_scores.csv` 列方向：`P1_vina_affinity` / `P3_gnina_affinity` **越低越好**；`P0_CNNscore` / `P2_CNNaffinity` / `P4_RTM_vina` / `P5_RTM_gnina` **越高越好**。

### 02 阶段二（42 MB）

| 路径 | 内容 |
|------|------|
| `docking_9dkb_URAT1/` | URAT1 结果：打分 CSV（1,583 行）、summary、done_ids、**3,166 个姿态** |
| `docking_7alv_NLRP3/` | NLRP3 结果：打分 CSV（1,583 行）、summary、**3,166 个姿态** |
| `receptors/` | 9DKB 与 7ALV 已准备受体 PDBQT + 质控 JSON |
| `ligands/ligands_p05_pdbqt/` | 1,583 个临床池配体 PDBQT + manifest |
| `pareto/` | 双靶合并主表（1,580 行 × 34 列）、shortlist、summary |
| `candidates/` `candidates_tau85/` | 提名候选（50 / 85 个）+ 多样性精选 + 提名依据 JSON |
| `cheminformatics/` | ADMET 与规则过滤结果 |
| `screening_pool/` | 上游 NLRP3 ML 打分与 p05 对接池 |
| `config/` `scripts/` `logs/` | P2 生产配置、运行脚本、服务器日志 |

### 03 阶段三（184 KB）

`pose_qc_dual.csv`（7 个化合物双靶综合分）、`pose_qc_table.csv`（16 行逐靶明细）、`residue_map.json`（关键残基定义）、`poses/{NLRP3,URAT1}/REP_*/ligand.pdb`（QC 用配体姿态）、`md_systems_built.txt`。

---

## 三、未包含的内容

以下数据体量大且属于中间过程，按需再单独导：

| 内容 | 位置 | 大小 |
|------|------|------|
| RTMScore 中间文件（分块 SDF、per-mol SDF） | `server_dock_maestro_prep/work/rtmscore_{vina,gnina}/` | 850 MB |
| 逐分子对接日志（19,679 个） | `server_dock_maestro_prep/work/logs/` | 78 MB |
| 阶段一配体库 SDF/PDBQT（对接输入，非产出） | `server_dock_maestro_prep/work/ligands_{sdf,pdbqt}/` | 102 MB |
| MD 轨迹与 CHARMM-GUI 建模文件 | `md_dual_target/systems/`、`charmm-gui-8417953108/` | 1.2 GB |

注：RTMScore 的**最终打分**已包含在 `01_.../scores/rtm_*_pose_scores.csv` 中，未包含的只是中间产物。
