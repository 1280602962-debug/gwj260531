# DualFourClass-Bench（Dual_Target_Docking）

JCIM **evaluation / benchmark** 课题：双靶对接的四状态实验比较评价。  
**不是**通用决策臂 / 新打分函数论文。  
贡献定位：八靶对上的评价设定审计 + 可复现评价协议 **DualFourClass-Bench**。它不是通用、代表性或 LIT-PCBA 规模的 benchmark suite。
本 DualFourClass 树是一项独立科学课题，不导入、不依赖仓库根目录下的其他项目。

## 一句话主张

同一套对接分数在不同实验状态比较下回答不同的评价问题。主指标是两条口袋匹配方向 AUROC 及其较弱方向摘要。在具有双端实验测量的化合物中，双端低活性对照可能掩盖对单靶活性化合物的识别不足。配体化学与口袋对应对照是解释这些差异所必需的。

## 快速入口（投稿用）

| 文档 | 用途 |
|------|------|
| **[`docs/MANUSCRIPT_JCIM_EN.md`](docs/MANUSCRIPT_JCIM_EN.md)** | ★ **组装后的英文主稿** |
| **[`docs/MANUSCRIPT_JCIM_ZH.md`](docs/MANUSCRIPT_JCIM_ZH.md)** | ★ **组装后的中文工作稿**（`python3 docs/assemble_manuscript_zh.py`） |
| **[`submission_pack/`](submission_pack/)** | ★ **投稿切片：文稿、SI、锁定表、图、核对报告** |
| **[`docs/STATISTICAL_LOCK_V1.md`](docs/STATISTICAL_LOCK_V1.md)** | ★ **八靶对 Table 2/3 估计量与 canonical CSV** |
| **[`docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`](docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md)** | 英文 SI（Tables S1–S13） |
| [`data/jcim_bench_v0/CLAIM_CEILING.md`](data/jcim_bench_v0/CLAIM_CEILING.md) | 可写 / 禁止写的 claim |
| [`docs/REPO_CLEANUP_NOTES.md`](docs/REPO_CLEANUP_NOTES.md) | 仓库删了什么、为何保留 |

章节工作稿在 `docs/*_DRAFT_ZH_JCIM_V1.md` 与 `docs/*_SECTION_JCIM_EN_V1.md`。组装后的主稿才是投稿文本。

## Canonical 数据

| 用途 | 文件 |
|------|------|
| Table 2（EGFR/HER2、AChE/BChE、PIK3CA/mTOR） | `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` |
| Table 2 / Table 3（其余五对） | `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv` |
| Table 3（前三对 Dual-versus-neither） | `data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv` |
| 主结果索引 | `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` |
| BindingDB 八对剩余（Table S11b） | `data/jcim_novelty_v0/tables/external_slice_summary_v1.csv` |
| 评价合约 | `data/jcim_novelty_v0/tables/DUALFOURCLASS_EVALUATION_CONTRACT_v1.json` |

历史三对 BindingDB 快照 `external_slice_summary_202608_contract_v1.csv` 不是 Table S11。PIK3CA/PIK3CB 已从主分析撤出，面板目录仅作撤出记录，不是 Table 2。

## 保留的数据目录

```text
data/
├── protocols/
├── public_pair_selection/
├── jcim_j0j1_v0/                   # 供给审计
├── jcim_bench_v0/
├── jcim_strengthen_t0t1_v0/        # 前三对 Table 2 CI
├── jcim_chembl_universe_v0/        # 五对生产分数与 Table 2 行
├── jcim_novelty_v0/                # 主结果索引、BindingDB、簇 bootstrap
├── jcim_independent_dock_v0/       # 独立 GNINA
├── jcim_multiseed_v0/
├── jcim_structure_robust_v0/
├── egfr_her2_panel120_v0/
├── ache_bche_panel_v0/
├── pik3ca_mtor_panel48_rdkit_v0/
└── pik3ca_pik3cb_panel_v0/         # withdrawn; not Table 2
```

## 本地运行

```bash
cd Dual_Target_Docking
python3 scripts/check_local_env.py
python3 scripts/primary/bootstrap_primary.py
python3 data/jcim_novelty_v0/scripts/validate_revision_v1.py
```

- 最小依赖：[`requirements-analysis.txt`](requirements-analysis.txt)
- 说明文档：[`docs/LOCAL_RUN.md`](docs/LOCAL_RUN.md)
- 该入口只核**零对接分析**，不重跑 Vina / RTM / GNINA

重建中英文主稿：

```bash
python3 docs/assemble_manuscript_en.py
python3 docs/assemble_manuscript_zh.py
```

## 状态

- 主分析：八个靶对，统一 θ = 6.0；PIK3CA/PIK3CB 不进入 Table 2
- BindingDB 八对独立来源剩余：没有任何靶对被包装为外部评价集
- 英文主稿 `docs/MANUSCRIPT_JCIM_EN.md`；中文工作稿 `docs/MANUSCRIPT_JCIM_ZH.md`
- Zenodo DOI：（发布后填这里；不要从当前分支签发）
