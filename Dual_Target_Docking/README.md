# DualFourClass-Bench（Dual_Target_Docking）

JCIM **evaluation / benchmark** 课题：双靶对接四类硬负评测。  
**不是**通用决策臂 / 新打分函数论文。  
与仓库内 `JNK1_Selectivity_Project/` **无关**。

## 一句话主张

双靶对接应建成 dual / A_only / B_only / neither 任务；主指标用**口袋匹配方向 AUROC**；必须报平凡基线与混淆对照。公开严格硬负极稀缺；K=4 冻结评测集上对接增量高度对靶依赖。

## 快速入口（投稿用）

| 文档 | 用途 |
|------|------|
| **[`docs/JCIM_P0_COMPLETION_GUIDE.md`](docs/JCIM_P0_COMPLETION_GUIDE.md)** | ★ **下一步：Zenodo + 英文稿怎么做完** |
| [`data/jcim_bench_v0/CLAIM_CEILING.md`](data/jcim_bench_v0/CLAIM_CEILING.md) | 可写 / 禁止写的 claim |
| [`data/jcim_strengthen_t0t1_v0/analysis/PRIMARY_METRIC_V2.md`](data/jcim_strengthen_t0t1_v0/analysis/PRIMARY_METRIC_V2.md) | 口袋匹配主表 |
| [`data/jcim_strengthen_t0t1_v0/analysis/A_GROUP_VERDICT.md`](data/jcim_strengthen_t0t1_v0/analysis/A_GROUP_VERDICT.md) | 混淆对照裁决 |
| [`data/jcim_strengthen_t0t1_v0/analysis/B_GROUP_VERDICT.md`](data/jcim_strengthen_t0t1_v0/analysis/B_GROUP_VERDICT.md) | E8 / enrichment / PM110 |
| [`docs/JCIM_STRENGTHENING_PLAN_V1.md`](docs/JCIM_STRENGTHENING_PLAN_V1.md) | 加厚规划（历史+路线） |
| [`docs/REPO_CLEANUP_NOTES.md`](docs/REPO_CLEANUP_NOTES.md) | 本仓库删了什么、为何保留 |

## 保留的数据目录

```text
data/
├── protocols/                      # PAIR_ROLES_APPROVED_JCIM.yaml 等
├── public_pair_selection/          # ChEMBL 缓存 + 供给审计原料
├── jcim_j0j1_v0/                   # 49 对供给审计
├── jcim_bench_v0/                  # 基准汇总 + CI 图 + CLAIM_CEILING
├── jcim_strengthen_t0t1_v0/        # 口袋匹配/混淆/E8/enrichment/PM110 分析
├── stage_m_v0/                     # Track B=Weak 测量审计（正文可引用）
├── egfr_her2_panel120_v0/          # EGFR 案例分数（θ=6）
├── egfr_her2_panel40_reprep_rdkit_v0/  # 统一 RDKit prep
├── egfr_her2_panel40_v0/           # 协议/cognate/表（姿态已从 git 去掉）
├── pik3ca_mtor_panel48_rdkit_v0/   # PM 主面板 RDKit
├── pik3ca_mtor_panel48_v0/         # LigPrep 对照 + cognate E8/E16（姿态已去掉）
├── pik3ca_mtor_panel110_rdkit_v0/  # PM 扩面
├── ache_bche_panel_v0/
└── pik3ca_pik3cb_panel_v0/
```

大姿态文件在本地  
`/mnt/d/CADD paper exercise/dual target docking/results/`  
上传 Zenodo 时按 `jcim_strengthen_t0t1_v0/POSE_UPLOAD_CHECKLIST.md` 打包。

## 复现主分析（零对接）

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
python3 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

## 状态

- 对接 + 混淆对照 + E/enrichment/PM110：**完成**  
- **当前缺口 = P0**：Zenodo DOI + 英文稿（见 P0 指南）  
- Zenodo DOI：（发布后填这里）
