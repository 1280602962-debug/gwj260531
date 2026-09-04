# DualFourClass-Bench（Dual_Target_Docking）

JCIM **evaluation / benchmark** 课题：双靶对接四类硬负评测。  
**不是**通用决策臂 / 新打分函数论文；**也不是**名为 D-DRAF 一类的“新型 Framework”。  
贡献定位：四靶对 *proof-of-concept formulation audit* + 可复现评价协议 **DualFourClass-Bench**。它不是通用、代表性或 LIT-PCBA 规模的 benchmark suite。
本 DualFourClass 树是一项独立科学课题，不导入、不依赖仓库根目录下的其他项目。

## 一句话主张

双靶对接应建成 dual / A_only / B_only / neither 四状态任务，主指标是两条口袋匹配方向 AUROC 及其最差方向判别摘要（不是四分类器或新 scoring function）。在当前支架分组评价中，表观判别高度依赖靶对，把 docking 加到 ECFP4 后只产生很小的 CV AUROC 增量改善。建立的是评价体系，不是新算法或 comprehensive suite。

## 快速入口（投稿用）

| 文档 | 用途 |
|------|------|
| **[`docs/MANUSCRIPT_JCIM_EN.md`](docs/MANUSCRIPT_JCIM_EN.md)** | ★ **组装后的英文主稿**（投稿以这份为准） |
| **[`docs/JCIM_PROJECT_AUDIT_AND_ACTION_PLAN_2026-08-27.md`](docs/JCIM_PROJECT_AUDIT_AND_ACTION_PLAN_2026-08-27.md)** | ★ **当前权威投稿审计、P0/P1 行动清单与停止规则** |
| **[`docs/MANUSCRIPT_JCIM_ZH.md`](docs/MANUSCRIPT_JCIM_ZH.md)** | ★ **组装后的中文工作稿**（`python3 docs/assemble_manuscript_zh.py`） |
| **[`docs/REVISION_PHASE1_CLOUD_VS_LOCAL.md`](docs/REVISION_PHASE1_CLOUD_VS_LOCAL.md)** | ★ 本轮修订：云端已完成 vs 必须本地 |
| [`docs/JCIM_P0_COMPLETION_GUIDE.md`](docs/JCIM_P0_COMPLETION_GUIDE.md) | 历史归档/排版指南；不能替代当前科学补强清单 |
| **[`docs/JCIM_PREWRITING_CHECKLIST_V1.md`](docs/JCIM_PREWRITING_CHECKLIST_V1.md)** | ★ 写作前注意事项 + 逐项核对 |
| [`docs/RESULTS_SECTION_JCIM_EN_V1.md`](docs/RESULTS_SECTION_JCIM_EN_V1.md) | ★ JCIM 英文 Results 规范源文件 |
| [`docs/METHODS_DRAFT_ZH_JCIM_V1.md`](docs/METHODS_DRAFT_ZH_JCIM_V1.md) | 中文 Methods 工作稿 |
| [`docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`](docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md) | 中文 SI 工作稿（盒子/cognate/敏感性；仅已有数据） |
| [`docs/RESULTS_DRAFT_ZH_JCIM_V1.md`](docs/RESULTS_DRAFT_ZH_JCIM_V1.md) | 中文 Results 对齐稿 |
| [`docs/DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md`](docs/DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md) | 中文 Discussion 局限（ChEMBL max-only 等） |
| [`data/jcim_bench_v0/CLAIM_CEILING.md`](data/jcim_bench_v0/CLAIM_CEILING.md) | 可写 / 禁止写的 claim |
| [`docs/POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](docs/POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md) | Intro/Abstract：“framework”怎么写才不像包装 |
| [`docs/STATISTICAL_LOCK_V1.md`](docs/STATISTICAL_LOCK_V1.md) | ★ **投稿统计锁：Table 2/3 估计量与 Table 2 CI 唯一来源** |
| [`docs/SUBMISSION_BRANCH_POLICY_V1.md`](docs/SUBMISSION_BRANCH_POLICY_V1.md) | ★ 冻结/ChatGPT 历史快照 vs 投稿整合分支 |
| [`docs/MANUSCRIPT_LOCK_INVENTORY_V1.md`](docs/MANUSCRIPT_LOCK_INVENTORY_V1.md) | 资产来源与 CI/estimand 冲突审计 |
| [`data/jcim_strengthen_t0t1_v0/analysis/A_GROUP_VERDICT.md`](data/jcim_strengthen_t0t1_v0/analysis/A_GROUP_VERDICT.md) | 混淆对照裁决 |
| [`data/jcim_strengthen_t0t1_v0/analysis/B_GROUP_VERDICT.md`](data/jcim_strengthen_t0t1_v0/analysis/B_GROUP_VERDICT.md) | E8 / enrichment / PM110 |
| [`docs/JCIM_STRENGTHENING_PLAN_V1.md`](docs/JCIM_STRENGTHENING_PLAN_V1.md) | 加厚规划（历史+路线） |
| **[`docs/JCIM_SUPPLEMENTARY_EXPERIMENTS_PLAN_V2.md`](docs/JCIM_SUPPLEMENTARY_EXPERIMENTS_PLAN_V2.md)** | 历史补实验方案；已生成的是同批 ChEMBL unused-pool 内部敏感性，不是外部验证 |
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

## 本地运行

最小本地复现入口：

```bash
cd Dual_Target_Docking
python3 scripts/check_local_env.py
bash scripts/run_local_repro.sh
```

- 最小依赖：[`requirements-analysis.txt`](requirements-analysis.txt)
- 说明文档：[`docs/LOCAL_RUN.md`](docs/LOCAL_RUN.md)
- 该入口只重算**零对接分析**，不重跑 Vina / RTM / GNINA

重建中英文主稿：

```bash
python3 docs/assemble_manuscript_en.py
python3 docs/assemble_manuscript_zh.py
```

文献阻断 / 时间分割 / assay 机器审计（仍为零对接）：

```bash
bash scripts/run_phase1_revision.sh
```

## 复现主分析（零对接）

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
python3 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

## 状态

- 对接 + 混淆对照 + E/enrichment/PM110 + A4 max→median + B5 两对受体替换：**完成**  
- 章节稿已按 formulation 主线重构；英文主稿 `docs/MANUSCRIPT_JCIM_EN.md`；中文工作稿 `docs/MANUSCRIPT_JCIM_ZH.md`  
- **提交前科学性缺口**：MCL1 topology-aware pose QC 或证据降级、主结果 cluster uncertainty 的最终核对、K=4/无可评估外部 docking 的严格 claim ceiling；详见 `docs/JCIM_PROJECT_AUDIT_AND_ACTION_PLAN_2026-08-27.md`。
- **提交合规缺口**：排版、版本化 Zenodo DOI 与从原始记录到主表的干净环境复现。
- Zenodo DOI：（发布后填这里）
