# Dual-Target Docking

独立课题：单分子双靶配体的 **passenger / moiety-resolved** 对接评测与开放基准。  
**与 `JNK1_Selectivity_Project/` 无关，请勿混入该目录。**

## 现行主张（一句话）

整分子独立对接会把第二药效团算进单口袋分数；用 **moiety 计分 + 分靶校准 + 短板双靶决策** 改善 dual-vs-single 排序，并发布 Dual-VSDS-Moiety。

## 目录结构

```
Dual_Target_Docking/
├── README.md
├── docs/
│   ├── NMI_SUBMISSION_PLAN_MOIETY.md      # ★ 现行投稿 / 执行规划
│   ├── PROJECT_MASTER_PLAN.md             # ★ 课题总览（精简）
│   ├── EGFR_HER2_DIAGNOSTIC_DEMO.md       # ★ 第一张诊断表怎么跑
│   ├── PUBLIC_TARGET_PAIR_SELECTION_REPORT.md
│   ├── NMI_REFERENCE_PAPER_PLAYBOOK.md     # 高分文对标（流程参考）
│   ├── DUAL_TARGET_SCORING_IMPLEMENTATION.md  # 打分组件（校准/短板等基线）
│   ├── DUAL_MULTI_TARGET_DOCKING_SURVEY.md
│   ├── REFERENCES_AND_MOLECULES.md
│   ├── DUAL_TARGET_COCRYSTAL_CATALOG_NOTES.md
│   ├── DATA_COLLECTION_FOR_SYNTHETIC_CHEMISTS.md
│   └── DUAL_TARGET_PAPERS_*.md            # 文献编目
├── data/
│   ├── public_pair_selection/             # 冻结靶点对 + ChEMBL 四类
│   ├── literature/
│   ├── schema/
│   └── dual_target_structures/
└── scripts/
    ├── audit_public_target_pairs.py
    └── export_egfr_her2_fourclass.py
```

## 快速入口

| 文档 | 用途 |
|------|------|
| [**红队清单 + 现在最该干什么**](docs/CRITIQUE_AND_NEXT_STEPS.md) | 关键质疑与最小实验决策树（先做这个） |
| [**NMI / 执行规划**](docs/NMI_SUBMISSION_PLAN_MOIETY.md) | 主张、实验骨架、Go/No-Go |
| [**课题总览**](docs/PROJECT_MASTER_PLAN.md) | 是什么 / 做什么 / 不做什么 |
| [**EGFR/HER2 诊断协议**](docs/EGFR_HER2_DIAGNOSTIC_DEMO.md) | whole-mol vs moiety 怎么跑 |
| [**公开靶点对**](docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) | 三对冻结说明 |
| [四类分子 CSV](data/public_pair_selection/egfr_her2_fourclass_chembl_ids.csv) | Dual / A-only / B-only 名单 |
| [文献与分子](docs/REFERENCES_AND_MOLECULES.md) | DOI + 共晶种子 |
| [高分对标](docs/NMI_REFERENCE_PAPER_PLAYBOOK.md) | VSDS-VD 等流程怎么学 |

## 不做

新对接采样器；药物协同 / DTI GNN；PROTAC 三元主线；把细胞表型写成双靶结合证明。
