# Dual-Target Docking

独立课题目录：双靶/多靶小分子对接方法调研、打分融合设计、共晶编目与问答总结。  
**与 `JNK1_Selectivity_Project/` 无关，请勿混入该目录。**

## 目录结构

```
Dual_Target_Docking/
├── README.md
├── docs/
│   ├── NMI_DUAL_COMPATIBILITY_PLAN.md           # ★ 总方案（问答+分析+路线图）
│   ├── PROJECT_MASTER_PLAN.md                   # ★ 课题总规划（是什么/为什么/怎么做/多远）
│   ├── NMI_REFERENCE_PAPER_PLAYBOOK.md          # ★ NMI/高分文：最可参考流程对标
│   ├── researchstudio_audit/                    # ★ ResearchStudio 思路/可完成性审计
│   ├── ars_audit/                               # ★ academic-research-skills (ARS) 课题分析
│   ├── DUAL_MULTI_TARGET_DOCKING_SURVEY.md      # 文献调研
│   ├── DUAL_TARGET_SCORING_IMPLEMENTATION.md    # 打分实现蓝图
│   ├── DUAL_TARGET_COCRYSTAL_CATALOG_NOTES.md   # 共晶目录说明
│   ├── DUAL_TARGET_DOCKING_QA_SUMMARY.md        # 问答总结
│   └── REFERENCES_AND_MOLECULES.md              # 文献链接 + 双靶分子总表
└── data/
    ├── schema/                                  # CSV 表头模板
    │   ├── activity_pairs.schema.csv
    │   ├── docking_runs.schema.csv
    │   └── metrics_report.schema.csv
    └── dual_target_structures/
        ├── README.md
        └── dual_target_cocrystal_catalog.csv     # 共晶种子目录
```

## 快速入口

| 文档 | 用途 |
|------|------|
| [**课题总规划**](docs/PROJECT_MASTER_PLAN.md) | 是什么 / 为什么 / 怎么做 / 创新点 / 离落地多远 |
| [**数据收集（合成向）**](docs/DATA_COLLECTION_FOR_SYNTHETIC_CHEMISTS.md) | 大白话：要收什么表、怎么填、两周怎么干 |
| [**JMC级双靶文献100篇**](docs/DUAL_TARGET_PAPERS_JMC100.md) | 含DOI链接的文献池（优先新文） |
| [**摘要复核分类100篇**](docs/DUAL_TARGET_PAPERS_JMC100_INDEX.md) | 合成+活性确认，并标 fused/linked 等 |
| [**NMI 总方案**](docs/NMI_DUAL_COMPATIBILITY_PLAN.md) | 全部思考/分析/数据与对接方案/阶段路线 |
| [**NMI 对标文章**](docs/NMI_REFERENCE_PAPER_PLAYBOOK.md) | 流程最像的 NMI/高分文 + 章节骨架怎么抄 |
| [**ResearchStudio 审计**](docs/researchstudio_audit/RESEARCHSTUDIO_AUDIT.md) | 思路漏洞、scoop、可完成性（Idea skills） |
| [**ARS 课题分析**](docs/ars_audit/ARS_TOPIC_ANALYSIS.md) | FINER RQ、方法蓝图、DA 攻击面、AI 失败模式 |
| [**K-Dense skills 映射**](docs/KDENSE_SKILLS_MAPPING.md) | scientific-agent-skills 对本课题能帮什么 |
| [文献与分子总表](docs/REFERENCES_AND_MOLECULES.md) | 全部参考文献链接、文章简介、双靶分子/PDB 信息 |
| [问答总结](docs/DUAL_TARGET_DOCKING_QA_SUMMARY.md) | 问题诊断、创新点、指标、PROTAC 边界 |
| [文献调研](docs/DUAL_MULTI_TARGET_DOCKING_SURVEY.md) | 方法与评价综述 |
| [实现蓝图](docs/DUAL_TARGET_SCORING_IMPLEMENTATION.md) | 校准 + softmin + 硬负样本融合 |
| [共晶说明](docs/DUAL_TARGET_COCRYSTAL_CATALOG_NOTES.md) | Tier / Morphy 分类要点 |

## 核心结论（一句话）

可发表增量在 **任务级校准与短板敏感融合**，不在再做一个通用 docking sampler；PROTAC 公开数据只验证 linked/双功能分支。
