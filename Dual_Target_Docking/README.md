# Dual-Target Docking

独立课题：**架构无关的双靶对接决策尺子**（Dual-VSDS-Decision）。  
**与 `JNK1_Selectivity_Project/` 无关，请勿混入该目录。**

## 现行主张（一句话）

双靶分子结构种类很多，不假设乘客/moiety 机制。课题要解决的是：两端独立对接后的**朴素融合无法可靠区分 Dual 与 A-only/B-only**；交付一套可复现的决策协议（姿态 QC → 可选重打分 → 校准 → 短板/门控）与公开四类标签基准。

## 目录结构

```
Dual_Target_Docking/
├── README.md
├── docs/
│   ├── PROJECT_MASTER_PLAN.md                  # ★ 课题总览（现行）
│   ├── NMI_SUBMISSION_PLAN_DECISION_RULER.md   # ★ 现行投稿主张
│   ├── CRITIQUE_AND_NEXT_STEPS.md              # ★ 红队 + 下一步
│   ├── EGFR_HER2_DIAGNOSTIC_DEMO.md
│   ├── NMI_SUBMISSION_PLAN_MOIETY.md           # 已废弃（乘客主线）
│   └── ...
├── data/
│   ├── public_pair_selection/
│   ├── diag_egfr_her2/                         # 40 分子诊断面板
│   └── ...
└── scripts/
```

## 快速入口

| 文档 | 用途 |
|------|------|
| [**课题总览**](docs/PROJECT_MASTER_PLAN.md) | 到底解决什么问题 |
| [**红队 + 下一步**](docs/CRITIQUE_AND_NEXT_STEPS.md) | 冻结方向与最小闭环 |
| [**决策尺子投稿规划**](docs/NMI_SUBMISSION_PLAN_DECISION_RULER.md) | Claim ladder / Go-No-Go |
| [**EGFR/HER2 诊断协议**](docs/EGFR_HER2_DIAGNOSTIC_DEMO.md) | 第一张诊断表怎么跑 |
| [40 面板](data/diag_egfr_her2/panel_v0_40.csv) | Dual / A_only / B_only / neither |
| [公开靶点对](docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) | 三对冻结说明 |

## 不做

乘客/moiety 封面故事；新对接采样器；药物协同 / DTI GNN；PROTAC 三元主线；把细胞表型写成双靶结合证明。
