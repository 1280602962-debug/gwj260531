# 开源对接栈

生产排序读出是 **Π\* = P2（gnina CNNaffinity）**，见 [`PROTOCOL_SELECTION_RESULT.md`](PROTOCOL_SELECTION_RESULT.md) 与 `config/docking_production_p2.yaml`。  
AutoDock Vina 用于协议比较中的 P1，不是临床库主表。

| 步骤 | 工具 |
|------|------|
| 受体 | gemmi + Open Babel → PDBQT，pH 7.4 |
| 配体 | RDKit ETKDG + Meeko |
| 生产对接 | **gnina**，`cnn_scoring=rescore`，`score_mode=cnnaff`，exh=32 |
| 协议对照 | Vina 1.2.5（P1）；RTMScore 仅敏感性（P4/P5，不生产） |

临床库一键：`bash scripts/run_funnel_p2.sh`（[`LOCAL_AGENT_TASKS.md`](LOCAL_AGENT_TASKS.md)）。

| 靶点 | PDB | 中心 (Å) | 盒边长 (Å) |
|------|-----|----------|------------|
| URAT1 | 9DKB | 99.97, 102.97, 105.70 | 22³ |
| NLRP3 | 7ALV | 16.76, 35.45, 125.71 | 20³ |

`dock_score` 存为越低越好（CNNaffinity 取负）。只做池内百分位，不与实验 \(K_i\)、也不与历史 Glide 分混比。  
P2 自对接：集合内可出现近原生姿，但 CNNaffinity Top-1 可能失败——见 [`REDOCK_SMOKE_ANALYSIS.md`](REDOCK_SMOKE_ANALYSIS.md)。
