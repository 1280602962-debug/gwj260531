# URAT1 / NLRP3 痛风双节点 — 临床药物重定位

面向高尿酸血症/痛风的 **URAT1（代谢）+ NLRP3（炎症）** 双节点临床药物重定位。现有 P2（gnina CNNaffinity）结果作为冻结基线；新版路线将姿态生成、姿态质控、排序与偏倚校正拆分，并通过靶点特异验证和多证据门控重新确认双靶候选。

**新版项目计划：** [`docs/PROJECT_REPLAN_MOLECULAR_DIVERSITY.md`](docs/PROJECT_REPLAN_MOLECULAR_DIVERSITY.md)（§12：哪些必须本地）  
**C1 本地战役（可测双靶候选）：** [`docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md`](docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md) · [`config/campaign_c1.yaml`](config/campaign_c1.yaml) · [`config/docking_c1.yaml`](config/docking_c1.yaml)  
**C2 课题重建路线（转运循环阻断，预登记）：** [`docs/PROJECT_ROUTE_C2_TRANSPORT_CYCLE.md`](docs/PROJECT_ROUTE_C2_TRANSPORT_CYCLE.md) · [`config/campaign_c2.yaml`](config/campaign_c2.yaml)  
**现有写作入口：** [`docs/MANUSCRIPT.md`](docs/MANUSCRIPT.md)（Results 仍是冻结 P2；C1 数字未出）  
**引言 / Methods / Results：** [`docs/INTRO_DRAFT_CN.md`](docs/INTRO_DRAFT_CN.md) · [`docs/METHODS_DRAFT_CN.md`](docs/METHODS_DRAFT_CN.md) · [`docs/RESULTS_DRAFT_CN.md`](docs/RESULTS_DRAFT_CN.md)  
**协议结果：** [`docs/PROTOCOL_SELECTION_RESULT.md`](docs/PROTOCOL_SELECTION_RESULT.md)  
**冻结 P2 漏斗（已归档，不要当 C1）：** [`docs/LOCAL_AGENT_TASKS.md`](docs/LOCAL_AGENT_TASKS.md) · `bash scripts/run_funnel_p2.sh`

**首投：** *Molecular Diversity*（拒稿后转 *JCAMD*）。不声称已验证双靶抑制剂。

## 科学定位

> 不把单一对接分解释为实验活性；NLRP3 分类只提供生物学缩库证据，双靶候选由靶点特异结构验证、偏倚校正、相互作用、稳定性和药物化学证据共同提名；Pareto 前沿不等于跟进名单。

## 三套数据（禁止混用）

| 数据集 | 用途 |
|--------|------|
| 临床库 8319 / P≥0.5 池 1588 / P2 完整案例 **1580** | 主筛选 |
| TrueDecoy / RandomDecoy | 只选 Π\* |
| 8973 distill | 仅 URAT1 回顾（可选 SI） |

## 文档

| 文档 | 内容 |
|------|------|
| [新版项目计划](docs/PROJECT_REPLAN_MOLECULAR_DIVERSITY.md) | Molecular Diversity 投稿路线、协议重构、补算包与候选提名规则 |
| [C1 本地战役](docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md) | 可测双靶候选的预注册执行书：Rank/Acid 双轨、必须本地的 gnina/MD、停止规则 |
| [C2 转运循环阻断路线](docs/PROJECT_ROUTE_C2_TRANSPORT_CYCLE.md) | 课题重建：以"与转运循环不相容"替代占据型打分排序；机制匹配负控制、T1/T2/T3 可观测量、预登记通过与失败标准 |
| [文稿指南](docs/MANUSCRIPT.md) | 现有结构、claim、lead |
| [工作流](docs/WORKFLOW.md) | 现行命令 |
| [差异化](docs/DIFFERENTIATION_VS_PLK1_NLRP3.md) | vs 激酶–NLRP3 文 |
| [数据事实](docs/DATA_FACT_CHECK.md) | 投稿前数字 |
| [开源对接](docs/OPEN_SOURCE_DOCKING.md) | P2 栈 |
| [重对接烟雾](docs/REDOCK_SMOKE_ANALYSIS.md) | pose vs 排序 |

## 快速开始

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
pip install -r requirements.txt
python3 scripts/00_prepare_data.py
python3 scripts/02_train_asymmetric_models.py --no-oat-transfer
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto
JOBS=8 bash scripts/run_funnel_p2.sh
```

TrueDecoy 重建与对接池说明见 `data/benchmarks/urat1_true_decoy/`。

## 许可

MIT License — 见 [LICENSE](LICENSE)。
