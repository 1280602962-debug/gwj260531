# URAT1 / NLRP3 痛风双节点 — 临床药物重定位

面向高尿酸血症/痛风的 **URAT1（代谢）+ NLRP3（炎症）** 双节点计算漏斗：先在 TrueDecoy/RandomDecoy 上锁定开源对接协议 **Π\* = P2（gnina CNNaffinity）**，再做 NLRP3 ML 缩库 → 双靶对接 → Pareto 审计 → 化学提名。

**写作入口：** [`docs/MANUSCRIPT.md`](docs/MANUSCRIPT.md)  
**引言 / Methods / Results：** [`docs/INTRO_DRAFT_CN.md`](docs/INTRO_DRAFT_CN.md) · [`docs/METHODS_DRAFT_CN.md`](docs/METHODS_DRAFT_CN.md) · [`docs/RESULTS_DRAFT_CN.md`](docs/RESULTS_DRAFT_CN.md)  
**协议结果：** [`docs/PROTOCOL_SELECTION_RESULT.md`](docs/PROTOCOL_SELECTION_RESULT.md)  
**本地漏斗：** [`docs/LOCAL_AGENT_TASKS.md`](docs/LOCAL_AGENT_TASKS.md) · `bash scripts/run_funnel_p2.sh`

**首投：** *Molecular Diversity*（拒稿后转 *JCAMD*）。不声称已验证双靶抑制剂。

## 科学定位

> 选择对接协议，而不是默认打分；NLRP3 分类只缩库，URAT1 由结构排序；Pareto 前沿不等于跟进名单。

## 三套数据（禁止混用）

| 数据集 | 用途 |
|--------|------|
| 临床库 8319 / P≥0.5 池 1588 / P2 完整案例 **1580** | 主筛选 |
| TrueDecoy / RandomDecoy | 只选 Π\* |
| 8973 distill | 仅 URAT1 回顾（可选 SI） |

## 文档

| 文档 | 内容 |
|------|------|
| [文稿指南](docs/MANUSCRIPT.md) | 结构、claim、lead |
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
