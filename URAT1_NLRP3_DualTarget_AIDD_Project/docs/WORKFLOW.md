# 当前计算工作流

生产协议 **Π\* = P2**（gnina CNNaffinity，exh=32）是**冻结战役**。任务细节见 [`LOCAL_AGENT_TASKS.md`](LOCAL_AGENT_TASKS.md)。

新产品目标（可测双靶候选）走 **C1**：[`LOCAL_C1_CANDIDATE_CAMPAIGN.md`](LOCAL_C1_CANDIDATE_CAMPAIGN.md)。不要用 `run_funnel_p2.sh` 跑 C1。

```
TrueDecoy + RandomDecoy @ 9DKB → 锁定 P2
        ↓
临床库 8319 → NLRP3 ML P(active)≥0.5 → 1588 → P2 完整案例 **1580**
        ↓
P2 对接 9DKB + 7ALV → 池内百分位 → Pareto（审计）
        ↓
Module F：dual-dock 门控 + 类药/MW/骨架审计 → 提名
        ↓
少数分子 MD（URAT1：膜+脂双层；NLRP3：水盒子）
```

```bash
# 漏斗已归档；复现审计/提名：
python3 scripts/14_candidate_nomination.py --tau 90 --mw-min 200 --mw-max 550
python3 scripts/archive_p2_export.py   # 从 docking_export_20260820 重建 slim 表与 SI
```

写作：[`MANUSCRIPT.md`](MANUSCRIPT.md) · [`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md) · [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)
