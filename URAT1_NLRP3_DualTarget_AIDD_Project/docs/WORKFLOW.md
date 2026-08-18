# 当前计算工作流

生产协议 **Π\* = P2**（gnina CNNaffinity，exh=32）。任务细节见 [`LOCAL_AGENT_TASKS.md`](LOCAL_AGENT_TASKS.md)。

```
TrueDecoy + RandomDecoy @ 9DKB → 锁定 P2
        ↓
临床库 8319 → NLRP3 ML P(active)≥0.5 → ~1588
        ↓
P2 对接 9DKB + 7ALV → 池内百分位 → Pareto（审计）
        ↓
Module F：dual-dock 门控 + 类药/MW/骨架审计 → 提名
        ↓
少数分子 MD（URAT1：膜+脂双层；NLRP3：水盒子）
```

```bash
# 生产漏斗（本地已有 gnina 时）
JOBS=8 bash scripts/run_funnel_p2.sh

python3 scripts/10_admet_druglikeness.py
python3 scripts/11_chemical_space_novelty.py
python3 scripts/13_pareto_robustness.py
python3 scripts/14_candidate_nomination.py --tau 90 --mw-min 200 --mw-max 550

python3 scripts/select_md_candidates.py --n-novel 2 --n-controls 2
python3 scripts/export_md_ready_candidates.py
```

跟进名单读 `results/candidates/nominated_shortlist_diverse.csv`，不要用裸 `pareto_shortlist.csv`。  
分数只在同一引擎、同一靶点、同一池内比百分位；不与历史 Glide 分混表。
