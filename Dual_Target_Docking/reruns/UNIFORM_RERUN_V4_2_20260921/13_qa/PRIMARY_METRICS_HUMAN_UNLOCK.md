# PRIMARY_METRICS_HUMAN_UNLOCK

Authorized: yes
Scope: formal V4.2 PRIMARY metrics only
Run: UNIFORM_RERUN_V4_2_20260921

Unlock mechanism:
- Source `scripts/rerun_v4_2/formal_metrics_lib.py` keeps `METRICS_EXECUTION_UNLOCKED = False`.
- Only `scripts/rerun_v4_2/run_primary_metrics_v42.py` and `scripts/rerun_v4_2/post_primary_metrics_qa_v42.py` may set the module flag True in-process.
- Other QA guards stay in place.

This round computes only primary deltas:
- M1−M0
- M1b−M1
- M2−M0
- M3−M0

Secondary remain secondary and are not computed:
- M3−M1b
- M3−M1

Forbidden:
- old scripts/analysis/*
- old M1/M3 masters
- current_score_master.csv
- old AUROC / old sensitivity
- population change, theta change, imputation, auto-fix
- sensitivity after PASS
