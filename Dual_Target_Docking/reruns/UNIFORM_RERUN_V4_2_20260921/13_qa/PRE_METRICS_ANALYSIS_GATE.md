# PRE_METRICS_ANALYSIS_GATE

PRE_METRICS_ANALYSIS_GATE = PASS
METRICS_IMPLEMENTATION_GATE = PASS

AUROC was not computed.
Formal metrics implementation exists at scripts/rerun_v4_2/formal_metrics_lib.py and is locked.
Even if PRE_METRICS_ANALYSIS_GATE=PASS, metrics stay locked until METRICS_IMPLEMENTATION_GATE=PASS
and a separate human unlock of METRICS_EXECUTION_UNLOCKED.

implementation_issues=[]
blocking_regression=[]

No docking, no scoring, no 8070 raw parse, no old AUROC reuse, no auto-fix.
