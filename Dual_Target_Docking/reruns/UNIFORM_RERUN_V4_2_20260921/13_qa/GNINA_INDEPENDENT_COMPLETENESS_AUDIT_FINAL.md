# GNINA INDEPENDENT COMPLETENESS AUDIT FINAL

Source: `10_gnina_independent_docking/GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv`
Generated from the `status` column. Not hand-edited.

- expected: **1592**
- accounted: **1592**
- SUCCESS: **1586**
- TIMEOUT: **6**
- AUROC: **NOT_COMPUTED**

This is the unique formal V4.2 M3 completeness record.

`GNINA_INDEPENDENT_COMPLETENESS_AUDIT_PRE_TOPOLOGY_REPAIR.md` reports 1590 SUCCESS / 2 TIMEOUT. That file is HISTORICAL_PRE_TOPOLOGY_REPAIR_ONLY and is not formal M3 authority.

Why 1590/2 became 1586/6:

- 1590/2 is the pre-topology-repair intermediate state.
- A frozen affected set of 29 ligands / 61 physical jobs then received a topology representation repair (see `00_protocol/M3_MACROCYCLE_INPUT_REPAIR_ADDENDUM.yaml`).
- 4 repaired jobs timed out.
- Those 4 plus the original 2 timeouts give the final 6 TIMEOUT jobs.
- Final formal M3 is therefore 1586 SUCCESS / 6 TIMEOUT.

n_unique_physical_jobs = 1592
n_physical_job_seed is not used for M3 (M3 is seed42 only).
