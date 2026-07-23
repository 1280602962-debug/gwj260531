# Dual-VSDS CSV schemas

Header templates for the Dual-Target Docking pipeline. Example rows are placeholders only.

| File | Purpose |
|------|---------|
| `activity_pairs.schema.csv` | Paired dual/A-only/B-only/inactive labels + design_type |
| `docking_runs.schema.csv` | Per-target docking scores, PoseBusters, optional RTMScore |
| `metrics_report.schema.csv` | Fusion method evaluation report |

See `docs/NMI_DUAL_COMPATIBILITY_PLAN.md` for field semantics and claim boundaries.
