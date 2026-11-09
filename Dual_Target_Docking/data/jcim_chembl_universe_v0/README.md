# ChEMBL universe pair census (v0)

Post-hoc **dump-level** enumeration of dual-target four-state supply among human `SINGLE PROTEIN` targets in ChEMBL 37.

This is **not** Table 2. Current pair inclusion is `scripts/analysis/analysis_config.py` `PRIMARY_PAIRS` and `data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv`. Current Track-B protocol lock: `tables/track_b_local_run_v1.yaml` (`do_not_now` is a historical design snapshot; see `docs/TRACK_B_FINAL_EXECUTION_STATUS.md`).

Deleted `analysis/*.md` files (`DOCKING_PLAN_V1.md`, `TIER1_DOCKING_ROSTER_V1.md`, `FEASIBLE_PAIR_LADDER_V1.md`, and the rest of that folder) are **git history only / historical / unavailable**. They are not current instructions.

| Path | Role |
|------|------|
| `tables/` | Pair counts, eligibility audit, Track-B panels |
| `tables/pair_eligibility_audit_s14_v1.csv` | Current pair inclusion/exclusion audit |
| `tables/track_b_local_run_v1.yaml` | Locked Track-B production yaml |
| `tables/track_b_panels/` | 110-ligand CSVs (strict 6.5/5.5 construction, seed 20260729) |
| `tables/track_b_panel_summary_v1.csv` | Pool sizes and sampled class counts |
| `data/jcim_chembl_universe_v0/scripts/prep_track_b_ligands_v1.py` | Recovered ETKDGv3 seed 20260727 + meeko 0.7.1 (method snapshot) |
| `tables/track_b_ligand_prep_status_v1.csv` | Recovered 550-row prep status; exact ligand PDBQT not deposited |
| `tables/receptor_freeze_v1.csv` | Locked PDB / Vina cognate / declined alternative per receptor |
| `local_track_b_v0/` | Production scores and poses; see that README |
| `cache/` | Local dump (gitignored) |

Primary analysis labels use θ=6.0. Strict 6.5/5.5 is panel construction only (`docs/LABEL_PROVENANCE.md`).
