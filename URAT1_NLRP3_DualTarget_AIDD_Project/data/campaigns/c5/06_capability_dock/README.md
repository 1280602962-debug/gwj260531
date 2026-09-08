# URAT1 capability docking products

Local D1/D2 outputs. How to run: `docs/URAT1_CAPABILITY_DOCK.md`.

`d1_job_inventory.csv` is the 27-cell matrix. `d1_capability_summary.csv` / `d2_capability_summary.csv` are the scored capability tables.

**Deposited 2026-09-08 (this branch):** D1 20 new + 7 reused GNINA cells; D2 9 native self-docks. Pose SDFs under `d1/**/seed*` and `d2/**/seed*` are tracked; stdout/log side files stay untracked. D1 and D2 RMSD tables must not be pooled.

Engine: gnina v1.3.2, `exhaustiveness=32`, `num_modes=9`, `cnn_scoring=rescore`, CPU `--no_gpu`. Crystal refs for D2 use Open Babel bond perception (CCD five-letter residue names break bare RDKit PDB import).
