# T0 skips

- T0.7 median/confidence≥8/Homo sapiens: mols_*.json stores single float (max pChEMBL) per ChEMBL ID; no per-assay median or confidence fields available locally.
- GNINA 9-pose fairness rescore (2026-08-21): not run in this cloud environment (no `gnina` binary; frozen K=4 panel pose coordinates for modes 2–9 were not in the repo at the time). **Superseded 2026-08-24**: the user ran the rescore locally and pushed real results (`data/jcim_bench_v0/analysis/GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`). Pocket-matched summary_min moves by only −0.04 to +0.08 across K=4 and never exceeds Vina.
