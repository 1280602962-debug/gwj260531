# GNINA 9-pose fairness — STATUS: SUPERSEDED (now done, 2026-08-24)

**This skip note is historical.** The user ran the fair best-of-9 GNINA rescore
locally (outside this cloud environment) and pushed real results on
2026-08-24. See [`GNINA_BEST9_STATUS.md`](GNINA_BEST9_STATUS.md) (worst-pocket
comparison) and
[`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`](GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md)
(true directional pocket-matched comparison, comparable to Table 2 Vina).
Headline: moving GNINA from mode-1 to best-of-9 changes pocket-matched
summary_min by −0.04 to +0.08 across the K=4 pairs and does not exceed Vina on
any pair; PIK3CA/mTOR remains the strongest pair. The remainder of this file
is kept for provenance (why this environment could not run it as of 2026-08-21).

---

**Original note (2026-08-21):** rescore all Vina modes (typically 9) with GNINA CNN so RTM (best-of-9) and GNINA are comparable.

## Why this was not executed here

1. **No `gnina` binary** in this cloud environment (`command -v gnina` fails). The recorded path from the original run (`/mnt/d/CADD paper exercise/gnina/bin/gnina`, v1.3.2) is a local Windows/WSL install, not present on this VM.
2. **Main-panel 9-pose coordinates are not in the repo.** Existing `scores_gnina_long.csv` files contain **mode_01 only**. Frozen K=4 panel directories (`ache_bche_panel_v0`, `pik3ca_pik3cb_panel_v0`, `pik3ca_mtor_panel48_rdkit_v0`, `egfr_her2_panel120_v0`) have **zero** `mode_*.pdbqt`. Holdout docking wrote multi-model `out.pdbqt`, but without GNINA that still cannot close the Table 2 channel-asymmetry item.
3. Fabricating a 9-pose GNINA table is forbidden (`CLAIM_CEILING.md`; Methods critique #5).

## What remains true in the manuscript

- GNINA results already reported are **mode_01 `--cnn_scoring rescore --minimize`**.
- Limitations must keep the sentence that RTM is best-of-9 while GNINA is mode-1; this check did **not** repair that asymmetry.
- To run later (local machine with gnina + archived Vina `out.pdbqt` / `mode_01`–`mode_09`): split models, rescore each, take CNN best-of-K, recompute pocket-matched `summary_min`. Do not mix newly rescored GNINA with old mode_01 numbers in the same primary cell.

See `data/jcim_bench_v0/analysis/GNINA_STATUS.md` and `data/jcim_strengthen_t0t1_v0/analysis/T0_SKIPS.md`.
