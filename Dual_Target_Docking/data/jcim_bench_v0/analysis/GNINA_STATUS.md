# GNINA rescore — STATUS: DONE (best-of-9 CNN minimize/rescore)

**Updated:** 2026-08-24 after fair all-mode GNINA rescore

**9-pose fairness:** DONE as of 2026-08-24 (previously NOT RUN in this cloud environment on 2026-08-21; see `GNINA_NINE_POSE_SKIP_V1.md` for the superseded history). The user ran the rescore locally and pushed real results.

## Binary
`/mnt/d/CADD paper exercise/gnina/bin/gnina` (v1.3.2), CPU `--no_gpu`

## Protocol
- Input: **all** Vina `mode_01`…`mode_09` PDBQT → SDF (Open Babel)
- `gnina --cnn_scoring rescore --minimize --seed 20260727 --cpu 1`
- Per ligand–target: take **max CNNscore** over up to 9 modes
- The worst-pocket comparison (`GNINA_BEST9_STATUS.md`) uses `min(score_A, score_B)` on per-end best-of-9 CNN for both contrasts
- The **true** directional pocket-matched GNINA (D vs A_only on pocket B; D vs B_only on pocket A) is in `GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`
- mode_01-only tables retained as `scores_gnina_*_mode01_backup.csv`

## Packs
| Pack | Jobs | ok |
|------|-----:|---:|
| AChE/BChE | 1702 | 1695 |
| PIK3CA/PIK3CB | 1787 | 1787 |
| PM48 RDKit | 864 | 864 |
| EGFR/HER2 panel120 | 1978 | 1978 |
| PM110 RDKit | 2068 | 2068 |

## Tables
- Per-pack: `tables/scores_gnina_long.csv`, `tables/scores_gnina_best.csv`
- Worst-pocket comparison: `jcim_bench_v0/tables/gnina_mode01_vs_best9_*.csv` ([`GNINA_BEST9_STATUS.md`](GNINA_BEST9_STATUS.md))
- True pocket-matched comparison (K=4 + PM48/PM110 stability): `jcim_bench_v0/tables/gnina_pocket_matched_mode01_vs_best9_*.csv` ([`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`](GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md))

## Claim update
RTMScore and GNINA now share the same pose coverage (best-of-9 over the same Vina modes).
Directional pocket-matched GNINA (best-of-9) does not change which pair ranks best and
does not exceed Vina on any K=4 pair (see `GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`).
