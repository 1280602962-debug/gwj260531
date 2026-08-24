# GNINA rescore — STATUS: DONE (best-of-9 CNN minimize/rescore)

**Updated:** 2026-08-24 after fair all-mode GNINA rescore

**9-pose fairness (2026-08-21):** NOT RUN — see `GNINA_NINE_POSE_SKIP_V1.md`. Existing tables below remain mode_01 only.

## Binary
`/mnt/d/CADD paper exercise/gnina/bin/gnina` (v1.3.2), CPU `--no_gpu`

## Protocol
- Input: **all** Vina `mode_01`…`mode_09` PDBQT → SDF (Open Babel)
- `gnina --cnn_scoring rescore --minimize --seed 20260727 --cpu 1`
- Per ligand–target: take **max CNNscore** over up to 9 modes
- Pocket-matched arm uses `min(score_A, score_B)` on per-end best-of-9 CNN
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
- Comparison: `jcim_bench_v0/tables/gnina_mode01_vs_best9_*.csv`
- Full report: [`GNINA_BEST9_STATUS.md`](GNINA_BEST9_STATUS.md)

## Claim update
RTMScore and GNINA now share the same pose coverage (best-of-9 over the same Vina modes).
