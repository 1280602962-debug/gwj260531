# All-pairs stack local run — complete (2026-09-07)

Destination: same 8-row article (`PROJECT_IDENTITY_LOCK_V1.md`).
Does not restock Table 2 / retitle.

## Completed locally

| Step | Result |
|------|--------|
| A GNINA independent JAK1/TYK2 | `tables/scores_gnina_independent_jak1_tyk2_v1.csv` (210 ok / 2 timeout skip / 8 parse fail) |
| B five-seed Vina | `tables/multiseed/scores_vina_mode1_seed{20260727,20260811–14}.csv` |
| C RTM best-of-9 | `tables/scores_rtm_best9_v1.csv` (pockets via ProDy; no `-gen_pocket`) |
| D GNINA CNN best-of-9 | `tables/scores_gnina_cnn_best9_v1.csv` (9742/9742) |
| E holdout Vina | `allpairs_stack/holdout/tables/holdout_scores_vina_mode1_v1.csv` (593 ok / 3 skip) |

Poses/logs are gitignored; score tables and pocket PDBs are in this commit.
