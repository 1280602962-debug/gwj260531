# Track-B final execution status

Date: 2026-09-18  
Locked yaml (historical design snapshot; **do not edit `do_not_now`**): `data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml`

The yaml `do_not_now` list is the 2026-09-04 freeze of what that local-run contract deferred. Later work completed some of those items as **sensitivity / independent / rescoring** channels. Completing them does not rewrite the yaml.

Deleted `analysis/DOCKING_PLAN_V1.md` is **git history only**. It is not a current instruction. Current Track-B pack README: `data/jcim_chembl_universe_v0/local_track_b_v0/README.md`.

## yaml `do_not_now` versus later execution

| yaml `do_not_now` item | Later status | Where |
|------------------------|--------------|-------|
| `f2_f10_protonation_sensitivity` | still not done | remains deferred; not a current gap for writing freeze |
| `pik3cb_2y3a_e32_cognate_retest` | still not done | not a primary pair |
| `gnina_independent_search` | done as **INDEPENDENT GNINA** on JAK1/TYK2 (Track-B) plus EGFR/HER2 and PIK3CA/mTOR | `results/canonical/computational_robustness.csv`; not eight-pair GNINA |
| `rtm_or_gnina_cnn_rescore` | done as **SAME-POSE RESCORING** on deposited Vina poses | `scores_rtm_best9_v1.csv`; `scores_gnina_cnn_best9_v1.csv` |
| `five_seed_vina` | done as five-seed Vina on the five Track-B pairs | `local_track_b_v0/tables/multiseed/scores_vina_mode1_seed{seed}.csv` |
| `ctsk_ctss_vina` | still not done | pair excluded (covalent); correctly out of primary |
| `pik3ca_pik3cb_redock` | still not done | not a primary pair |

## Protocol levels (current)

See `docs/PROTOCOL_LEVELS.md`. Track-B production Vina is PRIMARY for the five ordinary pairs (E=8). Items later executed from `do_not_now` are not promoted into PRIMARY.
