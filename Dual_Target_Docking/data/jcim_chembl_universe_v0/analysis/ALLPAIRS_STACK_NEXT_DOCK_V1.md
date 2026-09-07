# All-pairs stack — next docking queue (aligned to PR tip 2026-09-05)

Canonical locks on branch `cursor/chembl-exhaustive-pair-census-0b1a`:

- Identity: `analysis/PROJECT_IDENTITY_LOCK_V1.md` → **8 main-table rows**
- Local docking: `analysis/LOCAL_RECOMPUTE_PACK_V1.md` + `tables/five_pair_local_recompute_v1.yaml`
- Holdout IDs (frozen): `local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_*.csv`

Primary production Vina for the five new pairs is **done**. Zero-dock + dump-gated
+ BindingDB count-only are **done on cloud**. Do **not** start local jobs until
unrelated machine docking finishes and a human says go.

## Queue (submit in this order)

| Step | What | Scope | ~Jobs | Official script |
|------|------|-------|------:|-----------------|
| **A** | GNINA independent search | **JAK1/TYK2 only** | 220 | `dock_track_b_gnina_independent_v1.py` |
| **B** | Five-seed Vina | all 5 pairs | 4400 | `dock_track_b_fiveseed_v1.py` (reuse 20260727) |
| **C** | RTM best-of-9 | all 5 | ~550 | `rescore_track_b_rtm_v1.py` |
| **D** | GNINA CNN rescore (`rescore --minimize`) | all 5 | ~9900 modes | `rescore_track_b_gnina_cnn_v1.py` |
| **E** | Holdout Vina | frozen IDs | ~596 | dock frozen `HO*` panels — **do not re-draw** |

Seeds: primary **20260727** (reuse); extra **20260811–14**; holdout draw **20260731** (already drawn).

Holdout drawn sizes: 20/20/20 ×4; JAK1/JAK2 = **20/20/18** (Murcko cap 3). Report 18 B-only; do not relax cap.

## Explicitly not in this queue

- Restock Table 2 / retitle article (wait until stack complete)
- PIK3CA/PIK3CB redock; 2Y3A production; CTSK/CTSS Vina
- Independent GNINA on F2/F10, JAK1/JAK2, or either PPAR pair
- Crystal swaps / PM110 / E=16 on the five new pairs
- Re-draw holdout; BindingDB hard-dock; LigPrep; seed 42
- F2/F10 protonation (still deferred)

## Protocol notes (must keep)

- Same Meeko ligands across Vina seeds (ETKDG 20260727) — **do not re-embed**
- Timeout 600 s → skip; `TORSDOF ≥ 25` skip (GNINA independent same as Vina)
- Known prior timeout: `J1TYK2_092` @ 3LXP
- 9V8H: keep PG08-NL peptide; 3LXP: JH1 not JH2; 5U3Q: agonist 1 / altloc A
- Primary manuscript numbers stay seed 20260727; five-seed is sensitivity only
- Table 2 CIs = non-stratified bootstrap (not class-preserving Track B CIs)

## Submit when free

```bash
cd Dual_Target_Docking/data/jcim_chembl_universe_v0
python3 scripts/dock_track_b_gnina_independent_v1.py --workers 6 --timeout 600
python3 scripts/dock_track_b_fiveseed_v1.py --workers 8 --timeout 600
python3 scripts/rescore_track_b_rtm_v1.py
python3 scripts/rescore_track_b_gnina_cnn_v1.py --workers 6
# then dock frozen holdout panels (IDs already in five_pair_dump_gated_v1/)
```
