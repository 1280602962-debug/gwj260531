# Local recompute pack (five post-census pairs)

Machine-readable twin: `tables/five_pair_local_recompute_v1.yaml`.

Cloud is running the **zero-dock** stack on the already-scored production
Vina mode-1 file (`scripts/analyze_five_pair_stack_v1.py`). Poses are
gitignored; this machine has no Vina / RTM / GNINA. You submit the jobs
below on the same local machine that already finished Layer-3 + production
Vina.

This pack does **not** restock Table 2 and does **not** change the title.
Destination identity: `PROJECT_IDENTITY_LOCK_V1.md` (8-row main table after
withdrawing PIK3CA/PIK3CB). “Track B” remains a folder name.

Seed stays **20260727** (dock / ETKDG) and **20260729** (panel shuffle).
Holdout draw seed is **20260731**. No LigPrep. No seed 42.

---

## 0. Do now (local) vs do not do now

**Do now**

1. Five-seed Vina (reuse 20260727; redock 20260811–14).
2. RTM `rtmscore_model1` best-of-9 on production 9-mode poses.
3. GNINA 1.3.2 CNN `--cnn_scoring rescore --minimize` best-of-9 on the same poses.
4. Independent GNINA search on **JAK1/TYK2 only**.
5. *(done on cloud)* Dump-gated tables and BindingDB/PubChem count-only.
   Holdout IDs are frozen. Do not re-draw. JAK1/JAK2 drawn panel is
   20/20/18 (Murcko cap). See `local_track_b_v0/analysis/FIVE_PAIR_DUMP_GATED_V1.md`.

**Do not do now**

- Restock Table 2 or retitle the article.
- Re-dock PIK3CA/PIK3CB; 2Y3A production swap.
- CTSK/CTSS ordinary Vina.
- Crystal swaps or PM110 / E=16 on these five pairs.
- Independent GNINA on F2/F10, JAK1/JAK2, or either PPAR pair.
- Dock a holdout panel before the 20260731 member lists are frozen.
- F2/F10 protonation enum; LigPrep; seed 42.

---

## A. Five-seed Vina

Same receptors, boxes, exhaustiveness 8, `num_modes=9`, `energy_range=3`,
frozen Meeko ligands (ETKDG seed 20260727 — do **not** re-embed per Vina seed).

| Seed | Action |
|------|--------|
| 20260727 | Reuse `local_track_b_v0/tables/scores_vina_mode1_v1.csv`. Do not redock. |
| 20260811–20260814 | New docks into `local_track_b_v0/multiseed/{seed}/`. |

```bash
# from Dual_Target_Docking/data/jcim_chembl_universe_v0/
python3 scripts/dock_track_b_fiveseed_v1.py --workers 8 --timeout 600
```

Expected new jobs: **4 × 1,100 = 4,400**. Record skip/fail the same way as
production (`timeout`, `skip_torsdof>=25`). Do not silently drop failures.

Primary manuscript numbers remain seed 20260727. Five-seed is a sensitivity
(Table S54 style): directional `summary_min` and Dual-vs-neither (`vina_mean`).
Do not replace Table 2 with a multi-seed mean.

---

## B. RTM and GNINA CNN (same Vina poses)

If `local_track_b_v0/poses/` is missing, regenerate production poses first:

```bash
python3 scripts/prep_track_b_ligands_v1.py --workers 8   # if PDBQT cache missing
python3 scripts/dock_track_b_production_v1.py --workers 8 --timeout 600
```

Then:

```bash
python3 scripts/rescore_track_b_rtm_v1.py
python3 scripts/rescore_track_b_gnina_cnn_v1.py --workers 6
```

RTM: `rtmscore_model1`, best-of-9. GNINA: `--cnn_scoring rescore --minimize`,
best-of-9 CNNscore (fallback CNNaffinity). Defaults point at the same local
binaries used for K=4 (`/home/gwj/miniconda3/...`, `/mnt/d/CADD paper exercise/gnina`).
Override with `--rtm-python`, `--gnina`, `--obabel` if those paths moved.

---

## C. Independent GNINA — JAK1/TYK2 only

Original rule: independent GNINA search only on formulation-gap pairs
(EGFR/HER2, PIK3CA/mTOR). The only new pair that qualifies is **JAK1/TYK2**.

```bash
python3 scripts/dock_track_b_gnina_independent_v1.py --workers 6 --timeout 600
```

220 jobs (110 × 2). Seed 20260727, E=8, nine modes. Readout = mode-1
`minimizedAffinity`. Frozen Meeko ligands, same 6N7A / 3LXP receptors and boxes.
Per-job timeout 600 s and `TORSDOF ≥ 25` skip match production Vina
(`J1TYK2_092` already timed out on 3LXP at 600 s). Timeout → recorded skip,
partial output deleted. This is not a multi-engine bake-off.

---

## D. Dump-gated and BindingDB (done)

Already run against the frozen ChEMBL 37 dump (tarball SHA-256
`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).

- Tables: `local_track_b_v0/tables/five_pair_dump_gated_v1/`
- Verdict: `local_track_b_v0/analysis/FIVE_PAIR_DUMP_GATED_V1.md`
- BindingDB/PubChem: `local_track_b_v0/tables/five_pair_crossdb_v1/`
  and `local_track_b_v0/analysis/FIVE_PAIR_CROSSDB_V1.md`

Holdout IDs are frozen (`HOLDOUT_SEED=20260731`). JAK1/JAK2 is 20/20/18
because Murcko cap 3 blocked two leftover B-only ligands. Do not relax
the cap. Do not re-draw. Do not hard-dock BindingDB.

To reproduce only:

```bash
python3 scripts/analyze_five_pair_dump_gated_v1.py --sqlite PATH_TO_chembl_37.db
python3 scripts/five_pair_bindingdb_pubchem_count_v1.py --sqlite PATH_TO_chembl_37.db
```

---

## E. After local scores exist (back to cloud or local)

Do **not** paste class-preserving Track B CIs into Table 2. Table 2
`summary_min` CIs are ligand-level **non-stratified** resamples of the
dual+A+B pool (`unified_threshold_sensitivity_v2.csv` estimand). The
zero-dock script already computed that for production Vina. Repeat the
same estimand on five-seed / RTM / GNINA score tables when they land.

Independent GNINA analysis is JAK1/TYK2 only, same formulation-gap
contrast as EGFR/HER2.

---

## F. Claim ceiling (this pack)

Allowed: submit the local docking/rescoring jobs above; keep production
Vina 20260727 as the primary five-pair readout; leftover holdout IDs
are already drawn (JAK1/JAK2 = 20/20/18).

Forbidden: LigPrep or seed 42; expanding independent GNINA beyond
JAK1/TYK2; treating JAK1/JAK2 leftover as ineligible; CTSK Vina;
PIK3CB re-dock; restocking Table 2 from these scores before the
all-pairs stack (including local channels) is actually complete.
