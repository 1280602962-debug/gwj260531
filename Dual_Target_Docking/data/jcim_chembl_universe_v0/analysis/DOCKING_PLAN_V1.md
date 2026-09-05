# Docking plan (locked 2026-09-04) — local Vina pack

Machine-readable twin: `tables/track_b_local_run_v1.yaml`.

This is the pack to run **now** on a local machine that has AutoDock Vina 1.2.7. Cloud finished panel extraction and Meeko ligand prep status; cloud has **no Vina**.

**Do now:** receptor PDBQT + cognate boxes → Layer-3 cognate RMSD → production Vina on five pairs.

**Do not do now** (agreed, deferred):

1. F2/F10 open-source protonation sensitivity.
2. 2Y3A / GD9 cognate retest at exhaustiveness 32.

Do not start those two until this pack is finished and a later written go-ahead exists.

This file froze the **local Vina pack** (five new pairs, seed 20260727 / 20260729). It does not by itself restock Table 2 from census ranks, and it does not rewrite 2026-07-23 as an eight-pair freeze.

Destination article identity is now `PROJECT_IDENTITY_LOCK_V1.md`: same paper, 8 main-table rows after withdrawing PIK3CA/PIK3CB, five new pairs get the already-written all-pairs analysis stack. Seed stays **20260727** (dock / ETKDG) and **20260729** (panel shuffle). Do not switch to 42.

---

## 0. Locked decisions (do not reopen while running)

1. **CTSK/CTSS is out of ordinary noncovalent Vina.** Both holos are reversible-covalent (4X6H Cys25–3XT 1.83 Å; 9GJ2 Cys25–KH0 1.78 Å). 4X6H/9GJ2 stay as structure records only.
2. **PIK3CA/PIK3CB already-docked numbers are kept** as a declared special case, like EGFR/HER2. Table 2 0.500 [0.350, 0.650] = human-PIK3CB labels scored in a **mouse p110δ** pocket (2WXF / O35904). Do **not** re-dock. Do **not** replace 2WXF with 2Y3A/4BFR in this pack. 2Y3A E32 is a later cognate-only check, not a production swap.
3. **New Vina is five pairs only.** Extension / replication panel. Count **three new systems** (coagulation, JAK, PPAR), not five pairs.
4. **Prep is RDKit ETKDGv3 + meeko.** No LigPrep. Protonation / tautomer enum is a later F2/F10-only sensitivity, not a protocol change for this run.

---

## A. Already docked — do not re-run

| Pair | Role | Receptors | Notes |
|------|------|-----------|-------|
| PIK3CA/mTOR | ordinary thick primary | 4L23 / X6K, 4JT6 / X6K | Keep. ATP not FRB. |
| AChE/BChE | ordinary thick primary | 4EY7 / E20, 4BDS / THA | Keep. |
| EGFR/HER2 | special case: supply-limited | 3POZ / 3RCD | Keep EH110. min HN = 7. |
| PIK3CA/PIK3CB | special case: receptor-identity | 4L23 + **2WXF = mouse p110δ** | Keep numbers. Not isoform control. |

---

## B. What to dock now (Track B)

Five pairs, eight unique receptors (JAK1 and PPARA reused). Panels already extracted (strict 6.5/5.5, small-mol filter, seed 20260729, 32/32/32/14 = 110).

| # | Pair | Panel CSV | min HN | Receptor A | Receptor B | Declared site |
|---|------|-----------|-------:|------------|------------|---------------|
| 1 | F2/F10 | `tables/track_b_panels/panel_F2_F10_v1.csv` | 108 | **4UDW / N6L** | **2JKH / BI7** | S1 / S1–S4 |
| 2 | JAK1/TYK2 | `tables/track_b_panels/panel_JAK1_TYK2_v1.csv` | 91 | **6N7A / KEV** | **3LXP / IZA** | JH1 ATP both ends |
| 3 | JAK1/JAK2 | `tables/track_b_panels/panel_JAK1_JAK2_v1.csv` | 53 | **6N7A / KEV** | **8BXH / C87** | JH1 ATP both ends |
| 4 | PPARG/PPARA | `tables/track_b_panels/panel_PPARG_PPARA_v1.csv` | 82 | **9V8H / BRL** | **6LXA / EPA** | LBD; 9V8H is ternary (+PG08-NL) |
| 5 | PPARA/PPARD | `tables/track_b_panels/panel_PPARA_PPARD_v1.csv` | 82 | **6LXA / EPA** | **5U3Q / 7UJ** | LBD; 7UJ is agonist 1, not PEG |

Volume: 5 × 110 × 2 = **1,100** production jobs + 8 × 9 cognate poses. JAK1 (6N7A) and PPARA (6LXA) are prepared once and reused.

---

## C. Explicitly not in this pack

| Item | Why |
|------|-----|
| CTSK/CTSS Vina | Covalent holos. |
| CREBBP/BRD4, GPCR, SLC6, F2/PRSS1, JAK3/TYK2, OPRM1/OPRK1 | Failed a gate or undeclared site. |
| PIK3CA/PIK3CB re-dock / 2Y3A production | Special case kept; 2Y3A E32 is deferred and cognate-only. |
| F2/F10 protonation enum | Deferred. |
| GNINA independent search, RTM / GNINA CNN rescore, five-seed | **Authorized now** as the local recompute pack (`LOCAL_RECOMPUTE_PACK_V1.md`). Independent GNINA is JAK1/TYK2 only. |
| LigPrep; seed 42 | Forbidden protocol change. |

---

## D. Frozen search settings (same as Methods)

| Parameter | Value |
|-----------|--------|
| Engine | AutoDock Vina **1.2.7**, default `vina` function |
| Ligand prep | largest organic fragment → RDKit AddHs → **ETKDGv3 seed 20260727** → MMFF ≤ 200 → meeko 0.7.1 |
| Dock seed | **20260727** |
| `num_modes` | 9 |
| `energy_range` | 3 kcal mol⁻¹ |
| Exhaustiveness | **8** (fallback **16** only if Layer 3 fails the &lt; 2 Å gate at E=8) |
| Box | cognate heavy-atom AABB + **5 Å / axis**, any edge &lt; **20 Å** raised to 20 |
| Cognate gate | best-of-9 heavy-atom RMSD **&lt; 2.0 Å** (pose-generation check, not top-1 recovery) |
| Skip | `TORSDOF ≥ 25` (same as frozen panel dockers) |
| Receptor prep | water + cognate removed; `mk_prepare_receptor.py`; default altloc A |

**9V8H:** keep the **PG08-NL peptide** in the receptor. Remove water and **BRL** only. Do not pretend this is a binary PPARγ–rosiglitazone LBD.

**3LXP:** IZA is in **JH1**, not JH2. Do not swap to a JH2 crystal.

RMSD is not identity proof (2WXF passed at 0.405 Å). Layer 2 already signed these eight PDBs; Layer 3 only tests pose generation.

---

## E. Local directory (create on the Vina machine)

```text
Dual_Target_Docking/data/jcim_chembl_universe_v0/local_track_b_v0/
  receptors/{4UDW,2JKH,6N7A,3LXP,8BXH,9V8H,6LXA,5U3Q}_receptor.pdbqt
  cognates/{4UDW_N6L,2JKH_BI7,6N7A_KEV,3LXP_IZA,8BXH_C87,9V8H_BRL,6LXA_EPA,5U3Q_7UJ}.{pdb,pdbqt,sdf}
  boxes/{4UDW,2JKH,6N7A,3LXP,8BXH,9V8H,6LXA,5U3Q}_box.json
  ligands_sdf/          # from prep_track_b_ligands_v1.py
  ligands_pdbqt/
  cognate_qc/
  poses/{PDB}/{panel_id}/mode_01.pdbqt … mode_09.pdbqt
  logs/
  tables/job_status.csv
  tables/layer3_cognate_rmsd_v1.csv
  tables/scores_vina_mode1_v1.csv
```

Do not commit large pose trees unless a later upload checklist says so.

---

## F. Execution order (stop if a gate fails)

Work from repo root `Dual_Target_Docking/`. Paths below are relative to `data/jcim_chembl_universe_v0/`.

### F1. Regenerate ligand PDBQT (cloud binaries are gitignored)

```bash
python3 scripts/prep_track_b_ligands_v1.py --workers 8
```

Expect **550/550 ok** (`tables/track_b_ligand_prep_status_v1.csv`). Copy or point production jobs at `cache/track_b_ligands/pdbqt/{panel_id}.pdbqt`. Same seed 20260727 as Methods. No LigPrep.

### F2. Fetch the eight PDBs and build receptor / cognate / box

For each row in the table in §B:

1. Download the asymmetric-unit PDB from RCSB (do not use a biological assembly that duplicates the pocket).
2. Confirm polymer UniProt via SIFTS matches the lock table (`tables/receptor_freeze_v1.csv`).
3. Extract protein ATOM/TER (plus **PG08-NL peptide on 9V8H**).
4. Extract the first instance of the listed cognate CCD (`N6L`, `BI7`, `KEV`, `IZA`, `C87`, `BRL`, `EPA`, `7UJ`).
5. Box = cognate heavy-atom AABB + 5 Å/axis, min edge 20 Å. Write `boxes/{PDB}_box.json` with `center_*` and `size_*`.
6. `mk_prepare_receptor.py` → `{PDB}_receptor.pdbqt`.
7. Prepare the cognate ligand the same way as panel ligands (RDKit if needed + meeko), or from the crystal coordinates for Layer 3 only.

Reuse 6N7A and 6LXA; do not rebuild them per pair.

### F3. Layer 3 — cognate best-of-9 RMSD (hard stop)

For each of the eight receptors, Vina-redock the cognate into that receptor’s box:

```text
vina --receptor receptors/{PDB}_receptor.pdbqt
     --ligand   cognates/{PDB}_{CCD}.pdbqt
     --center_x … --center_y … --center_z …
     --size_x … --size_y … --size_z …
     --exhaustiveness 8 --num_modes 9 --energy_range 3
     --seed 20260727 --cpu 8
     --out cognate_qc/{PDB}_{CCD}_out.pdbqt
```

Gate: **best-of-9 heavy-atom RMSD &lt; 2.0 Å**. If E=8 fails, one prespecified fallback at **E=16**. If both fail, **stop that receptor** — do not production-dock it and “check later.”

Write `local_track_b_v0/tables/layer3_cognate_rmsd_v1.csv` with PDB, CCD, exhaustiveness, top-1 RMSD, best-of-9 RMSD, pass/fail.

Production Vina starts only after **all eight** pass.

### F4. Production Vina

Both ends of each of the five panels. Same box, seed, modes, energy range. Exhaustiveness = the value that passed Layer 3 for that receptor (8 unless the fallback was required).

≈ 1,100 jobs. Record every job:

`N_attempted / N_successful / N_failed` and a reason (`timeout`, `skip_torsdof`, parse error). Do not silently drop failures.

Primary readout: Vina **mode-1** energy. Convert to `S = −E` for AUROC (higher is better), same as Table 2.

Suggested layout mirrors `data/ache_bche_panel_v0/scripts/dock_panel.py` (config file + `seed = 20260727`).

### F5. Scores and AUROC (after docking, still this pack)

On unified θ = 6.0 labels **and** on the panel’s strict 6.5/5.5 class column (already in the CSVs):

- Dual vs A-only → pocket **B**
- Dual vs B-only → pocket **A**
- `summary_min = min(AUROC_D/A, AUROC_D/B)`
- Dual vs neither (`vina_mean`) as a formulation contrast only; do not treat JAK1/TYK2 neither n=14 as a thick neither law

Count **systems**: JAK two pairs = one system; PPAR two pairs = one system.

Do **not** replace Table 2 with these numbers. Do **not** pick more census pairs after seeing AUROCs.

---

## G. Cloud already finished (do not redo unless files are missing)

| Step | Status | Artifact |
|------|--------|----------|
| Four-state panels | done 2026-09-04 | `tables/track_b_panels/panel_*_v1.csv` |
| Small-mol min HN vs QC | 108 / 91 / 53 / 82 / 82 | `tables/track_b_panel_summary_v1.csv` |
| Meeko prep status | 550/550 ok on cloud | `tables/track_b_ligand_prep_status_v1.csv` (re-run F1 locally) |

---

## H. Deferred — written so they are not forgotten, not authorised now

### H1. F2/F10 protonation sensitivity (later)

Open-source protonation enum on the F2/F10 panel only (e.g. Dimorphite-DL). Re-prep + re-dock that pair. Compare directional AUROC and S1-cation poses to the RDKit default. **Not LigPrep. Not a K=4 redo. Not applied to JAK/PPAR in the first pass.**

### H2. 2Y3A cognate at E=32 (later)

Mouse Pik3cb (Q8BTI9), cognate GD9, exhaustiveness **32**, best-of-9 RMSD only. Purpose: nail whether the historical ~3.85 Å E8/E16 fail stands at higher search. **Not** a production re-dock of PIK3CA/PIK3CB. **Not** a licence to call 2Y3A the Table 2 receptor.

---

## I. How to report (when scores exist)

- Estimand unchanged: pocket-matched directional AUROC; `summary_min` is a worst-direction summary, not a new score.
- 9V8H: Methods sentence on PG08-NL.
- JAK ChEMBL labels are target-level, not JH1 vs JH2. Docking tests the **declared JH1 ATP pocket**.
- Do not pool PIK3CA/PIK3CB 0.500 into a kinase-isoform success story.
- Do not pool CTSK/CTSS (not docked).

## J. Claim ceiling (this pack)

Allowed: run this five-pair local Vina pack; keep PIK3CA/PIK3CB as a receptor-identity special case; list H1/H2 as deferred.

Forbidden: LigPrep or seed 42; CTSK/CTSS Vina; 2Y3A/2WXF production swap in this run; replacing Table 2; calling nine pairs a complete dual-target suite; starting H1/H2 before this pack is done.
