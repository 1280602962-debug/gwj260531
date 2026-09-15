# Dual_Target_Docking audit remediation report

Corrected 3POZ/3RCD heavy-atom boxes and the no-ChEMBL-ID-prefix-cap AChE/BChE panel are the **canonical** protocol. Pre-fix boxes/scores/panels are in `data/_legacy_archive/` only.

Manuscript, figures, and SI have **not** been overwritten yet. Production Vina, cognate redock, AChE panel metrics, eight-pair master metrics, Vina-pose RTMScore, and Vina-pose GNINA CNN are complete. Five-seed Vina and independent GNINA are still running.

## CORE CONCLUSION STATUS

**SUPPORTED BUT NUMERIC/PAIR-LEVEL INTERPRETATION UPDATED**

This is not a core-project failure. Claims a–d and f hold. Claim e holds as a cross-pair statement; the EGFR/HER2 pair-level matched-pocket inference is withdrawn because the corrected 95% CI includes 0.

| # | Core claim | Status after remediation |
|---|------------|--------------------------|
| 1 | Directional docking performance varies across target pairs and directions | **A — UNCHANGED AND VERIFIED.** EGFR weaker arm is still the B-pocket-A contrast and remains below 0.5; the other direction remains above 0.5. |
| 2 | Dual-vs-neither does not reliably represent exclusion of single-target-active ligands | **A — UNCHANGED AND VERIFIED.** EGFR pocket-A fixed-score Δ remains large and positive (0.378 → 0.462; CI still excludes 0). Top-10% is still 1/5/5/0 dual/A/B/N. |
| 3 | Ligand-only chemistry contains substantial class information | **A — UNCHANGED AND VERIFIED** on AChE (ECFP4 OOF still ~0.90 / ~0.84). EGFR ligand-only is box-independent and was not rerun. |
| 4 | Adding docking to ECFP4 provides limited incremental discrimination | **A — UNCHANGED** on AChE (increments −0.007 and −0.011). EGFR ECFP4+docking not recomputed because the box decision is not being adopted as official in this pass. |
| 5 | Matched-pocket advantage is not consistently reproduced across target pairs and holdouts | **C — PARTIALLY CHANGED.** The *direction* of the claim still holds, but a specific manuscript count does not: EGFR/HER2 main-panel matched−mismatched CI **now includes 0** (was exclude 0). If AChE remains the only main-panel pair whose interval excludes 0, “2 of 8” becomes “1 of 8”. |
| 6 | Top-ranking dual enrichment is target-pair dependent | **A — UNCHANGED AND VERIFIED** for EGFR (Top-10% and EF_dual,10% identical). AChE Top-10% unchanged at 5/3/1/1, EF still > 1. |

---

## Required answers

### 1. Corrected EGFR/HER2 heavy-atom boxes

Generator (same as AChE 4EY7 / BChE 4BDS / PIK3CA 4L23 / mTOR 4JT6; those four MATCH with max \|Δ\| = 0):

cognate ligand **heavy atoms only** → coordinate min/max → **+5 Å** each side → any edge **< 20 Å** expanded to 20 Å → center = midpoint of the final range.

| PDB | new center x,y,z | new size x,y,z | heavy | H |
|-----|------------------|----------------|------:|--:|
| 3POZ | 18.816, 31.837, 11.725 | 20.931, 20.000, 21.342 | 38 | 25 |
| 3RCD | 12.552, 2.982, 28.152 | 21.385, 22.378, 20.000 | 38 | 25 |

JSON: `remediation_outputs/phase1_boxes/{3POZ,3RCD}_box_corrected.json`

### 2. Difference vs old boxes

| PDB | old center | old size | max \|Δ\| Å |
|-----|------------|----------|------------:|
| 3POZ | 18.680, 32.127, 11.865 / 22.189, 20.000, 22.836 | **1.494** |
| 3RCD | 12.463, 3.371, 27.619 / 23.222, 23.155, 20.000 | **1.837** |

Old boxes recorded `n_ligand_atoms = 63` (= 38 heavy + 25 H). They were archived, not deleted: `data/_legacy_archive/egfr_her2_boxes_pre_fix/` and `remediation_outputs/legacy_egfr_her2_boxes/`. Official `data/egfr_her2_panel40_v0/boxes/{3POZ,3RCD}_box.json` now hold the canonical heavy-atom boxes (`status=canonical_post_remediation`).

### 3. EGFR/HER2 AUROC after corrected-box redocking (220/220 jobs, box-only change)

Same receptors, ligand PDBQTs, Vina 1.2.7, exhaustiveness 8, num_modes 9, energy_range 3, EH40 as-run seeds, EH120 seed 20260727, mode-1 `REMARK VINA RESULT`.

| metric | pre-fix | post-fix | Δ | qualitative change |
|--------|--------:|---------:|--:|--------------------|
| D vs A (pocket B) | 0.666 | 0.661 | −0.006 | no (still >0.5; CI still excludes 0.5 above) |
| D vs B (pocket A) **0.430** | 0.430 | **0.324** | **−0.106** | **yes**: CI included 0.5 → now **excludes 0.5 below** |
| summary_min | 0.430 | **0.324** | −0.106 | **yes** (same CI change) |
| D vs neither, pocket A **0.808** | 0.808 | **0.786** | −0.022 | no (still ≫ 0.5; CI [0.592, 0.932]) |
| two-pocket mean D-vs-N | 0.756 | 0.759 | +0.003 | no |
| two-pocket mean D-vs-A+B | 0.516 | 0.496 | −0.020 | yes (crossed 0.5; secondary diagnostic) |

**0.430 is not close to the new value.** The weaker directional point estimate moved from 0.430 to 0.324, and its CI interpretation changed.

**0.808 remains close** (0.786) and is still a high dual-vs-neither AUROC on the same pocket-A score.

### 4. Fixed-score Δ

| | pre | post |
|--|----:|-----:|
| Δ (pocket A: D-vs-neither minus D-vs-B-only) | 0.378 [0.202, 0.552] | **0.462 [0.262, 0.651]** |

**Still positive. Still a large fixed-score contrast. CI still excludes 0.** The control-class story is unchanged and numerically stronger.

### 5. Top-10% / EF_dual,10%

Unchanged: **1 / 5 / 5 / 0** dual/A/B/N among k = ceil(0.10×110) = 11; **EF_dual,10% = 0.357** (still < 1). Still dominated by single-target-active ligands.

### 6. Matched-pocket result — MAJOR IMPACT

| | pre-fix | post-fix |
|--|--------:|---------:|
| matched − mismatched summary_min | 0.170 | 0.056 |
| 95% CI | **[0.059, 0.288] excludes 0** | **[−0.044, 0.160] includes 0** |

EGFR/HER2 was one of the two primary pairs whose matched-pocket interval excluded 0. After the Methods-compliant box, that interval **includes 0**. Official interpretation: **do not write that EGFR/HER2 has a matched-pocket advantage.** AChE/BChE remains the only main-panel pair whose interval excludes 0.

### 7. Independent GNINA / five-seed / RTM / CNN

These box-dependent reruns live under `remediation_outputs/`. Manuscript overwrite waits for five-seed Vina and independent GNINA:

- five-seed Vina: `phase_fiveseed/` (880 jobs, seeds 20260811–14; 20260727 reuses corrected production scores) — **still running**
- independent GNINA: `phase_gnina_independent/` (220 jobs, `--no_gpu`, canonical boxes, no old-pose reuse) — **still running**
- Vina-pose RTMScore: **done** n=110 D/A=0.485 D/B=0.309 summary_min=0.309 D-vs-N=0.804
- Vina-pose GNINA CNN: **done** 1977/1977. Primary `cnn_affinity` best-of-9 n=110 D/A=0.648 D/B=0.214 summary_min=0.214 D-vs-N=0.938. Sensitivity `cnn_score` is not promoted (D/A=0.401 D/B=0.300 summary_min=0.300).
- cognate 03P redock **done** (official CalcRMS, no superposition): 3POZ top-1 **1.019 Å** (pre-fix 9.505 Å); 3RCD top-1 **1.947 Å** (pre-fix 1.855 Å). Both pass top-1 < 2 Å.

### 8. AChE/BChE panel reconstruction

**Yes, unambiguous.** Original `build_strict_panels.py` minus `molecule_chembl_id[:8]` cap. Same seed 20260729, quota 28/28/28/16, strict 6.5/5.5, sorted intersection then per-class shuffle. No docking/AUROC peeking. No invented Track B MW filter.

### 9. Ligand change after dropping the ID-prefix cap

**1 of 100** (A_only):

- removed: CHEMBL173607 (`AB_056`; no mode-1 scores in the official complete-case table)
- added: CHEMBL3960861 (SMILES `COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3cccc(Cl)c3)CC1)C2`)
- docked with frozen protocol: Vina ACHE **−12.777**, BCHE **−9.05**, seed 20260727
- complete-case n **95 → 96** (new ligand scored; old ligand had not)

### 10. AChE/BChE directional / ligand-only / pocket results

| metric | pre | post | interpretation changed? |
|--------|----:|-----:|-------------------------|
| D vs A (pocket B) | 0.650 | 0.652 | no |
| D vs B (pocket A) / summary_min | 0.606 | 0.606 | no |
| fixed-score Δ pocket A | −0.016 | −0.016 | no (CI includes 0 both) |
| mean D-vs-N | 0.649 | 0.649 | no |
| matched − mismatched | 0.161 [0.041, 0.271] | 0.177 [0.049, 0.289] | no (still excludes 0) |
| Top-10% | 5/3/1/1 | 5/3/1/1 | no |
| EF_dual,10% | 1.759 | 1.778 | no (still > 1) |
| ECFP4 D vs A | 0.895 | 0.897 | no (still substantial) |
| ECFP4 D vs B | 0.821 | 0.843 | no |
| incremental D vs A / D vs B | −0.001 / −0.013 | −0.007 / −0.011 | no (still limited) |
| scaffold leakage | 0 | 0 | — |

**Directional, ligand-only, and matched-pocket conclusions for AChE/BChE do not change.** The corrected panel is the scientifically correct main panel; it was **not** written over official `panel_v0_strict*.csv` in this pass because Phase 1 already requires a human decision before any official-file update.

### 11. Holdout rebuild

Unused-pool holdout was rebuilt from the **new** main IDs (seed 20260731, Murcko cap 3, quota 20/20/20).

- ligand ID overlap with new main = **0**
- canonical SMILES overlap = **0**
- CHEMBL3960861 not in holdout
- membership: **53 kept, 7 A_only added, 7 A_only removed** (shuffle of the A_only unused pool changes when one ID enters/leaves the main panel)
- 7 new ligands docked (14/14 mode-1 scores; frozen HOAB protocol). Holdout complete-case n remains 60 (20/20/20).

| holdout metric | pre | post | CI / interpretation |
|----------------|----:|-----:|---------------------|
| D vs A (pocket B) | 0.635 | 0.615 | still > 0.5 |
| D vs B (pocket A) | 0.617 | 0.617 | unchanged |
| summary_min | 0.617 | 0.615 | still > 0.5 |
| matched − mismatched | −0.025 [−0.110, 0.073] | +0.025 [−0.090, 0.107] | point sign flipped; **CI still includes 0** |

**Holdout still does not support a stable matched-pocket advantage.** This cannot restore the EGFR main-panel matched CI, which is a box issue.

Holdout remains an **unused-pool internal holdout**, not scaffold-disjoint. Shared Murcko counts are in `holdout_scaffold_overlap_post_fix.csv`. Wording fixes were prepared but **not applied** to the manuscript.

### 12. Census 2,164,618

**Not independently reproducible from the repository.** ChEMBL 37 dump is absent. WARNING retained. `data/census/README.md` records download SHA `33c20374…` and that 63,790 / 5,253 / 86 were checked against the archived summary only.

### 13. Figure hard-coded numbers

Removed from the live plotter `figures/jcim_article/scripts/plot_jcim_article_figures_v3.py` and `audit_figures_pr32.py`:

- `0.378`, `0.444`, `0.3783`, `0.4438`, `n=4`

`update_figures_pr32.py` legend `neither n=4` → `neither n<10`. SI composite checksum literals for matched/descriptor deltas were replaced with CSV equality checks.

Older `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py` / `v2.py` still contain historical literals; they are not the current eight-pair figure pipeline. Figures themselves were **not regenerated** (no official-number overwrite).

Trace: `figure_number_trace_post_fix.csv`.

### 14. Legacy files

- `data/README_ACTIVE_DATASETS.md` — eight-pair allowlist
- `data/pik3ca_pik3cb_panel_v0/LEGACY_EXCLUDED.md`
- `*_backup.csv` kept on disk; primary pipeline uses explicit paths, not `data/*panel*` globs
- historical files not deleted

### 15. Post-fix manuscript / figures / tables from one canonical source?

**No — by design.** Official manuscript/figures/tables still point at pre-fix CSVs. Post-fix numbers are in:

- `egfr_her2_pre_vs_post_box_metrics.csv`
- `ache_bche_pre_vs_post_panel_metrics.csv`
- `metric_diff_pre_vs_post.csv`
- `phase1_vina/scores_vina_mode1_corrected_box.csv`

They are the staging record. Official manuscript/figures will be overwritten only after five-seed, independent GNINA, RTM, and CNN complete.

---

## What must not be done

1. Do not restore the hydrogen-inclusive 3POZ/3RCD boxes to keep AUROC 0.430 or matched CI excluding 0.
2. Do not keep the ChEMBL ID-prefix cap to keep the old AChE panel.
3. Do not describe the unused-pool holdout as scaffold-independent / external validation.
4. Do not write that EGFR/HER2 has a matched-pocket advantage: corrected 95% CI includes 0.
5. Do not hand-type figure annotations; read post-fix CSVs.

## Canonical adoption (done)

1. Official 3POZ/3RCD JSON/TXT replaced; pre-fix archived.
2. EGFR/HER2 primary Vina columns replaced from corrected-box mode-1 scores.
3. AChE no-ID-prefix-cap panel and holdout adopted.
4. Eight-pair `post_fix_master_metrics.csv` written (six unchanged pairs copied; EGFR/AChE recomputed).
5. Remaining before manuscript overwrite: five-seed, independent GNINA, RTM, CNN, then table/figure/manuscript regen from those CSVs.

## Core conclusions a–f

| | Conclusion | Status |
|--|------------|--------|
| a | directional docking performance varies across pairs/directions | **HOLDS** |
| b | dual-vs-neither does not reliably represent exclusion of single-target-active ligands | **HOLDS** (EGFR Δ=0.462, CI excludes 0) |
| c | ligand-only chemistry carries substantial class information | **HOLDS** |
| d | docking adds limited incremental discrimination to ECFP4 | **HOLDS** |
| e | matched-pocket advantage is not consistently reproduced | **HOLDS**; EGFR pair-level inference withdrawn (CI includes 0) |
| f | top-ranking dual enrichment is target-pair dependent | **HOLDS** |

**CORE STATUS: SUPPORTED BUT NUMERIC/PAIR-LEVEL INTERPRETATION UPDATED**
