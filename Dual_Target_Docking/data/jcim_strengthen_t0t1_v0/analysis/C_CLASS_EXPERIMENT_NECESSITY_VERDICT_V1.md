# C-class experiment necessity — fact-check (2026-08-24)

> **Update:** A4 (full-panel max vs median) and B5 (PIK3CA/PIK3CB receptor swap) are **done**. See `data/jcim_novelty_v0/analysis/A4_B5_STATISTICAL_AUDIT_V1.md`. Do not keep writing max pChEMBL as uncomputed, or receptor swap as a unidirectional PIK3CA collapse. The same 4JPS/5DXT crystals raise PIK3CA/PIK3CB `summary_min`. Remaining C-class items below are still optional / not required for JCIM submission.

Question: were the three remaining “C-class” experiments (max vs median pChEMBL; 1000 unused-pool panels; PLIF / rotamer / PDBFixer+Reduce) actually doable as previously described, and which still need to be run before this JCIM article?

**Short answer:** none of the three is required to keep the current claims honest. Two of the three were oversold. If anything is still worth running locally, it is **side-chain/rotamer (and optional PLIF) on materials already in git**, not a 1000-panel campaign and not PDBFixer redocking.

This note is an internal necessity verdict. Do **not** promote the 27-ligand median diagnostic below to a Supporting Information table.

---

## 1. What the previous advice got wrong

| Previous statement | Verdict | Correction |
|---|---|---|
| Frozen `mols_*.json` cannot compute median; live ChEMBL API can, **without redocking** | **Mostly true** | Cache is max-only (confirmed: dict of `chembl_id → float`). API is reachable. Relabeling reuses frozen scores. Class membership **can flip**, so n per class changes. |
| “1000 independent panels from the unused strict pool” | **False as advertised** | After main+holdout, leftover strict hard-negs cannot even supply **two** extra non-overlapping 30/30/30 panels on any pair; PIK3CA/PIK3CB and EGFR supply **zero**. See §3. |
| EGFR unused strict `B_only` = 0 | **True for strict 6.5/5.5** | All 7 strict HER2-selectives are already in EH110. Under the **actual** EGFR construction rule (θ = 6.0) leftover `B_only` = **14**, still too few for another 20/20/20. |
| Unused remaining ≈ AChE 639/141/30, PIK3CB 554/8/19, PM 1503/37/39 | **True** (strict, after main+holdout) | Recomputed from `mols_*.json` ∩ panel CSVs. IDs match. |
| ProLIF can read frozen **K=4 production** `mode_01.pdbqt` from this git repo | **False** | Main-panel pose trees are not in git (`ache_bche_panel_v0/poses` etc. absent). User almost certainly has them locally (GNINA best-of-9 was run locally; `POSE_UPLOAD_CHECKLIST.md` points at a local results tree). |
| Holdout `out.pdbqt` can support holdout-only PLIF | **True** | HOAB 120/120, HOPM 120/120, HOAP 118/120. Missing only `HOAP_028` (boron, already excluded from AUROC). |
| Crystal-swap poses are unavailable | **False / understated** | `jcim_structure_robust_v0/poses/{4JPS,5DXT,4JSX}` has **48/48** `out.pdbqt` per crystal. This is the right PLIF set for the PIK3CA collapse. |
| Rotamer on 4L23 vs 4JPS/5DXT is doable from deposited PDBs, zero docking | **True** | `4L23_protein.pdb` plus `4JPS_protein.pdb` / `5DXT_protein.pdb` are in the repo. Cα superposition is already done; side-chain is the missing layer. |
| PDBFixer+Reduce is a script re-run | **False** | New receptor prep ⇒ new docking. Not required if Limitations stay explicit. |
| BindingDB can expand the unused ChEMBL pool then resample | **Not done, and not a free lunch** | Count-level Table S12 exists. No InChIKey merge, no pChEMBL-aligned labels, no docking of BindingDB-only ligands. |

---

## 2. Experiment A — max vs median pChEMBL

### Doability

| Step | Docking? | Blocked? |
|---|---|---|
| Re-fetch assay-level `pchembl_value` from ChEMBL API | No | No (API 200 in this environment) |
| Re-label frozen ligands; reuse Vina/RTM/GNINA scores | No | No |
| Report class-flip counts **and** AUROC on the relabeled set | No | Must not pretend n is unchanged |
| confidence ≥ 8 / Homo sapiens | No | Assay endpoint, not every activity row; extra API work, still zero docking |

`T0_SKIPS.md` is right that **the cache** cannot do this. It is wrong if read as “the experiment is impossible.”

### Diagnostic (not a paper result)

Live ChEMBL pull, 2026-08-24, first-N ligands of each class (PM110: 5 dual + 5 A_only + 5 B_only; AChE/BChE: 4+4+4). 27 ligands × 2 ends = 54 ligand–target rows. All 54 API maxima **exactly match** the frozen cache (data integrity of `mols_*.json` is fine).

- 16 / 54 rows have max ≠ median (mostly well-assayed PM duals; many AChE rows are n = 1 so max = median).
- Strict 6.5/5.5: **2 / 27** class flips in this sample (CHEMBL521851 dual → gray, mTOR median 6.335; CHEMBL5281758 A_only → gray, AChE median 6.38 vs max 6.52).
- Unified **θ = 6.0** (Table 2 rule): **0 / 27** flips in this sample.

Because this is first-N, not the full panels, **do not** write “median does not change Table 2.” It only shows: (i) inflation is real on multi-assay duals; (ii) the primary θ = 6.0 rule is less brittle than strict 6.5 in this slice; (iii) a full SI table is still optional reviewer-proofing, not a claim-blocker.

Raw rows: `data/jcim_strengthen_t0t1_v0/tables/max_vs_median_diagnostic_sample_v1.csv`.

### Still needed?

**No, not for scientific integrity.** Keep Limitations item 2. Optional local SI if a reviewer insists; then run **all frozen-panel IDs**, report flips, and recompute pocket-matched AUROC on the relabeled set. Do not stop at ΔAUROC.

---

## 3. Experiment B — “1000 independent unused-pool panels”

### Leftover strict pool after main panel **and** holdout

Recomputed from `mols_*.json` with the same 6.5/5.5 rule as `build_holdout_candidate_pool_v1.py`. Used IDs = frozen panel CSV ∪ holdout CSV (PM uses PM110, which already covers PM48).

| Pair | Strict pool dual/A/B | After main | After main+holdout | Extra non-overlap 30/30/30 | Extra non-overlap 20/20/20 |
|---|---|---|---|---:|---:|
| PIK3CA/mTOR | 1552 / 80 / 81 | 1523 / 57 / 59 | **1503 / 37 / 39** | **1** | 1 |
| AChE/BChE | 687 / 189 / 78 | 659 / 161 / 50 | **639 / 141 / 30** | **1** | 1 |
| PIK3CA/PIK3CB | 602 / 56 / 67 | 574 / 28 / 39 | **554 / 8 / 19** | **0** | **0** |
| EGFR/HER2 (strict) | 951 / 39 / 7 | 925 / 22 / **0** | same (no holdout) | **0** | **0** |

Hard-negative supply, not dual supply, is the wall. The entire strict A_only universe is 80 (PM), 56 (PIK3CB), 39 (EGFR). You cannot draw 1000 **non-overlapping** 30-ligand A_only sets from 80 molecules.

### θ = 6.0 leftover (not what Methods currently promises)

If one instead resampled under the **reporting** rule θ = 6.0, leftover after main+holdout would allow ~5–6 extra non-overlapping 30/30/30 panels on the three thick pairs, and still **0** on EGFR (`B_only` leftover 14). Still not 1000. Methods 2.3 currently says **严格供给池**, so the strict table is the one that matters.

### Substitutes that are **not** the experiment

| Substitute | What it actually is | Allowed? |
|---|---|---|
| With-replacement bootstrap of already-docked ligands, B = 1000 | Near-duplicate of the existing ligand bootstrap (B = 2000) | Only if disclosed as such; **do not** call it unused-pool panel resampling |
| Dock leftover unused hard-negs, then draw 1000 **overlapping** 30/30/30 subsets | A new docking campaign (PM 76 HN, AChE 171, PIK3CB 27, EGFR-θ 183, × 2 pockets) plus a highly correlated subset distribution | Honest, expensive, still not “1000 independent panels” |
| BindingDB merge → new pool → 1000 panels | Requires identity matching, activity alignment, **and** docking | Not started; Table S12 is count-level only |

Holdout (one unused-pool draw per thick pair) is already the right, honest answer to “does the protocol only work on the construction sample?” Results 3.9 should keep that framing.

### Still needed?

**Do not run, and do not promise, 1000 independent unused-pool panels.** Manuscript language that treats leftover docking as sufficient to unlock that distribution is overstated and is corrected in Methods 2.3 / Results 3.8.

---

## 4. Experiment C — PLIF / rotamer / PDBFixer+Reduce

### What is actually in git

| Material | In this git repo? | Use |
|---|---|---|
| K=4 / PM48-RDKit / PM110 production `mode_01.pdbqt` | **No** | Main-panel PLIF only if local pose archive exists |
| Holdout `HOAB/HOAP/HOPM/poses/*/*/out.pdbqt` | **Yes** (358/360; missing HOAP_028) | Holdout wrong-pocket PLIF, scoring-free, no new docking |
| Crystal-swap PM48 poses on 4JPS / 5DXT / 4JSX | **Yes** (48 × 3) | Best PLIF set for the PIK3CA AUROC collapse |
| Receptor `*_protein.pdb` for 4L23, 4JPS, 5DXT, 4JT6, 4JSX | **Yes** | Side-chain RMSD / rotamer vs Cα (Table S10 already exists) |
| PDBFixer+Reduce-prepared receptors | **No** | Would be a new docking campaign |

Contact-count on holdout mode-1 poses is already done (Table S11). Residue-level PLIF would go **beyond** burial counts; it is not a re-run of contact-count.

### Still needed?

| Item | Needed for claim integrity? | If the user wants one more mechanism probe |
|---|---|---|
| Side-chain / rotamer on 4L23 vs 4JPS vs 5DXT pocket residues (n = 20 list in `POCKET_MECHANISM_VERDICT_V1.md`) | No (Cα already shows conservation is not sufficient) | **Highest value / lowest cost.** Zero docking. Directly addresses the open “why did 5DXT collapse at 0.343 Å pocket Cα” hole. |
| ProLIF on crystal-swap PM48 poses | No | Second. Poses are in git. Do not mix LigPrep-era PM48 poses. |
| ProLIF on holdout poses | No | Third. Poses are in git. May or may not beat contact-count; report as residue-level, not “mechanism solved.” |
| ProLIF on main K=4 | No | Only if local poses are confirmed. Git is not enough. |
| PDBFixer + Reduce + redock | No | Follow-up paper / reviewer revision at most. Keep Limitations item 6. |

---

## 5. What still needs to be done (priority)

### Science — not required before this submission

Open Discussion items that must **stay open** unless a real new experiment closes them:

1. Holdout wrong-pocket ≥ matched (survives potency/size matching; contact-count does not match Vina magnitude; main-panel direction is the opposite).
2. PIK3CA crystal-swap collapse (Cα consistent in direction, not quantitative; 5DXT local Cα 0.343 Å still collapsed).
3. PIK3CA/mTOR Δ vs descriptor baseline CI still includes 0.
4. BindingDB `equal_only` vs ChEMBL pChEMBL alignment is convention, not per-record verified.

Optional local extras, in order:

1. **Rotamer / pocket side-chain RMSD** (repo-complete, zero docking).
2. **max vs median on all frozen-panel IDs** (API + relabel, zero docking) — only to close a predictable reviewer comment.
3. **PLIF on crystal-swap and/or holdout poses already in git.**
4. Nothing that is sold as “1000 independent panels.”
5. No PDBFixer+Reduce campaign for this article.

### Submission engineering — these **are** still blocking

From `docs/JCIM_PREWRITING_CHECKLIST_V1.md`, unchanged by the C-class fact-check:

- No Abstract / Conclusions / Cover letter
- Methods / Intro / Discussion still Chinese drafts; English Results exists
- No Zenodo DOI / Data and Software Availability
- Figs 1, 2, 4, 5 missing; forest ~180 dpi
- No TOC graphic, reference library, keywords
- Main-panel poses still need a Zenodo upload (`POSE_UPLOAD_CHECKLIST.md`); they are not in git

Those are writing and data-release tasks, not new docking.

---

## 6. Claim-ceiling reminders if any optional experiment is later run

- Do not write that median “confirmed” Table 2 from the 27-ligand diagnostic.
- Do not write that overlapping subset draws of already-docked ligands are unused-pool panel resampling.
- Do not write PLIF/rotamer “solved” the holdout paradox or the crystal-swap collapse.
- Do not mix LigPrep poses into any PLIF.
- Do not dock BindingDB-only ligands and call them ChEMBL-rule hard-negs without identity and activity alignment.
