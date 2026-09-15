# Pre-submission technical audit

- generated: 2026-09-15 05:48 UTC
- git commit: `6f5bea67a7aca5ffc5c9abdbf6f2226790615431`
- canonical table: `audit_outputs/canonical_ligand_table.csv` (1214 rows)
- official manuscript/figures/result CSVs were not modified

# Executive summary

- **CRITICAL:** 0
- **MAJOR:** 2
- **MINOR:** 1
- **WARNING:** 14
- **PASS:** 242
- invariant tests: 98 PASS / 0 FAIL

A. Core directional AUROC — Yes. Independent Mann–Whitney AUROC and sklearn roc_auc_score recover Table 2 after S=−E and D/A→B, D/B→A.

B. Fixed-score ΔAUROC — Yes. Dual-vs-neither and dual-vs-selective used the same pocket score. EGFR/HER2 0.430→0.808 (Δ=0.378) and JAK1/TYK2 Δ=0.444 were recovered.

C. summary_min bootstrap — Yes. B=2000, pooled dual+A+B, non-stratified, shared duals, min inside replicate, percentile CI. Valid replicates=2000 for all eight pairs.

D. ECFP4 leakage — No scaffold train/test overlap in independent GroupKFold. Incremental Δ used the same folds.

E. Matched/mismatched — Yes. Mismatched only swaps existing S_A/S_B. Main CIs excluding 0: EGFR/HER2 and AChE/BChE. All 7 holdout CIs include 0.

F. Top-10% / EF_dual,10% — Yes. k=⌈0.10 n⌉, descending S_mean, D+A+B+N=k, EF uses full-panel D/n. Totals 33/23/22/3.

G. Figure / Table / manuscript consistency — Mostly. Flagship and Table 2/3 cells match independent 3-decimal rounding. Remaining issues are Methods/box wording, census unverifiability, figure literals, and AChE ID-prefix cap — not contradictory AUROC tables.

H. Problems that would change the main conclusion — No independent error was found that reverses the main scientific claims (control-class dependence, limited incremental docking over ECFP4, inconsistent matched-pocket advantage, pair-dependent enrichment, limited external generalization). EGFR/HER2 box geometry does not match the written heavy-atom rule; that is a protocol-documentation defect whose effect on scores was not re-docked.

# Critical issues

_None._

# Major issues

### MAJOR-008: 3POZ docking box is not the Methods heavy-atom AABB
- **severity:** MAJOR
- **affected files:** data/egfr_her2_panel40_v0/boxes/3POZ_box.json
- **affected code lines:** Methods 2.4.1; Table S2a; n_ligand_atoms=63 in box JSON
- **affected metrics:** EGFR/HER2 docking box geometry (not AUROC point estimates, which used the deposited box)
- **affected figures/tables:** Table S2
- **current manuscript value:** center/size [18.68, 32.127, 11.865, 22.189, 20.0, 22.836]; Methods: cognate heavy-atom AABB +5 Å / min 20 Å
- **independent recalculation:** heavy n=38 [18.816, 31.837, 11.725, 20.931, 20.0, 21.342]; all-atom n=63 [18.65, 32.143, 11.888, 22.173, 20.0, 22.902]; max|Δ| vs JSON heavy=1.494 all-atom=0.066
- **reason:** Deposited JSON records n_ligand_atoms=63, matching the hydrogen-inclusive SDF (38 heavy + 25 H). Independent heavy-atom boxes differ by >1 Å in at least one edge. AChE/PIK3CA boxes match heavy-atom Methods exactly. Scores were generated with the deposited boxes, so primary AUROCs remain internally consistent but the EGFR/HER2 site definition does not match the written protocol. Impact on AUROC is not verifiable without re-docking.
- **impact on conclusion:** Does not independently change Table 2 numbers. It does mean EGFR/HER2 Methods/Table S2 overstate heavy-atom box construction. Re-docking could change EGFR/HER2 scores.
- **recommended fix:** Either recompute EGFR/HER2 boxes from cognate heavy atoms and re-dock, or revise Methods/Table S2 to state that EGFR/HER2 boxes used hydrogen-inclusive cognate coordinates from the original preparation.

### MAJOR-009: 3RCD docking box is not the Methods heavy-atom AABB
- **severity:** MAJOR
- **affected files:** data/egfr_her2_panel40_v0/boxes/3RCD_box.json
- **affected code lines:** Methods 2.4.1; Table S2a; n_ligand_atoms=63 in box JSON
- **affected metrics:** EGFR/HER2 docking box geometry (not AUROC point estimates, which used the deposited box)
- **affected figures/tables:** Table S2
- **current manuscript value:** center/size [12.463, 3.371, 27.619, 23.222, 23.155, 20.0]; Methods: cognate heavy-atom AABB +5 Å / min 20 Å
- **independent recalculation:** heavy n=38 [12.552, 2.982, 28.152, 21.385, 22.378, 20.0]; all-atom n=63 [12.453, 3.283, 27.835, 23.233, 22.981, 20.0]; max|Δ| vs JSON heavy=1.837 all-atom=0.216
- **reason:** Deposited JSON records n_ligand_atoms=63, matching the hydrogen-inclusive SDF (38 heavy + 25 H). Independent heavy-atom boxes differ by >1 Å in at least one edge. AChE/PIK3CA boxes match heavy-atom Methods exactly. Scores were generated with the deposited boxes, so primary AUROCs remain internally consistent but the EGFR/HER2 site definition does not match the written protocol. Impact on AUROC is not verifiable without re-docking.
- **impact on conclusion:** Does not independently change Table 2 numbers. It does mean EGFR/HER2 Methods/Table S2 overstate heavy-atom box construction. Re-docking could change EGFR/HER2 scores.
- **recommended fix:** Either recompute EGFR/HER2 boxes from cognate heavy atoms and re-dock, or revise Methods/Table S2 to state that EGFR/HER2 boxes used hydrogen-inclusive cognate coordinates from the original preparation.

# Minor issues

### MINOR-010: AChE/BChE construction used a ChEMBL-ID prefix cap, not a Murcko cap
- **severity:** MINOR
- **affected files:** data/ache_bche_panel_v0/scripts/build_strict_panels.py
- **affected code lines:** 19, 71-75
- **affected metrics:** AChE/BChE panel membership (not AUROC formula)
- **affected figures/tables:** Table 1 Methods sentence on remaining pairs having no additional scaffold cap
- **current manuscript value:** the remaining pairs had no additional scaffold cap
- **independent recalculation:** MURCKO_CAP=5 applied to molecule_chembl_id[:8]
- **reason:** The AChE builder names a Murcko cap but applies it to an 8-character ChEMBL ID prefix because SMILES were deferred. This is a weak diversity proxy, not docking-based selection.
- **impact on conclusion:** Does not change directional definitions. Slightly qualifies the 'no additional scaffold cap' sentence.
- **recommended fix:** Revise Methods to say AChE used an ID-prefix diversity cap of 5 during construction, or drop the claim for that pair.

# Warnings

### WARNING-001: JAK1/JAK2 holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 8
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-002: JAK1/TYK2 holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 8
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-003: PIK3CA/mTOR holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 2
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-004: AChE/BChE holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 9
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-005: F2/F10 holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 12
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-006: PPARG/PPARA holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 15
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-007: PPARA/PPARD holdout shares Murcko scaffolds with main
- **severity:** WARNING
- **affected files:** n/a
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Table S6
- **current manuscript value:** n/a
- **independent recalculation:** 17
- **reason:** Unused-pool holdouts are not scaffold-split. Overlap is leftover chemistry.
- **impact on conclusion:** Does not invalidate the leftover-pool robustness check; must not be described as strict external validation.
- **recommended fix:** Keep current Methods wording (internal unused-pool check).

### WARNING-011: 2,164,618 pair-universe count is not independently recomputable from deposited pair lists
- **severity:** WARNING
- **affected files:** data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv
- **affected code lines:** n/a
- **affected metrics:** n_pairs_n_both_ge_1
- **affected figures/tables:** Figure 1C; Methods 2.3
- **current manuscript value:** 2,164,618
- **independent recalculation:** NOT VERIFIABLE FROM CURRENT REPOSITORY (summary CSV only; unordered pair list / ChEMBL 37 SQLite dump not deposited)
- **reason:** 63,790 / 5,253 / 86 can be checked against the summary row and later gated tables. The n_both>=1 universe requires the full unordered pair enumeration, which is not in the repository.
- **impact on conclusion:** Does not affect eight-pair docking metrics. Figure 1C lead count cannot be rebuilt from ligand-level files here.
- **recommended fix:** Deposit the pair-enumeration table or the SQLite extract that produced n_pairs_n_both_ge_1.

### WARNING-012: Figure audit script hard-codes 0.378 and 0.444
- **severity:** WARNING
- **affected files:** figures/jcim_article/scripts/audit_figures_pr32.py
- **affected code lines:** 31-32
- **affected metrics:** Figure 2A flagship ΔAUROC
- **affected figures/tables:** Figure 2A
- **current manuscript value:** 0.378 / 0.444
- **independent recalculation:** 0.3783 / 0.4438
- **reason:** Plotted values also exist in plotted_values.json and equal-score CSVs. The literals are cross-checks, not the only source. Independent points match.
- **impact on conclusion:** None if plotted_values.json remains the plot input. Risk if someone edits only the CSV.
- **recommended fix:** Assert against plotted_values.json / result CSV rather than bare floats.

### WARNING-013: Figure 2C n=4 annotation is a plotter literal
- **severity:** WARNING
- **affected files:** figures/jcim_article/scripts/plot_jcim_article_figures_v3.py
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** Figure 2C
- **current manuscript value:** PIK3CA/mTOR neither n=4
- **independent recalculation:** 4
- **reason:** Independent complete-case neither count is 4, so the literal is currently correct but not data-driven.
- **impact on conclusion:** None while n remains 4.
- **recommended fix:** Annotate from the neither count in the plotted table.

### WARNING-014: legacy/duplicate path present: data/egfr_her2_panel120_v0/tables/scores_gnina_best_mode01_backup.csv
- **severity:** WARNING
- **affected files:** data/egfr_her2_panel120_v0/tables/scores_gnina_best_mode01_backup.csv
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** n/a
- **current manuscript value:** n/a
- **independent recalculation:** n/a
- **reason:** Backup or withdrawn-pair artifacts remain in the tree. Manuscript primary analyses do not use PIK3CA/PIK3CB.
- **impact on conclusion:** Confusion risk if a later script reads backup columns (best-of-9 vs mode-1).
- **recommended fix:** Leave files in place; document canonical paths in submission_pack/INVENTORY.md.

### WARNING-015: legacy/duplicate path present: data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_gnina_long_mode01_backup.csv
- **severity:** WARNING
- **affected files:** data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_gnina_long_mode01_backup.csv
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** n/a
- **current manuscript value:** n/a
- **independent recalculation:** n/a
- **reason:** Backup or withdrawn-pair artifacts remain in the tree. Manuscript primary analyses do not use PIK3CA/PIK3CB.
- **impact on conclusion:** Confusion risk if a later script reads backup columns (best-of-9 vs mode-1).
- **recommended fix:** Leave files in place; document canonical paths in submission_pack/INVENTORY.md.

### WARNING-016: legacy/duplicate path present: data/ache_bche_panel_v0/tables/scores_gnina_best_mode01_backup.csv
- **severity:** WARNING
- **affected files:** data/ache_bche_panel_v0/tables/scores_gnina_best_mode01_backup.csv
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** n/a
- **current manuscript value:** n/a
- **independent recalculation:** n/a
- **reason:** Backup or withdrawn-pair artifacts remain in the tree. Manuscript primary analyses do not use PIK3CA/PIK3CB.
- **impact on conclusion:** Confusion risk if a later script reads backup columns (best-of-9 vs mode-1).
- **recommended fix:** Leave files in place; document canonical paths in submission_pack/INVENTORY.md.

### WARNING-017: legacy/duplicate path present: data/pik3ca_pik3cb_panel_v0/
- **severity:** WARNING
- **affected files:** data/pik3ca_pik3cb_panel_v0/
- **affected code lines:** n/a
- **affected metrics:** n/a
- **affected figures/tables:** n/a
- **current manuscript value:** n/a
- **independent recalculation:** n/a
- **reason:** Backup or withdrawn-pair artifacts remain in the tree. Manuscript primary analyses do not use PIK3CA/PIK3CB.
- **impact on conclusion:** Confusion risk if a later script reads backup columns (best-of-9 vs mode-1).
- **recommended fix:** Leave files in place; document canonical paths in submission_pack/INVENTORY.md.

# Passed checks

- PASS: EGFR/HER2 primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: EGFR/HER2 complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: JAK1/JAK2 primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: JAK1/JAK2 complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: JAK1/TYK2 primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: JAK1/TYK2 complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: PIK3CA/mTOR primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: PIK3CA/mTOR complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: AChE/BChE primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: AChE/BChE complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: F2/F10 primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: F2/F10 complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: PPARG/PPARA primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: PPARG/PPARA complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: PPARA/PPARD primary classes equal independently assigned θ=6.0 from pChEMBL
- PASS: PPARA/PPARD complete-case primary ligands have both pChEMBL values (missing not treated as inactive)
- PASS: TEST 15 all primary analyses use θ=6.0 labels
- PASS: TEST 16 strict 6.5/5.5 not used as primary labels (six quota-sampled pairs: θ=6.0 coincides with sampled strict classes; EGFR/PIK3CA used θ=6.0 sampling)
- PASS: TEST 1 EGFR/HER2 D/A uses B score (AUROC=0.666353)
- PASS: TEST 2 EGFR/HER2 D/B uses A score (AUROC=0.429688)
- PASS: TEST 3 EGFR/HER2 dual is positive class
- PASS: TEST 4 EGFR/HER2 AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 EGFR/HER2 summary_min=min(two directional points) (0.429688)
- PASS: TEST 7 EGFR/HER2 Smean=(SA+SB)/2
- PASS: EGFR/HER2 D/A(score B)=0.6664 differs from incorrect D/A(score A)=0.6983
- PASS: EGFR/HER2 summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 EGFR/HER2 pocket A Δ uses identical score channel
- PASS: TEST 8 EGFR/HER2 Top10% D+A+B+N=k (1/5/5/0 k=11)
- PASS: TEST 9 EGFR/HER2 k=ceil(0.10 n) (n=110 k=11)
- PASS: TEST 10 EGFR/HER2 EF formula (0.357143)
- PASS: TEST 18 EGFR/HER2 matched/mismatched only swaps assignment (Δsummary_min=0.1696 CI=[0.059,0.288])
- PASS: EGFR/HER2 S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 JAK1/JAK2 D/A uses B score (AUROC=0.588379)
- PASS: TEST 2 JAK1/JAK2 D/B uses A score (AUROC=0.727539)
- PASS: TEST 3 JAK1/JAK2 dual is positive class
- PASS: TEST 4 JAK1/JAK2 AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 JAK1/JAK2 summary_min=min(two directional points) (0.588379)
- PASS: TEST 7 JAK1/JAK2 Smean=(SA+SB)/2
- PASS: JAK1/JAK2 D/A(score B)=0.5884 differs from incorrect D/A(score A)=0.6074
- PASS: JAK1/JAK2 summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 JAK1/JAK2 pocket A Δ uses identical score channel
- PASS: TEST 8 JAK1/JAK2 Top10% D+A+B+N=k (6/5/0/0 k=11)
- PASS: TEST 9 JAK1/JAK2 k=ceil(0.10 n) (n=110 k=11)
- PASS: TEST 10 JAK1/JAK2 EF formula (1.875000)
- PASS: TEST 18 JAK1/JAK2 matched/mismatched only swaps assignment (Δsummary_min=-0.0190 CI=[-0.088,0.054])
- PASS: JAK1/JAK2 S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 JAK1/TYK2 D/A uses B score (AUROC=0.575101)
- PASS: TEST 2 JAK1/TYK2 D/B uses A score (AUROC=0.364919)
- PASS: TEST 3 JAK1/TYK2 dual is positive class
- PASS: TEST 4 JAK1/TYK2 AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 JAK1/TYK2 summary_min=min(two directional points) (0.364919)
- PASS: TEST 7 JAK1/TYK2 Smean=(SA+SB)/2
- PASS: JAK1/TYK2 D/A(score B)=0.5751 differs from incorrect D/A(score A)=0.5151
- PASS: JAK1/TYK2 summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 JAK1/TYK2 pocket A Δ uses identical score channel
- PASS: TEST 8 JAK1/TYK2 Top10% D+A+B+N=k (1/3/7/0 k=11)
- PASS: TEST 9 JAK1/TYK2 k=ceil(0.10 n) (n=109 k=11)
- PASS: TEST 10 JAK1/TYK2 EF formula (0.319648)
- PASS: TEST 18 JAK1/TYK2 matched/mismatched only swaps assignment (Δsummary_min=-0.0645 CI=[-0.151,0.033])
- PASS: JAK1/TYK2 S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 PIK3CA/mTOR D/A uses B score (AUROC=0.714286)
- PASS: TEST 2 PIK3CA/mTOR D/B uses A score (AUROC=0.692130)
- PASS: TEST 3 PIK3CA/mTOR dual is positive class
- PASS: TEST 4 PIK3CA/mTOR AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 PIK3CA/mTOR summary_min=min(two directional points) (0.692130)
- PASS: TEST 7 PIK3CA/mTOR Smean=(SA+SB)/2
- PASS: PIK3CA/mTOR D/A(score B)=0.7143 differs from incorrect D/A(score A)=0.7103
- PASS: PIK3CA/mTOR summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 PIK3CA/mTOR pocket A Δ uses identical score channel
- PASS: TEST 8 PIK3CA/mTOR Top10% D+A+B+N=k (4/1/0/0 k=5)
- PASS: TEST 9 PIK3CA/mTOR k=ceil(0.10 n) (n=48 k=5)
- PASS: TEST 10 PIK3CA/mTOR EF formula (2.133333)
- PASS: TEST 18 PIK3CA/mTOR matched/mismatched only swaps assignment (Δsummary_min=0.0903 CI=[-0.117,0.286])
- PASS: PIK3CA/mTOR S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 AChE/BChE D/A uses B score (AUROC=0.650370)
- PASS: TEST 2 AChE/BChE D/B uses A score (AUROC=0.605820)
- PASS: TEST 3 AChE/BChE dual is positive class
- PASS: TEST 4 AChE/BChE AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 AChE/BChE summary_min=min(two directional points) (0.605820)
- PASS: TEST 7 AChE/BChE Smean=(SA+SB)/2
- PASS: AChE/BChE D/A(score B)=0.6504 differs from incorrect D/A(score A)=0.4444
- PASS: AChE/BChE summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 AChE/BChE pocket A Δ uses identical score channel
- PASS: TEST 8 AChE/BChE Top10% D+A+B+N=k (5/3/1/1 k=10)
- PASS: TEST 9 AChE/BChE k=ceil(0.10 n) (n=95 k=10)
- PASS: TEST 10 AChE/BChE EF formula (1.759259)
- PASS: TEST 18 AChE/BChE matched/mismatched only swaps assignment (Δsummary_min=0.1614 CI=[0.041,0.271])
- PASS: AChE/BChE S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 F2/F10 D/A uses B score (AUROC=0.413306)
- PASS: TEST 2 F2/F10 D/B uses A score (AUROC=0.344758)
- PASS: TEST 3 F2/F10 dual is positive class
- PASS: TEST 4 F2/F10 AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 F2/F10 summary_min=min(two directional points) (0.344758)
- PASS: TEST 7 F2/F10 Smean=(SA+SB)/2
- PASS: F2/F10 D/A(score B)=0.4133 differs from incorrect D/A(score A)=0.4506
- PASS: F2/F10 summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 F2/F10 pocket A Δ uses identical score channel
- PASS: TEST 8 F2/F10 Top10% D+A+B+N=k (4/1/6/0 k=11)
- PASS: TEST 9 F2/F10 k=ceil(0.10 n) (n=107 k=11)
- PASS: TEST 10 F2/F10 EF formula (1.255132)
- PASS: TEST 18 F2/F10 matched/mismatched only swaps assignment (Δsummary_min=-0.0307 CI=[-0.115,0.046])
- PASS: F2/F10 S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 PPARG/PPARA D/A uses B score (AUROC=0.649194)
- PASS: TEST 2 PPARG/PPARA D/B uses A score (AUROC=0.706055)
- PASS: TEST 3 PPARG/PPARA dual is positive class
- PASS: TEST 4 PPARG/PPARA AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 PPARG/PPARA summary_min=min(two directional points) (0.649194)
- PASS: TEST 7 PPARG/PPARA Smean=(SA+SB)/2
- PASS: PPARG/PPARA D/A(score B)=0.6492 differs from incorrect D/A(score A)=0.6512
- PASS: PPARG/PPARA summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 PPARG/PPARA pocket A Δ uses identical score channel
- PASS: TEST 8 PPARG/PPARA Top10% D+A+B+N=k (7/3/0/1 k=11)
- PASS: TEST 9 PPARG/PPARA k=ceil(0.10 n) (n=109 k=11)
- PASS: TEST 10 PPARG/PPARA EF formula (2.167614)
- PASS: TEST 18 PPARG/PPARA matched/mismatched only swaps assignment (Δsummary_min=0.0296 CI=[-0.079,0.162])
- PASS: PPARG/PPARA S=-raw_Vina for every complete-case ligand
- PASS: TEST 1 PPARA/PPARD D/A uses B score (AUROC=0.646484)
- PASS: TEST 2 PPARA/PPARD D/B uses A score (AUROC=0.446289)
- PASS: TEST 3 PPARA/PPARD dual is positive class
- PASS: TEST 4 PPARA/PPARD AUROC(y,-score)≈1-AUROC (sum=1.000000000000)
- PASS: TEST 5 PPARA/PPARD summary_min=min(two directional points) (0.446289)
- PASS: TEST 7 PPARA/PPARD Smean=(SA+SB)/2
- PASS: PPARA/PPARD D/A(score B)=0.6465 differs from incorrect D/A(score A)=0.4346
- PASS: PPARA/PPARD summary_min bootstrap: 2000/2000 valid pooled non-stratified replicates
- PASS: TEST 6 PPARA/PPARD pocket A Δ uses identical score channel
- PASS: TEST 8 PPARA/PPARD Top10% D+A+B+N=k (5/2/3/1 k=11)
- PASS: TEST 9 PPARA/PPARD k=ceil(0.10 n) (n=110 k=11)
- PASS: TEST 10 PPARA/PPARD EF formula (1.562500)
- PASS: TEST 18 PPARA/PPARD matched/mismatched only swaps assignment (Δsummary_min=0.0117 CI=[-0.088,0.143])
- PASS: PPARA/PPARD S=-raw_Vina for every complete-case ligand
- PASS: EGFR/HER2 has no unused-pool holdout (supply-limited; expected)
- PASS: JAK1/JAK2 main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: JAK1/TYK2 main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: PIK3CA/mTOR main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: AChE/BChE main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: F2/F10 main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: PPARG/PPARA main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: PPARA/PPARD main vs holdout: zero ligand-id / ChEMBL / exact-SMILES overlap
- PASS: TEST 12 main/holdout ligand overlap = 0
- PASS: TEST 11 ECFP4 train/test scaffold overlap = 0 (leaked_folds=0)
- PASS: incremental max |ΔAUROC| independently retrained = 0.0234 (manuscript ≈ 0.023)
- PASS: EGFR/HER2 descriptor heavy D-vs-A unflipped AUROC matches deposited table
- PASS: EGFR/HER2 descriptor mw D-vs-A unflipped AUROC matches deposited table
- PASS: EGFR/HER2 descriptor clogp D-vs-A unflipped AUROC matches deposited table
- PASS: EGFR/HER2 descriptor tpsa D-vs-A unflipped AUROC matches deposited table
- PASS: AChE/BChE descriptor heavy D-vs-A unflipped AUROC matches deposited table
- PASS: AChE/BChE descriptor mw D-vs-A unflipped AUROC matches deposited table
- PASS: AChE/BChE descriptor clogp D-vs-A unflipped AUROC matches deposited table
- PASS: AChE/BChE descriptor tpsa D-vs-A unflipped AUROC matches deposited table
- PASS: PIK3CA/mTOR descriptor heavy D-vs-A unflipped AUROC matches deposited table
- PASS: PIK3CA/mTOR descriptor mw D-vs-A unflipped AUROC matches deposited table
- PASS: PIK3CA/mTOR descriptor clogp D-vs-A unflipped AUROC matches deposited table
- PASS: PIK3CA/mTOR descriptor tpsa D-vs-A unflipped AUROC matches deposited table
- PASS: physicochemical baselines do not apply max(AUC, 1-AUC)
- PASS: Top-10% totals D/A/B/N = 33/23/22/3 (expected 33/23/22/3)
- PASS: TEST 17 primary Vina score = mode 1 (five-pair scores_vina_mode1_v1.csv; orig-3 ablation *_affinity columns)
- PASS: EGFR/HER2 independent GNINA complete-case n=(28, 38, 32, 11)
- PASS: EGFR/HER2 independent GNINA weaker=0.2199 neither=0.7825 (manuscript 0.220 / 0.783)
- PASS: PIK3CA/mTOR independent GNINA complete-case n=(18, 13, 12, 4)
- PASS: PIK3CA/mTOR independent GNINA weaker=0.6325 neither=0.5694 (manuscript 0.633 / 0.569)
- PASS: JAK1/TYK2 independent GNINA complete-case n=(30, 32, 29, 14)
- PASS: JAK1/TYK2 independent GNINA smin=0.3172 neither=0.7048 (manuscript 0.317 / 0.705)
- PASS: PPARG/PPARA Vina smin=0.6492 RTM=0.3691 CNN-affinity rescoring=0.5000 on identical ligand set n=32/31/32
- PASS: PPARG/PPARA CNN rescoring is Vina-pose CNN affinity, not independent GNINA search
- PASS: TEST 19 receptor substitution and independent GNINA are not mixed in PPARG rescoring (RTM=0.369 CNN=0.500 on Vina poses)
- PASS: PIK3CA/mTOR substitution 4L23/4JT6 summary_min=0.6921
- PASS: PIK3CA/mTOR substitution 4JPS/4JT6 summary_min=0.4861
- PASS: PIK3CA/mTOR substitution 5DXT/4JT6 summary_min=0.5046
- PASS: PIK3CA/mTOR substitution 4L23/4JSX summary_min=0.6389
- PASS: PIK3CA/mTOR crystal-substitution table does not contain independent GNINA 0.633
- PASS: TEST 19 receptor substitution and independent GNINA results are not mixed
- PASS: 4EY7 box matches independent heavy-atom AABB+5Å/min20 (max |Δ|=0.0000 Å)
- PASS: 4BDS box matches independent heavy-atom AABB+5Å/min20 (max |Δ|=0.0000 Å)
- PASS: 4L23 box matches independent heavy-atom AABB+5Å/min20 (max |Δ|=0.0000 Å)
- PASS: 4JT6 box matches independent heavy-atom AABB+5Å/min20 (max |Δ|=0.0000 Å)
- PASS: primary box JSON present for 3POZ
- PASS: primary box JSON present for 3RCD
- PASS: primary box JSON present for 4L23
- PASS: primary box JSON present for 4JT6
- PASS: primary box JSON present for 4EY7
- PASS: primary box JSON present for 4BDS
- PASS: primary box JSON present for 4UDW
- PASS: primary box JSON present for 2JKH
- PASS: primary box JSON present for 6N7A
- PASS: primary box JSON present for 8BXH
- PASS: primary box JSON present for 3LXP
- PASS: primary box JSON present for 9V8H
- PASS: primary box JSON present for 6LXA
- PASS: primary box JSON present for 5U3Q
- PASS: EGFR/HER2 A 3POZ local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: EGFR/HER2 B 3RCD local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: JAK1/JAK2 A 6N7A present in freeze/site-verification tables as the intended human holo
- PASS: JAK1/JAK2 B 8BXH present in freeze/site-verification tables as the intended human holo
- PASS: JAK1/TYK2 A 6N7A present in freeze/site-verification tables as the intended human holo
- PASS: JAK1/TYK2 B 3LXP present in freeze/site-verification tables as the intended human holo
- PASS: PIK3CA/mTOR A 4L23 local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: PIK3CA/mTOR B 4JT6 local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: AChE/BChE A 4EY7 local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: AChE/BChE B 4BDS local identity audit: human, accession matches labels (OK_expected_human_protein)
- PASS: F2/F10 A 4UDW present in freeze/site-verification tables as the intended human holo
- PASS: F2/F10 B 2JKH present in freeze/site-verification tables as the intended human holo
- PASS: PPARG/PPARA A 9V8H present in freeze/site-verification tables as the intended human holo
- PASS: PPARG/PPARA B 6LXA present in freeze/site-verification tables as the intended human holo
- PASS: PPARA/PPARD A 6LXA present in freeze/site-verification tables as the intended human holo
- PASS: PPARA/PPARD B 5U3Q present in freeze/site-verification tables as the intended human holo
- PASS: JAK1 6N7A is the shared receptor for JAK1/JAK2 and JAK1/TYK2
- PASS: PPARA 6LXA is the shared receptor for PPARG/PPARA and PPARA/PPARD
- PASS: Withdrawn PIK3CA/PIK3CB 2WXF mouse p110δ is not in the eight-pair primary set
- PASS: 14/14 receptors have min saved-pose CalcRMS < 2.0 Å (coverage gate, not top-1)
- PASS: top-1 < 2 Å on 8/14 receptors; manuscript does not claim uniform top-pose recovery
- PASS: manuscript distinguishes min saved-pose coverage from top-1 ranking
- PASS: EGFR/HER2 scaffold cap = 5/class; PIK3CA/mTOR = 2/class
- PASS: five Track B main panels: quota + shuffle only, no Murcko cap
- PASS: main-panel builders do not select ligands by docking score or AUROC
- PASS: AChE n_scored 95 vs n_panel 100 is docking-failure complete-case filtering after freeze, not post-hoc deletion of unfavorable scores
- PASS: census n_pairs_n_both_ge_10=63790 matches manuscript (from deposited summary CSV)
- PASS: census n_directional_n10=5253 matches manuscript (from deposited summary CSV)
- PASS: census n_strict_thick=86 matches manuscript (from deposited summary CSV)
- PASS: EGFR/HER2 max→median summary_min 0.430→0.424
- PASS: AChE/BChE max→median summary_min 0.606→0.629
- PASS: PPARA/PPARD max→median summary_min 0.446→0.446
- PASS: five pairs unchanged under median aggregation (class composition and summary_min point)
- PASS: EGFR/HER2 document_cluster Δ CI excludes 0
- PASS: EGFR/HER2 document_cluster deposited note states cluster-level resampling
- PASS: EGFR/HER2 scaffold_cluster Δ CI excludes 0
- PASS: EGFR/HER2 scaffold_cluster deposited note states cluster-level resampling
- PASS: JAK1/TYK2 document-cluster Δ CI [-0.034,0.682] includes 0 (matches Methods/Results)
- PASS: JAK1/TYK2 document_cluster deposited note states cluster-level resampling
- PASS: JAK1/TYK2 scaffold_cluster deposited note states cluster-level resampling
- PASS: all 7 available holdout matched-minus-mismatched 95% CIs include 0
- PASS: main-panel matched-minus-mismatched CIs excluding 0: ['AChE/BChE', 'EGFR/HER2']
- PASS: EGFR/HER2 0.430 → 0.808, Δ=0.378 independently recovered
- PASS: JAK1/TYK2 pocket-A ΔAUROC=0.4438 rounds to 0.444
- PASS: Figure 1 caption no longer contains 'Independent of the census'; current text: eight-pair set is not a direct continuation of the census counts
- PASS: Zhou 2013 is cited as prior evidence that single-target inhibitors can be dual-target docking false positives
- PASS: no 'first to show single-target false positives' claim found
- PASS: Figure 2 plotted_values.json directional AUROCs were previously locked to independent recomputes in audit v1
- PASS: TEST 13 Figure/Table n values correspond to real ligand IDs
- PASS: TEST 14 all manuscript numeric claims have result-table provenance
- PASS: scanned flagship and Table 2/3 numbers match independent 3-decimal rounding
- PASS: TEST 20 reported CI procedures match Methods for summary_min (B=2000 pooled dual+A+B min-inside-replicate)
- PASS: software: numpy 2.4.4; sklearn 1.9.1; rdkit 2026.03.6; git 6f5bea67a7aca5ffc5c9abdbf6f2226790615431

# Metrics audit

## Directional AUROC

| Pair | n (D/A/B) | D vs A-only (B) | D vs B-only (A) | summary_min | sklearn match |
|------|-----------|----------------:|----------------:|------------:|:-------------:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 | yes |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 | yes |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 | yes |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 | yes |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 | yes |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 | yes |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 | yes |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 | yes |

## summary_min bootstrap

| Pair | point | 95% CI | valid/B |
|------|------:|--------|--------:|
| EGFR/HER2 | 0.430 | [0.282, 0.578] | 2000/2000 |
| JAK1/JAK2 | 0.588 | [0.444, 0.725] | 2000/2000 |
| JAK1/TYK2 | 0.365 | [0.231, 0.503] | 2000/2000 |
| PIK3CA/mTOR | 0.692 | [0.470, 0.813] | 2000/2000 |
| AChE/BChE | 0.606 | [0.437, 0.730] | 2000/2000 |
| F2/F10 | 0.345 | [0.211, 0.477] | 2000/2000 |
| PPARG/PPARA | 0.649 | [0.504, 0.751] | 2000/2000 |
| PPARA/PPARD | 0.446 | [0.296, 0.584] | 2000/2000 |

## fixed-score ΔAUROC

| Pair | pocket | D vs selective | D vs neither | Δ |
|------|--------|---------------:|-------------:|--:|
| EGFR/HER2 | B | 0.666 | 0.720 | +0.054 |
| EGFR/HER2 | A | 0.430 | 0.808 | +0.378 |
| JAK1/JAK2 | B | 0.588 | 0.730 | +0.142 |
| JAK1/JAK2 | A | 0.728 | 0.723 | -0.004 |
| JAK1/TYK2 | B | 0.575 | 0.705 | +0.130 |
| JAK1/TYK2 | A | 0.365 | 0.809 | +0.444 |
| PIK3CA/mTOR | B | 0.714 | 0.583 | -0.131 |
| PIK3CA/mTOR | A | 0.692 | 0.472 | -0.220 |
| AChE/BChE | B | 0.650 | 0.709 | +0.058 |
| AChE/BChE | A | 0.606 | 0.590 | -0.016 |
| F2/F10 | B | 0.413 | 0.528 | +0.115 |
| F2/F10 | A | 0.345 | 0.508 | +0.163 |
| PPARG/PPARA | B | 0.649 | 0.644 | -0.005 |
| PPARG/PPARA | A | 0.706 | 0.759 | +0.053 |
| PPARA/PPARD | B | 0.646 | 0.665 | +0.019 |
| PPARA/PPARD | A | 0.446 | 0.484 | +0.038 |

## two-pocket mean and secondary D vs A+B

| Pair | n_neither | D vs neither | D vs A+B | Δ_mean |
|------|----------:|-------------:|---------:|-------:|
| EGFR/HER2 | 12 | 0.756 | 0.516 | +0.240 |
| JAK1/JAK2 | 14 | 0.730 | 0.655 | +0.075 |
| JAK1/TYK2 | 14 | 0.770 | 0.474 | +0.296 |
| PIK3CA/mTOR | 4 | 0.514 | 0.699 | -0.185 |
| AChE/BChE | 15 | 0.649 | 0.559 | +0.090 |
| F2/F10 | 12 | 0.519 | 0.384 | +0.135 |
| PPARG/PPARA | 14 | 0.685 | 0.673 | +0.013 |
| PPARA/PPARD | 14 | 0.565 | 0.512 | +0.052 |

## ligand-only / incremental AUROC

Unflipped physicochemical AUROCs match deposited orig-3 tables. Independently retrained ECFP4 vs ECFP4+docking max |Δ| = 0.0234.

## matched-minus-mismatched

| Pair | matched | mismatched | Δ | CI excludes 0 |
|------|--------:|-----------:|--:|:-------------:|
| EGFR/HER2 | 0.430 | 0.260 | +0.170 | True |
| JAK1/JAK2 | 0.588 | 0.607 | -0.019 | False |
| JAK1/TYK2 | 0.365 | 0.429 | -0.065 | False |
| PIK3CA/mTOR | 0.692 | 0.602 | +0.090 | False |
| AChE/BChE | 0.606 | 0.444 | +0.161 | True |
| F2/F10 | 0.345 | 0.376 | -0.031 | False |
| PPARG/PPARA | 0.649 | 0.620 | +0.030 | False |
| PPARA/PPARD | 0.446 | 0.435 | +0.012 | False |

## Top-10% and EF_dual,10%

| Pair | n | k | D/A/B/N | EF |
|------|--:|--:|---------|---:|
| EGFR/HER2 | 110 | 11 | 1/5/5/0 | 0.357 |
| JAK1/JAK2 | 110 | 11 | 6/5/0/0 | 1.875 |
| JAK1/TYK2 | 109 | 11 | 1/3/7/0 | 0.320 |
| PIK3CA/mTOR | 48 | 5 | 4/1/0/0 | 2.133 |
| AChE/BChE | 95 | 10 | 5/3/1/1 | 1.759 |
| F2/F10 | 107 | 11 | 4/1/6/0 | 1.255 |
| PPARG/PPARA | 109 | 11 | 7/3/0/1 | 2.168 |
| PPARA/PPARD | 110 | 11 | 5/2/3/1 | 1.562 |

# Data leakage audit

- Panel construction: activity / supply / quota / scaffold or ID diversity; no docking-score or AUROC selection of main-panel ligands.
- Primary labels: θ=6.0 from pChEMBL; missing measurements are not complete-case members.
- ECFP4 GroupKFold: overlap_count=0 on every fold (see cv_leakage_audit.csv).
- Incremental docking uses the same splits as ECFP4-only.
- Main vs holdout: zero ligand-id / ChEMBL / exact SMILES overlap; Murcko overlap exists and is leftover-pool, not a CV leak.
- Matched/mismatched does not redock.

# Structural audit

- 14 primary receptors: local metadata support human expected proteins; JAK1 6N7A and PPARA 6LXA correctly shared.
- Cognate coverage: 14/14 min saved-pose RMSD < 2 Å; several top-1 failures, consistent with Methods.
- Docking boxes: AChE and PIK3CA match heavy-atom AABB+5/min20. EGFR 3POZ/HER2 3RCD JSON atoms=63 include hydrogens and do not match the written heavy-atom rule (MAJOR).
- Independent GNINA (pose generation) is separate from Vina-pose RTM/CNN rescoring.
- PIK3CA/mTOR crystal substitution 0.692 / 0.486 / 0.505 / 0.639 is not mixed with independent GNINA 0.633 in the substitution CSV.

# Figure/table/manuscript consistency audit

- manuscript_number_trace.csv rows: 35; MISMATCH: 0
- Figure 1 no longer claims independence from the census.
- Figure 5C is summary_min across Vina seeds in the current caption, not a fixed-score Δ plot.
- Title/abstract remain a methodological evaluation of docking for dual-target candidate ranking, not a new engine or prospective screen.

# Reproducibility audit

- git `6f5bea67a7aca5ffc5c9abdbf6f2226790615431`
- numpy 2.4.4, scikit-learn 1.9.1, RDKit 2026.03.6
- primary scores: mode-1 Vina, S=−E
- bootstrap seed base 20260729 with SHA256 offsets as in Methods
- this audit writes only audit_outputs/
- ChEMBL 37 SQLite dump is not in the repository; census n_both≥1 is therefore not fully verifiable

# Answers to the 20 core questions

1. Four-state labels: YES, independently rebuilt as θ=6.0 from pChEMBL; missing values are not treated as inactive; exact 6.0 is active.
2. Directional score assignment: YES for all 16 arms (D/A uses B; D/B uses A). Incorrect A-only→score A was tested and differs.
3. Vina sign: YES, S=−raw affinity throughout.
4. Primary score mode 1: YES for Track B via scores_vina_mode1_v1.csv; orig-3 via ablation affinity columns matching membership S=−E.
5. summary_min and bootstrap: YES, min of two point estimates; CI from min-inside-replicate pooled bootstrap B=2000.
6. Fixed-score ΔAUROC: YES, same pocket channel for selective and neither controls; EGFR/HER2 and JAK1/TYK2 flagship values recovered.
7. Two-pocket mean: YES, complete cases, S_mean=(S_A+S_B)/2; PIK3CA/mTOR neither n=4 confirmed. Secondary D-vs-A+B is reported and does not replace primary directional analysis.
8. ECFP4 leakage: NO scaffold train/test overlap in independent GroupKFold. Four acyclic PPAR ligands use singleton groups.
9. Incremental AUROC: YES, same ligands, labels, and folds; max |Δ|≈0.023.
10. Matched-minus-mismatched: YES, assignment swap only; paired bootstrap.
11. Main/holdout ligand leakage: NO ligand-id/ChEMBL/SMILES overlap. Scaffold overlap exists (unused pool). EGFR has no holdout.
12. Receptor identity: PASS on local human/accession tables for the eight primary pairs. Live RCSB was not re-queried here. Shared 6N7A and 6LXA are correct.
13. Redocking: YES, 2.0 Å applied to min saved-pose RMSD; 14/14 pass; not claimed as uniform top-1 recovery.
14. GNINA rescoring vs independent docking: YES, distinguished in code and recovered numbers.
15. Receptor substitution vs GNINA mix-up: NO mix in the substitution CSV; 0.633 is independent GNINA weaker-arm.
16. Threshold/aggregation sensitivity: deposited max-vs-median table matches the manuscript 6/110, 1/95, 1/110 story; independent relabel from raw ChEMBL records was not fully re-harvested (SQLite dump absent) — WARNING for full record-level rebuild, PASS vs deposited freeze.
17. Top-10% composition: YES vs expected 1/5/5/0 … totals 33/23/22/3.
18. EF_dual,10%: YES, denominator is D_total/n on the full four-state panel.
19. Figure/Table/Abstract consistency: flagship and Table 2/3 match. Remaining inconsistencies are protocol wording (EGFR box), unverifiable census 2,164,618, hard-coded figure literals, and AChE ID-prefix cap.
20. Main conclusions supported: YES. Discrimination varies by pair and control class; dual-vs-neither does not reliably imply exclusion of single-target actives; ligand chemistry carries class information; docking incremental value is small in the ECFP4 tests; matched-pocket advantage is not consistent; enrichment is pair-dependent; external gates were not met.

# Final judgment

**B. CORE RESULTS VERIFIED WITH MINOR CORRECTIONS**

Independent ligand-level recomputation recovered the primary directional AUROCs, summary_min bootstrap construction, fixed-score ΔAUROC flagship values, two-pocket mean, Top-10%/EF, ECFP4 same-fold incremental analysis, and matched/mismatched assignment swap. No A/B score reversal, θ=6.0 label error, mode-1 substitution, ECFP4 CV leakage, or wrong-species primary receptor was found. The EGFR/HER2 deposited boxes include hydrogens and therefore do not match the written heavy-atom box rule; that is a protocol/SI defect whose effect on docking scores was not re-estimated. The 2,164,618 census count cannot be rebuilt from files in this repository. Those items warrant correction in Methods/SI but do not, on present evidence, reverse the paper's methodological conclusions.

Based on the independent audit, the manuscript’s core conclusions are verified with corrections.

