# EH40_23 CASE REPORT

## 1. Identity & panel label
- **EH40_23** / CHEMBL3633938 / class **A_only**
- Panel pChEMBL: EGFR **6.35**, HER2 **5.57** (threshold 6.0)
- SMILES: `N#CC(C#N)=CNc1ccc2ncnc(Nc3ccc(OCc4cccc(F)c4)c(Cl)c3)c2c1`
- Architecture in freeze: unknown (not reassigned here)

## 2. ChEMBL audit → **label_ok**
- Human EGFR (CHEMBL203) / HER2 (CHEMBL1824) API pull: see `tables/eh40_23_chembl_activities.csv`, `tables/label_audit.md`
- API sparse (1 activity/end in this pull) but pchembl_max pattern consistent with A_only when both ends present; **not** reclassifying as dual
- Anchors TAK-285 / lapatinib: dual labels remain consistent with panel + API max≥6 both ends (see label_audit.md)

## 3. Warhead flags
| ligand | quinazoline | malononitrile_like | acrylamide | warhead_risk |
|--------|-------------|--------------------|------------|--------------|
| EH40_23 | True | True | False | **high** |
| EH40_01 | False | False | False | low |
| EH40_02 | True | False | False | low |

Note: dicyano/malononitrile-like warhead; electrophile risk; noncovalent docking does not model covalent addition | CYS nearby: 3POZ mode6 cys_min=3.699A; 非共价对接可能高估/语义错位；协议未建模共价

## 4. Pose vs TAK-285 / lapatinib (both ends)
- RTM-best modes: 23→3POZ **mode_06**, 3RCD **mode_02**; 01→2/3; 02→5/1
- Hinge H-bond proxy: **yes on both ends** (MET793 / MET801)
- Pose closer to **lapatinib** (MCS RMSD ~0.1–0.25 Å; chem Tanimoto 0.53) than TAK-285 (0.32 chem; MCS ~0.27–0.33 Å)
- Details: `tables/interaction_summary.csv`, `pose_similarity.csv`, `pose_compare_verdict.md`
- PDBs: `poses_pdb/` (RDKit export; Open Babel not available)

## 5. Why clash gate failed
- Expected: both ends are **geometrically clean** classical type-I poses (n_clash_lt_2.2Å = 0; hinge satisfied)
- Clash/PB-like gate cannot down-rank a **chemotype homolog** that docks like lapatinib

## 6. Hard-negative class
- **Primary: type_II_chemotype_homolog**
- Secondary note: type_IV_covalent_mismatch (secondary: malononitrile-like + EGFR CYS797 nearby on 3POZ)
- Not type_I_score_artifact (unlike EH40_18 which RTM suppressed); not type_III_label_noise given audit label_ok

## 7. Implication for protocol v0.1
- Geometric clash gate: **retain as negative result** for this failure mode (do not retune thresholds to force-drop EH40_23)
- Optional **chemotype warning** (anilinoquinazoline / lapatinib-like TKI + optional malononitrile flag) as a separate diagnostic layer — **suggestion only; do not rescore/brush ranks in v0 freeze**

## 8. Next step recommendation
Keep EH40_23 as the canonical **type_II chemotype-homolog hard negative** in the EGFR/HER2 diagnostic panel, and design v0.1 evaluation around reporting “warning flags” separately from gated scores rather than tightening clash cutoffs.
