# Historical pre-audit EGFR/HER2 docking boxes

These JSON/TXT files are the **deposited pre-remediation** 3POZ and 3RCD boxes.

They record `n_ligand_atoms = 63` (38 heavy atoms + 25 hydrogens) and match an
all-atom (hydrogen-inclusive) AABB more closely than the Methods heavy-atom rule.

**historical pre-audit boxes; not used in post-remediation primary analysis.**

Corrected heavy-atom boxes (AABB of cognate heavy atoms + 5 Å/side, min edge 20 Å)
are computed in `../phase1_boxes/` and compared in:

- `../egfr_3POZ_box_old_vs_corrected.csv`
- `../her2_3RCD_box_old_vs_corrected.csv`

Do not delete these files. They document the boxes that generated the pre-fix
Table 2 EGFR/HER2 scores.
