# EGFR/HER2 cognate pose recovery (local)

Cloud inventory: `data/jcim_novelty_v0/tables/cognate_artifact_inventory_v1.csv`.

Git currently has EGFR/HER2 receptors (`3POZ_receptor.pdbqt`, `3RCD_receptor.pdbqt`) and a crystal PDB named `3POZ_cocrystal_03P.pdb`. It does **not** have the nine-mode cognate Vina outputs needed to recompute top-1 / top-3 / best-of-9 RMSD the way 4EY7, 4BDS, and 2WXF were recomputed.

## Search order on the local calculation disk

Typical local root from the README:

`/mnt/d/CADD paper exercise/dual target docking/results/`

Look for:

- Vina 9-mode cognate PDBQT for 3POZ (ligand 03P) and 3RCD (ligand TAK-285)
- crystal ligand coordinates and atom mapping
- receptor/box JSON used for the production gate
- logs with exhaustiveness and seed 20260727

If found, copy into `data/egfr_her2_panel40_v0/cognate_qc/` (or panel120 equivalent), add the pair to `cognate_rank_qc_v1.py`, and recompute ranked RMSD.

## Status (2026-08-26)

Historical nine-mode production PDBQTs were not recovered. Reconstructed QC is deposited under `data/egfr_her2_panel40_v0/cognate_qc/` and re-audited with `cognate_rank_qc_v1.py`. Manuscript Table S3 uses the topology-checked CalcRMS values (3POZ 9.505 / 6.227 / 0.760 Å; 3RCD 1.855 / 1.394 / 1.394 Å). Do not replace reconstructed QC with a later redock unless the SI notes the rebuild.

## If the files are gone

Re-redock under the frozen protocol (Meeko receptor, cognate box, Vina 1.2.7, nine poses, seed 20260727, exhaustiveness 8). Label the result **reconstructed QC**, not an original production artifact. Do not replace the historical pass/fail sentence unless the reconstructed RMSD is deposited and the SI says it was rebuilt.
