# END_TO_END_DATAFLOW_AUDIT

Read-only map. No correctness judgment in this file except listing conversion risks.

## Channel 1 — primary 14-receptor method ablation

raw receptor PDB (`01_receptors/{pdb}_protein_raw.pdb`)
→ PDBFixer (`*_protein_fixer*.pdb`) → optional PTR→TYR (`*_protein_ptr_fallback.pdb`)
→ prepared PDB (`*_protein_prepared.pdb`) → Meeko (`*_receptor.pdbqt`)
script: `scripts/rerun_v4_2/phase2_prepare_receptors.py`
software: PDBFixer 1.12 / OpenMM 8.6.1 / Meeko 0.7.1
may change: hydrogens (Meeko), atom types, charges; PTR→TYR only 6N7A/8BXH/3LXP
must not change: deposited sequence except frozen PTR fallback; retained/removed components

cognate crystal JSON + AABB pad 5 Å min edge 20 Å → `03_boxes/{pdb}_box.json`
script: historical V4.1 box reused / Phase 2 does not rewrite boxes
may change: none if frozen; risk if alt box reused primary coordinates

frozen ligand SDF `02_ligands/sdf/{eid}.sdf` (reused from V4.1, 785)
→ Meeko PDBQT `02_ligands/pdbqt/{eid}.pdbqt`
may change: AD4 types, Gasteiger charges, explicit H, macrocycle CGn/Gn glue (CGn=real atom, Gn=glue)

Vina five-seed `06_vina_fiveseed/jobs/{pair}__{pdb}__{eid}__seed{s}/`
script: `phase7_fiveseed_production.py` vina 1.2.7 exh 16 modes 9 energy_range 6 cpu 1
pair-expanded seed records: 8070
unique physical job-seed: 7960
may change: coordinates (search); must not change ligand graph

Vina pose → old GNINA M1 `09_secondary_scoring/gnina_rescore/`
script: `ablation_phase4_gnina_rescore.py` score_only no_gpu scoring vina cnn_scoring rescore
risk: CG0 not valid AutoDock type (historical)

Vina pose → M2 `09_secondary_scoring/rtmscore/`
script: `ablation_phase578_rtmscore.py` frozen SDF graph + SMILES IDX coords, model1, 10 Å pockets

frozen/converted ligand → old M3 `10_gnina_independent_docking/`
script: `ablation_phase10_12_gnina.py` CG0→C drop G0; seed 42 exh 16 modes 9
risk: ring-closure bond not restored (CGn/Gn topology)

## Channel 2 — alternative receptor experiment (separate freeze)

authority: `ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml`
8 alts: 4L2Y 4JT5 4EY6 1P0M 3SHC 2Y5F 6KAX 5U46
own prep + own box + alt redock `07_alt_cognate_redocking` + alt production `08_alt_vina_seed42`
A2/B1 replaces A only; A1/B2 replaces B only; kept side reuses primary score
NOT the same channel as analysis_freeze allowed_pdb 4JPS/5DXT/4JSX

## Channel 3 — V2 representation repair (NOT final evidence in this audit)

`gnina_rescore_representation_v2/` and `jobs_topology_repair_v2/` are active writes.
This audit does not treat them as official M1/M3.
