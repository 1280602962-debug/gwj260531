# RMSD definition — cognate redock v0

- Reference: `tables/4L23_cocrystal_X6K.pdb`, `tables/4JT6_cocrystal_X6K.pdb`
- Atoms: heavy atoms only
- Pose atom map: meeko `REMARK SMILES IDX` pairs `(smiles_atom_1based, pdbqt_atom_1based)` onto PI-103 SMILES template
- Symmetry: min CalcRMS over template graph automorphisms
- No superposition (docking/prepared-complex frame)
- Gate: both ends best_of_9 < 2.0 Å
