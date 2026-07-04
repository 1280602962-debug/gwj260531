# Extended MD Analysis Summary — Compound 2231 vs JNK1/2/3

## Simulation context

- Single-replica 200 ns MD per isoform (JNK1/3ELJ, JNK2/3E7O, JNK3/3TTI).
- Production statistics use frames from **50–200 ns** (first 50 ns discarded as equilibration).
- Ligand heavy atoms were subject to a **2.0 kcal/mol/Å² positional restraint on ligand heavy atoms** during production; ligand RMSD values therefore reflect restrained motion relative to the starting pose, not unrestrained binding-mode sampling.
- MM-GBSA (MMPBSA.py, igb=8, frames 15001–20000) provides **relative ranking** across isoforms within this setup; values are **not** absolute binding affinities.

## 1. Production-window RMSD statistics

Full table: `09_production_rmsd_percentiles.csv`. Figure: `10_ligand_rmsd_violin.png`.

### JNK1 (3ELJ)
- **Protein RMSD**: mean 3.32 Å, median 3.30 Å, p5–p95 = 2.82–3.85 Å
- **Ligand RMSD**: mean 0.63 Å, median 0.57 Å, p5–p95 = 0.44–1.09 Å
- **Complex RMSD**: mean 11.14 Å, median 12.58 Å, p5–p95 = 6.74–12.81 Å

### JNK2 (3E7O)
- **Protein RMSD**: mean 2.58 Å, median 2.56 Å, p5–p95 = 1.98–3.31 Å
- **Ligand RMSD**: mean 1.77 Å, median 1.74 Å, p5–p95 = 1.20–2.26 Å
- **Complex RMSD**: mean 3.24 Å, median 3.20 Å, p5–p95 = 2.76–3.84 Å

### JNK3 (3TTI)
- **Protein RMSD**: mean 3.75 Å, median 3.80 Å, p5–p95 = 2.37–5.13 Å
- **Ligand RMSD**: mean 1.10 Å, median 1.08 Å, p5–p95 = 0.77–1.46 Å
- **Complex RMSD**: mean 3.96 Å, median 4.00 Å, p5–p95 = 2.76–5.14 Å

## 2. Ligand vs protein RMSD correlation

- **JNK1**: Pearson r = -0.232 (p = 3.22e-182, n = 15000)
- **JNK2**: Pearson r = 0.378 (p = 0.00e+00, n = 15000)
- **JNK3**: Pearson r = 0.330 (p = 0.00e+00, n = 15000)

Scatter plots: `11_rmsd_correlation.png`.

## 3. Flexible residues (RMSF > 3 Å)

- **JNK1** (13 residues): 167 (12.71 Å), 168 (11.31 Å), 169 (9.53 Å), 327 (7.58 Å), 170 (6.89 Å), 171 (5.91 Å), 328 (5.21 Å), 172 (4.42 Å)
- **JNK2** (4 residues): 274 (4.20 Å), 275 (3.55 Å), 27 (3.06 Å), 328 (3.00 Å)
- **JNK3** (16 residues): 323 (15.27 Å), 322 (13.66 Å), 321 (12.17 Å), 320 (10.72 Å), 319 (9.55 Å), 167 (9.28 Å), 168 (7.84 Å), 318 (7.49 Å)

Full list: `12_rmsf_above_3A.csv`. Plot: `12_rmsf_highlights.png`.

## 4. MM-GBSA per-residue decomposition (top 5 favorable protein residues)

### JNK1 — ligand total ΔG = -12.48 kcal/mol
- ASN 108: -3.46 kcal/mol
- LEU 104: -2.35 kcal/mol
- ILE  26: -2.15 kcal/mol
- VAL 152: -1.86 kcal/mol
- LEU 162: -1.86 kcal/mol

### JNK2 — ligand total ΔG = -8.07 kcal/mol
- VAL 110: -1.80 kcal/mol
- MET 113: -1.72 kcal/mol
- VAL 151: -1.10 kcal/mol
- SER 153: -1.01 kcal/mol
- ALA 105: -0.69 kcal/mol

### JNK3 — ligand total ΔG = -7.69 kcal/mol
- ILE  25: -2.70 kcal/mol
- ASN 107: -1.40 kcal/mol
- ALA 106: -1.15 kcal/mol
- VAL  33: -0.87 kcal/mol
- GLN 110: -0.76 kcal/mol

Figure: `13_decomp_top5.png`.

## 5. Persistent H-bonds (top 3 by occupancy)

### JNK1
- MOL_348@O2 ↔ ASN_108@ND2: 68.4%
- MOL_348@N4 ↔ MET_105@N: 6.8%
- ILE_26@HD11 ↔ MOL_348@N2: 4.5%

### JNK2
- SER_153@OG ↔ MOL_343@N2: 7.8%
- MOL_343@O2 ↔ ASN_106@ND2: 7.1%
- SER_147@O ↔ MOL_343@N2: 2.3%

### JNK3
- MOL_342@O2 ↔ GLN_30@NE2: 18.5%
- GLN_110@OE1 ↔ MOL_342@N2: 7.1%
- ALA_46@HB3 ↔ MOL_342@N2: 1.3%

## 6. MM-GBSA energy components

All systems reported an **internal potential inconsistency warning** in `gbsa.dat`; interpret component magnitudes cautiously.

- **JNK1**: VDW -44.01, EEL -10.47, EGB 28.23, ESURF -5.65, **DELTA TOTAL -31.90** kcal/mol
- **JNK2**: VDW -28.10, EEL -6.93, EGB 22.00, ESURF -3.79, **DELTA TOTAL -16.84** kcal/mol
- **JNK3**: VDW -29.28, EEL -8.19, EGB 25.23, ESURF -3.90, **DELTA TOTAL -16.16** kcal/mol

Table: `14_mm_gbsa_components.csv`. Figure: `14_mm_gbsa_components.png`.

## Interpretation notes

- Differences in ligand RMSD across isoforms may partly reflect dissimilar starting poses and protein backbone fluctuations rather than isoform selectivity.
- A single trajectory per target does not support statistical selectivity claims; replicates and unrestrained ligand simulations would be needed for stronger conclusions.
- MM-GBSA ranking here is best viewed as a **hypothesis-generating** comparison consistent with the docking pose, not experimental validation.
