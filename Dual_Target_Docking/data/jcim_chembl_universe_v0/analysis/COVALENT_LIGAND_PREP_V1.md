# Covalent-capable cognates — Vina ligand prep (4X6H, 9GJ2)

Applies only to **CTSK/CTSS**. The other eight new receptors are ordinary noncovalent holos.

The benchmark uses one noncovalent engine (AutoDock Vina). That is a **pre-reaction noncovalent approximation**, declared here, not a claim that the crystal complexes are noncovalent.

## Rule

| Step | 4X6H (CTSK) | 9GJ2 (CTSS) |
|------|-------------|-------------|
| Receptor PDB | keep **4X6H** | keep **9GJ2** |
| Box | heavy-atom AABB of the crystallographic occupant (**3XT** or **I37**) + 5 Å/axis, min edge 20 Å | heavy-atom AABB of **KH0** + 5 Å/axis, min edge 20 Å |
| Free ligand for Meeko / Vina | **I37** (pre-reaction nitrile), never 3XT | reconstructed **pre-reaction α-ketoamide 13b**, never KH0 as extracted |
| Layer-3 pose-gold | **I37** coordinates in 4X6H | **no** valid noncovalent CCD gold; do not report RMSD vs KH0 as cognate recovery |

If a ligand is cut from the covalent complex, the warhead does **not** snap back to the free electrophile. Scoring that adduct as a free ligand mixes chemical mechanisms into the same Vina table as thrombin/JAK/PPAR.

## 4X6H chemistry (RCSB CCD, 2026-09-04)

Entry contains both states, which matches Boríšek et al.: reversible tight-binding glycinonitriles with covalent and noncovalent Cys25 occupancy.

| CCD | Role | Name (RCSB) | Formula | MW |
|-----|------|-------------|---------|----:|
| **I37** | pre-reaction nitrile (dock this) | 4-amino-N-{1-[(cyanomethyl)carbamoyl]cyclohexyl}-3-fluorobenzamide | C16 H19 F N4 O2 | 318.346 |
| **3XT** | reacted 2-iminoethyl form (do not dock) | 4-amino-3-fluoro-N-(1-{[(2Z)-2-iminoethyl]carbamoyl}cyclohexyl)benzamide | C16 H21 F N4 O2 | 320.362 |

Paper: Boríšek J et al. *J Med Chem* 2015;58:6928–6937. DOI 10.1021/acs.jmedchem.5b00746.

## 9GJ2 chemistry

- Protein/ligand identity: human cathepsin S (P25774), KH0 = ketoamide **13b**, 1.15 Å. **Keep the receptor.**
- KH0 formula on RCSB: C31 H41 N5 O7 (MW 595.687) — bound representation, not the free α-ketoamide.
- Mechanism: α-ketoamide 13b is a reversible-covalent warhead (thiohemiketal at catalytic Cys). The same 13b forms a covalent thiohemiketal with CatL Cys25 (Steigerwald et al. analog; do not treat 9GJ2 as ordinary noncovalent).
- Citation: RCSB primary citation is still **To Be Published** (authors Falke S., Karnicar K., Usenik A., Lindic N., Sekirnik A., Reinke P.Y.A., Guenther S., Turk D., Meents A.). **No DOI.** Cite `PDB 9GJ2`, not a journal.

Reconstruct 13b from the deposition / CCD connectivity by restoring the α-ketoamide carbonyl (pre-reaction), not by deleting the Cys–ligand bond and keeping the thiohemiketal carbons.

## What to write in Methods (exact claim)

> 4X6H and 9GJ2 were retained as receptor structures because the crystallographic ligands occupy the catalytic papain-fold pockets. Both inhibitors are reversible-covalent (nitrile / α-ketoamide). Docking used the corresponding pre-reaction noncovalent ligand representations (4X6H: CCD I37; 9GJ2: reconstructed α-ketoamide 13b) rather than the covalently reacted CCD (3XT / KH0). Vina scores on this pair are a noncovalent approximation and are not chemically equivalent to the thrombin, JAK, or PPAR rows.

## Forbidden

- Extract 3XT or KH0, protonate, Meeko, and call that the cognate.
- Fail or replace 4X6H / 9GJ2 solely because they are covalent-capable.
- Report Layer-3 RMSD of docked I37 against 3XT, or of reconstructed 13b against KH0, as ordinary cognate recovery.
- Cite a journal article for 9GJ2 while RCSB still says To Be Published.
