# Pose compare verdict — EH40_23 vs TAK-285 / lapatinib

## Modes (RTM best)
| ligand | 3POZ mode | 3RCD mode |
|--------|-----------|-----------|
| EH40_01 TAK-285 | 2 | 3 |
| EH40_02 LAPATINIB | 5 | 1 |
| EH40_23 | 6 | 2 |

## Does EH40_23 show hinge + backpocket on both ends?
Yes (distance proxies):
target  hinge_hbond  hinge_min_N_to_bb_A hydrophobic_tail_in_backpocket cys_covalent_site_nearby  cys_min_A clash_flag
  3POZ yes (MET793)                3.056                            yes                      yes      3.699         no
  3RCD yes (MET801)                3.300                            yes                       no      5.552         no
- Both ends: hinge H-bond proxy to MET793/MET801 bb **yes** (~3.05–3.30 Å)
- Hydrophobic tail / far-from-hinge atoms present
- Clash gate: **no** clashes <2.2 Å (explains why geometric gate failed)

## More similar to which anchor?
- **Chemistry (ECFP4 Tanimoto)**: vs EH40_02 lapatinib **0.5287** >> vs EH40_01 **0.3232**
- **Pose MCS RMSD** (shared anilinoquinazoline-like scaffold): vs lapatinib **~0.10 Å (3POZ) / 0.25 Å (3RCD)**; vs TAK-285 ~0.27 / 0.33 Å
- **Pharmacophore (hingeN/Cl/F-aryl)**: nearly identical to lapatinib (pharm RMSD ~0.14–0.21 Å)

→ Binding-mode and chemotype both closer to **lapatinib** than TAK-285.

## Chemotype-driven hard negative?
**Supported.** EH40_23 is an anilinoquinazoline TKI homolog of lapatinib that adopts the same classical type-I hinge+backpocket pose on EGFR and HER2 in noncovalent docking/RTM, while panel biology remains A_only (HER2 pChEMBL<6). High dual scores are therefore expected from shared kinase chemotype × shared pocket geometry, not from validated dual pharmacology.
