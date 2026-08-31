# PLIF / pocket-contact verdict v1

Method: geometric heavy-atom proximity ≤ 4.5 Å to the 20 frozen
PIK3CA pocket residues (ProLIF/RDKit conversion segfaults on these PDBs;
this is the SOP-allowed equivalent occupancy snapshot).

Top occupancy-shift residues (4JPS/5DXT vs 4L23): Met772, Leu807, Gln859, Thr856, Cys838, Glu849, Phe930, Asp933

Allowed claim:
> The performance shift coincided with altered interaction patterns at residues
> Met772, Leu807, Gln859, Thr856, Cys838, Glu849, Phe930, Asp933, providing a structural **hypothesis** for the receptor sensitivity.

Forbidden: residue X caused the AUROC change; PLIF explains the opposite
PIK3CA/PIK3CB shift.
