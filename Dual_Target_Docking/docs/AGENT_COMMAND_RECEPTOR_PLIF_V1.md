# Local follow-up: residue-level interaction snapshot for receptor swaps

Optional **P1** support for Results 3.4 / 3.6. Not required to establish that
receptor realization changes `summary_min`. Do **not** run MD.

Git currently has **no full main-panel pose archive**, so this analysis must
use the local pose workspace (`/mnt/d/...` or the original Vina `_out.pdbqt`
directories), not cloud-invented contacts.

## Question

What interaction-pattern changes **accompany** the PIK3CA swap

- PIK3CA/mTOR: 4L23 → 4JPS / 5DXT (`summary_min` 0.692 → 0.486 / 0.505)
- PIK3CA/PIK3CB: the same crystals **raise** `summary_min` (0.500 → 0.691 / 0.685)

## Scope

Pocket A only (PIK3CA). Frozen pocket B scores stay unused.

Use the 20 PIK3CA pocket residues already listed in
`data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/POCKET_MECHANISM_VERDICT_V1.md`
(Met772, Trp780, Ile800, Lys802, Leu807, Asp810, Leu814, Tyr836, Cys838,
Ile848, Glu849, Val850, Val851, Ser854, Thr856, Gln859, Met922, Phe930,
Ile932, Asp933). Report **5–10** residues with the largest occupancy shift;
do not fingerprint the whole protein.

## Method (ProLIF or equivalent)

For each of 4L23 / 4JPS / 5DXT, on the **same** PM48 ligands that entered
Table S30:

- H-bond, hydrophobic contact, aromatic/π, hinge-adjacent contacts
- Binary 0/1 per ligand–residue (or occupancy %)
- Heatmap: rows = representative dual / A-only / B-only ligands;
  columns = selected residues; panels = receptor

`contact_count` (Table S11) is **not** a substitute: it is a scoring-independent
heavy-atom count on wrong-pocket holdout poses.

## Claim freeze

Allowed:

> The performance shift coincided with altered interaction patterns at residues
> X/Y/Z, providing a structural **hypothesis** for the receptor sensitivity.

Forbidden:

> Residue X caused the AUROC change.
> PLIF explains the opposite PIK3CA/PIK3CB shift.

The opposite-direction pair already shows that a PIK3CA-only residue story
cannot be a complete mechanism.
