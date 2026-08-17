# Wrong-pocket mechanism verdict v1 — why the holdout wrong-pocket control is not worse than pocket-matched

> Script: `scripts/wrong_pocket_contact_v1.py`; raw output: `wrong_pocket_contact_v1_output.txt`.
> Inputs: already-docked holdout mode-1 poses (`HOAB/HOAP/HOPM/poses/*/*/out.pdbqt`) and frozen receptor
> PDBQTs. **No new docking; no rescoring.**
> Question: `HOLDOUT_VERDICT.md` showed `wrong_pocket_control_vina` ≥ `pocket_matched_vina` on all three
> holdout pairs (e.g. PM: 0.788 vs 0.765; AChE/BChE: 0.643 vs 0.618; PIK3CB: 0.520 vs 0.425). Is this a
> Vina-scoring artifact, or does it show up in the raw docked geometry, independent of the scoring function?

## Method

`wrong_pocket_control_vina` compares dual vs A_only using the **pocket-A** Vina score (the pocket both
classes are, by label, potent at) instead of pocket-B, and symmetrically for dual vs B_only using
pocket-B. If dual ligands are simply larger/more "generically dockable" molecules, this same-pocket
comparison would show separation from Vina score alone even with no pocket-specific selectivity signal.

To test this independent of the scoring function, we computed a crude, scoring-free geometric proxy
directly from the same mode-1 docked poses already committed to the repo:

- **contact_count** = number of ligand heavy atoms in the mode-1 pose with at least one receptor heavy
  atom within 4.0 Å (a coarse burial/steric-contact count, not a validated PLIF).
- Computed separately for pocket A and pocket B, for every holdout ligand in HOAB, HOPM, and HOAP
  (every ligand is already docked into both pockets to build the four-class panel, so no new docking was
  needed).
- AUROC of `contact_count` alone (dual vs A_only in pocket A; dual vs B_only in pocket B) mirrors exactly
  the `wrong_pocket_control_vina` comparison, but uses only 3D geometry, not the Vina energy function.

## Results

Mean heavy-atom count by class (pooled across both pockets, all three pairs combined trend is consistent
within each pair; see per-pair numbers in the script output):

| pair | class | mean n_heavy |
|---|---|---:|
| AChE/BChE (HOAB) | dual | 34.8 |
| AChE/BChE (HOAB) | A_only | 33.8 |
| AChE/BChE (HOAB) | B_only | 29.6 |

Own-pocket geometric contact-count AUROC (mirrors `wrong_pocket_control_vina`; scoring-free):

| prefix | pair | D vs A_only, pocket A contact_count AUROC | D vs B_only, pocket B contact_count AUROC |
|---|---|---:|---:|
| HOAB | AChE/BChE | **0.581** | **0.706** |
| HOPM | PIK3CA/mTOR | **0.552** | **0.698** |
| HOAP | PIK3CA/PIK3CB | **0.622** | **0.714** |

Mean contact_count by pocket × class on HOAB (representative; full numbers for all three pairs in the raw
output file):

| pocket | class | n | mean contact_count |
|---|---|---:|---:|
| A (AChE) | dual | 20 | 32.15 |
| A (AChE) | A_only | 20 | 30.10 |
| B (BChE) | dual | 20 | 28.10 |
| B (BChE) | B_only | 20 | 22.90 |

## Interpretation

1. **The wrong-pocket-control-is-not-worse pattern reproduces at the level of raw docked geometry, without
   any scoring function.** A simple heavy-atom contact count computed directly from the mode-1 poses
   separates dual from A_only/B_only in their *own* pocket with AUROC 0.55–0.71 on **all three** holdout
   pairs — consistent in direction and magnitude with the Vina-score-based `wrong_pocket_control_vina`
   results. This is stronger evidence than the 2D heavy-atom/TPSA covariate analysis alone (§3.4), because
   it is measured directly on the 3D docked pose, not on the free ligand.
2. **This is consistent with, and mechanistically grounds, the confounding narrative already in §3.4**:
   dual ligands in these panels are on average larger (more heavy atoms) and simply bury more surface area
   in whichever pocket they are docked into, inflating both the "correct" and the "wrong" pocket score in
   the same direction. The pocket-matched signal reported in the main tables is therefore a mixture of a
   real directional component and this ligand-size/burial component; the wrong-pocket control isolates the
   latter and shows it is non-trivial on its own.
3. **This does not mean docking carries zero pocket-specific information in general.** Using the identical
   `wrong_pocket_control_vina` definition, the frozen main panels show pocket-matched summary_min clearly
   above wrong-pocket summary_min for all four pairs (Supporting Information Table S6; e.g. AChE/BChE
   0.606 vs 0.444, PM 0.692 vs 0.602). On the unused-pool holdout the gap reverses (wrong-pocket ≥
   matched). We do not have a resolved explanation for this main-panel-vs-holdout contrast (candidates:
   smaller holdout n=20/class giving noisier point estimates; curated main panels controlling chemotype
   diversity in a way the raw unused pool does not); it is reported as an open discrepancy, not smoothed
   over. What the contact-count analysis does establish is that the size/burial confound is real, large,
   and reproducible at the pose level, and is sufficient by itself to produce the holdout's wrong-pocket
   pattern without invoking scoring-function noise.

## Claim implication

The manuscript may state, with this evidence: *the wrong-pocket-control result on the unused-pool holdout
is not an artifact of the Vina scoring function; a scoring-free geometric contact count computed directly
from the same docked poses reproduces the same pattern (AUROC 0.55–0.71 across all three holdout pairs),
indicating that a ligand-size/burial confound — not pocket-specific chemistry — accounts for a material
share of the apparent separation.* Do **not** claim that pocket-matched docking carries no independent
information at all; the main-panel pocket-matched-vs-wrong-pocket gap (Table S6) still shows a residual
advantage for pocket matching that this contact-count proxy does not erase.
