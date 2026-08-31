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
- AUROC of `contact_count` alone (dual vs A_only in pocket A; dual vs B_only in pocket B) uses the same
  class contrast as `wrong_pocket_control_vina`, but only 3D geometry, not the Vina energy. This is a
  parallel control, not a test that the two AUROCs must match in magnitude.

## Results

Mean heavy-atom count by class from `holdout_ligand_scores_v1.csv` (RDKit heavy; not PDBQT pose counts):

| pair | class | mean n_heavy (`holdout_ligand_scores_v1.csv`) |
|---|---|---:|
| AChE/BChE (HOAB) | dual | 35.1 |
| AChE/BChE (HOAB) | A_only | 34.0 |
| AChE/BChE (HOAB) | B_only | 29.5 |
| PIK3CA/mTOR (HOPM) | dual | 33.5 |
| PIK3CA/mTOR (HOPM) | A_only | 32.3 |
| PIK3CA/mTOR (HOPM) | B_only | 31.0 |
| PIK3CA/PIK3CB (HOAP) | dual | 34.5 |
| PIK3CA/PIK3CB (HOAP) | A_only | 31.6 |
| PIK3CA/PIK3CB (HOAP) | B_only | 28.3 |

Own-pocket geometric contact-count AUROC (same comparison as `wrong_pocket_control_vina`, scoring-free; **not** a magnitude match):

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

1. **Scoring-free geometry shows a real size/burial confound, mainly on the B arm — not a magnitude match to Vina.** Own-pocket contact-count AUROC is clearly above chance on dual vs B_only (0.698–0.714) on all three pairs, and closer to chance on dual vs A_only (0.552–0.622). That split tracks ligand size: dual vs B_only mean heavy-atom gaps are large (AChE/BChE 35.1 vs 29.5; PM 33.5 vs 31.0; PIK3CB 34.5 vs 28.3), whereas dual vs A_only gaps are small (35.1 vs 34.0; 33.5 vs 32.3; 34.5 vs 31.6). This is pose-level evidence of a burial confound independent of the Vina energy function. It is **not** "consistent in magnitude" with `wrong_pocket_control_vina`:

   | pair | Vina wrong-pocket summary_min (D/A, D/B) | contact_count (A / B) | contact min |
   |---|---:|---:|---:|
   | AChE/BChE | 0.643 (0.643 / 0.653) | 0.581 / 0.706 | **0.581** |
   | PIK3CA/mTOR | 0.788 (0.788 / 0.858) | 0.552 / 0.698 | **0.552** |
   | PIK3CA/PIK3CB | 0.520 (0.640 / 0.520) | 0.622 / 0.714 | **0.622** |

   The PM mismatch is the clearest: Vina wrong-pocket 0.788 vs contact min 0.552. On PIK3CB, contact min (0.622) even **exceeds** Vina wrong-pocket (0.520). Contact count therefore cannot be said to reproduce or explain the Vina holdout pattern.

2. **This is consistent with, but does not close, the confounding narrative in §3.4.** Dual ligands in these holdout draws are larger than B_only ligands and bury more surface in pocket B. The A-arm size gap is too small to carry the Vina A-arm wrong-pocket numbers, especially on PM. The pocket-matched signal in the main tables remains a mixture of any directional component and this ligand-size/burial component; contact count isolates only part of the latter.

3. **This does not mean docking carries zero pocket-specific information, nor that the holdout wrong-pocket result is fully explained.** Using the identical `wrong_pocket_control_vina` definition, the frozen main panels show pocket-matched summary_min above wrong-pocket summary_min for all four pairs (Table S6; e.g. AChE/BChE 0.606 vs 0.444, PM 0.692 vs 0.602). On the unused-pool holdout the gap reverses (wrong-pocket ≥ matched). We do not have a resolved explanation for this contrast. Report it as an open discrepancy. Do not claim that contact count is sufficient by itself to produce the holdout Vina wrong-pocket pattern.

Mean heavy-atom counts are from `holdout_ligand_scores_v1.csv` (RDKit heavy), not from PDBQT pose atom counts.

## Claim implication

The manuscript may state: *a scoring-free geometric contact count computed from the same docked poses shows a real ligand-size/burial confound, especially on the B arm (AUROC 0.698–0.714); it does not reproduce Vina wrong-pocket in magnitude (PM: Vina 0.788 vs contact min 0.552), so it does not explain the holdout wrong-pocket result.* Do **not** claim that the wrong-pocket-control-is-not-worse pattern "reproduces at the geometry level" in both direction and magnitude, and do **not** claim pocket-matched docking carries no independent information — the main-panel pocket-matched-versus-wrong-pocket gap (Table S6) is still present.
