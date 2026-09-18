# Historical leftover files (not current)

These files still store the **pre-freeze ECFP incremental max |Δ| = 0.0112** (JAK1/JAK2 D vs A) or the matching pre-freeze ligand-ML / descriptor copies. They are kept as evidence of the deposited tree. **Do not write from them.**

Current sources:

| Claim | File |
|---|---|
| Figure 3A/B, SI Table S5 ECFP rows, `fig3B_max_abs` | `results/canonical/ecfp4_incremental_information.csv` |
| Figure S1A best-descriptor points | `results/canonical/descriptor_baselines.csv` |
| Submission-pack tables | `submission_pack/tables/canonical/` (copy of freeze canonical) |

Freeze value: max |Δ| = **0.0234** at PPARA/PPARD D vs B. Manuscript/SI three-decimal claim is **0.023**.

## Leftover files

| Path | What it still contains | Why it is wrong for writing |
|---|---|---|
| `data/jcim_novelty_v0/tables/incremental_information_v1.csv` | JAK1/JAK2 D vs A Δ = −0.0112; EGFR D vs A ECFP4 = 0.7992 | Pre-freeze RDKit/sklearn fold IDs. Figure 3 used this until the freeze redraw. |
| `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv` | Same 0.0112 / 0.7992 stack | Duplicate of the leftover incremental table. |
| `data/jcim_strengthen_t0t1_v0/tables/ligand_ml_baseline_scaffold_cv_v1.csv` | `auroc_ml` EGFR D vs A = 0.7992 | Pre-freeze ECFP4 OOF; Vina rank AUROC column is still Table 2 and is unused by the current plotter. |
| `data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv` | JAK1/TYK2 best cLogP `summary_min` = 0.5806 | Freeze `descriptor_baselines.csv` is 0.5796 on rdkit 2026.03.6. |
| `docs/FINAL_PRE_SUBMISSION_AUDIT.md` | Audit dated 2026-09-17 marked ECFP max \|Δ\| 0.0112 as PASS | Predates the freeze rebuild. Not a current number lock. |

Each leftover CSV has a sidecar `*.HISTORICAL_NOT_CURRENT.md` in the same directory.

`close_publication_from_canonical.py` no longer overwrites these four CSVs. The official plotter `update_figures_pr32.py` reads freeze canonical files only for Figure 3 and Figure S1A.
