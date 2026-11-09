# PR #35 item-by-item (no whole-branch merge)

Source: PR #35 `cursor/unify-cognate-rmsd-qc-0b1a`, unique commit `bb3f134f7f8895f2f361b58fdb485530b677d1c4`  
Base of that PR: `cursor/jcim-submission-pack-0b1a` (pre-activity-eligible, pre-PR #38).  
Current freeze baseline: `17435413410290960da3437de640509705f00b51`.

Raw evidence remains on the PR #35 branch and in git history. This round copies nothing from that commit into canonical RMSD.

| Item | What PR #35 did | Adopt? | Reason |
|---|---|---|---|
| One 14-receptor chemically mapped CalcRMS table as the primary RMSD source | Added `primary_cognate_rmsd_calcrrms_v1.csv` | **Principle yes, file no** | Current tree already has `data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv` (14 rows, CalcRMS). That is the freeze RMSD source. |
| EGFR 3POZ / HER2 3RCD numbers 9.505 / 0.760 and 1.855 / 1.394 | Reconstructed-QC poses labeled as such | **No** | Current table is `canonical_heavy_atom_redock` after corrected boxes: EGFR top-1=best=1.019 Å, HER2 1.947 Å. Replacing with 9.505 would reintroduce a number the manuscript now treats as historical. |
| Non-EGFR/HER2 CalcRMS values (AChE 0.339, mTOR best 0.445, JAK2 best 0.807, …) | Same metric | **Already present** | Current 14-row table matches those best-Å values. No copy needed. |
| `assemble_primary_cognate_rmsd_v1.py` | Assembler for the PR #35 table | **Keep as evidence only** | Not in the current tree. Do not restore it to overwrite corrected-box EGFR RMSD. |
| `math.isfinite` on fixed-membership scores | `parse_finite_energy`; 0 extra exclusions; 545 members unchanged | **Yes, as shared helper** | `analysis_config.parse_finite` is what master/build and score parsers use. Multiseed script path is absent here; do not restore the old audit runner. |
| `multiseed_fixed_membership_exclusions_v1.csv` | Empty extra-exclusion log | **No** | Membership for current Table 2 is activity-eligible complete-case on `current_score_master.csv`, not the 545-row five-seed intersection. |
| Checksum manifest / `validate_revision_v1` / five-round audit docs | Packaging QC | **No** | Replaced by independent freeze verification. Old PASS is not acceptance. |
| Manuscript/SI wording that one CalcRMS table covers all 14 receptors | Language | **Already on PR #38** | Current SI Fig S3 caption already distinguishes historical 9.505 from current 1.019. |
| Dropping internal “最终认可 / 槽位” Chinese audit asides | Language | **Already on later PRs** | Do not merge the old manuscript files. |

## Handling recommendation for PR #35 itself

- Do **not** merge PR #35 into this freeze or into main.
- Do **not** close or delete the PR/branch until the user confirms.
- Treat it as a closed scientific decision record: unification adopted; reconstructed-QC EGFR/HER2 Å rejected for current SI.
