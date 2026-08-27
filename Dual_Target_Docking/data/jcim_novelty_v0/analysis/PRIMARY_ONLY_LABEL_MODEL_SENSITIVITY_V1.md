# Adjudicated primary-only label and model sensitivity v1

## Scope

This is a decision-targeted primary-source sensitivity analysis, not a claim that every ChEMBL record has been paper-level verified. Explicitly adjudicated review-derived, target-mapping-error, and cellular/pathway-surrogate rows are removed. Unreviewed high-confidence rows remain unchanged. A missing arm after removal is `unresolved_missing_arm`, never inactive.

- Threshold: pActivity >= 6.0.
- Bootstrap: ligand resampling, B=2000, seed base=20260828.
- Explicitly excluded activity rows: 11.
- Class changes versus API-max: 1.
- Unresolved ligands after exclusion: 1.

## Changed or unresolved labels

| Pair | Ligand | API-max | Adjudicated primary-only | A max | B max |
|---|---|---:|---:|---:|---:|
| EGFR/HER2 | EH120_059 | A_only | unresolved_missing_arm | 7.54 | None |

## Model sensitivity

| Pair | Model | D/A | D/B | summary_min [95% CI] | delta vs API-max | delta vs frozen | unresolved |
|---|---|---:|---:|---:|---:|---:|---:|
| EGFR/HER2 | vina | 0.6622 | 0.4297 | 0.4297 [0.2824, 0.5815] | 0.0 | 0.0 | 1 |
| EGFR/HER2 | rtmscore | 0.5164 | 0.3527 | 0.3527 [0.2109, 0.5011] | 0.0 | 0.0 | 1 |
| AChE/BChE | vina | 0.6504 | 0.6058 | 0.6058 [0.4484, 0.7348] | 0.0 | 0.0 | 0 |
| AChE/BChE | rtmscore | 0.5185 | 0.545 | 0.5185 [0.3413, 0.6296] | 0.0 | 0.0 | 0 |
| PIK3CA/PIK3CB | vina | 0.6905 | 0.5 | 0.5 [0.3495, 0.6454] | 0.0 | 0.0 | 0 |
| PIK3CA/PIK3CB | rtmscore | 0.705 | 0.5421 | 0.5421 [0.3763, 0.6952] | 0.0 | 0.0 | 0 |
| PIK3CA/mTOR | vina | 0.7143 | 0.6921 | 0.6921 [0.4721, 0.7976] | 0.0 | 0.0 | 0 |
| PIK3CA/mTOR | rtmscore | 0.6151 | 0.6574 | 0.6151 [0.3889, 0.758] | 0.0 | 0.0 | 0 |

## Interpretation guardrail

The `adjudicated_primary_only` scenario is suitable as a transparent SI sensitivity table. It cannot yet support the sentence 'all labels were reconstructed exclusively from primary papers'. That stronger statement requires paper-level tier assignment for every record that can determine a ligand's maximum on either target.
