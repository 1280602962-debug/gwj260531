# Novel review panel — keep labels v1

Agent structural triage of `novel_review_panel` (50 compounds). Values: `yes` / `no` / `unsure`.

## Counts

| keep | n |
|------|---|
| yes | 21 |
| no | 21 |
| unsure | 8 |

## By stratum

| stratum | yes | no | unsure |
|---------|-----|----|--------|
| erg_top15 | 8 | 4 | 3 |
| erg_median15 | 8 | 5 | 2 |
| tc_core_boundary10 | 4 | 4 | 2 |
| hinge0_erg_only10 | 1 | 8 | 1 |

## Decision rules used

1. Prefer clear ATP-site hinge chemotypes (anilinopyrimidine, pyrrolopyrimidine, imidazopyridine, pyridopyrimidinone, isoquinolinone, aminopyridine, etc.).
2. Single monoacrylamide with solvent/Cys-plausible vector.
3. Reject PAINS, amidoxime/amidine acryloyl warheads, hydrazine-N acrylamides, obvious peptidomimetics without hinge, MCR/poly-CN junk.
4. `hinge_hits=0` is a soft prior only — several FNs kept when structure shows hinge (e.g. CHEMBL3410050, CHEMBL333483).
5. `unsure` = atypical but not obvious junk; do **not** expand Novel on these alone.

## Calibration implications (for Step B)

- **ErG-only / hinge0 stratum**: 1 yes / 8 no / 1 unsure → supports tightening Novel to require hinge-like feature OR stronger structural filter beyond ErG alone.
- **High ErG without hinge** still produced several `no` → do not trust ErG rank alone.
- **Boundary Tc≈0.22**: mix of yes/no; amidino-warhead and hydantoin rejected even at boundary.
- Suggested provisional gate after this panel: Novel keep if `(hinge_like>=1 OR manual_yes) AND not (amidine_acryloyl OR PAINS OR hydrazine_acryloyl)`.

## File

- Labels: `L7/novel_review_panel_keep_v1.csv`
