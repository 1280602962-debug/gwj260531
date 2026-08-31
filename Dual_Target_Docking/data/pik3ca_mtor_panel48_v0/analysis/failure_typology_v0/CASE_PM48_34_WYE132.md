# CASE PM48_34 — WYE-132 (rescued B_only)

## Identity
- **PM48_34** / CHEMBL601661 / **B_only** / WYE-132
- pChEMBL: PIK3CA **5.93** / mTOR **9.72**

## Ranks
| arm | rank |
|-----|------|
| vina_mean | **10** |
| rtm_min_z | **40** |

## Pose / score
| target | RTM mode | hinge | occ vs X6K | RTM z |
|--------|----------|-------|------------|-------|
| 4L23 | 6 | yes | 0.55 | **−1.38** |
| 4JT6 | 1 | yes | 0.68 | +0.28 |

PIK3CA shortfall end collapses under RTM → min fusion demotes correctly.

## Typology
**T1 score_artifact_rescued** — analogue of EGFR EH40_18.

## Notes
4JT6 job has **8 valid modes** (phantom mode9); recorded in Limitations. Does not change rescue conclusion.

## Protocol implication
Keep min / min_z as necessary component; cite as positive control that RTM can fix Vina B-only pollution.
