# M5 — Cleaned pre-registration arm list

## Primary competition arms (≤4)

| # | id | role | Why kept |
|---|----|------|----------|
| 1 | `vina_mean` | docking baseline | Frozen protocol score; Stage M shows PIK3CA/mTOR signal vs volume |
| 2 | `heavy_atoms` | trivial control | EGFR: volume ≥ docking on `min(D/A,D/B)`; mandatory negative control |
| 3 | `rtm_min` | ≤1 rescoring | Avoid panel-wide z; if used under LOTO, calibrate only on train folds |
| 4 | `size_deconfounded_vina` | ≤1 mechanism placeholder | Formula defined to address EGFR volume confound; scores optional later |

File: `data/protocols/CANDIDATE_ARMS_V1_STAGE_M.yaml`

## Demoted / excluded from primary competition

| id | reason |
|----|--------|
| `rank_consensus_*` | EXPLORATION_DERIVED — selected after seeing panel40 CI |
| `rtm_min_z` | S1 No-Go; prep-sensitive (see M4); not a frozen universal primary |
| shortfall / clash / flags-gated | diagnostic only; not into primary score |
| pooled Dual-vs-A∪B as sole headline | cancels directional signals (M1) |

## Relation to Stage M numbers

- EGFR/HER2: vina_mean D/A≈0.69, D/B≈0.31; `heavy_atoms` beats docking on `summary_min` → M3 No-Go on this pair.
- PIK3CA/mTOR: docking arms ~0.60–0.69 directional, volume ~0.46 → M3 Go on this pair.
- Therefore Track B cannot claim a universal arm win from current Exploration pairs alone.

## Gate

**M5 = Go** if this list is frozen ≤4 primary arms with demotions documented (this file + YAML).
