# BindingDB external-slice feasibility (final SI packaging)

**Status:** completed negative feasibility result. **No docking.**

Pre-frozen filters (literature → structure → ECFP4 < 0.70 → class/source gates)
yielded **0** target pairs eligible for primary external evaluation.
`packaged_as_external_evaluation = 0` for all audited pairs.

## Funnel (machine-readable)

- Wide: `data/jcim_novelty_v0/tables/bindingdb_external_feasibility_flow_v1.csv`
- Long: `data/jcim_novelty_v0/tables/bindingdb_external_feasibility_funnel_long_v1.csv`
- Per-layer source: `data/jcim_novelty_v0/tables/external_candidate_flow.csv`
- Gate summary: `data/jcim_novelty_v0/tables/external_slice_summary_v1.csv`

## Condensed flow

```
raw BindingDB (θ=6.0 paired)
        ↓
literature exclusion
        ↓
structure exclusion
        ↓
ECFP4 max-sim < 0.70
        ↓
class / source concentration gates
        ↓
0 eligible pairs  →  no docking
```

## Claim ceiling

- Allowed: strict database-external bidirectional hard negatives are scarce under this freeze.
- Forbidden: relaxing thresholds after seeing counts; calling remaining ligands an external set;
  presenting BindingDB as external validation.

| pair | after ECFP n | dual/A/B/neither | gate | docked |
|---|---:|---|---|---|
| EGFR/HER2 | 216 | 180/10/20/6 | insufficient | no |
| AChE/BChE | 85 | 4/8/14/59 | insufficient | no |
| PIK3CA/PIK3CB | 112 | 9/0/3/100 | insufficient | no |
| PIK3CA/mTOR | 98 | 91/4/1/2 | insufficient | no |
| MCL1/Bcl-xL | 3 | 1/0/2/0 | insufficient | no |
