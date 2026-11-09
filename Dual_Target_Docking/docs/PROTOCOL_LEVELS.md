# Protocol levels

Current analysis uses four named levels. Do not collapse them in captions or SI.

## PRIMARY

Production AutoDock Vina 1.2.7, current labels (θ=6.0 four-state), activity-eligible complete-case.

- PIK3CA/mTOR: exhaustiveness **E=16**
- All other primary pairs and holdouts: **E=8**
- EGFR/HER2 boxes: corrected cognate heavy-atom (`*_box_corrected.json`)
- Seeds: production `20260727`; Table 2 is this seed

## SENSITIVITY

Same ligands / labels / score definition unless noted. Not Table 2.

- PIK3CA/mTOR **E=8** on PM48
- PIK3CA/mTOR **PM110** (E=16)
- Receptor substitution **4JPS / 5DXT / 4JSX** on PM48
- Five-seed Vina (`20260727`, `20260811–14`), including AChE available-case and fixed-membership
- Label-threshold, max/median aggregation, cluster bootstrap, non-stratified bootstrap, unused-pool holdout

## INDEPENDENT GNINA

Independent pose generation (not CNN rescore of Vina poses). Three pairs only:

- EGFR/HER2
- PIK3CA/mTOR
- JAK1/TYK2

There is no eight-pair GNINA campaign.

## SAME-POSE RESCORING

RTMScore and GNINA CNN applied to already generated Vina poses. Not independent docking.

## Track-B yaml

`data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml` `do_not_now` remains a historical design snapshot. See `docs/TRACK_B_FINAL_EXECUTION_STATUS.md`.
