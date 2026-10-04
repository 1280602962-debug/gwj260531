# STAGE2_IMPL_ACCEPTANCE

Status after the implementation repair (artificial data + read-only prep).

INPUT_POPULATION_CHECKS=PASS
SECOND_STAGE_IMPLEMENTATION_READY=YES
SECOND_STAGE_COMPUTE_EXECUTED=NO

This file records checks of the **actual module functions** and writers. It does not reuse the earlier 27/27 prep table as implementation evidence.

## What was run

- `jcim_stage2_impl_check.py` on artificial tables only: **59/59 PASS**
- Temporary output directory `/tmp/jcim_stage2_impl_*` (not `results/jcim_stage2`)
- `jcim_stage2_audit.py --write-report` read-only prep: **32/32 PASS**
- Independent raw-CSV recount (no stage2 filter import): **938 / 934 / 928**
- `jcim_stage2_compute.py` with no flag: refused, `COMPUTE_EXECUTED=NO`
- `--dry-run`: no write
- `import jcim_stage2_lib` / `import jcim_stage2_compute`: no real compute
- `results/jcim_stage2` was **not** created
- `--execute-compute` was **not** run

## Artificial module evidence

| area | observed |
|---|---|
| B sign | dual-vs-neither minus dual-vs-single is positive when neither is easier, negative when single is easier |
| B empty class | neither empty → `delta_negative=NA`, `EMPTY_NEITHER_CLASS`; n_single kept |
| B pairing | one dual index recorded per replicate; single and neither draws differ; CI from replicate delta array |
| C zero | identical pocket scores → `delta_pocket=0` |
| C swap | flipping correct/wrong flips the point and CI bounds |
| C extra_drop | computed from member sets, not hardcoded |
| D test fold | single-class test fold still receives OOF predictions; unit not blocked |
| D train fold | single-class train fold → `TRAIN_SINGLE_CLASS`, no partial OOF table |
| D scaler | `module_d` probe: scaler mean equals that fold's train M0 values; train/test IDs disjoint |
| D members | three models and supporting M0 share the same ligand IDs |
| D1 | layer-1 descriptors; acyclic Murcko `""`; pair+id join prevents cross-pair SMILES overwrite |
| E | TIMEOUT listed, not imputed; incomplete seeds → reason, no IQR; 5 finite seeds → IQR written |
| F | TIMEOUT dual dropped from common set; unreplaced AUC delta is 0 when A scores are identical; affected delta changes when B alt scores change |
| I | M1/M2/M3 read PRIMARY paired points; M1b requires identical dual/A_only/B_only members across M1b−M0, M1−M0, and M1b−M1; a third-group-only extra dual (`ID_EXTRA`, M1b+M1 valid, M0 missing) yields STOP and no delta |
| D scaler | `module_d` `scaler_probe` mean equals the fold's actual train M0 values; train/test IDs disjoint |
| I states | CONCORDANT_IMPROVE, CONCORDANT_WORSEN, DISCORDANT, NO_CHANGE |
| CSV | first NA row plus later delta row keeps all schema fields; empty table keeps header; B–I path reread matches fixed schemas |
| Manifest | E/D/H/I bootstrap_B is null, not a false B=10000 claim |

## Formal inputs unchanged vs d365c901

PRIMARY tables, population, score masters, folds, mapping, and the four freeze YAML files have empty git diff against `d365c90183473c566162552fff1483006519195a`.

## Not done

- No real AUROC / bootstrap / chemistry models / joint ranks on official ligands
- No G/J/K
- No live `--execute-compute`
