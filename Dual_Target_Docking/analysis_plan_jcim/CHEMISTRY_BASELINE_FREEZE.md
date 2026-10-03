# CHEMISTRY_BASELINE_FREEZE

Scope: M0 formal directional tasks only (`dual vs A_only` and `dual vs B_only`).
Rule type for model parameters: TYPE B (`NEW_PRE_SPECIFIED_BEFORE_ANALYSIS`).
Folds: TYPE A existing file.

This file does not train models and does not compute descriptor distributions as scientific results.

## Structure source

Unique SMILES authority:

- `reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv`
- column `canonical_smiles`

Forbidden:

- substituting `isomeric_smiles` after a canonical parse failure;
- hand-editing SMILES;
- automatic repair;
- silent deletion of unparsable ligands.

If RDKit cannot parse `canonical_smiles` for a ligand used by a `pair × arm`, that **unit** is:

- `data_readiness = BLOCKED`
- `blocking_reason = CHEMISTRY_SMILES_PARSE_FAILURE`

Other pairs/arms are not automatically blocked.

## Independent observations

Use `parent_alias_registry.csv`. Drop alias rows.

AB_040 is `DUPLICATE_PARENT_ALIAS` of AB_046. It is not an independent observation. AB_046 `fold_id=1` is the kept fold assignment.

## D1. Composition description (no model)

Eight descriptors, all recomputed from the same RDKit Mol built from `canonical_smiles`:

- MolWt
- MolLogP
- TPSA
- NumHDonors
- NumHAcceptors
- NumRotatableBonds
- RingCount
- FormalCharge = `sum(atom.GetFormalCharge() for atom in mol.GetAtoms())`

Do not mix in `pair_ligand_mapping.formal_charge`.

Also fixed:

- Bemis–Murcko scaffold from the same Mol;
- ECFP4 nearest-neighbor similarity, defined below.

Do not add descriptors. Do not drop descriptors after significance tests.

### ECFP4 nearest neighbor

For each pair and each directional arm:

- dual vs A_only, or dual vs B_only;
- same pair only;
- for each ligand, maximum Tanimoto to ligands in the **other** class of that arm;
- exclude aliases;
- do not search across pairs.

This remains a description. It is not a third chemistry baseline AUROC unless a later freeze adds that analysis.

ECFP parameters for this similarity are the same as D3.

## D2. Physicochemical logistic regression

Features: the eight descriptors above.

Preprocessing: `StandardScaler` fitted on the **training fold only**, then applied to that fold's test ligands. Fitting a scaler on the full data before CV is forbidden.

Model: L2 logistic regression with the parameters in the shared block below.

## D3. ECFP4 logistic regression

RDKit Morgan fingerprint:

- radius = 2
- nBits = 2048
- useChirality = false
- useFeatures = false
- binary fingerprint

No other fingerprint, radius, or bit length.

## Shared logistic-regression parameters (TYPE B)

- penalty = L2
- C = 1.0
- solver = liblinear
- class_weight = None
- max_iter = 1000
- random_state = 271828

No hyperparameter search.

Forbidden learners: random forest, XGBoost, LightGBM, SVM, neural networks.

These conventional defaults were fixed before chemistry performance was computed. They are not historical preregistration.

## Folds (TYPE A)

Unique file: `results/canonical/model_fold_assignments.csv`.

Fixed: 5 folds. Do not resplit, randomly split, or shuffle into new groups.

Arms in that file are the formal directional arms (`D_vs_A`, `D_vs_B`).

## Fold execution rules

Stage 1 / later readiness counts, for each `pair × arm × fold`:

- training dual count
- training negative count
- test dual count
- test negative count

If any **training** fold for a `pair × arm` contains only one class:

- that `pair × arm` is `BLOCKED` with `blocking_reason = TRAIN_SINGLE_CLASS`;
- do not resplit;
- do not skip the fold silently and train on the rest.

If a **test** fold contains only one class:

- `evidence_status = TEST_SINGLE_CLASS_WARNING`;
- do not block the `pair × arm`;
- do not report that fold's AUROC.

When later permitted, the formal chemistry AUROC is computed from the concatenation of all five test-fold out-of-fold predictions. Training AUROC is not a formal result.

If the merged OOF labeled set still has only one class:

- that `pair × arm` is `BLOCKED` with `blocking_reason = OOF_SINGLE_CLASS`.

Module-level readiness:

- if at least one `pair × arm` is executable, the chemistry module is not globally BLOCKED;
- blocked units are listed in notes.

## D4. Docking incremental model

Only one increment:

- ECFP4 versus ECFP4 + the single relevant M0 score.

`dual vs A_only` adds `M0_score_B` only.

`dual vs B_only` adds `M0_score_A` only.

Forbidden: adding both pocket scores, `min(score_A, score_B)`, or M1/M1b/M2/M3 scores.

The docking score is standardized with training-fold mean and SD, then applied to the matching test fold. Full-data standardization is forbidden.

Scientific question only: after ECFP4 already encodes 2-D chemistry, does the M0 docking score add discrimination?

## Software provenance

Stage 4 records Python, RDKit, and scikit-learn versions actually importable in the current environment. Versions are provenance. Do not reinstall or upgrade to change parse success.

## Stage 4 chemistry readiness (audit, not a rule change)

Recorded in this environment:

- Python 3.11.15
- RDKit 2025.09.5
- scikit-learn 1.9.0

`canonical_smiles` parse: 808 / 808 `MolFromSmiles` success. No `CHEMISTRY_SMILES_PARSE_FAILURE` unit.

Alias drop: 1 row (`AB_040`) removed from the 929-row fold file; 928 rows remain.

Training folds: every `pair × arm × fold` has both classes. `TRAIN_SINGLE_CLASS` = 0 units.

Test folds: no single-class test fold. `TEST_SINGLE_CLASS_WARNING` = 0.

OOF merged class counts are two-class for every `pair × arm`. `OOF_SINGLE_CLASS` = 0.

Module `data_readiness = READY_ZERO_DOCKING`.

Note (not a blocker): some fold OOF n values are 1 below the formal population (example: F2/F10 dual 31 in folds vs 32 in the population). Stage 1 does not resplit and does not treat this as `TRAIN_SINGLE_CLASS`. A later compute phase must join folds to the formal directional population and report any unmatched ligand; it must not invent a new split.
