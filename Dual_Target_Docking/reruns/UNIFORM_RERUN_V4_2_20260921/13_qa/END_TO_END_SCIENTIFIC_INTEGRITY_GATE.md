# END_TO_END_SCIENTIFIC_INTEGRITY_GATE

END_TO_END_SCIENTIFIC_INTEGRITY_GATE = PASS

New M1/M3 task-final files were audited. Historical old M1/M3 remain documented as non-candidates.
Official analysis candidates: M0, NEW M1/M1b, M2, NEW M3.

V2 running: False ; new M1 master exists: True ; new M3 master exists: True

## Blocking issues
- none

## Per-module

### Receptor preparation
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Ligand preparation
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Box
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Vina redocking
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Vina production
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### M0
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### M1/M1b old
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = PARTIAL
SCIENTIFIC_QA_PASS = NO — 61 INPUT_FORMAT_FAILURE; mixed representation; not the formal candidate

### M1/M1b new
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### M2
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### M3 old
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = NO — 29 ligands INVALID_INPUT_REPRESENTATION; not the formal candidate

### M3 new
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Alternative receptor
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Identity/reuse
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Labels/activity
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Parsing
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Missingness
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

### Cross-software interfaces
EXECUTION_COMPLETE = YES
PARSING_COMPLETE = YES
SCIENTIFIC_QA_PASS = YES

## NEW M1 summary
- master 1592 status {'SUCCESS': 1569, 'ATOM_MAPPING_FAILURE': 17, 'MISSING_NO_VINA_POSE': 6}
- long 14097; SUCCESS poses 14097; mapping-fail poses 153; accounted 14250/14250
- M1b<M1 0; mode1 mismatch 0; raw-log mismatch 0

## NEW M3 summary
- master 1592 status {'SUCCESS': 1586, 'TIMEOUT': 6}
- repaired 61/61; parse disagreements 0; AB_001 TIMEOUT kept 2

## Historical old channels (not formal candidates)
- CRITICAL retained as provenance: old M3 CG0→C/drop G0 does not restore frozen-SDF ring closure (29 ligands / 61 jobs).
- MAJOR retained as provenance: old M1 61 jobs INPUT_FORMAT_FAILURE (CG0).
- Constraint: within-method physical reuse only (same receptor + global ligand). Methods are not required to share scores.
- Phase 7 levels: 8070 pair-expanded seed records AND 7960 unique physical job-seed combinations. Both layers were checked.

No AUROC / summary_min / ΔAUC / Top-K / EF computed.
No original files overwritten.
