# PRE_AUROC audit run V2

Task: PRE_AUROC_METHOD_ABLATION_INTEGRITY_AND_LIGAND_REPRESENTATION_AUDIT_V2

AUROC / summary_min / ΔAUC / Top-K / EF / Jaccard: not computed.

Original masters, `gnina_rescore/`, Phase 7, and canonical tables were not overwritten.

## Phase 0 — authority

- Current M0/M1/M2/M3 runners read `ablation_config.py`, not `analysis_freeze.yaml` method fields.
- `gnina.independent_docking_in_v4_1: false` and old `allowed_pdb` are STALE_PROVENANCE_ONLY.
- Meeko 0.7.1 `mk_prepare_ligand.py --rigid_macrocycles` confirmed before addenda.
- No STOP.

## Phase 1 — counts

| Level | Result |
|---|---|
| pair-expanded | 1614 = 1608 SUCCESS + 6 TIMEOUT |
| physical jobs | 1592 = 1586 with Vina pose + 6 MISSING_NO_VINA_POSE |
| physical poses | 1562 × 9 + 24 × 8 = 14250 |

## Phase 2 — ligand representation

- 785/785 frozen SDF and PDBQT present.
- Special types actually observed: CG0 and G0 only, 29 ligands.
- Every CG0 is in REMARK SMILES IDX (real carbon).
- No G0 is in SMILES IDX (glue only).

## Phase 3 — M3 topology

- 29/29 affected ligands: M3_INPUT_TOPOLOGY_FAIL.
- Frozen SDF has the ring-closure bond; historical CG0→C / drop G0 does not restore it in the torsion tree.
- Cartesian distances (~1.51–1.53 Å) recorded as notes only, not used as proof.
- 61 physical M3 jobs are INVALID_INPUT_REPRESENTATION, not TIMEOUT.

## Phase 4 — original 61 M1 failures

- 61/61 GNINA_INPUT_TYPE_FAILURE (`CG0` not a valid AutoDock type).
- All already n_retries=1; same-input retry forbidden.
- Not treated as the final M1 missing set.

## Phase 5 / 11 addenda

- `00_protocol/GNINA_RESCORE_REPRESENTATION_ADDENDUM.yaml` written before rescoring.
- `00_protocol/M3_MACROCYCLE_INPUT_REPAIR_ADDENDUM.yaml` written before M3 redock.
- Affected set frozen from phases 2–3 only.

## Phase 10 — M3 reparse

- 1590 SUCCESS: finite CNNscore from `out.pdbqt` MODEL 1.
- 2 TIMEOUT kept (AB_001 = NPHIQFCGKIKGSA-YQOHNZFASA-N × 4EY7/4BDS).
- That ligand is not in the special-type set: TIMEOUT remains missing (no rescue).
- 14/14 receptors spot-checked.
- Original master preserved.

## Phase 14–17 (partial)

- M0_LOCKED = YES (pair-expanded 1614; score = −vina_mode1_affinity; 6 TIMEOUT exact).
- M2_LOCKED = YES (physical 1586 SUCCESS, long table 14250, model1.pth, 14 × 10 Å pockets, mode1 master = long).
- Identity: 808 / 807 / 785 / 1614 / 1592. Shared physical jobs 22 (6N7A=13, 6LXA=9). M0 reuse conflicts 0.
- AB_040 alias only; not an official M0 observation.
- EH120_059 and AB_087: activity_eligible=0, class empty, not mapped to neither.

## In progress

- Unified M1/M1b representation rebuild + score_only of all 14250 poses (8 workers; resumable; first jobs SUCCESS, RMS=0, M1b≥M1).
- 29 corrected M3 ligands prepared with Meeko 0.7.1 `--rigid_macrocycles` (no CGn/Gn).
- 61 pre-frozen FAIL physical jobs redocking into `jobs_topology_repair_v2/` (concurrency 2, seed 42, timeout 2400s). Old invalid M3 outputs kept.

## Not yet

- M1 verified tables complete
- M3 topology-verified master
- METHOD_MISSINGNESS_AUDIT
- PRE_AUROC GATE_V2
