# P0 Computational Package Status (C1–C5)

Locked narrative: **Option A** (`confirmed_motivation.md`, `confirmed_contribution.md`).

| ID | Deliverable | Status | Path |
|----|-------------|--------|------|
| **C1** | Chemotype novelty audit (690, 2157) | **DONE** | `results/chemotype_novelty/` |
| **C2** | Multi-seed redock protocol | Protocol only — needs licensed docking | `docs/protocols/C2_C3_pose_md_replica_protocol.md` |
| **C3** | MD replica mini-panel protocol | Protocol only — needs MD compute | same |
| **C4** | Pre-registered IC50/SI analysis | **LOCKED script** (waiting for assay) | `scripts/c4_preregistered_ic50_analysis.py`, `results/assay/`, `results/assay_analysis/` |
| **C5** | Selectivity-method autopsy table | **DONE** | `results/selectivity_autopsy/` |

## C1 headline (ECFP4)

| Compound | maxTc vs literature refs | Nearest ref | maxTc vs ChEMBL JNK pool | Nearest ChEMBL |
|----------|--------------------------|-------------|--------------------------|----------------|
| 690 | 0.23 | Q63 | 0.27 | CHEMBL101035 |
| 2157 | 0.23 | Q63 | 0.27 | CHEMBL1761572 |

Interpretation: **ECFP4-distant** from curated JNK set and from E1/CC-90001/SP600125; still discuss hinge-binder pharmacophore risk in Discussion (fingerprint ≠ binding-mode novelty).

## C5 headline

| Method | Verdict | Used for purchase? |
|--------|---------|--------------------|
| Δsel_dock direction | FAIL (<55%) | **NO** |
| Gly87 occupancy | FAIL (non-discriminative) | **NO** |
| ML selective classifier F1 | FAIL (=0) | **NO** |
| ML family p_family≥6.0 | PASS as recall only | YES (activity gate only) |

## Next actions

1. Fill `results/assay/ic50_raw.csv` when wet-lab returns → re-run `python3 scripts/c4_preregistered_ic50_analysis.py`.
2. Execute C2/C3 under **licensed** Schrödinger **or** open-tool re-dock/MD per protocol + `SOFTWARE_LICENSE_NOTE.md`.
3. PaperSpine writing can proceed from confirmed Option A + C1/C5 tables.
