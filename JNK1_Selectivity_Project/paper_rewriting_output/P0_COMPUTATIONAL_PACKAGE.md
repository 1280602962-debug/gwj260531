# P0/P1 Computational Package Status

Locked narrative: **Option A** (`confirmed_motivation.md`, `confirmed_contribution.md`).  
PaperSpine writing aids: `citation_bank.md`, `section_blueprints.md`.

| ID | Deliverable | Status | Path |
|----|-------------|--------|------|
| **C1** | Chemotype novelty audit (690, 2157) | **DONE** | `results/chemotype_novelty/` |
| **C2** | Multi-seed pose consensus (Vina open fallback) | **DONE** | `results/pose_consensus/` |
| **C3** | MD replica mini-panel | Protocol only | `docs/protocols/C2_C3_pose_md_replica_protocol.md` |
| **C4** | Pre-registered IC50/SI analysis | **LOCKED** (await assay) | `results/assay/`, `results/assay_analysis/` |
| **C5** | Selectivity-method autopsy | **DONE** | `results/selectivity_autopsy/` |
| **C7** | PAINS / physchem risk | **DONE** | `results/purchase_risk/` |
| **C11** | 2231 vs 690/2157 opportunity cost | **DONE** | `results/c11_2231_comparison/` |

## C1 headline

| Compound | maxTc vs literature refs | maxTc vs ChEMBL JNK | Interpretation |
|----------|--------------------------|---------------------|----------------|
| 690 | 0.23 (Q63) | 0.27 | ECFP4-distant |
| 2157 | 0.23 (Q63) | 0.27 | ECFP4-distant |

## C2 headline (Vina, 3 seeds, exhaustiveness 16)

| Compound | Isoform | Mean Vina | Pairwise RMSD mean | Consensus (pairs ≤2 Å ≥66%) |
|----------|---------|-----------|--------------------|-----------------------------|
| 690 | JNK1/2/3 | −8.27 / −8.35 / −8.93 | 0.99 / 0.62 / 1.22 | **PASS / PASS / PASS** |
| 2157 | JNK1/2/3 | −6.57 / −6.22 / −6.49 | 0.27 / 1.73 / 1.41 | **PASS / FAIL / PASS** |

Note: Open-source receptor prep (meeko, `-a` bad-res ignore); **does not replace Glide ranks**. Soften JNK2 pose claims for 2157.

## C5 headline

Δsel / Gly87 / ML selective F1 → **FAIL**; purchase **decoupled**; family ML recall gate only.

## C7 headline

690, 2157, 2231: **no PAINS**; Lipinski OK (E1 has 1 violation — control).

## C11 headline

2231: strongest MD bias + best score_JNK1, but **pose_grade C / pass_md_overall 否** → not purchased for RQ-A; document as opportunity cost for RQ-B.

## Next actions

1. Wet-lab → fill `results/assay/ic50_raw.csv` → re-run C4.  
2. C3 MD replicas when compute available.  
3. Draft manuscript from `section_blueprints.md` (Intro / Methods / RQ-C Results ready now).  
4. Methods: use license-safe Glide wording **or** cite Vina C2 as pose-consensus evidence.
