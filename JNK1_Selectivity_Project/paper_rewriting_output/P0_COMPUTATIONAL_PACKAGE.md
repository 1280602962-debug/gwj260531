# P0/P1 Computational Package Status

Locked narrative: **Option A**.  
**Purchase (updated):** **690 + 2231** (was 690+2157).  
Wait-window plan: `WAIT_WINDOW_PLAN_690_2231_zh.md`.

| ID | Deliverable | Status | Path |
|----|-------------|--------|------|
| **C1** | Chemotype novelty (690, **2231**) | **DONE** | `results/chemotype_novelty/` |
| **C2** | Vina multi-seed (690, **2231**) | **DONE** | `results/pose_consensus/` |
| **C3** | Unrestrained MD replicas | **NEXT (Week 2–3)** | protocol ready |
| **C4** | IC50/SI lock **v2 = 690+2231** | **LOCKED** | `results/assay_analysis/` |
| **C5** | Selectivity autopsy | **DONE** | `results/selectivity_autopsy/` |
| **C7** | PAINS | **DONE** | `results/purchase_risk/` |
| **C11** | Purchase rationale 690+2231 | **DONE** | `results/c11_2231_comparison/` |

## C1 headline

| Compound | maxTc vs literature refs | maxTc vs ChEMBL JNK |
|----------|--------------------------|---------------------|
| 690 | 0.23 | 0.27 |
| 2231 | 0.22 | 0.29 |

## C2 headline (Vina, 3 seeds)

| Compound | Isoform consensus |
|----------|-------------------|
| 690 | JNK1/2/3 all **PASS** |
| 2231 | JNK1/JNK3 **PASS**; **JNK2 FAIL** (pairwise RMSD mean 2.31 Å) |

→ Reinforces archived grade-C / off-isoform pose risk; **C3 unrestrained MD replicas are mandatory** before over-claiming 2231 bias.

## C4 lock

- Version: `c4_v2_locked_2026-07-16_purchase_690_2231`
- Primary: ≥1 of {690,2231} any isoform IC50 ≤ 10 µM
- Secondary: SI≥3 vs both JNK2 and JNK3 (hypothesis focus on **2231**)

## Month priority

1. **C3** unrestrained MD replicas for 690+2231  
2. PaperSpine Intro/Methods/RQ-C expansion  
3. C8′ / C9 / C10 as capacity allows  
4. Do not idle-wait for compound arrival  
