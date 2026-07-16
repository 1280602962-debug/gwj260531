# Confirmed Core Contribution (Option A)

**Status:** LOCKED  
**Depends on IC50?** Strengthens RQ-A if active; **RQ-C stands without IC50**.

---

## Core contribution (one sentence)

We deliver an end-to-end commercial-library pipeline (ML → Glide → ADMET → MD → purchase) that enriches **pose-credible JNK-family binder candidates**, and we show on a literature JNK benchmark that **common computational isoform-selectivity predictors (Δsel docking, Gly87 occupancy, ML selectivity labels) fail as purchase criteria**. Wet-lab purchases (**690** as family/pose-QC anchor; **2231** as the strongest MD JNK1-bias hypothesis despite grade-C overall MD flags) test family enrichment primarily and isoform preference only secondarily under pre-registered SI rules.

---

## Evidence already in archive

| Claim element | Evidence |
|---------------|----------|
| Funnel | ~4979 docked → 157 (F1∧F2) → 25 ADMET → 16 MD → purchase **690+2231** |
| Δsel docking direction | VSW single-PDB accuracy **43%** (3/7); ensemble archive **29%** (2/7) |
| Gly87 occupancy | 5/5 benchmarks `occ_JNK1=True`; distances 0.59–1.18 Å — non-discriminative |
| ML selectivity classifier | Positive n=8; test **F1 = 0** |
| ML family gate | High recall, decoy FPR 95.3% @ p_family≥6.0; enrichment via ranking (EF1%=9.2) |
| MD QC ≠ selectivity | E1 hinge mis-rank; SP600125 active with low hinge — documented |
| Purchase decoupling | Shortlist chosen for family activity / pose credibility / bias hypothesis, **not** pass_selectivity |

## Evidence pending

| Element | Status |
|---------|--------|
| JNK1/2/3 IC50 for 690, 2231, E1, CC-90001 | Wet-lab (compounds on order) |
| C2 Gnina layer + redock | Open-source; in wait-window plan |
| C3 unrestrained MD replicas | Highest-priority compute |
| C1/C4/C5/C7/C11 for 690+2231 | Largely done |

---

## Allowed strong claims

- Computational isoform-selectivity filters tested here are **insufficient** for JNK purchase decisions.  
- Pipeline produces a **documented, pose-QC’d shortlist** suitable for family-activity assays.  
- MD pass criteria are **pose/stability QC**, not isoform adjudication.

## Forbidden / must-soften claims

- “Discovered JNK1-selective inhibitors” (unless SI rule met *and* powered — currently not).  
- “MD confirmed selectivity.”  
- “Kinase-selective / kinome-clean.”  
- Hit-rate statistics from n=2 new molecules.

---

## Contribution check (pre-writing)

| Gate | Status |
|------|--------|
| Motivation ↔ contribution aligned | YES (Option A) |
| Selective-discovery not core | YES |
| Negative method result primary | YES |
| Wet-lab framed as enrichment/calibration | YES |
