# Scoop check: Dual-VSDS-Decision

**Date:** 2026-07-24  
**Skills:** ResearchStudio `scoop_check` + `paper_search`  
**Full narrative:** [`../docs/RESEARCH_DIRECTION_REFREEZE.md`](../docs/RESEARCH_DIRECTION_REFREEZE.md)

## Inputs

- **Problem:** Dual-end docking + naive score fusion fails Dual vs A-only/B-only ranking.
- **Novelty:** Architecture-agnostic decision protocol (pose QC → rescoring → calibration → weak-end gating) + four-class benchmark; NOT passenger; NOT new sampler.

## Verdict

**Level 3.5 / 5 — partially crowded, not fully scooped.**

### DELTA

Prior dual-target SBVS fuses scores to enrich known duals (often vs decoys). The open gap is a **reproducible four-class hard-negative decision problem** and a protocol that audits poses and penalizes one-ended evidence across architectures.

## Top threats

| Paper | Overlap | Surviving gap |
|-------|---------|---------------|
| Pérez-Castillo 2017 fusion | 4 | Not Dual vs A/B-only; no modern QC/calibration shortfall |
| Zhou 2013 dual kinase docking | 3.5 | Diagnoses FP; no corrective decision ruler |
| Jaiteh 2018 prospective dual VS | 3.5 | Discovery case, not reusable protocol |
| Fromer 2024 Pareto MolPAL | 3 | Acquisition optimization ≠ docking decision audit |
| Gu 2025 VSDS-VD / RTMScore | 2.5 | Single-target VS components |

## Direction ranking

1. **Architecture-agnostic decision ruler** — adopt  
2. Pure four-class benchmark resource — fallback  
3. Passenger/moiety — reject (conflicts with covering all dual architectures)

## Journal tier

Default **JCIM / J. Cheminform. / Digital Discovery**. NMI only if C3+C4+C5 all strong.
