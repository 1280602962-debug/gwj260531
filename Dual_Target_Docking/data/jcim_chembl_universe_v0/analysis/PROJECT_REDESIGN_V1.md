# Objective redesign of the DualFourClass project

Written from the evidence in this directory, not from sunk cost. Nothing here changes Table 2, K = 4, or any manuscript file; it is a proposal.

## 1. The question the paper should actually answer

Not “is our four-pair benchmark good enough,” but:

> **In public data, how many dual-target pairs can support a four-state docking evaluation at all — and what does that ceiling do to any benchmark built on them?**

That question is now answered with numbers at every gate (`FEASIBLE_PAIR_LADDER_V1.md`): 10.9 M possible pairs → 63,790 with paired data → 5,253 four-state → 86 thick → 26 scientific → 19 with human holo structures → **17** after keeping only drug-like small molecules → **12** independent target systems.

This is the contribution nobody has published. It is a resource, it is falsifiable, and it does not depend on how many pairs were docked.

## 2. Direction change

| Now | Proposed |
|-----|----------|
| Four-pair benchmark suite (“DualFourClass-Bench”) | **Feasibility landscape** of public dual-target docking evaluation + a small verified demonstration set |
| n = 4 is the headline | The ladder (17 / 12 systems) is the headline; docked pairs are the worked examples |
| “Bench” implies reusable suite | Drop the suite claim; publish the ladder, the label rules, and the identity/QC protocol as the reusable part |
| Pair-dependence asserted from 4 pairs | Pair-dependence stated as *observed on the verified pairs*, with the ladder explaining why more pairs are hard to get |

The framing change costs no compute and removes the single most likely reviewer objection.

## 3. Which docking results survive

Receptor identity was verified for every receptor actually used (`RECEPTOR_IDENTITY_AUDIT_V1.md`).

| Pair | Receptors verified | Supply gate | Keep? |
|------|--------------------|-------------|-------|
| **PIK3CA/mTOR** | 4L23 / 4JT6 both human, correct accession, PI-103 cognate on both ends | thick (min HN 71 small-molecule) | **Keep as primary** |
| **AChE/BChE** | 4EY7 / 4BDS both human, correct accession | thick (65) | **Keep as primary** |
| **EGFR/HER2** | 3POZ / 3RCD both human, correct accession | **not** thick (7) | **Keep as declared supply-limited case**, never as a thick panel |
| **PIK3CA/PIK3CB** | **2WXF is mouse p110δ (O35904), not human PIK3CB** | fails after small-molecule filter (56 → 40) and has no human structure | **Withdraw from primary table**; report as receptor-identity failure case |

So three of four docked pairs survive unchanged. All the surrounding work — five-seed sensitivity, unified RDKit prep, prep-sensitivity analysis, ligand-only baselines, assay max-versus-median, threshold grids, failed-ligand stress tests, mixed-library enrichment, claim ceilings — survives, because none of it depends on the withdrawn pair.

The withdrawn pair is not wasted: it becomes a concrete demonstration that cognate-RMSD QC **cannot** catch a wrong-protein receptor, which is exactly the failure mode a benchmark protocol must guard against.

## 4. What to add, in priority order

1. **Zero-compute, mandatory.** Fold the ladder, the receptor-identity audit, and the ligand-identity QC into Methods/SI. Correct the Table S30 swap sentence. Restate K as three verified docked pairs plus one documented failure.
2. **Cheap, high value.** Publish the identity-QC protocol as the reusable artifact: accession + organism + entity-coverage check **before** cognate redocking, plus the drug-like small-molecule label filter. This is the part other groups can adopt.
3. **Optional expansion, only if you want more pairs.** The homogeneous soluble subset from the ladder: F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD, CTSK/CTSS, plus the two kept pairs — 8 pairs across 6 independent systems, one pocket-prep protocol (kinase ATP / protease / hydrolase / NR LBD), roughly 1,300–1,800 Vina jobs and ~12 receptor freezes. Exclude F2/PRSS1 (trypsin is an antitarget).
4. **Do not do now.** The 6 GPCR pairs, 2 SLC6 pairs, and CREBBP/BRD4. Membrane constructs and the CREBBP domain choice add prep confounds that would weaken, not strengthen, the measurement claim. List them as the declared frontier with their H3 evidence.

## 5. Honest statement of what the paper then claims

- Four-state (dual / A-only / B-only / neither) is the right formulation, and worst-direction summary is the right primary readout.
- Public data admit at most 17 such pairs over 12 independent systems, and only 2 of those are genuinely cross-family; the rest are paralog selectivity. Public dual-target docking evaluation is therefore **structurally biased toward selectivity**, not polypharmacology.
- On the verified pairs, apparent discrimination is strongly pair-dependent and docking adds little over ECFP4.
- Receptor identity, not search effort, is the dominant silent failure mode; cognate RMSD does not detect it.

## 6. What must not be claimed

- No “comprehensive” or “general” dual-target benchmark.
- No treating 17, 19, or 86 as docked or as independent replicates.
- No PIK3CA/PIK3CB isoform-control result until a defensible human or explicitly-declared surrogate receptor exists.
- No expansion of K without a written protocol amendment.
