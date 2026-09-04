# Project redesign — dual-target docking benchmark, rebuilt in the right order

Supersedes the first draft of this file. Written from the census/audit evidence plus the 2024–2026 benchmark literature. Nothing here changes Table 2 or K = 4; it is a plan.

## 1. The ordering mistake this fixes

The original project ran: hand-list candidate pairs → check ChEMBL supply → freeze receptors → dock → discover afterwards that one receptor was the wrong protein.

The correct order, which this directory now implements, is: **census → gates → identity audit → site verification → cognate QC → dock.** Under that order, PIK3CA/PIK3CB is excluded at G3 (human PIK3CB has zero PDB entries) and never becomes a docked pair. It is therefore not "a main result we corrected" — it is a pair the protocol rejects prospectively, and it belongs in the paper only as the retrospective proof that the protocol is necessary.

## 2. It is still a benchmark paper

What the field currently calls a dual-target benchmark:

- **DualDiff** (NeurIPS 2024, arXiv 2410.20688) — 12,917 target pairs over 438 targets, pairs derived from **drug-synergy annotations**, pockets from PDB **or AlphaFold** located with **P2Rank**, one "reference ligand" per target, evaluation by Vina score. There is no paired experimental activity requirement, no four-state label, no verified holo receptor per end, and no hard-negative class.
- Recent reuse of DrugComb/NCI-ALMANAC synergy tensors for target-pair benchmarking follows the same pattern.

What the field now demands of any benchmark:

- **PLINDER** (2024/2025) — leakage-aware splits over protein, pocket, protein–ligand-interaction and ligand similarity; PoseBusters physical-validity checks; novelty-stratified test sets.
- **LIT-PCBA leakage audit** (arXiv 2507.21404, 2025) — documents severe leakage, redundancy and analog bias in a de-facto standard benchmark and explicitly recommends building new sets with stricter overlap, scaffold and query controls.
- **TopU-LBVS** (NeurIPS 2026 D&B, under review) — property-matched hard decoys at fixed 1:40, with a paired random-decoy control showing EF@1% falls roughly four-fold when decoys get hard.
- **Systematic Investigation of Dual-Target-Directed Ligands** (2024, PMC11215722) — independent evidence that real DTDL programmes are driven by binding-pocket similarity and shared chemical space, i.e. the public dual space really is paralog-biased.

Verified URLs for every item above: `LITERATURE_2024_2026_DUAL_BENCHMARKS_V1.md`. Track B must not start until Layer-2 site verification is logged (`SITE_VERIFICATION_CHECKLIST_V1.md`).

The gap is obvious: **nobody has established the experimental-and-structural feasibility layer that dual-target benchmarks assume.** DualDiff asserts 12,917 usable pairs. This work shows that under paired experimental four-state labels with verified human holo receptors, the number is **17**, of which **8** are dockable with conventional rigid-receptor protocols, spanning **6** independent target systems.

So the paper stays a benchmark paper. It changes from *"here is our four-pair suite"* to *"here is what a dual-target docking benchmark can be built from, the protocol that keeps it honest, and a verified demonstration set."* That is a benchmark-construction contribution, not a retreat.

## 3. Contributions, in order of strength

1. **Feasibility ladder** (`FEASIBLE_PAIR_LADDER_V1.md`): 10.9 M possible pairs → 63,790 paired → 5,253 four-state → 86 thick → 26 scientific → 19 with human holo → 17 after ligand identity → 12 independent systems → 8 conventional-pocket dockable. Fully reproducible from a frozen ChEMBL 37 dump plus RCSB.
2. **Three-layer receptor protocol** (`TIER1_DOCKING_ROSTER_V1.md`): identity → site → cognate RMSD, with three demonstrated failure modes that each pass the layer below: a mouse p110δ receptor passing cognate RMSD at 0.405 Å; presented 9-mer peptides passing an accession match; zymogen-derived proteases whose light chain hijacks the entity length check.
3. **Ligand identity QC**: four-state labels restricted to drug-like small molecules, with Bemis–Murcko scaffold diversity reported per pair (two pairs fall below the thick gate; no surviving pair is a single congeneric series).
4. **Verified demonstration docking** on Tier-1 pairs, with the existing audit depth (unified RDKit/meeko prep, cognate gate, five-seed sensitivity, receptor-realization swaps, failure typology, ligand-only baselines, threshold grids).
5. **Structural bias finding**: only 2 of 17 feasible pairs are genuinely cross-family. Public dual-target evaluation measures **within-family selectivity**, not cross-pathway polypharmacology. This is a claim about the data, and it constrains every benchmark built on public data, including DualDiff-style ones.

## 4. What survives from the current project

| Item | Status |
|------|--------|
| Four-state formulation + worst-direction directional AUROC as primary | **Keep** — core methodological claim, independent of n |
| PIK3CA/mTOR docking (4L23 / 4JT6) | **Keep** — receptors verified human, correct accession, PI-103 on both ends |
| AChE/BChE docking (4EY7 / 4BDS) | **Keep** — receptors verified human, correct accession |
| EGFR/HER2 docking (3POZ / 3RCD) | **Keep as declared supply-limited case** — receptors verified; min HN = 7, never a thick panel |
| Five-seed sensitivity, unified RDKit prep, prep-sensitivity, failure typology, assay-context audit, ligand-only baselines, threshold grid, mixed-library enrichment, BindingDB supply freeze, claim-ceiling discipline | **Keep all** — none depends on the excluded pair |
| Receptor-realization swaps on PIK3CA/mTOR (4JPS, 5DXT, 4JSX) | **Keep** — all verified human, correct accession |
| PIK3CA/PIK3CB panel and its swap arm (2WXF frozen) | **Remove from results**; reuse as the protocol-validation case study |
| "DualFourClass-Bench" as a comprehensive/general suite | **Rename and rescope** — a verified demonstration set inside a mapped landscape |

## 5. Two execution tracks

**Track A — no new compute, required.** Reorder Methods to census-first. Move the ladder and the three-layer protocol into the main text. Restate the docked set as the Tier-1 pairs that pass all gates. Convert PIK3CA/PIK3CB into the protocol-validation subsection. Correct the Table S30 swap sentence that held 2WXF frozen.

**Track B — new compute, user's local machine, optional.** Dock the 6 new Tier-1 pairs (F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD, CTSK/CTSS): 10 receptor freezes, roughly 220 Vina jobs per pair at 110 ligands × 2 ends, ~1,300–1,500 jobs total, plus per-pair cognate QC, five-seed check and failure typology to match existing depth. Result: 8 docked pairs over 6 independent systems, with the JAK and PPAR triples supporting a within-system replication check that n = 3 cannot give.

Track B is what turns "pair-dependence observed on three pairs" into "pair-dependence quantified across six independent systems". It is the difference between a case series and a benchmark. It is not required for Track A to be publishable, and it must not be started before per-pair Layer-2 site verification is logged (`SITE_VERIFICATION_CHECKLIST_V1.md`).

## 6. Paper outline

1. **Introduction** — dual-target docking needs a four-state formulation; existing dual-target benchmarks are built on synergy annotations and predicted pockets; nobody has audited how many pairs public data can actually support.
2. **Methods** — 2.1 census and gates G1–G4; 2.2 structural feasibility G3; 2.3 ligand identity QC G4; 2.4 three-layer receptor protocol (identity → site → cognate); 2.5 docking and scoring; 2.6 statistics on independent systems, not on pairs.
3. **Results** — 3.1 the ladder (headline); 3.2 what the 86 thick pairs actually are (qHTS/CYP/metal/paralog); 3.3 Tier-1 demonstration docking, directional AUROC by pair and by system; 3.4 protocol validation: the wrong-species receptor, the presented-peptide trap, the protease chain trap; 3.5 sensitivity analyses.
4. **Discussion** — public dual-target space is paralog-biased; implications for synergy-derived, predicted-pocket benchmarks; what a future cross-pathway benchmark would require (new experimental pairing data, not more compute).
5. **Conclusions** — bounded.

## 7. Claim ceiling for the new framing

Allowed: the ladder numbers; "8 dockable pairs over 6 independent systems under these gates"; naming DualDiff-style benchmarks as unaudited on this axis; the three demonstrated receptor failure modes.

Forbidden: calling Tier 1 comprehensive or representative; treating 8 pairs as 8 independent replicates (use the 6 systems); claiming any of the 17 is docked before it is; presenting PIK3CA/PIK3CB as a corrected result rather than an excluded pair; expanding K without a written amendment.
