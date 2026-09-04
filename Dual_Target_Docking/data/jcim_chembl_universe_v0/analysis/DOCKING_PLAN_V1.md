# Docking plan (locked 2026-09-04)

Written amendment after the user’s decision:

1. **CTSK/CTSS is out of the ordinary noncovalent Vina campaign.** Both holos are reversible-covalent (4X6H Cys25–3XT 1.83 Å; 9GJ2 Cys25–KH0 1.78 Å). A unified AutoDock Vina benchmark cannot absorb that mechanism without mixing chemical tasks. Do not covalent-dock them in this paper. Receptor freeze for 4X6H/9GJ2 remains as a structure record; they are **not** Track B production pairs.
2. **PIK3CA/PIK3CB docking already done is kept**, in the same spirit as EGFR/HER2: a **declared special case**, not a thick ordinary pair and not an isoform-control result. Do **not** re-dock it. Do **not** replace 2WXF. Report Table 2 0.500 [0.350, 0.650] as human-PIK3CB labels scored in a **mouse p110δ** pocket (O35904), i.e. a receptor-identity protocol case.
3. **New Vina is five pairs only**, all ordinary noncovalent declared pockets, all co-tested four-state thick after the small-molecule filter.

Table 2 / K = 4 stay frozen. This plan does not unfreeze K and does not restock Table 2 from census ranks.

## A. What is already docked (do not re-run)

| Pair | Role in the paper | Receptors | Notes |
|------|-------------------|-----------|-------|
| PIK3CA/mTOR | ordinary thick primary | 4L23 / X6K, 4JT6 / X6K | Keep. ATP not FRB. |
| AChE/BChE | ordinary thick primary | 4EY7 / E20, 4BDS / THA | Keep. |
| EGFR/HER2 | **special case: supply-limited** | 3POZ / 3RCD | Keep existing EH110. min HN = 7. Never call it a thick panel. |
| PIK3CA/PIK3CB | **special case: receptor-identity failure** | 4L23 + **2WXF = mouse p110δ** | Keep existing numbers. Forbidden: “PI3Kα/β isoform control”; re-docking on 2Y3A/4BFR or a homology model without a new written exception. |

Primary directional table remains these four rows, with the two special-case labels above. That is the EGFR/HER2 precedent applied to PIK3CA/PIK3CB.

## B. What gets newly docked (Track B)

Five pairs, eight unique new receptors (JAK1 and PPARA reused). All passed co-tested four-state thick gates (`min` strict A-only/B-only small-molecule ≥ 50).

| # | Pair | min HN (small-mol) | Receptor A | Receptor B | Declared site |
|---|------|-------------------:|------------|------------|---------------|
| 1 | F2/F10 | 108 | 4UDW / N6L | 2JKH / BI7 | S1 / S1–S4 |
| 2 | JAK1/TYK2 | 91 | 6N7A / KEV | 3LXP / IZA | JH1 ATP both ends |
| 3 | JAK1/JAK2 | 53 | 6N7A / KEV | 8BXH / C87 | JH1 ATP both ends |
| 4 | PPARG/PPARA | 82 | 9V8H / BRL | 6LXA / EPA | LBD; 9V8H is ternary (+PG08-NL) |
| 5 | PPARA/PPARD | 82 | 6LXA / EPA | 5U3Q / 7UJ | LBD |

These are an **extension / replication panel**, not a replacement for Table 2. Do not pick extra pairs after seeing AUROCs.

Independent new systems = **3** (coagulation; JAK family; PPAR family), not 5.

## C. What is explicitly not docked

| Pair | Why |
|------|-----|
| **CTSK/CTSS** | Reversible-covalent warheads. Out of ordinary Vina. No covalent docking campaign in this paper. |
| CREBBP/BRD4 | G5: domain (HAT vs bromodomain) undeclared. |
| GPCRs / SLC6 | G5: membrane conformational state. |
| F2/PRSS1 | Trypsin antitarget, not a designed dual. |
| JAK3/TYK2, OPRM1/OPRK1 | Failed small-molecule thick gate (48 and 46). |

## D. Execution order (nothing after a skipped gate)

Cloud has RDKit + meeko + the ChEMBL 37 SQLite. Cloud **does not have Vina**. Steps 1–2 can run here; steps 3–5 need a local Vina machine.

1. **Extract four-state panels** from ChEMBL 37, co-tested only (`set(map_A) ∩ set(map_B)`). Same endpoints and max pChEMBL as the census. Labels: strict 6.5/5.5 for AChE-style thick pairs; record θ = 6.0 counts in parallel. Sample with the same class-quota + deterministic shuffle used on the frozen panels (seed **20260729**; target depth **110** ligands: dual / A-only / B-only, and neither **only if** `strict_neither` ≥ 20 after small-molecule filter). JAK1/TYK2 strict neither after small-mol = 25 (neither arm included); do not force a neither arm on pairs that only look thick directionally.

   **Done 2026-09-04** (`scripts/extract_track_b_panels_v1.py`). All five pairs sampled 32/32/32/14 = 110. Small-molecule min HN matches the QC table (108 / 91 / 53 / 82 / 82). CSVs: `tables/track_b_panels/panel_*_v1.csv`.
2. **RDKit ETKDG + meeko PDBQT.** No LigPrep. Same as Methods (ETKDGv3 seed **20260727**, largest fragment, MMFF 200).

   **Done 2026-09-04** (`scripts/prep_track_b_ligands_v1.py`). 550/550 ligands wrote PDBQT (`tables/track_b_ligand_prep_status_v1.csv`). Binaries live in gitignored `cache/track_b_ligands/`. Re-run the script on a Vina machine to regenerate them.
3. **Layer 3 cognate best-of-9 RMSD** on the eight new receptors, box = cognate heavy-atom AABB + 5 Å/axis, min edge 20 Å. Fail the receptor if no near-crystal pose is generated (RMSD is not identity proof; 2WXF passed).
4. **Production Vina** on both ends of each of the five pairs.
5. Optional, only if matching frozen depth: five-seed Vina, failure typology (N_attempted / N_successful / N_failed).

Approximate new Vina volume: 5 pairs × ~110 ligands × 2 ends ≈ **1,100** production jobs, plus 8 × 9 cognate poses. Not 6 pairs / 1,300–1,500 (that count included CTSK/CTSS).

## E. How to report new AUROCs

- Estimand unchanged: pocket-matched directional AUROC (dual vs A-only uses pocket B; dual vs B-only uses pocket A); `summary_min` as worst-direction summary.
- Count **systems**, not pairs: JAK1/TYK2 and JAK1/JAK2 are one JAK system; PPARG/PPARA and PPARA/PPARD are one PPAR system.
- 9V8H: Methods sentence on PG08-NL; do not over-interpret PPAR subtype conformation from one ternary crystal.
- ChEMBL JAK labels are target-level, not JH1 vs JH2. Docking tests discrimination **at the declared JH1 ATP pocket**, not proof that every ChEMBL active binds there.
- Do not pool CTSK/CTSS into a mean AUROC with the five ordinary pairs (they are not being docked).
- Do not pool PIK3CA/PIK3CB 0.500 into a “kinase isoform success” narrative.

## F. Claim ceiling (this amendment)

Allowed: five new ordinary noncovalent pairs as a predeclared extension; PIK3CA/PIK3CB retained as a special case analogous to EGFR/HER2; CTSK/CTSS listed among census-evaluable pairs but excluded from Vina for mechanism.

Forbidden: calling CTSK/CTSS a completed Vina pair; covalent docking as if it were the same experiment; replacing Table 2 with extension AUROCs; presenting PIK3CA/PIK3CB as human p110β isoform control; docking CREBBP/BRD4 or GPCRs under this plan.
