# Kinase-family representative audit (V2)

RUN_ID = UNIFORM_RERUN_V4_2_20260921  
Document = KINASE_FAMILY_REPRESENTATIVE_AUDIT_V2  
Status = DESIGN_ONLY; awaiting human confirmation  
Docking / AUROC / summary_min / Top-K / EF = not run, not read

Pocket residue sets are defined once from the primary cognate ligand at 4.0 / 4.5 / 5.0 Å. Protocol froze `BOX_PAD_A=5.0` for docking boxes, not a pocket-residue cutoff. No cutoff was chosen after seeing RMSDs. Pocket-perturbation bins were pre-declared: NEAR < 1.5 Å, INTERMEDIATE ≤ 3.0 Å, FAR > 3.0 Å (pocket CA RMSD).

## Family 1 pairs evaluated

| Pair | Both ends have experimental holo alt? | Fewest extra variables? | One-at-a-time possible? | Comparability | Pocket |
|---|---|---|---|---|---|
| EGFR/HER2 | YES (1XKK; 3PP0) | NO | YES | CLEAN-NEAR / MODERATE-FAR | 1XKK NEAR all cutoffs; 3PP0 FAR all cutoffs |
| JAK1/JAK2 | YES (4EI4; 8BM2) | YES | YES | CLEAN / CLEAN | 4EI4 NEAR all cutoffs; 8BM2 NEAR at 4.0/4.5, INTERMEDIATE at 5.0 |
| JAK1/TYK2 | INCOMPLETE_V2 | UNKNOWN | YES in principle | 4EI4 CLEAN; TYK2 alt not re-verified | V1 3NZ0/4GVJ ~2.1–2.3 Å not recomputed |
| PIK3CA/mTOR | YES (4L2Y; 4JT5/4JSX/4JSV) | YES | YES | CLEAN / CLEAN | 4L2Y and all three mTOR alts NEAR at all cutoffs |

## EGFR/HER2

Primary: 3POZ (EGFR, 03P/TAK-285, chain A, 1.50 Å) / 3RCD (HER2, 03P/TAK-285, chain A, 3.21 Å).

**1XKK vs 3POZ**

- Same UniProt P00533, kinase domain, ATP site, no partner.
- 1XKK ligand FMM (lapatinib/GW572016).
- No engineered mutations. 1XKK has N-terminal expression-tag residues; 3POZ does not declare tags.
- Motif CA RMSD: αC 0.579 Å; DFG 0.092 Å; A-loop 0.523 Å.
- Pocket CA RMSD 0.557 / 0.565 / 0.822 Å at 4.0 / 4.5 / 5.0.
- Tier: CLEAN. Pocket: NEAR.
- This substitution can be read as mainly holo conformation (plus a different cognate ligand).

**3PP0 vs 3RCD — do not treat as LOW**

Verified facts:

- 3PP0 mutations from mmCIF `entity.pdbx_mutation` and `struct_ref_seq_dif`: **M706A, Q711L, M712L**. Not in the unified HER2 4.5 Å pocket set (726–1004).
- 3RCD declares **Mutation(s): 0**. Residues 702–712 are labeled expression tag. Deposited 340-aa `entity_poly` sequences of 3RCD and 3PP0 are **identical**. Residues 706/711/712 in 3RCD are ALA/LEU/LEU — the same amino acids as the 3PP0 “mutations”.
- **Conflict (not silently resolved):** construct ranges 706–1009 vs 713–1028; mutation vs expression-tag labels. Sequence identity is the chemical fact.
- Copy number 2 vs 4. Resolution 2.25 vs 3.21. Ligand 03Q (SYR127063) vs 03P (TAK-285). 3PP0 paper (Aertgeerts 2011 JBC) describes a kinase dimer and αC flexibility; Lys753–Glu770 salt bridge absent.
- Motif CA RMSD 3RCD-A vs 3PP0-B: αC 8.61 Å; DFG 0.355 Å; A-loop 3.707 Å; kinase core 720–980 = 15.259 Å.
- Pocket CA RMSD 6.468 / 7.303 / 7.500 Å. Completeness 1.0.
- Tier: **MODERATE** (copy/dimer-context + metadata conflict). Pocket: **FAR**.

HER2 kinase holos remain scarce (V1 hard-pass alts = 5). 3PP0 is still the only verified non-insertion non-covalent alternative in this V2 set. Scarcity does not make FAR into NEAR.

## JAK1/JAK2

Primary: 6N7A (KEV, AB, PTR 1034/1035) / 8BXH (C87/momelotinib, A, PTR 1007/1008).

**4EI4 vs 6N7A**

- Both JH1, UniProt 854–1154, Homo sapiens, ATP site.
- Both PTR at 1034 and 1035 on both chains. No engineered mutations.
- Pocket CA RMSD 0.724 / 0.706 / 0.698 Å. A-loop 1021–1050 CA RMSD 1.622 Å.
- Tier: CLEAN. Pocket: NEAR at all three cutoffs.
- Residual: primary extract is AB; freeze whether 4EI4 uses A only or AB before any prep. V4.2 PTR→TYR fallback, if later applied, must be the same documented rule — not a new decision.

**8BM2 vs 8BXH**

- Same 2024 J. Med. Chem. paper. Both JH1 840–1132. Both PTR 1007/1008.
- His-tag differs by one Gly (MGHHHHHH vs MHHHHHH). Not a pocket mutation.
- Copy number 2 vs 1.
- Pocket CA RMSD 0.628 / 0.589 / **1.636** Å. Bin is NEAR at 4.0/4.5 and INTERMEDIATE at 5.0 (**CUTOFF_SENSITIVE**).
- Tier: CLEAN. Pocket: report both bins; do not pick 4.5 Å because it is more favorable.

JAK1/TYK2 is not eligible to win: 3LXP carries five engineered mutations (C936A, Q969A, E971A, K972A, C1142A) plus PTR1054; V2 did not re-extract 3NZ0/4GVJ.

## PIK3CA/mTOR

Primary: 4L23 (p110α 1–1068 + p85 318–615, X6K/PI-103, p85 REMOVE) / 4JT6 (mTORδN–mLST8, X6K/PI-103, extra mTOR and mLST8 REMOVE; auth A = label C).

**4L2Y vs 4L23 — same-construct ligand-induced**

- Same paper (Zhao 2014 ACS Med. Chem. Lett.). Identical p110 and p85 sequences. Same 1–1068 / 318–615.
- Ligand XXK (compound 9d) vs X6K (PI-103), same ATP/catalytic site.
- Pocket CA RMSD 0.250 / 0.278 / 0.294 Å. Completeness 1.0. No mutations.
- Tier: CLEAN. Pocket: NEAR at all cutoffs.
- **Not in V1 comparison.** HISTORICALLY_EXPOSED_CANDIDATE = NO.

**4JPS vs 4L23 — not the same thing as 4L2Y**

- p85 307–593 vs 318–615. p110 mutations **M232K, L233K** (not in the unified pocket set 772–933).
- Pocket CA RMSD 0.541 / 0.611 / 0.583 Å (NEAR).
- Tier: **MODERATE**. This is same-site ligand-bound p110+p85, but extra construct/mutation variables.
- HISTORICALLY_EXPOSED_CANDIDATE = YES (manifest + `analysis_freeze.yaml` allowed_pdb). Not excluded for that reason.

**5XGH vs 4L23**

- p110 8–1055, p85 322–598. No declared mutations. p85 present.
- Pocket CA RMSD 0.591 / 0.569 / 0.523 Å (NEAR).
- Tier: **MODERATE** (trimmed construct).

**5DXT vs 4L23**

- p110 107–1068. **No p85.** Ligand 5H5/GDC-0326.
- Pocket CA RMSD 4.931 / 5.784 / 5.752 Å (FAR).
- Tier: **CONFOUNDED**. This is construct/component/state change, not a same-construct ligand-induced conformation. V1 INTERMEDIATE was wrong.

**mTOR 4JT5 / 4JSX / 4JSV vs 4JT6**

All three are the Yang 2013 Nature mTORδN–mLST8 series.

| Alt | Ligand | Pocket 4.5 Å | Tier | Exposure | Perturbation class |
|---|---|---|---|---|---|
| 4JT5 | P2X (pp242, inhibitor) | 0.722 NEAR | CLEAN | V1 discussion | same-construct ligand-induced |
| 4JSX | 17G (Torin2, inhibitor) | 0.354 NEAR | CLEAN | manifest + freeze | same-construct ligand-induced |
| 4JSV | ADP+Mg (nucleotide) | 0.295 NEAR | CLEAN | V1 selected optional | same-construct; ligand-class change |

4JSX is not excluded because it is historically named. 4JT5 is preferred for the B-side only because it stays in the inhibitor ligand class of 4JT6; that is a construct/ligand-class comparison, not a performance comparison. If the human freeze wants the smallest pocket RMSD among CLEAN inhibitor holos, 4JSX is equally eligible (TIE with 4JT5 on comparability).

## Who is the kinase-family representative?

Selection order used (no performance, no resolution-only, no “avoid old names”):

1. Both ends need a qualified experimental holo alternative.  
2. Fewest non-conformation variables.  
3. Both ends support one-receptor-at-a-time.  
4. CLEAN > MODERATE > CONFOUNDED.  
5. If still tied: metadata completeness, pocket completeness, no extra repair.

**EGFR/HER2 is not the cleanest.** JAK1/JAK2 is cleaner than EGFR/HER2: both JAK arms are CLEAN with matched JH1 construct and matched PTR; HER2 3PP0 is MODERATE-FAR with copy/dimer-context nuisance.

**JAK1/JAK2 vs PIK3CA/mTOR = TIE on the formal rules.**

Reasons they are not split further without inventing a score:

- Both pairs have CLEAN alternatives on both ends and can run A2/B1 then A1/B2.
- PIK3CA 4L2Y is the tightest same-paper / identical-sequence ligand-induced pair in the whole kinase set, and it is NEAR at all three cutoffs. It is not historically exposed.
- mTOR 4JT5/4JSX/4JSV are NEAR at all three cutoffs under the already-frozen REMOVE policy.
- JAK2 8BM2 is CUTOFF_SENSITIVE (NEAR vs INTERMEDIATE). JAK2 copy number is 2 vs 1. JAK PTR requires applying the already-documented V4.2 fallback if those receptors are later prepared.
- PIK3CA/mTOR crystals contain p85 or mLST8 that V4.2 already removes; that is matched if the same policy is applied, but it is still a component-policy dependency.

Do not force a unique Family 1 winner. Human freeze must pick **one**:

- **Option K1:** PIK3CA/mTOR with A2=4L2Y, B2=4JT5 (4JSX equally CLEAN if preferred).
- **Option K2:** JAK1/JAK2 with A2=4EI4, B2=8BM2 (copy A; report 5.0 Å INTERMEDIATE).

EGFR/HER2 remains available as an optional FAR-arm demonstration, not as the Family 1 representative.

TYK2 is not in the TIE. V2 verification of TYK2 alternatives was not completed; 3LXP already carries five engineered mutations.
