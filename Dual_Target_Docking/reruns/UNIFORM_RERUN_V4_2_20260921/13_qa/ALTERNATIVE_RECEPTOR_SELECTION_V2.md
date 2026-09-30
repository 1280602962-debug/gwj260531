# Alternative-receptor selection V2 (fact verification and ablation freeze)

RUN_ID = UNIFORM_RERUN_V4_2_20260921  
Document = ALTERNATIVE_RECEPTOR_SELECTION_V2  
Status = DESIGN_FROZEN_AWAITING_HUMAN_CONFIRMATION  
Docking = NOT STARTED  
AUROC / summary_min / Top-K / EF / LEGACY or V4.2 performance = NOT READ, NOT USED

This document does not overwrite V1. V1 files remain in `13_qa/` and `00_protocol/`.

This task does not search for a receptor that improves docking. It only verifies identities, splits comparability from pocket perturbation, and freezes a one-receptor-at-a-time ablation design.

No receptor preparation, redocking, production docking, or AUROC was run.

---

## Method and sources

Independently extracted from RCSB mmCIF, wwPDB validation XML, and primary papers for the planned primary/alternative set (32 structures). Pocket residues were defined **once per target** from the primary cognate ligand at 4.0 / 4.5 / 5.0 Å. All alternatives of a target use that same residue set. Sequence mapping + Kabsch pocket CA RMSD. Mutation/PTM from `_entity.pdbx_mutation`, `_struct_ref_seq_dif`, and `_pdbx_struct_mod_residue` — not from V1 comments.

Protocol freeze has `BOX_PAD_A=5.0` for docking boxes. It does **not** freeze a pocket-residue cutoff. All three distances are reported. Bins were pre-declared in `scripts/rerun_v4_2/alt_receptor_fact_verify_v2.py`: NEAR < 1.5 Å, INTERMEDIATE ≤ 3.0 Å, FAR > 3.0 Å.

V1 RCSB availability inventory (2318 rows) was **not re-queried**. V2 availability file carries that inventory only as per-target counts. No new PDB search.

Performance files were not opened. `alternative_receptor_manifest.csv` was opened for ID list only. Score CSVs named there were not opened.

---

## 1. V1 factual errors

Listed only where V2 independent sources disagree with V1 text or tables. “UNKNOWN in V1” is not counted as an error if V1 marked it unknown.

1. **Single LOW / INTERMEDIATE / HIGH axis.** V1 mixed “how many extra variables” with “how large is the pocket change”. HER2 3PP0, PIK3CA 5DXT, and PPARG 7AWC are the clearest casualties.

2. **HER2 3PP0 called LOW.** Unified-pocket CA RMSD vs 3RCD is 6.47 / 7.30 / 7.50 Å (FAR). αC CA RMSD = 8.61 Å; kinase-core 720–980 CA RMSD = 15.26 Å. This is not a low-perturbation conformation-only swap.

3. **3PP0 mutations not identified.** V1 wrote “surface engineered mutations”. mmCIF lists **M706A, Q711L, M712L**. They are N-terminal and **not** in the unified 4.5 Å pocket set. V1 did not state the residues or check pocket membership.

4. **3RCD vs 3PP0 mutation/construct conflict not disclosed.** Deposited 340-aa sequences are identical. 3RCD labels 702–712 as expression tag (ALA/LEU/LEU at 706/711/712). 3PP0 labels those residues engineered mutations. V1 treated 3PP0 as mutated and 3RCD as clean without listing the conflict.

5. **4L2Y omitted.** Same paper and identical p110/p85 sequences as 4L23. Pocket CA RMSD 0.25–0.29 Å. This is the cleanest PIK3CA alternative. V1 never compared it.

6. **PIK3CA 4JPS called LOW.** 4JPS has **M232K, L233K** and p85 307–593 vs 4L23 318–615. Pocket is NEAR, but comparability is MODERATE, not CLEAN/LOW.

7. **PIK3CA 5DXT called INTERMEDIATE.** No p85; construct 107–1068; pocket CA RMSD 5.78 Å (FAR). This is CONFOUNDED (component + construct + state), not a modest intermediate.

8. **PPARG 7AWC called INTERMEDIATE.** Peptide PG08-NL is absent; pocket CA RMSD 5.62 Å (FAR). CONFOUNDED, not INTERMEDIATE.

9. **BChE “same engineered-mutation list” stated without residues.** The conclusion happens to be correct, but V1 did not list them. Both 4BDS and 1P0M have **N17Q, N455Q, N481Q, N486Q** (auth 17/455/481/486; UniProt 45/483/509/514). 4BDS `entity.pdbx_mutation` is only `YES`; residue identities come from `struct_ref_seq_dif`.

10. **F10 “both have a mutation” treated as matching.** V1 did not identify the residue. Both 2JKH and 2Y5F have the **same** mutation: **R372E (auth A:150)**. The identity match is now verified; V1 did not verify it.

11. **JAK PTR not compared.** V1 called 4EI4 and 8BM2 LOW without checking phosphorylation. Both primary and alternative JAK1 structures have PTR 1034/1035; both JAK2 structures have PTR 1007/1008. That part is actually matched — V1 did not establish it.

12. **Pocket cutoff 4.5 Å presented as if frozen.** Protocol froze box pad, not pocket distance.

13. **Validation reports not used.** Several primaries have material quality issues that V1 did not record (3RCD clashscore 11.83, rota outliers 14.05%; 4JT series clashscore 14–19).

14. **V1 Family 1 pick used “avoid a second kinase pair”.** That is family-balance language, not the V2 priority order (qualified holo → fewest extra variables → one-at-a-time → CLEAN > MODERATE → completeness). Under V2 rules EGFR/HER2 loses to JAK1/JAK2 and is not uniquely better than PIK3CA/mTOR.

15. **V1 avoided 4JSX as the named mTOR alt to reduce historical overlap.** Historical exposure is not an exclusion rule. 4JSX remains CLEAN-NEAR and is eligible.

16. **First-pass ligand inference for 5U46.** V2 automated ligand pick initially selected B7G (detergent). The orthosteric agonist is **7T1 (GW501516)**. V1 freeze YAML already had 7T1. Pocket RMSD used the primary 7UJ residue set and is unaffected.

V1 items that **replicated**:

- 1XKK, 4EY6, 6KAX, 5U46 identities and same-campaign relationship.
- 4EY6 sequence/construct identity with 4EY7.
- Hirudin retained on 3SHC.
- 5U3Q copy A is not an alternative receptor.
- PPARG has no peptide-retained CLEAN alternative.
- Predicted structures excluded.

---

## 2. Former LOW receptors after splitting the two axes

| Target | PDB | V1 | V2 comparability | V2 pocket (4.0/4.5/5.0) | Reclassification |
|---|---|---|---|---|---|
| EGFR | 1XKK | LOW | CLEAN | NEAR / NEAR / NEAR | Keep as clean-near |
| HER2 | 3PP0 | LOW | MODERATE | FAR / FAR / FAR | **Not LOW** |
| JAK1 | 4EI4 | LOW | CLEAN | NEAR / NEAR / NEAR | Keep as clean-near |
| JAK2 | 8BM2 | LOW | CLEAN | NEAR / NEAR / **INTERMEDIATE** | Cutoff-sensitive |
| PIK3CA | 4JPS | LOW | MODERATE | NEAR / NEAR / NEAR | **Not CLEAN** |
| PIK3CA | 5XGH | LOW | MODERATE | NEAR / NEAR / NEAR | **Not CLEAN** |
| PIK3CA | 4L2Y | not in V1 | CLEAN | NEAR / NEAR / NEAR | New clean-near |
| PIK3CA | 5DXT | INTERMEDIATE | CONFOUNDED | FAR / FAR / FAR | **Not INTERMEDIATE** |
| mTOR | 4JT5, 4JSX, 4JSV | LOW | CLEAN | NEAR / NEAR / NEAR | Keep as clean-near |
| AChE | 4EY6 | LOW | CLEAN | NEAR / NEAR / NEAR | Keep |
| BChE | 1P0M | LOW | CLEAN | NEAR / NEAR / NEAR | Keep; mutations now listed |
| F2 | 3SHC | LOW | MODERATE | NEAR / NEAR / NEAR | Peptide-length nuisance |
| F10 | 2Y5F | LOW | CLEAN | NEAR / NEAR / NEAR | Keep; R372E matched |
| PPARA | 6KAX | LOW | CLEAN | NEAR / NEAR / NEAR | Keep |
| PPARD | 5U46 | LOW | CLEAN | NEAR / NEAR / NEAR | Keep |
| PPARG | 7AWC | INTERMEDIATE | CONFOUNDED | FAR / FAR / FAR | **Not INTERMEDIATE** |

True “mainly change holo conformation” substitutions (CLEAN + NEAR at least at 4.0 and 4.5 Å):

1XKK, 4EI4, 8BM2 (with 5.0 Å caveat), 4L2Y, 4JT5, 4JSX, 4JSV, 4EY6, 1P0M, 2Y5F, 6KAX, 5U46.

3SHC is NEAR but MODERATE (hirudin/light length). 3PP0 is FAR. 5DXT and 7AWC are CONFOUNDED-FAR.

---

## 3. EGFR/HER2 vs JAK1/JAK2

**JAK1/JAK2 is cleaner.** This is not a TIE between those two pairs.

- JAK1 4EI4: same 854–1154 construct, same PTR 1034/1035, pocket NEAR at all cutoffs.
- JAK2 8BM2: same 840–1132 construct, same PTR 1007/1008, same paper as 8BXH; pocket NEAR at 4.0/4.5, INTERMEDIATE at 5.0; copy 2 vs 1.
- EGFR 1XKK is CLEAN-NEAR, but HER2 3PP0 is MODERATE-FAR (αC 8.61 Å, copy 2 vs 4, mutation/tag labeling conflict).

Family 1 as a whole: **JAK1/JAK2 vs PIK3CA/mTOR is a TIE.** See `KINASE_FAMILY_REPRESENTATIVE_AUDIT_V2.md`. Human must pick Option K1 (PIK3CA/mTOR) or Option K2 (JAK1/JAK2). EGFR/HER2 is not the representative.

JAK1/TYK2 cannot win: TYK2 alternatives were not re-verified in V2; 3LXP already has five engineered mutations.

---

## 4. PIK3CA / mTOR perturbation classes

No performance used.

| PDB | vs primary | Class | Tier | Pocket |
|---|---|---|---|---|
| 4L2Y | 4L23 | same-construct ligand-induced conformation (same paper, identical sequences, XXK vs X6K) | CLEAN | NEAR |
| 4JPS | 4L23 | same partner class + distal mutations (M232K/L233K) + p85 boundary change | MODERATE | NEAR |
| 5XGH | 4L23 | construct-trimmed p110+p85 | MODERATE | NEAR |
| 5DXT | 4L23 | construct/component/state change (no p85, 107–1068) | CONFOUNDED | FAR |
| 4JT5 | 4JT6 | same-construct ligand-induced (pp242 vs PI-103) | CLEAN | NEAR |
| 4JSX | 4JT6 | same-construct ligand-induced (Torin2 vs PI-103) | CLEAN | NEAR |
| 4JSV | 4JT6 | same-construct nucleotide vs inhibitor (ADP+Mg) | CLEAN | NEAR |

4JPS / 4JSX / 5DXT are historically exposed. They were not excluded for that reason. 4L2Y is preferred for PIK3CA because it is CLEAN, not because it is new.

---

## 5. AChE / BChE mutation and construct background

**Yes, they match.**

AChE 4EY6 vs 4EY7: same Cheung 2012 paper; `entity_poly` sequences identical (544 aa); same auth 2–543 / UniProt 33–574; same six disulfides; same glycosylation sites; dimer AB; no mutations. CLEAN-NEAR.

BChE 1P0M vs 4BDS: both have the four glycosylation-site mutations **N17Q, N455Q, N481Q, N486Q**. Sequences identical (531 aa). Same 1–529 / 29–557 construct. Mutations are not in the unified gorge pocket set. Residual: 1P0M has 120 unobserved atoms vs 5 on 4BDS; `entity.pdbx_mutation` is a residue list on 1P0M and `YES` on 4BDS (conflict of annotation completeness, not of sequence). CLEAN-NEAR.

---

## 6. F2 / F10 partner and mutation match

**F2 — partner comparable, construct not identical.**

Both 4UDW and 3SHC retain hirudin variant-2 with sulfotyrosine and the thrombin light chain. That matches the V4.2 RETAIN policy. Hirudin is shorter on 4UDW (11 aa, I:555–565) than on 3SHC (13 aa, I:53–65). Light-chain N-terminus is also longer on 3SHC. No engineered mutations on either thrombin. Pocket CA RMSD ~0.20 Å (NEAR). Completeness metric is internally conflicted because of the 147-insertion numbering. Tier: **MODERATE**. Partner identity matches; peptide length does not.

**F10 — mutation identity matches.**

Both 2JKH and 2Y5F have **R372E (auth A:150)**. It is not in the unified 4.5 Å pocket set. Light chain present both (RETAIN). Residual: heavy C-terminus 475 vs 468; light start 126 vs 127; 2JKH has a CA ion that 2Y5F does not list. Tier: **CLEAN**. V1’s “both mutated” claim is now residue-verified.

---

## 7. PPAR representative

**Yes. PPARA/PPARD remains cleaner than PPARG/PPARA.**

PPARA 6KAX: same 2020 iScience campaign and 200–468 LBD as 6LXA; no peptide; pocket 0.13–0.15 Å; CLEAN-NEAR.

PPARD 5U46: same 2017 PNAS LBD series as 5U3Q; ligand 7T1 not B7G; copy A of 5U46 vs frozen primary copy B of 5U3Q; pocket 0.49–0.52 Å; CLEAN-NEAR. Do not use 5U3Q copy A.

PPARG 7AWC: loses PG08-NL (V4.2 RETAIN on 9V8H) and is FAR (5.37–5.62 Å). CONFOUNDED. PPARG/PPARA is not selected to complete the family.

---

## 8. Recommended 3–4 representative pairs

Family slots after re-evaluation:

| Family | Pair | Status |
|---|---|---|
| 1 kinase / ATP-site | **TIE: PIK3CA/mTOR (K1) or JAK1/JAK2 (K2)** | human picks one |
| 2 cholinesterase gorge | **AChE/BChE** | selected |
| 3 S1 protease | **F2/F10** | selected |
| 4 PPAR LBD | **PPARA/PPARD** | selected |

That is three fixed pairs plus one kinase TIE (four after the human pick). EGFR/HER2 and PPARG/PPARA are not recommended as representatives.

---

## 9. Ablation designs (A1/B1, A2/B1, A1/B2, optional A2/B2)

A2/B2 is not a first experiment. It is allowed only after both single substitutions are accounted.

### Family 1 option K1 — PIK3CA/mTOR (lean if human wants cutoff-stable identical-sequence arms)

- A1/B1: **4L23 / 4JT6**
- A2/B1: **4L2Y / 4JT6**
- A1/B2: **4L23 / 4JT5** (4JSX is an equally CLEAN inhibitor TIE)
- optional A2/B2: 4L2Y / 4JT5
- No CLEAN-NEAR plus CLEAN-FAR ladder on PIK3CA: 5DXT is FAR but CONFOUNDED. Do not use it as A3.

### Family 1 option K2 — JAK1/JAK2

- A1/B1: **6N7A / 8BXH**
- A2/B1: **4EI4 / 8BXH**
- A1/B2: **6N7A / 8BM2** (copy A)
- optional A2/B2: 4EI4 / 8BM2
- Freeze 4EI4 copy (A vs AB) before prep. Apply the already-documented PTR fallback uniformly if those receptors are later prepared.

### Family 2 — AChE/BChE

- A1/B1: **4EY7 / 4BDS**
- A2/B1: **4EY6 / 4BDS**
- A1/B2: **4EY7 / 1P0M**
- optional A2/B2: 4EY6 / 1P0M
- 4EY7 extract remains AB; box-defining ligand remains E20 A:604.

### Family 3 — F2/F10

- A1/B1: **4UDW / 2JKH**
- A2/B1: **3SHC / 2JKH**
- A1/B2: **4UDW / 2Y5F**
- optional A2/B2: 3SHC / 2Y5F
- Hirudin I and FXa light L remain RETAIN.

### Family 4 — PPARA/PPARD

- A1/B1: **6LXA / 5U3Q** (5U3Q copy B)
- A2/B1: **6KAX / 5U3Q**
- A1/B2: **6LXA / 5U46** (5U46 copy A)
- optional A2/B2: 6KAX / 5U46

### Optional later FAR demonstration (not a family representative)

EGFR/HER2 if the human wants a large kinase conformational contrast:

- A1/B1: 3POZ / 3RCD
- A2/B1: 1XKK / 3RCD (CLEAN-NEAR)
- A1/B2: 3POZ / 3PP0 (MODERATE-FAR)
- Do not interpret the B-side as conformation-only.

---

## 10. Intended variable vs residual nuisance

| Swap | Intended unique variable | Residual nuisance |
|---|---|---|
| 4L23 → 4L2Y | holo conformation / cognate ligand (XXK vs X6K) | none material; p85 still REMOVE |
| 4JT6 → 4JT5 | holo conformation / cognate ligand (P2X vs X6K) | 3.45–3.60 Å series quality; mLST8 1–323 vs 1–326 annotation |
| 4JT6 → 4JSX | same | historically exposed (metadata only); same series quality |
| 6N7A → 4EI4 | holo conformation / cognate ligand | copy A vs AB if not frozen; PTR fallback if prepared |
| 8BXH → 8BM2 | holo conformation / cognate ligand | 1 vs 2 copies; 5.0 Å INTERMEDIATE; His-tag Gly |
| 4EY7 → 4EY6 | holo conformation / gorge ligand (GNT vs E20) | none material |
| 4BDS → 1P0M | holo conformation / gorge ligand (CHT vs THA) | more missing atoms on 1P0M; annotation completeness |
| 4UDW → 3SHC | holo conformation / catalytic ligand | hirudin and light-chain length |
| 2JKH → 2Y5F | holo conformation / cation-inhibitor ligand | C-term/light start; CA ion on 2JKH only |
| 6LXA → 6KAX | holo conformation / fatty acid (PLM vs EPA) | none material |
| 5U3Q-B → 5U46-A | holo conformation / agonist (7T1 vs 7UJ) | which crystal copy is used |
| 3POZ → 1XKK | holo conformation / ATP-site ligand | N-term tag on 1XKK |
| 3RCD → 3PP0 | large kinase conformation (αC / A-loop) | copy/dimer context; mutation vs tag labels; 03Q vs 03P; 3RCD quality |

---

## 11. What can be called conformation vs representation

**Receptor conformation sensitivity** (CLEAN, extra variables limited to cognate ligand and holo geometry):

- 4L23 vs 4L2Y
- 4JT6 vs 4JT5 or 4JSX (4JSV also, with nucleotide-class caveat)
- 6N7A vs 4EI4
- 8BXH vs 8BM2 (copy nuisance remains)
- 4EY7 vs 4EY6
- 4BDS vs 1P0M
- 2JKH vs 2Y5F
- 6LXA vs 6KAX
- 5U3Q-B vs 5U46-A
- 3POZ vs 1XKK

**Broader receptor-representation sensitivity** (construct, partner, mutation, or activation-state confounders):

- 3RCD vs 3PP0 (FAR + copy/dimer + annotation conflict)
- 4L23 vs 4JPS (distal mutations + p85 boundary)
- 4L23 vs 5XGH (trimmed construct)
- 4L23 vs 5DXT (no p85 + FAR)
- 4UDW vs 3SHC (peptide-length; still NEAR)
- 9V8H vs 7AWC (peptide removed + FAR)

Do not write these two mechanisms as one.

---

## 12. Did old performance information affect any recommendation?

**No V2 recommendation used docking performance.**

Opened: receptor registry, component policy, `analysis_freeze.yaml` identity/component/alternative-ID fields, V1 structure metadata, RCSB/mmCIF/validation/papers.

Not opened: any AUROC table, summary_min, Top-K/EF, LEGACY or V4.2 scores, Phase 5 RMSD as a selector, seed-sensitivity metrics, GNINA/RTMScore, deposited alternative score CSVs named in the manifest.

Residual listing-overlap risk (disclosed, not used as a selector):

- IDs 4JPS, 5DXT, 4JSX appear in `analysis_freeze.yaml` `allowed_pdb` and in `alternative_receptor_manifest.csv`.
- V1 discussed most of the other alternatives. Those are marked `HISTORICALLY_EXPOSED_CANDIDATE=YES`.
- 4L2Y is the main CLEAN candidate that is **not** historically exposed.
- 4JSX was **not** dropped because it is an old name. It remains CLEAN-NEAR.
- V1’s choice of EGFR/HER2 “to avoid a second kinase pair” is rejected as the Family 1 rule. That rejection is from the V2 priority order, not from performance.

If any later reader remembers an old AUROC for 4JPS/4JSX/5DXT, that memory is a conflict-of-information risk. The V2 tables do not contain those numbers and were not ranked by them.

---

## Pocket definition freeze (candidates, not a silent pick)

Unified pocket = protein residues with any heavy atom ≤ *d* Å from any heavy atom of the **primary** cognate ligand, using the frozen primary auth copy.

| Target | Primary ligand | n residues 4.0 / 4.5 / 5.0 |
|---|---|---|
| EGFR | 03P | 21 / 25 / 27 |
| HER2 | 03P | 19 / 23 / 26 |
| JAK1 | KEV | 19 / 20 / 22 |
| JAK2 | C87 | 16 / 19 / 25 |
| PIK3CA | X6K | 13 / 18 / 21 |
| mTOR | X6K | 11 / 17 / 17 |
| AChE | E20 | 15 / 19 / 21 |
| BChE | THA | 11 / 15 / 16 |
| F2 | N6L | 20 / 22 / 22 |
| F10 | BI7 | 17 / 19 / 21 |
| PPARG | BRL | 16 / 22 / 23 |
| PPARA | EPA | 19 / 24 / 26 |
| PPARD | 7UJ | 21 / 25 / 28 |

Human may freeze one *d* later. Do not freeze 4.5 Å because it makes 8BM2 NEAR.

---

## Stop

V2 files are written. No receptor prep, no redock, no Vina/GNINA/RTMScore, no AUROC, no further PDB search. Candidate set and fact tables are for human confirmation.
