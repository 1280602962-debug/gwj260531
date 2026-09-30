# V4.2 alternative-receptor selection (structure-only freeze)

RUN_ID = UNIFORM_RERUN_V4_2_20260921  
Status = DESIGN_ONLY  
Docking = NOT STARTED  
AUROC / summary_min / Top-K / EF / LEGACY or V4.2 performance = NOT READ, NOT USED

This note freezes a receptor-sensitivity plan from experimental holo structures only. The 14 primary receptors are not modified. Predicted structures are excluded. Production alternative docking is forbidden until this file and `00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE.yaml` are manually confirmed.

## Method (what was and was not used)

Used:

- V4.2 frozen primary identities (PDB, UniProt, auth copy, cognate CCD, component policy)
- RCSB Search + GraphQL for Homo sapiens experimental entries of the 14 UniProt accessions
- Domain coverage against UniProt functional-domain ranges
- Non-additive cocrystal ligands as the holo requirement
- Partner/component comparison against the V4.2 extract policy
- Sequence-mapped pocket superposition on a construct-stratified representative subset (pocket = residues ≤4.5 Å from the primary cognate ligand)

Not used:

- Any AUROC, summary_min, Top-K, EF, or other screening metric
- V4.2 or LEGACY docking scores
- Phase 5 cognate redocking RMSD
- “Looks easier to prepare” as a reason to change PDB
- The prior manifest phrases `none_better_under_same_rule` / `later_JH1_entries_worse_reso` as decision rules

A prior V4.1/V4.2 `alternative_receptor_manifest.csv` already listed 4JPS, 5DXT, 4JSX and several declined IDs. Those IDs were re-opened as structures and re-classified. Overlap is marked `PRIOR_MANIFEST_OVERLAP`. That listing is a conflict-of-information risk, not a performance input; deposited score files named in that manifest were not opened.

Hard availability filter (all 14 targets): human + experimental + same UniProt + domain coverage ≥0.70 + ≥1 non-additive ligand. This produced 2318 human experimental records and 1718 non-primary hard-pass holos. Individual pocket superposition was done on a representative subset, not on every hard-pass entry.

## 1. Targets with a genuine LOW_PERTURBATION alternative

LOW means: same UniProt, same functional domain, same orthosteric/canonical site, same or highly similar construct, same partner/component background as the V4.2 primary extract, no declared key pocket mutation that defines a different biological state, experimental holo.

| Target | Primary | Frozen LOW alternative | Why LOW |
|---|---|---|---|
| EGFR | 3POZ | **1XKK** | Kinase-only, ATP site, lapatinib (FMM). Seq-align pocket CA RMSD 0.565 Å vs 3POZ. No extra partner. DOI 10.1158/0008-5472.CAN-04-1168 |
| HER2 | 3RCD | **3PP0** | One of the few human HER2 kinase holos that is not an exon-20 insertion or covalent complex. Ligand 03Q, same TAK-related chemical series family as 3POZ/3RCD (JBC 2011). Pocket complete. Surface engineered mutations recorded; not used as a docking-quality argument |
| JAK1 | 6N7A | 4EI4 / 4IVC / 5E1E | JH1-only ATP-site holos; pocket CA RMSD 0.68–0.97 Å |
| JAK2 | 8BXH | 8BM2 / 4IVA | JH1-only ATP-site holos; pocket CA RMSD 0.58 / 0.86 Å |
| TYK2 | 3LXP | 3NZ0 / 4GVJ | JH1 ATP-site holos; pocket CA RMSD ~2.1–2.3 Å. Same domain, not JH2 |
| PIK3CA | 4L23 | **4JPS** (also 5XGH) | p110α + p85 present, ATP/catalytic site. 4JPS pocket bb RMSD 0.609 Å. 5XGH same partner, 0.577 Å, lower resolution (2.97 Å) |
| mTOR | 4JT6 | **4JSV** (also 4JSX, 4JSP, 4JT5) | Same Nature 2013 mTORδN–mLST8 crystallization series as 4JT6; only the bound nucleotide/inhibitor changes. Pocket bb RMSD 0.28–0.72 Å |
| AChE | 4EY7 | **4EY6** (also 4EY5, 4M0E) | Same recombinant human AChE construct and Cheung 2012 paper as 4EY7; galantamine vs donepezil in the catalytic gorge. Pocket bb RMSD 0.317 Å |
| BChE | 4BDS | **1P0M** | Human BChE catalytic gorge; same engineered-mutation list as 4BDS. Pocket bb RMSD 0.177 Å |
| F2 | 4UDW | **3SHC** (also 2ZFF, 1H8D, 5AFY) | Catalytic site + **hirudin retained**, matching V4.2 RETAIN policy. Pocket bb RMSD 0.16–0.28 Å |
| F10 | 2JKH | **2Y5F** (also 2Y5G) | Same “Factor Xa – cation inhibitor complex” series as 2JKH. Seq-align pocket CA RMSD 0.141 Å |
| PPARA | 6LXA | **6KAX** (also 6KB3) | Same human PPARα LBD crystallization campaign as 6LXA; no extra peptide. Pocket CA RMSD 0.147 / 0.470 Å |
| PPARD | 5U3Q (copy B) | **5U46** (also 2ZNP after sequence remap) | Same PNAS 2017 LBD series as 5U3Q; GW501516 vs 7UJ. Pocket bb RMSD 0.505 Å. Copy A of 5U3Q is **not** an alternative receptor |

## 2. Targets with only INTERMEDIATE or HIGH alternatives as the realistic next step

**PPARG (9V8H)**  
Primary retains peptide PG08-NL plus rosiglitazone (BRL). The closest binary experimental holo is **7AWC** (BRL, no peptide): INTERMEDIATE_PERTURBATION (component background changed). RXR heterodimer holos (1FM6, 3DZY) are HIGH. There is no LOW alternative that keeps the peptide background.

**EGFR additional holos** (1M17, 2ITY, 4HJO, 5UGC)  
Same kinase/ATP site, but sequence-aligned pocket CA RMSD remains 7–12 Å, consistent with a different kinase conformational state. They are INTERMEDIATE receptor-state comparators, not LOW conformation-only swaps.

**PIK3CA 5DXT**  
p110α + GDC-0326 without the p85 niSH2 present in 4L23: INTERMEDIATE (component).

**AChE 1B41**  
Fasciculin-2 complex: HIGH (partner).

**HER2 8VB5 / 8U8X / 7JXH / 7PCD**  
Exon-20 insertion or covalent kinase complexes: HIGH. They must not be mixed with 3PP0 in a LOW interpretation.

**mTOR cryo-EM mTORC1/C2**  
Same UniProt and kinase domain coverage, but Raptor/Rictor/Rag/DEPTOR backgrounds: HIGH. Not interchangeable with the 4JT6 mTORδN–mLST8 crystal series.

## 3. What is not suitable as receptor-conformation sensitivity

- AlphaFold or any predicted model
- 5U3Q copy A (same experiment as primary copy B; V4.2 already froze A as LEGACY_ONLY)
- EGFR / HER2 extracellular-domain structures (wrong domain)
- TYK2 / JAK JH2 pseudokinase structures (wrong domain; already removed by coverage)
- FRB-only rapamycin–FKBP mTOR fragments without the kinase site used by 4JT6
- Apo-only entries (fail holo rule)
- Mixing LOW and HIGH into one “alternative receptor” average

HER2 is **suitable** only through 3PP0. The kinase-domain holo set is tiny (5 non-primary hard-pass entries). That scarcity is a structural fact, not a performance fact.

## 4. Cross-family pairs recommended for the formal one-at-a-time experiment

One pair per structural family, so the set is representative rather than kinase-heavy:

| Family | Pair | Status |
|---|---|---|
| Kinase / ATP-site RTK | **EGFR/HER2** | SELECTED |
| AChE/BChE catalytic gorge | **AChE/BChE** | SELECTED |
| Serine protease | **F2/F10** | SELECTED |
| PPAR nuclear-receptor LBD | **PPARA/PPARD** | SELECTED |

Evaluated, suitable for LOW, **not selected** (would over-represent kinases):

- JAK1/JAK2 and JAK1/TYK2: both arms have LOW JH1 holos
- PIK3CA/mTOR: both arms have unusually clean same-series LOW holos (4JPS / 4JSV)

If a later amendment wants a second kinase-adjacent pair, PIK3CA/mTOR is the cleanest construct-matched add-on. It is **not** in the primary four-pair freeze.

PPARG/PPARA is **not** selected because PPARG has no LOW alternative under the V4.2 peptide-retained background.

## 5. Why each selected pair is selected

**EGFR/HER2**  
Represents the protein-kinase ATP-site family. EGFR has many human kinase holos; 1XKK is a kinase-only ATP-site holo with low pocket CA RMSD. HER2 almost has no other fair kinase holo except 3PP0. Selecting this pair tests conformation sensitivity on an RTK pair without touching JAK or PI3K.

**AChE/BChE**  
Represents the cholinesterase gorge family. 4EY6 is the same recombinant human construct and paper as 4EY7. 1P0M matches the 4BDS engineered-mutation background. Partner backgrounds stay empty/monomer as in V4.2.

**F2/F10**  
Represents S1 serine proteases. 3SHC keeps hirudin, which V4.2 explicitly RETAINs on 4UDW. 2Y5F is the same Factor Xa cation-inhibitor series as 2JKH, so the light-chain/cation background is held as constant as the crystal series allows.

**PPARA/PPARD**  
Represents NR LBDs. Both primaries are peptide-free LBD holos (unlike 9V8H). 6KAX and 5U46 come from the same experimental campaigns as 6LXA and 5U3Q. Helix-12 motion is recorded as a known NR property, but these two alternatives are the lowest construct perturbation available.

## 6. Why each unselected pair is unselected

| Pair | Why not in the primary four |
|---|---|
| JAK1/JAK2 | LOW alternatives exist; omitted to avoid a second kinase ATP-site pair |
| JAK1/TYK2 | Same as above; TYK2 JH1 holos exist (3NZ0) |
| PIK3CA/mTOR | LOW alternatives exist and are construct-clean; omitted from the primary four for family balance. Optional add-on only |
| PPARG/PPARA | PPARG lacks a LOW alternative that keeps PG08-NL. 7AWC would force an INTERMEDIATE interpretation mixed with PPARA LOW |

No pair was dropped because a historical docking table was better or worse. Those tables were not read.

## 7. Could existing performance information have influenced selection?

**Declared risk: YES, a residual listing-overlap risk exists; performance files were not used.**

- I opened `alternative_receptor_manifest.csv`. It names 4JPS, 5DXT, 4JSX and declined IDs (4M0E, 5AFY, 7AWC, 8EXL, …) and contains result-flavored phrases such as `none_better_under_same_rule`. Those phrases were **not** adopted.
- I did **not** open the score CSVs listed in that manifest.
- I did **not** read Phase 5 RMSD audits, Phase 7 affinities, AUROC, summary_min, Top-K, or LEGACY comparison tables.
- Mitigation:  
  - EGFR alt is **1XKK**, not a recycled manifest ID  
  - mTOR alt is **4JSV**, not 4JSX (4JSX remains an equally LOW same-series structure)  
  - AChE alt is **4EY6**, not 4M0E (4M0E is also LOW and was re-classified as such)  
  - F2 alt is **3SHC**, not the previously declined fragment 5AFY (5AFY is also LOW by construct)  
  - PPARG 7AWC remains INTERMEDIATE and is **not** used to enter PPARG/PPARA  
  - PIK3CA 4JPS **is** a selected optional/evaluated LOW structure because it is the p85-matched crystal, and that overlap is explicitly marked

If this residual overlap is unacceptable, replace PIK3CA’s evaluated alt with 5XGH only, and keep 4JPS out of any later optional pair.

## Frozen one-receptor-at-a-time designs (no docking yet)

Notation: pair A/B, primary A1/B1.

### EGFR/HER2

- Primary: 3POZ / 3RCD  
- First swap: **1XKK / 3RCD** (A2/B1)  
- Second swap: **3POZ / 3PP0** (A1/B2)  
- A2/B2 (1XKK / 3PP0) only after both single swaps are accounted, and only if later authorized

### AChE/BChE

- Primary: 4EY7 / 4BDS  
- First: **4EY6 / 4BDS**  
- Second: **4EY7 / 1P0M**  
- A2/B2 later only

### F2/F10

- Primary: 4UDW / 2JKH  
- First: **3SHC / 2JKH**  
- Second: **4UDW / 2Y5F**  
- A2/B2 later only

### PPARA/PPARD

- Primary: 6LXA / 5U3Q  
- First: **6KAX / 5U3Q**  
- Second: **6LXA / 5U46**  
- A2/B2 later only  
- 5U3Q copy remains B

Preparation, if later authorized, must reuse the V4.2 recipe (PDBFixer heavy-only, Meeko 0.7.1, no Reduce2, no forbidden Meeko flags, no new box-hunting). One target is replaced per design. Boxes are redefined only from the alternative structure’s own cognate ligand by the frozen pad rule, not from Vina scores.

## What this freeze does not do

- Does not dock  
- Does not compute AUROC  
- Does not modify the 14 primaries  
- Does not add a PDB to the production 14  
- Does not authorize GNINA/RTMScore  
- Does not treat predicted models as receptors  
