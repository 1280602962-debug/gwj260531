# Section Blueprints — Option A (English manuscript)

**Target voice:** JCIM-style methods + negative result + prospective shortlist.  
**Word budget (flash tier):** ~3500–5000 main text (adjust later).

---

## Title (working)

Computational Isoform-Selectivity Filters Fail on JNK1/2/3 Near-Identical ATP Sites: An End-to-End Screening Pipeline for JNK-Family Binder Enrichment

*(Soft alt if IC50 succeeds: “…and Prospective Enzymatic Evaluation of Two Commercial-Library Candidates”)*

---

## Abstract (blueprint)

1. **Context:** JNK1 bias desirable (CC-90001/E1); ATP sites nearly identical.  
2. **Gap:** Cheap Δsel / residue heuristics / ML selective labels used as purchase criteria without hard negative tests.  
3. **What we did:** ML family gate → Glide → ADMET → MD QC → purchase; benchmark autopsy on literature panel.  
4. **Computational finding:** Δsel direction accuracy below usability; Gly87 non-discriminative; ML selective F1=0; purchase **decoupled**.  
5. **Wet-lab status:** Pre-registered JNK1/2/3 IC50 on **690, 2231** + E1/CC-90001 *(fill numbers)*.  
6. **Takeaway:** Family enrichment pipeline usable; isoform purchase filters not.

**Citations:** CIT-A2, CIT-B2, DAT-E1.

---

## 1 Introduction

| Paragraph | Content | Evidence / cites |
|-----------|---------|------------------|
| 1.1 Biology | JNK in fibrosis/inflammation; why JNK1 bias | CIT-A1–A3 |
| 1.2 Selectivity problem | Isoform vs kinome; ATP conservation | CIT-A5, CIT-B3 |
| 1.3 Prior computational approaches | Docking Δ, IFP, FEP, ML — mixed success on other kinases | CIT-B1–B4 |
| 1.4 This work | Pipeline + **negative benchmark as core**; wet-lab enrichment not selectivity discovery | confirmed_contribution |

**Forbidden here:** “We discovered selective JNK1 inhibitors.”

---

## 2 Methods

| Subsection | Must include | Source |
|------------|--------------|--------|
| 2.1 Data | ChEMBL curation; paired set; benchmark list | CIT-C1; DAT tables |
| 2.2 ML | XGBoost pActivity; p_family≥6.0 as **recall** gate; selective classifier attempted | CIT-C6; training_report |
| 2.3 Docking (Selection) | Historical Glide XP triage; PDB panel; state institutional Glide as triage only | project report; SOFTWARE_LICENSE_NOTE |
| 2.4 Selectivity metrics tested | Δsel_dock; Gly87; ML selective label — **not purchase gates** | DAT-E1 |
| 2.5 ADMET / shortlist | 157→25→16 | project report |
| 2.6 MD QC | RMSD + hinge; pass_md_overall | MD docs |
| 2.7 Purchase rule | **690** (grade A anchor) + **2231** (bias hypothesis, grade C risk) | C11 |
| 2.8 Open-source confirmation | Vina/Gnina multi-seed; unrestrained MD replicas | C2/C3 |
| 2.9 Enzymatic assays | JNK1/2/3 IC50 + C4 lock | C4 |

---

## 3 Results

### 3.1 Funnel (one figure)
4979 → 157 → 25 → 16 → 2 purchased.  
**Message:** selectivity tags not hard gates.

### 3.2 RQ-C: Selectivity-method autopsy (main table — C5)
Δsel FAIL; Gly87 FAIL; ML F1=0; family ML recall OK.  
**Explicit sentence:** “These filters were not used for purchase.”

### 3.3 Chemotype novelty (C1)
ECFP4 maxTc ~0.23 vs literature refs; distant from E1/CC-90001; caveat pharmacophore.

### 3.4 Pose consensus (C2)
Vina multi-seed RMSD consensus; Glide ranks historical. Soften if consensus fails.

### 3.5 Purchase rationale (690 + 2231)
690 = pose-QC / pan-leaning RQ-A anchor (grade A).  
2231 = strongest MD bias hypothesis for RQ-B despite grade C / overall MD fail.  
Selectivity filters **not** used as purchase gates.

### 3.6 Assay risk filters (C7)
PAINS/physchem clear/alert table.

### 3.7 Enzymatic IC50 (when available)
C4 table; RQ-A primary; RQ-B only if SI rule met.  
Controls: E1 direction; CC-90001 multi-isoform.

---

## 4 Discussion

| Point | Stance |
|-------|--------|
| Negative result value | Align with Kinase-Bench / FEP selectivity literature: cheap filters insufficient on JNK |
| Pipeline usefulness | Family enrichment still justified if RQ-A hits |
| Why not JMC-style selective discovery | No SAR/kinome/cell; n=2; 2231 not bought |
| MD ≠ selectivity | E1/SP600125 counterexamples |
| Limitations | No kinome; C2 open prep simplified; single-replica MD historically |
| Future | Optional 2231 buy; open MD replicas; kinome if budget |

---

## 5 Conclusions

1. Δsel / Gly87 / ML selective labels **fail** as JNK isoform purchase criteria on our benchmark.  
2. End-to-end pipeline yields pose-credible family shortlist.  
3. Wet-lab IC50 calibrates enrichment *(result-dependent)*.  
4. Isoform selectivity remains an experimental problem for JNK.

---

## Figures / Tables plan

| Item | Content |
|------|---------|
| Fig 1 | Funnel |
| Fig 2 | Δsel direction / confusion (benchmark) |
| Fig 3 | Gly87 distances (non-discriminative) |
| Fig 4 | Optional: 690/2157 poses |
| Table 1 | C5 autopsy summary |
| Table 2 | C1 novelty |
| Table 3 | Purchase vs 2231 (C11) |
| Table 4 | IC50 / SI (C4) |

---

## Contribution check (before drafting Results prose)

- [x] Motivation Option A locked  
- [x] C5 table exists  
- [x] C1 table exists  
- [ ] IC50 numbers (pending wet-lab)  
- [ ] License-safe docking Methods sentence chosen  
