# PAPER_OUTLINE_TRACK_A — DualFourClass diagnostic / benchmark note

**Working title:** Why pooled dual-target docking metrics look random — and why public data rarely supports a strict four-class test  

**Default venue:** Molecular Informatics / JCAMD (hybrid). JCIM Article only after approved K=4 docking (J2–J5).

## 1. Introduction
- Dual-target design vs docking decision task  
- Problem: mean/min fusion reported with pooled AUROC  
- Contribution preview: directional metric + baseline discipline + supply audit  

**Pointers:** `docs/JCIM_ROUTE_ASSESSMENT_V1.md` §0–3; Stage M verdict  

## 2. Task formalization (Fig1)
- Four classes; θ=6 vs strict 6.5/5.5; untested ≠ inactive  
- Algebra: pooled ≈ weighted average of D/A and D/B  

## 3. Measurement audit on exploration pairs (Fig2–3, Fig5)
- EGFR/HER2 unified RDKit EH110: docking ≤ cLogP  
- PIK3CA/mTOR: docking > volume (note LigPrep until J2)  
- Prep sensitivity (M4)  

**Data:** `track_a_starter_v0/tables/eh110_*`; `stage_m_v0/tables/*`  

## 4. Public-data supply ceiling (Fig4) — C2
- J0: 49 pairs audited; 3 non-metal strict-Y  
- EGFR B_only_strict = 7 as hard ceiling  

**Data:** `tables/j0_strict_label_supply.csv`  

## 5. Implications for evaluation practice
- Report directional AUROCs + trivial baselines  
- Choose pairs by supply, not only literature heat  
- Freeze ligand prep  

## 6. Limitations & outlook
- K=2 docked pairs today; K=4 requires user-approved J2–J5  
- No wet validation; no method win claimed  

## TODO before submission
- [ ] Draft figure PDFs from CSV seeds  
- [ ] Decide Mol. Inf. vs wait for JCIM K=4  
- [ ] If JCIM: wait for approval → J2 PM48 RDKit + J3 receptor freeze (do not self-start)
