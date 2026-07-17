# Target Journal Research (Non-OA / Hybrid Subscription Path)

> PaperSpine `target-journal-research` · ARS venue-fit · 2026-07-17  
> Constraint: **prefer subscription / hybrid with optional APC** (author may choose non-OA). Avoid Gold-OA-only venues.

## Official Requirements Snapshot

### 1) Journal of Molecular Modeling (Springer) — **Hybrid**

| Dimension | Requirement | Implication for Our Paper |
|---|---|---|
| Publishing model | Hybrid (subscription possible) | Non-OA path OK |
| Scope URL | https://link.springer.com/journal/894/aims-and-scope | **Updated hard gates** |
| Commercial tools | “Studies that use commercially available techniques … for docking, ADMET or MD are **discouraged**” | **Glide XP narrative is a desk-reject risk** |
| Pre-dock MD | Protein for docking must be equilibrated **≥500 ns** | Not in current plan |
| Ensemble docking | **≥3 conformations** (multi-crystal or ≥500 ns MD) for consensus poses | Current main funnel is **single-state 9DKB + 7ALV** |
| Pose MD | **3× ≥300 ns** or **1× ≥500 ns**; RMSD-only short MD discouraged | Planned **50–100 ns** fails |
| New methodology | Welcome if extensively validated | Our “funnel” is system-level, not new force-field |

**Fit verdict: ❌ Not recommended as primary target under current data budget.**  
Previous outline ranking JMM #1 is **outdated relative to current Aims & Scope**. Meeting JMM would require a **methods-heavy rebuild** (open-source docking emphasis, 500 ns-scale MD, multi-conformation consensus)—essentially a different paper.

### 2) Journal of Computer-Aided Molecular Design (Springer) — **Hybrid**

| Dimension | Requirement | Implication |
|---|---|---|
| Publishing model | Hybrid | Non-OA path OK |
| Scope | Theory + application of CADD; ML, docking, MD, chemoinformatics | Funnel + audit fits |
| Experimental priority | Priority to papers with experimental validation; novel methods may use **retrospective validation** | Our 8973 track + y-scramble/AD = retrospective defense |
| Novelty bar | Method/protocol clarity over disease marketing | Emphasize asymmetric design + Pareto≠nomination |
| Typical length | Full research article | Need complete Methods + SI |

**Fit verdict: ✅ Best primary non-OA fit** if claims stay methodological and MD ≥50–100 ns with redock is completed. Still expects strong reproducibility and honest limitations.

### 3) Chemical Biology & Drug Design (Wiley) — **Hybrid / subscription with OA option**

| Dimension | Requirement | Implication |
|---|---|---|
| Publishing model | Subscription journal with OA option (APC only if OA chosen) | Non-OA path OK |
| Scope | Chemical biology + drug design; mechanistic insight; concept papers | Disease dual-node story welcome |
| Experimental pressure | Often prefers ligand/biology insight; “discourage additional experiments” culture helps but pure in silico still scrutinized | Frame as **hypothesis generation** + comorbidity pharmacology of canagliflozin |
| Speed | Marketed as relatively fast first decision | Good secondary |

**Fit verdict: ✅ Strong secondary** if MD + interaction story is polished and canagliflozin mechanism boundary is crystal-clear.

### 4) Journal of Molecular Graphics and Modelling (Elsevier) — **Hybrid**

| Dimension | Notes | Implication |
|---|---|---|
| Model | Hybrid | Non-OA possible |
| Bar | Historically selective (~low acceptance); graphics/modeling focus | Possible tertiary; less ideal for “clinical library funnel” narrative than JCAMD |

**Fit verdict: ⚠️ Tertiary backup.**

### Venues to avoid for this package

| Venue | Why |
|---|---|
| JCIM / J Med Chem / Nat Commun | Height / wet-experiment expectation exceeds current package |
| Gold-OA-only mega-journals as first choice | Violates non-OA preference |
| Molecular Diversity | Often wants richer discovery packaging; hybrid but less protocol-centric than JCAMD |

## Venue Writing Contract (recommended: JCAMD)

- **Rewards**: reproducible protocol; retrospective enrichment; clear negative controls (colchicine, PAINS demotion); open discussion of claim boundaries.
- **Filters out**: “first dual inhibitor”; commercial black-box without open alternatives; metric dumps without contribution mapping.
- **Expected novelty**: **system/protocol novelty under asymmetric data**, not a new scoring function.
- **Expected evidence**: redock ≤2 Å; 8973 AUC/EF; modules A–F numbers; MD RMSD/contacts + qualitative MM-GBSA; no Ki claims.
- **Claim strength**: hypothesis-generating computational nomination.

## Recent Exemplars to Imitate (logic, not text)

| Paper | Why learn from it |
|---|---|
| Suo et al. Nat Commun 2025 (URAT1 cryo-EM + CHARMM-GUI POPC MD) | Gold-standard membrane MD Methods language (cite, do not try to match 1 μs×5 for JCAMD) |
| Eurycoma Nat Commun 2025 dual anti-gout | Complementary wet dual route — cite as **different paradigm** |
| JCAMD 2025–26 docking reproducibility / case studies | Tone: protocol + failure modes |
| PLK1/NLRP3 asymmetric VS (closest method cousin) | Explicit Differentiation paragraph required |

## Fit Assessment Summary

| Criterion | JMM | JCAMD | CBDD |
|---|---|---|---|
| Non-OA possible | ✅ | ✅ | ✅ |
| Match current MD budget (50–100 ns) | ❌ | ✅ | ✅ |
| Match Glide XP usage | ❌ discouraged | ⚠️ OK if open scripts + versions | ⚠️ OK |
| Match funnel narrative | ❌ prefers hardcore MD | ✅ | ✅ disease angle |
| Desk-reject risk now | **High** | Medium | Medium–High if overclaim |
| Recommended rank | Drop as #1 | **#1** | **#2** |
