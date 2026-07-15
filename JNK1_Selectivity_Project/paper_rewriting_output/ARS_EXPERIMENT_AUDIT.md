# ARS Codex — Experiment Selection & Computational Gap Audit

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan + validate (advisory; non-rerunnable computational campaign)
- Origin Date: 2026-07-15
- Verification Status: ANALYZED
- Version Label: ars_exp_audit_v1
- Package: Imbad0202/academic-research-skills-codex @ skills/academic-research-suite
- Project materials: JNK1_Selectivity_Project (local evidence authoritative)
- Wet-lab constraint: **only JNK1/JNK2/JNK3 enzymatic IC50** (no kinome panel, no cell assay assumed)

---

## 1. Research Question Under Test (clarified)

| Layer | Question | Can current purchase answer it? |
|-------|----------|----------------------------------|
| RQ-A | Does the pipeline enrich **JNK-family binders** (any isoform IC50 activity)? | **Yes — primary** (690 + 2157 + E1/CC-90001) |
| RQ-B | Do purchased molecules show **JNK1 isoform preference** (≥3–10× vs JNK2 and JNK3)? | **Weakly / poorly powered** — design not optimized for this |
| RQ-C | Are computational isoform-selectivity predictors trustworthy for purchase? | **Already answered computationally (negative)**; wet-lab is calibration, not required for this claim |
| RQ-D | Kinome selectivity (JNK vs p38/ERK/other)? | **No** under stated constraint |

**ARS plan verdict:** Purchase set is **acceptable for RQ-A**, **inadequate as a prospective test of RQ-B**, **sufficient calibration for RQ-C narrative**, **silent on RQ-D**.

---

## 2. Final Molecule Selection — Objective Check

### Locked IDs

| Role | Library ID | HIT / source | MD overall | Why selected (project) |
|------|------------|--------------|------------|------------------------|
| New #1 | **690** | HIT103871685 | Pass A | G1, Tier1′, pan-leaning activity anchor |
| New #2 | **2157** | HIT101201113 | Pass A | G1, MD bias #2 (hinge J1≫J3) |
| +ctrl | E1 | in hand | Pass B | Literature JNK1-biased |
| +ctrl | CC-90001 | in hand | Pass A | Near-pan / functional JNK1 bias |

ID cross-check vs `md_shortlist_final.csv` / Table 27: **PASS (no mis-ID)**.

### Decision matrix (ARS-style)

| Criterion | 690 | 2157 | Best alternative not bought (2231) |
|-----------|-----|------|-------------------------------------|
| Pose credibility (pass_md_overall) | ✅ | ✅ | ❌ grade C |
| Family-activity prior (score/MMGBSA/ADMET) | ✅ | ✅ | ✅ (stronger score) |
| Prospective JNK1-bias hypothesis strength | Low (pan hinge) | Medium-weak (Δsel_dock **−1.05**) | **Highest** (hinge 0.91/≈0) |
| Matches “selectivity discovery” claim | ❌ | ❌ | Still not proof, but better hypothesis test |
| Fits “family enrichment + assay calibration” | ✅ | ✅ | Optional exploratory |

### Verdict on “有没有问题”

| Judgment | Detail |
|----------|--------|
| **No fatal experimental-design error** for **family-binder enrichment + assay calibration** | Two grade-A G1 candidates + two literature positives is a coherent minimal panel |
| **Yes — strategic mismatch** if the scientific question is still “JNK1 selectivity” | You did **not** buy the strongest MD-bias molecule; 2157’s docking Δsel is **anti-JNK1** |
| **Yes — statistical power defect** | n_new = 2 → cannot distinguish enrichment from luck; multiple isoform comparisons uncorrected |
| **Yes — claim hazard** | Same hinge-binder chemotype space → reviewers will ask kinome questions you cannot answer under constraint |

**Bottom line (ARS):** Selection is **defensible for RQ-A**, **not defensible as the best set for RQ-B**. Proceed only if success criteria are pre-registered around family activity, not selectivity discovery.

---

## 3. Defect Register (objective)

### D1 — Design / sampling defects

1. **n=2 new molecules** — enrichment rate undefined; one hit looks like luck.
2. **No purchased negative anchor (G4)** — cannot show “decoy stays inactive.”
3. **No purchase of best bias hypothesis (2231)** — strongest asymmetry signal left untested.
4. **Both new molecules are G1** — G1 vs G2 chemotype-strategy comparison **abandoned**.
5. **Single-replica MD** used for ranking — stochastic pose risk unquantified.

### D2 — Endpoint / claim defects

6. **Isoform selectivity underpowered** — 2 compounds × 3 isoforms; no pre-specified SI threshold analysis plan in purchase docs.
7. **Kinome silent** (constraint) — common hinge scaffolds → polypharmacology prior not zero.
8. **MD cannot adjudicate selectivity** — already shown by E1 mis-rank / SP600125 low-hinge activity; must not interpret post-hoc MD–IC50 agreement as validation of MD selectivity.

### D3 — Computational evidence defects (already in archive)

9. Selectivity predictors failed (Δsel 43%; Gly87 null; ML F1=0) — good honesty, but means purchase was **not** selectivity-driven.
10. 2157 Δsel_dock negative — if it later shows JNK1 bias, that is **discordant** with docking, not a predicted success.
11. Extended MD (2231) used ligand restraints — limited generalizability even for the unbought molecule.

### Fallacy scan (11/11 checked — advisory)

| # | Fallacy | Status for current plan |
|---|---------|-------------------------|
| 1 | Simpson’s paradox | Watch: pooling isoforms without per-isoform reporting |
| 2 | Ecological fallacy | N/A |
| 3 | Restriction of range | **RISK**: shortlist already heavily filtered → IC50 hits overstate pipeline |
| 4 | Multiple comparisons | **RISK**: 2×3 IC50 matrix + SI ratios without correction |
| 5 | HARKing | **HIGH RISK** if selectivity claimed after seeing data |
| 6 | p-hacking / endpoint switching | Pre-register: primary = any isoform IC50 < X µM |
| 7 | Causation from correlation | Do not claim MD caused activity |
| 8 | Underpowered null | n=2 cannot prove “pipeline fails” either |
| 9 | Overfitting narrative | Avoid fitting story to whichever isoform looks better |
| 10 | Dichotomania | Report continuous IC50/SI, not only “selective/not” |
| 11 | Missing controls | Mitigated by E1/CC-90001; missing SP600125 pan tool & G4 |

---

## 4. Pre-registered success criteria (required before assay)

| Outcome | Pre-define now |
|---------|----------------|
| **Primary success (RQ-A)** | ≥1 of {690,2157} shows IC50 ≤ 10 µM (or lab’s validated LOD) on **any** of JNK1/2/3 |
| **Secondary (directional)** | SI = IC50_off / IC50_JNK1; “preference” only if SI ≥ 3 vs **both** JNK2 and JNK3 |
| **Assay validity** | E1: JNK1 preference direction roughly retained; CC-90001: measurable multi-isoform activity |
| **Null / negative** | Both new inactive + controls active → pipeline enrichment failure (publishable as such) |
| **Forbidden post-hoc** | Calling MD hinge asymmetry “confirmed selectivity” without meeting SI rule |

---

## 5. Under “only JNK1/2/3 enzyme assays” — computational experiments still needed

Constraint: **no kinome, no cell**. Computation must (i) strengthen interpretability of the 3-isoform IC50, (ii) reduce over-claim risk, (iii) fill gaps that wet-lab will not cover.

### Tier P0 — do before or with IC50 write-up (high value / feasible)

| # | Computational experiment | Why it matters under 3-isoform-only wet lab | Deliverable |
|---|--------------------------|-----------------------------------------------|-------------|
| **C1** | **Chemotype novelty audit** of 690 & 2157 vs ChEMBL/known JNK inhibitors (ECFP4 max Tc, Murcko, Bemis scaffold) | Answers “just another hinge binder?” without kinome data | Table: maxTc to SP600125/E1/CC-90001/CC-930 + nearest ChEMBL JNK ligands |
| **C2** | **Pose consensus check**: re-dock 690/2157 with **≥3 seeds** or Glide/alternative scoring; report RMSD cluster of top poses @ 3ELJ/3E7O/3TTI | Shows purchased poses are not single-run flukes | Pose RMSD matrix + 2D interaction consistency |
| **C3** | **MD replica mini-panel** on purchased set only: **≥2 independent seeds × 3 isoforms** for 690 and 2157 (shorter OK, e.g. 20–50 ns), **no ligand restraint** if affordable | Replaces reliance on single-replica ranking; supports pose credibility when claiming IC50–pose links | Mean±SD hinge occ / ligand RMSD; pass-rate stability |
| **C4** | **Pre-registered SI analysis notebook**: compute IC50, pIC50, SI_J2/J1, SI_J3/J1, with bootstrap CIs when duplicate wells exist | Prevents HARKing on isoform ratios | Locked analysis script before unblinding |
| **C5** | **Selectivity-method autopsy figure (already mostly done)**: one main-text table of Δsel/Gly87/ML failures on **literature benchmark**, explicitly decoupled from purchase | Makes RQ-C the publishable contribution even if RQ-B fails | Benchmark table + “not used for purchase” statement |

### Tier P1 — strongly recommended if compute available

| # | Experiment | Why |
|---|------------|-----|
| **C6** | **Off-isoform physics check without new biology**: MM-GBSA or short MD free-energy proxy **per isoform** for 690/2157 **after** IC50 known — only as **consistency check**, never as selectivity proof | Explains pan vs preference outcomes mechanistically |
| **C7** | **PAINS / aggregator / reactivity / logS risk re-check** at purchase SMILES (RDKit + public filters) | Rules out trivial assay artifacts before blaming target |
| **C8** | **Cross-docking / induced-fit lite** on JNK1 Ile106 vs JNK2/3 Leu sites for 2157 | Tests whether “J1≫J3 hinge” has a structural story or is MD noise |
| **C9** | **Retrospective enrichment simulation**: using docking/ML scores of purchased vs random library draws — approximate hit-rate prior | Frames n=2 result against expected baseline |

### Tier P2 — optional / do not pretend they replace kinome

| # | Experiment | Note |
|---|------------|------|
| **C10** | Ligand-based off-target inference (SEA/SwissTargetPrediction) for 690/2157 | **Hypothesis only**; cannot substitute kinome panel; report as “predicted risk” |
| **C11** | Buy-time 2231 **in silico only** full comparison table vs 2157 | Documents opportunity cost of not purchasing 2231 |
| **C12** | Longer unrestrained MD on actives **after** IC50 | Post-hoc mechanism; not selection justification |

### Explicitly **not** required / not sufficient under your constraint

- More Gly87 / Δsel filters on the library — already falsified.
- Claiming kinome cleanliness from docking to p38 without experiment.
- Using MD hinge asymmetry alone to declare selectivity success.

---

## 6. Minimal computational package to ship with the 3-isoform IC50 paper

If time is scarce, do **only P0: C1–C5**. That package:

1. Shows molecules are identified and somewhat novel (C1),  
2. Shows poses are stable to protocol noise (C2–C3),  
3. Locks statistics against HARKing (C4),  
4. Keeps the strongest scientific claim (predictor failure) independent of luck in n=2 (C5).

---

## 7. ARS experiment plan summary

## Experiment Overview

- **Title**: Prospective JNK1/2/3 IC50 of pipeline shortlist (690, 2157) with literature controls
- **Objective**: Test pipeline **family-binder enrichment** (primary); observe isoform SI only as secondary
- **Hypothesis (primary)**: ≥1 new compound active ≤10 µM on ≥1 JNK isoform
- **Hypothesis (secondary, weak prior)**: 2157 shows JNK1 preference; 690 is pan-like
- **Type**: wet-lab enzymatic assay + supporting computation (C1–C5)

## Analysis Plan

- **Primary metric**: IC50 (or %inh → IC50) per isoform  
- **Success threshold**: see §4  
- **Comparison**: E1, CC-90001 same plate; optional historical SP600125 if available in lab stocks  
- **Multiplicity**: report all 6 new-compound IC50s; SI exploratory with pre-stated rule  

## Monitoring / integrity

- Blind plate layout if possible  
- Duplicate wells  
- Do not reinterpret MD rankings after seeing IC50 without labeling post-hoc  

---

## 8. Final ARS objective scores (advisory)

| Dimension | Score | Note |
|-----------|-------|------|
| ID correctness | 5/5 | HIT↔library OK |
| Fit to RQ-A (family enrichment) | 4/5 | Good; n=2 limits |
| Fit to RQ-B (isoform selectivity) | 2/5 | Wrong best molecule omitted; docking discordant |
| Assay calibration design | 4/5 | E1+CC-90001 present; SP600125 optional missing |
| Over-claim risk if uncareful | High | Must soft-claim selectivity / kinome |
| Computations still owed under wet-lab constraint | **C1–C5 mandatory** | See §5 |

**One-line ARS verdict:**  
**Buy 690+2157 is scientifically OK for testing whether the pipeline yields JNK-family activity; it is not OK as a decisive test of JNK1 selectivity. Under JNK1/2/3-only biology, prioritize computational chemotype novelty, pose/MD replicas, locked SI analysis, and the already-negative selectivity-predictor benchmark — not more failed selectivity filters.**
