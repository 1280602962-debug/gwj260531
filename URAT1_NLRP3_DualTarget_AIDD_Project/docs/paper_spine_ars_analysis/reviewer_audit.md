# Reviewer-Aware Audit

> PaperSpine `reviewer-audit` + ARS `academic-paper-reviewer` personas · 2026-07-17  
> Manuscript maturity: **pre-MD package** (outline + partial draft)

## 0. Field / Persona Configuration (ARS Phase 0)

| Role | Persona | Focus |
|---|---|---|
| EIC | Hybrid CADD journal editor (JCAMD-like) | Scope fit, overclaim, desk-reject |
| R1 Methods | Computational chemist | Docking/MD reproducibility, redock, splits |
| R2 Domain | Gout / transporter–inflammasome biologist | Mechanism boundaries, canagliflozin/EGCG biology |
| R3 Perspective | Drug-repurposing / chemoinformatics | Clinical-library novelty vs “pipeline paper” |
| Devil’s Advocate | Hostile skeptic | “Docking of indirect drugs is meaningless” |

---

## 1. Reviewer Value Map

| Reviewer criterion | What reviewers/editors want | Our manuscript evidence | Current weakness | Revision action |
|---|---|---|---|---|
| Novelty | Clear delta vs closest prior | Asymmetry + 8973 track + Pareto≠F nomination; vs Eurycoma / PLK1-NLRP3 | Still easy to misread as generic ML→dock→MD | One “Unlike …” sentence in Abstract + Intro; Differentiation box |
| Significance | Why dual-node clinical reuse matters | Gout dual pathology; canagliflozin comorbidity narrative | No wet validation | Explicit “hypothesis for experimental follow-up” + assay plan |
| Technical soundness | Reproducible docking/ML | Glide XP + scripts; NLRP3 AUROC 0.89; y-scramble | 【待填】 numbers; MD pending; single conformation | Fill versions/grids; complete redock; report MD parameters |
| Evidence sufficiency | Enough to back Core contribution | Funnel counts, Pareto 6, module F | MD/redock missing → contribution incomplete | P0 MD + redock before submission |
| Clarity | One spine, consistent roles | Outline revised (cana lead) | CN draft still EGCG-as-MD-lead in places | Align MANUSCRIPT_DRAFT_CN with REVISED outline |
| Venue fit | Match aims & evidence bar | JCAMD/CBDD OK | **JMM Aims & Scope hard fail** if kept as #1 | Retarget journal list |

---

## 2. Reviewer Objection Register

| ID | Severity | What the reviewer may say | Preemptive fix | Status |
|---|---|---|---|---|
| O1 | CRITICAL | “Commercial Glide XP + short MD is desk-reject at JMM.” | Do not submit to JMM under current Aims; target JCAMD/CBDD; disclose open Vina/gnina scripts as reproducibility adjunct without mixing scores | Open → journal retarget |
| O2 | CRITICAL | “No wet assay → reject discovery claims.” | Never claim discovery; title/abstract as funnel + nomination; Discussion assay plan | Partially fixed in outline |
| O3 | MAJOR | “Canagliflozin acts via SGLT2, not URAT1 direct binding—docking is nonsense.” | Frame docking as **structural hypothesis for engagement tests**, co-exist with pathway-adjacent pharmacology; 7ALV as exploratory | Needs Discussion paragraph |
| O4 | MAJOR | “EGCG is a classic PAINS—why feature it?” | Feature as **blind recovery + demotion case**, not lead; MD optional SI | Outline fixed; draft must sync |
| O5 | MAJOR | “Pareto 6 then nominate someone not on the front—arbitrary.” | Methods: τ=90 dual gate + clean filters + ranking keys; report n at each step | Have script; need Table clarity |
| O6 | MAJOR | “URAT1 ML failed—why trust any ML in the funnel?” | Separate roles: NLRP3 ML for compression only; URAT1 ML excluded from ranking; 8973 docking enrichment | Present in Results R1–R3 |
| O7 | MAJOR | “Same as PLK1/NLRP3 asymmetric paper.” | Dedicated Differentiation section (transporter, clinical library, Pareto≠nomination, no fixed 0.5/0.5) | Doc exists; must enter manuscript |
| O8 | MAJOR | “Redock missing → poses untrusted.” | Report lesinurad@9DKB and 7ALV analog RMSD ≤2 Å | Pending |
| O9 | MINOR | “Inconsistent Glide vs Vina wording across drafts.” | Single score source for all reported numbers | Mostly fixed; audit leftovers |
| O10 | MINOR | “INNOVATION_POINTS lists S_trap / generative path not in paper.” | Quarantine aspirational TAPE-GATE claims from submission package | Pending doc hygiene |

---

## 3. Editorial Fit Map

| Question | Answer |
|---|---|
| Why would an editor send this out? | Clear methods story on **asymmetric dual-node clinical repurposing** with auditable nomination—not another phytochemical VS. |
| Desk-reject triggers | JMM hard MD gates; “discovered dual inhibitor”; canagliflozin = URAT1 inhibitor; unfinished MD; mixed docking engines in Results |
| Best editorial pitch (cover letter) | Protocol paper for gout dual-node reuse under data asymmetry; complementary to wet dual-inhibitor programs |
| Recommended decision if reviewed today (pre-MD) | **Major revision / return without review** until P0 MD+redock filled |
| Recommended decision if P0 complete + claim-bounded | **Revise & resubmit / accept with revisions** at JCAMD or CBDD |

---

## 4. ARS Devil’s Advocate (strongest counter-argument)

**Counter**: Dual-node docking of pathway-adjacent drugs (canagliflozin) and polyphenols (EGCG) confuses binding with pharmacology and misleads readers.

**Best rebuttal the paper can actually support**: The paper’s product is a **ranked, audited hypothesis list** plus a **negative methodological result** (URAT1 ML insufficient; PAINS can dominate Pareto). Docking is evidence of *possible* pocket occupancy for experimental disproof, not proof of clinical mechanism. Readers who want dual *inhibitors* should look to Eurycoma-style wet programs; readers who want *reuse of clinical space* need this funnel.

If this rebuttal is not in Discussion §1, the paper fails significance review.
