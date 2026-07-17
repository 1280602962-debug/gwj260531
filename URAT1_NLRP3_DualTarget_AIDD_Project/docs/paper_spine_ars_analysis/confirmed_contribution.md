# Confirmed Contribution

> PaperSpine V4 Contribution-First gate · 基于 **已有可投稿数据**（非 TAPE-GATE 愿景清单）  
> 日期：2026-07-17 · 工具链：PaperSpine + ARS-Codex（deep-research / reviewer）

## Core Contribution

| Field | Content |
|---|---|
| Main contribution statement | We establish a **reproducible, asymmetric computational funnel** on a ChEMBL clinical library (n=8,319) that (i) uses NLRP3 classification to compress chemical space, (ii) ranks dual-node structural hypotheses with Glide XP at URAT1 9DKB and NLRP3 7ALV plus Pareto non-domination, (iii) independently validates docking-led URAT1 ranking on an 8,973-compound retrospective track, and (iv) applies non-docking chemical/model audits that **downgrade PAINS-prone Pareto hits (EGCG)** while **nominating cleaner pathway-adjacent candidates (canagliflozin)**—without claiming experimental dual-pocket inhibition. |
| Contribution type | **new system / new analysis-or-benchmark**（可复现重定位流程 + 方法学审计），不是 new wet-drug discovery |
| One-sentence reviewer payoff | A transparent dual-node gout repurposing pipeline that **explains why URAT1 cannot be ML-ranked alone**, and that **actively filters false-promiscuous Pareto hits** rather than overselling them. |

## Why This Contribution Is Needed

| Field | Content |
|---|---|
| Field problem | Gout is managed on two largely separate axes (urate lowering via URAT1; inflammation via NLRP3/IL-1β), while dual-node chemistry is emerging (Eurycoma/Nat Commun 2025; NLRP3/URAT1-IN-1) but clinical-library **computational** dual-node funnels remain ad hoc. |
| Specific gap | No published **asymmetric** clinical-library funnel that (a) refuses URAT1 ML as main ranker because of documented benchmark failure, (b) keeps an **independent** large URAT1 docking retrospective track, and (c) separates **Pareto mathematics** from **chemically audited nomination**. |
| Concrete challenge | Data asymmetry: URAT1 regression recovers only 2/4 uricosuric benchmarks; NLRP3 assays are heterogeneous; transporter vs inflammasome pockets are not interchangeable; docking scores ≠ affinity. |
| Why prior work leaves it unresolved | **Eurycoma 2025**: phenotypic → de novo dual inhibitors (wet + synthesis), not clinical-library reuse. **PLK1/NLRP3-style asymmetric VS**: kinase+NLRP3, fixed-weight fusion, commercial library—not gout transporter + clinical drugs. **Single-target NLRP3 ML+dock papers**: no URAT1 node, no Pareto/audit split. **SGLT2–gout epidemiology**: pharmacology without a dual-structure computational nomination protocol. |

## How This Paper Responds

| Field | Content |
|---|---|
| Design response | Asymmetric funnel: NLRP3 ML gate → dual Glide XP → Pareto on (S_U, S_N) → modules A–F audit/nomination; 8973 URAT1-only retrospective track; MD on benchmarks + canagliflozin (not EGCG as lead). |
| Evidence required | (1) Quantified ML asymmetry; (2) funnel counts 8319→1588→1451→6; (3) 8973 enrichment; (4) Pareto + control-drug behavior; (5) PAINS/ADMET/y-scramble/AD/nomination numbers; (6) redock RMSD; (7) MD pose stability for 5 systems. |
| Evidence available | (1)–(5) largely available in repo; nomination of canagliflozin scripted; draft figures 2–4. |
| Evidence missing | **P0**: 5× MD (canagliflozin lead); redock RMSD; full Methods versions/parameters; unify 【待填】 Glide numbers. Without these, contribution softens to “funnel architecture + nomination protocol” without conformational discussion. |

## Claim Boundary

| Field | Content |
|---|---|
| Strong claims allowed | Asymmetric funnel is justified by data; Pareto front is mathematically defined; EGCG can be recovered then **downgraded** by explicit filters; canagliflozin is a **cleaner computational nominee**; 8973 supports docking-led URAT1 ranking. |
| Claims to soften or avoid | Dual-pocket direct inhibition; 1+1>2 synergy; clinical recommendation; canagliflozin as lesinurad-like URAT1 inhibitor; EGCG as developable lead; “first dual URAT1–NLRP3 inhibitor discovered computationally.” |
| Novelty risk | “Just another ML→dock→MD pipeline.” **Answer**: asymmetry + independent retrospective track + Pareto≠nomination split + explicit PAINS demotion. |
| Significance risk | “No wet assay → so what?” **Answer**: hypothesis-generation + reproducible nomination for T2DM–gout comorbidity adjunct discussion; cite Eurycoma as complementary wet route. |
