# DualFourClass-Bench: JCIM novelty gap, what was run here, what remains local

> Internal planning + evidence note (2026-08-24). Not a manuscript section.  
> Frozen-score analyses live in `data/jcim_novelty_v0/`. Claim bounds: `data/jcim_bench_v0/CLAIM_CEILING.md`.

---

## 0. Verdict in one page

The 2013 JCIM paper (Zhou, Li, and Hou; four kinase pairs; docking vs noninhibitors; structure dependence; high false-positive rate among predicted duals) already occupies the sentence “dual-target docking can be benchmarked, is structure-dependent, and has many false positives.” That sentence is **not** a 2026 contribution.

What *can* still be a JCIM contribution is a different sentence, and only if the data support it:

> Conventional dual-versus-inactive (or pooled) evaluation can make docking look competent at dual-target recognition, while a directional Dual vs A-only / Dual vs B-only hard-negative task on the **same scores** does not.

That comparison is now computed on the frozen K = 4 Vina scores (**no new docking**). The honest result is:

| Pair | Directional `summary_min` | Dual vs neither (`vina_mean`) | Dual vs all non-duals | Interpretation |
|---|---:|---:|---:|---|
| EGFR/HER2 | 0.430 [0.285, 0.579] | **0.756 [0.562, 0.920]** (n_neg = 12) | 0.551 [0.439, 0.666] | Conventional Dual-vs-inactive **looks useful**; directional worst arm is below chance. This is the novelty-proof case. |
| AChE/BChE | 0.606 [0.450, 0.750] | 0.649 [0.486, 0.805] | 0.579 | Small point-estimate gap; CIs overlap. Not a reversal. |
| PIK3CA/PIK3CB | 0.500 [0.346, 0.662] | 0.559 [0.364, 0.746] | 0.556 | Conventional is also near chance. Not a reversal. |
| PIK3CA/mTOR | 0.692 | 0.514 [0.222, 0.806] | 0.674 [0.513, 0.822] | Dual vs neither **n = 4, underpowered**. Do **not** read the negative Δ as “hard-negatives are easier.” |

**Do not write** “conventional evaluation systematically overestimates dual-target docking on all pairs.”  
**Do write** “on EGFR/HER2, a Zhou-like Dual-versus-neither readout (AUROC 0.756) would have supported docking-based dual recognition, while the directional worst arm is 0.430; Dual vs all non-duals already collapses to 0.55, showing that the extra difficulty is the selectives.”

That is field-level new knowledge relative to 2013: not that false positives exist, but that **the benchmark formulation itself decides whether docking appears to work**, and that the hard negatives are A-only/B-only selectives rather than inactives.

JCIM fit after A4 + B5: **still not a large-scale benchmark suite** (K = 4 remains a reviewer target). It is a **methodological evaluation paper** with a curated panel. Remaining risk is synthesis and claim freeze, not missing core experiments. Direct remaining reviewer surfaces: K = 4 generalizability; EGFR-dominant formulation contrast; `summary_min` as a custom aggregator (Table S26 already shows ranking invariance); PAB_034 timeout transparency (now in Table S27/S30).

---

## 1. What is actually new vs Zhou, Li, Hou 2013

**Zhou, S.; Li, Y.; Hou, T.** Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors. *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: [10.1021/ci400065e](https://doi.org/10.1021/ci400065e).

They did:

- 4 kinase pairs (CDK2–GSK3B, EGFR–Src, Lck–Src, Lck–VEGFR2);
- single-target inhibitor vs noninhibitor docking;
- dual-target identification;
- structure dependence;
- high false-positive rate;
- conclusion that docking helps but needs other VS methods.

They did **not**:

- define Dual / A-only / B-only / neither as an experimental four-state space;
- use directional Dual vs A-only (pocket B) and Dual vs B-only (pocket A);
- put ligand physicochemical / ECFP baselines on the same task;
- show that Dual-vs-neither AUROC and directional `summary_min` **diverge** on the same protocol;
- treat chemotype-matched selectives as a second hard-negative layer.

So the novelty is **not** “we also have four pairs.” It is **benchmark formulation + confounder attribution** on experimental selectives.

---

## 2. A / B / C items: here vs local

| ID | Item | Docking? | Status 2026-08-24 | Where |
|---|---|---|---|---|
| **A1** | Conventional vs DualFourClass | No (frozen scores; `neither` already docked on main panels) | **DONE** | `tables/formulation_*.csv` |
| **A2** | Chemotype-matched hard negatives (ECFP4 Tanimoto) | No | **DONE** (T ≥ 0.7 empty; used 0.3 / 0.4 / 0.5) | `tables/chemotype_matched_hardneg_v1.csv` |
| **A3** | Incremental information (physchem / ECFP / docking / combinations) | No | **DONE** | `tables/incremental_information_v1.csv` |
| **A4** | max vs median pChEMBL on **all** frozen-panel IDs | No | **DONE** | `tables/assay_max_vs_median_*_v1.csv`; agreement table; SI Table S29. 7/110, 1/95, 1/99, 0/48 class flips; pair-level `summary_min` insensitive. Do not mix EGFR frozen 0.430 with API-max 0.417. |
| **B5** | Second-pair receptor swap | **Yes (done)** | **DONE** | PIK3CA/PIK3CB into 4JPS/5DXT, 2WXF frozen. `summary_min` 0.500 → 0.691 / 0.685 (opposite of PM48). PAB_034: 100/99/1 timeout on original and both alts. |
| **B6** | Prospective-style mixed-library enrichment | No | **DONE** | `tables/mixed_library_enrichment_v1.csv` (EGFR Top-10: 9/10 hard-negatives; EF5 = 0.66, worse than random) |
| **C7** | Wet-lab prospective | Wet lab | **Cannot do here** | Optional; not required to keep a benchmark paper honest |

Do **not** add more docking engines. Do **not** chase 4 → 10 target pairs if label quality drops.

---

## 3. A1 — formulation comparison (primary novelty experiment)

Source: `scripts/benchmark_formulation_v1.py` on frozen pocket-matched Vina.

Definitions (same ligands, same scores):

- **DualFourClass directional:** Dual vs A-only in pocket B; Dual vs B-only in pocket A; `summary_min` = min of the two arms (Table 2).
- **Conventional Dual vs neither:** Dual vs experimental inactives, scored with `vina_mean` (Zhou-like dual vs noninhibitors) and with `vina_worst` (AND-like).
- **Dual vs all non-duals:** Dual vs (A-only + B-only + neither); still not directional.
- **Single-target analogue:** (dual + A-only) vs (B-only + neither) in pocket A, and the symmetric B contrast.

EGFR/HER2 is the only pair where Dual-vs-neither is both (i) apparently strong and (ii) clearly above the directional worst arm. The two CIs barely overlap (0.562–0.920 vs 0.285–0.579). This is **not** a paired Δ (different negative sets); do not call it a significance test.

AChE and PIK3CA/PIK3CB do **not** show a formulation reversal: Dual-vs-neither is only ~0.04–0.06 higher and remains compatible with chance.

PIK3CA/mTOR Dual-vs-neither is unusable (neither n = 4). Dual vs all non-duals on that pair is 0.674, close to `summary_min` 0.692 — i.e. when selectives are the bulk of the negatives, pooled scoring does not inflate PM the way Dual-vs-neither inflates EGFR.

**Manuscript use:** main-text Table 3 + SI Fig S4. One Results paragraph after Table 2. Do not rebuild the paper around a universal overestimation law.

---

## 4. A2 — chemotype-matched hard negatives

Tanimoto ≥ 0.7 nearest-neighbor A-only/B-only sets are **empty** on these panels. That is itself a finding: DualFourClass hard negatives are activity-hard, not scaffold-matched.

At modest cuts:

- PIK3CA/PIK3CB Dual vs A-only: unmatched 0.691 → T ≥ 0.3 (n_neg = 11) **0.503**; distant T < 0.3 (n_neg = 16) **0.819**. The “good” arm is largely chemotype.
- EGFR Dual vs A-only: 0.666 → T ≥ 0.4 (n_neg = 15) 0.579 → T ≥ 0.5 (n_neg = 7) 0.556.
- AChE Dual vs A-only: 0.650 → T ≥ 0.3 (n_neg = 7) 0.571; T ≥ 0.4 empty.
- PIK3CA/mTOR matched n is 1–6 at T ≥ 0.4: do not interpret those AUROCs.

**Allowed claim:** part of the apparent docking signal disappears when hard negatives are required to be chemically closer to duals; T ≥ 0.7 matching is not supported by current supply.  
**Forbidden:** “chemotype matching proves docking has no pocket signal” (n too small at high T).

---

## 5. A3 — incremental information

Scaffold `GroupKFold` logistic AUROC (not the rank-AUROC in Table 2):

On every directional arm, **ECFP4 + docking ≈ ECFP4** (Δ typically |Δ| ≤ 0.01; sometimes docking slightly lowers CV AUROC). Physchem + docking does not recover ECFP. Logistic docking AUROC is often below rank docking AUROC (EGFR Dual vs A-only: logistic 0.591 vs rank 0.666) — report both, do not swap them.

**Allowed:** docking adds little incremental discrimination beyond 2D chemical structure under scaffold-aware CV.  
**Forbidden:** “docking is redundant with chemistry on every possible dual-target pair.”

This is the answer to Reviewer 3 (“why is this a docking benchmark rather than a chemical-series benchmark?”): because the benchmark *shows* that docking rank AUROC is not independent of chemotype, and that is a finding about docking evaluation, not a reason to abandon the docking task.

---

## 6. B6 — mixed-library enrichment (already partly in Results 3.2)

EGFR/HER2, library n = 110, 28 duals, Top-10 by `vina_mean`: 1 dual / 5 A-only / 4 B-only / 0 neither; EF10 = 0.39; hard-neg fraction = 0.90. Same pattern at EF5. This is the VS-facing version of A1: pooled dual ranking **enriches selectives, not duals**.

AChE/PIK3CB show modest EF5 ≈ 2.1 with hard-neg fraction 0.4. PM Top-10 EF = 1.6 (small library). Pair-dependent, consistent with Table 2.

---

## 7. A4 — assay aggregation (DONE; Table S29)

Full-panel re-fetch completed. Native contrast is API-max vs API-median. Label agreement = 1 − n_flip / n_scored:

- EGFR/HER2: 7/110 (93.6%); API-max 0.417 → median 0.424; frozen Table 2 is 0.430 because of EH120_060.
- AChE/BChE: 1/95 (98.9%); 0.606 → 0.629.
- PIK3CA/PIK3CB: 1/99 (99.0%); 0.500 → 0.500.
- PIK3CA/mTOR: 0/48 (100%); 0.692 → 0.692.

Numeric max ≠ median is common; class flips at θ = 6.0 are not. Do not promote the 27-ligand diagnostic. Do not write max pChEMBL as an unresolved fatal threat.

Details: `data/jcim_novelty_v0/analysis/A4_B5_STATISTICAL_AUDIT_V1.md`. Frozen Vina scores were not recomputed. The 27-ligand diagnostic is superseded and must not be promoted.

---

## 8. B5 — second-pair receptor swap (DONE; Table S30)

PIK3CA/PIK3CB panel redocked into already-prepared 4JPS and 5DXT; 2WXF held frozen; exhaustiveness 8.

- Original: 0.500 [0.347, 0.648]
- 4JPS: 0.691 [0.516, 0.779] (Δ +0.191); weak arm switches to frozen D/A
- 5DXT: 0.685 [0.506, 0.768] (Δ +0.185)

Opposite of PIK3CA/mTOR (0.692 → 0.486 / 0.505). Phrase as receptor-realization effect, not robustness, not collapse.

PAB_034: 100 attempted / 99 successful / 1 timeout on original 4L23 and both alts. Docking timeout, not a label filter. Same 99-ligand set as Table 2.

Do **not** dock chimeric 3T8M. Extra PIK3CB crystals are optional later work, not required.

---

## 9. C7 — wet lab

Not required for a benchmark/evaluation article. Would change the paper class. Out of scope for this agent.

---

## 10. How to change the paper without over-claiming

**Title (optional, not applied in this commit):**  
*Benchmarking Docking-Based Dual-Target Recognition with Directional Selectivity Hard Negatives*

**Central questions (nested):**

1. How does benchmark formulation affect the apparent ability of docking to recognize dual-target ligands?  
2. Can docking distinguish dual-actives from directional experimental selectives on both arms?

Keep (2) as the prespecified primary endpoint. (1) is the novelty comparison with 2013. Dual-versus-neither is a **nonselectivity-controlled comparator**, not “the conventional benchmark.”

**Do not** rename DualFourClass-Bench a “comprehensive benchmark suite.”  
**Do** call it a curated four-state panel plus evaluation protocol.

---

## 11. Remaining JCIM risk after this increment

Still real:

- K = 4 (call it data-constrained, not comprehensive). Novelty evidence is **EGFR/HER2-dominant**.
- ECFP still strong: A3 helps, but reviewers can still say the labels are chemotype-structured.
- Two receptor-swap pairs share PIK3CA; opposite directions are not a universal law.

Closed: A4 full-panel max vs median (Table S29); B5 second-pair receptor swap (Table S30); Dual-vs-neither is no longer called “the conventional benchmark”; 0.756 vs 0.430 is descriptive; `summary_min` ranking is aggregation-insensitive (Table S26); docking failures are censused including PAB_034; identifier prefix deleted; ECFP incremental / CV / wrong-pocket / descriptor wording tightened.

Do **not** add more docking engines, MD, or extra target pairs. Remaining work is quantitative synthesis, figures, abstract, SI freeze, Zenodo, JCIM typesetting.

If a reviewer asks for 10 pairs or 10 engines, the correct reply is the supply audit (Figure 2) and this note: more engines do not test formulation.
