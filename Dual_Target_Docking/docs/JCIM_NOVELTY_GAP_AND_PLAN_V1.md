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

JCIM fit after this increment: **still not a large-scale benchmark suite** (K = 4 remains a reviewer target). It is closer to a **methodological evaluation paper** with a curated panel. Direct submission risk is lower than before A1, but **not low**, until A4 (full-panel median) and ideally B5 (second-pair receptor swap) are in.

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
| **A4** | max vs median pChEMBL on **all** frozen-panel IDs | No | **SCRIPT READY; LIVE FULL-PANEL FETCH BLOCKED** | ChEMBL `activity.json?molecule_chembl_id=` currently HTTP 500 / timeout. 27-ligand diagnostic already exists and must **not** be promoted to SI. Re-run `scripts/assay_aggregation_max_vs_median_v1.py` locally when the activity filter is up. |
| **B5** | Second-pair receptor swap | **Yes** | **LOCAL ONLY** | Vina is not in this cloud image. Cheapest protocol: dock PIK3CA/PIK3CB panel into already-prepared 4JPS / 5DXT (same PIK3CA prep as PM48 swap), keep 2WXF scores. Alternative: EGFR/HER2 alt crystals (needs new prep). Playbook below. |
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

## 7. A4 — assay aggregation (script ready; not a completed SI table)

`scripts/assay_aggregation_max_vs_median_v1.py` re-fetches assay-level `pchembl_value` for every scored frozen-panel ligand, relabels at θ = 6.0, and recomputes directional AUROC on frozen Vina scores.

This environment: ChEMBL `/activity.json?molecule_chembl_id=` returned HTTP 500 / timeouts on 2026-08-24 after a successful `status.json`. Do not fabricate a full-panel table from the 27-ligand diagnostic (`max_vs_median_diagnostic_sample_v1.csv`; 0/27 θ = 6.0 flips in that slice).

**Local command** (when the activity endpoint accepts molecule filters):

```bash
python3 Dual_Target_Docking/data/jcim_novelty_v0/scripts/assay_aggregation_max_vs_median_v1.py
```

Outputs: `assay_max_vs_median_{ligand,summary,auroc,flips}_v1.csv`. Promote to SI only after the full panel exists. Then one Limitations sentence can be replaced with “median relabel flipped N ligands; summary_min moved by Δ.”

---

## 8. B5 — second-pair receptor swap (local Vina)

This cloud image has **no `vina` binary**. Alternate PIK3CA receptors **4JPS** and **5DXT** are already prepared in `data/jcim_structure_robust_v0/` (boxes + `*_receptor.pdbqt`) from the PM48 swap. That is the cheapest second pair:

1. Prepare PIK3CA/PIK3CB panel ligand PDBQTs with the frozen RDKit ETKDG + meeko protocol (same as main panel).
2. Dock each ligand into 4JPS and into 5DXT (exhaustiveness = 8 to match the PIK3CB main panel; seed 20260727; n_modes = 9).
3. Keep frozen 2WXF PIK3CB scores.
4. Recompute pocket-matched `summary_min` (Dual vs B-only uses the alternate PIK3CA score; Dual vs A-only still uses 2WXF).

Do **not** dock chimeric 3T8M. 2Y3A / 4BFR exist as extra PIK3CB crystals in `pik3ca_pik3cb_panel_v0/receptors/` if a B-end swap is wanted later; they need their own cognate QC pass in the write-up.

EGFR/HER2 is the scientifically sharper second pair (high sequence identity, directional failure on one arm) but needs new crystal selection + prep. Prefer PIK3CA/PIK3CB first because the PIK3CA prep is already frozen.

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
- Receptor swap still one pair until B5.
- max pChEMBL until A4 completes.
- ECFP still strong: A3 helps, but reviewers can still say the labels are chemotype-structured.

Closed in the claim-hardening round (`scripts/claim_hardening_v1.py`; manuscript rewrite): Dual-vs-neither is no longer called “the conventional benchmark”; 0.756 vs 0.430 is descriptive; `summary_min` ranking is aggregation-insensitive; docking failures are censused; identifier prefix deleted; ECFP incremental / CV / wrong-pocket / descriptor wording tightened.

If a reviewer asks for 10 pairs or 10 engines, the correct reply is the supply audit (Figure 2) and this note: more engines do not test formulation.
