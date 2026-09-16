# Methods 2 + Results 3 review (user paste vs locked figures/tables)

Date: 2026-09-14  
Scope: the pasted English Sections 2–3, checked against typeset Figures 1–6 / S1–S4, Tables 1–3 / S1–S9, and locked directional CSVs.  
Line: retrospective four-state evaluation of docking scores for dual-target discrimination. Do not restore harvest-history language.

Display order (locked): EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD.

---

## 1. Figure and table order

First-cite sequence in Results should be:

| Section | Cite | What it actually shows |
|---|---|---|
| 3.1 | Figure 1A–C (once) | Four-state scheme + supply funnel |
| 3.2 | Figure 2A; Table S4 | Fixed-score dual-versus-neither minus dual-versus-STA, same pocket |
| 3.2 | Figure 2B; Table 2 | Directional AUROCs and `summary_min` |
| 3.2 | Figure 2C; Table 3 | Two-pocket **mean** dual-versus-neither |
| 3.2 (optional) | Figure S1 | EGFR operating points only |
| 3.3 | Figure 3; Table S5; Figure S2 | Ligand chemistry; Figure S2 = Vina CIs + best-descriptor **points**, not Δ CIs |
| 3.4 | Figure 4A; Table S6 | Matched vs mismatched `summary_min` |
| 3.4 | Figure 5A; Table S7 | Independent GNINA |
| 3.4 | Figure 5B; Table S7 | PIK3CA/mTOR alternative receptors |
| 3.4 | Figure 5C | Five-seed `summary_min` **distributions** (not in Table S7) |
| 3.4 | Figure S4; Table S2 | Cognate RMSD |
| 3.4 (optional) | Figure S3 | PIK3CA/mTOR panel size and exhaustiveness |
| 3.5 | Figure 6A; Table S3 | Threshold reassignment on fixed panels |
| 3.5 | Figure 4B; Table S6 | Holdout `summary_min` |
| 3.5 | Figure 6B; Table S4 | Scaffold/document cluster intervals for the two largest Δ |
| 3.6 | Figure 6C,D; Table S8 | BindingDB remainder after independence filters |

Fixes needed in the paste:

1. Methods 2.3 lists pairs as EGFR/HER2, PIK3CA/mTOR, AChE/BChE, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD. That is not Table 1 / figure order. Results 3.1 already uses the locked order. Use the locked order in both sections.
2. 3.1 cites “Figure 1C … (Figure 1)” twice. Cite Figure 1 once.
3. 3.2 names JAK1/TYK2 before EGFR/HER2 for the largest fixed-score Δ. Figure 2A and the abstract lead with EGFR/HER2. Report EGFR/HER2 first, then JAK1/TYK2.
4. 3.4 assigns five-seed distributions to Table S7. They are Figure 5C. Table S7 is independent GNINA, alternative receptors, and PPARG rescoring.
5. 3.4 says Figure 5C shows the EGFR fixed-score task difference. Figure 5C is `summary_min` by seed. The EGFR Δ remaining positive belongs in one sentence, not as a caption of 5C.
6. 3.3 cites Figure S2 for “the 95% interval for this difference”. Figure S2 does not draw those Δ CIs. Use Table S5.
7. Do not first-cite Figure 3C as if it were a separate main figure panel for TPSA only; keep Figure 3 as one figure.
8. The 2018 publication-year split is not in the typeset SI. Either drop 2.6.2 / 3.5, or one Methods clause that it is archived and not a typeset SI table.

---

## 2. Facts that do not match the locked branch

Keep ChEMBL 37 as the sole named source. Do not add REST harvest, dump-versus-API, or the 2026-08-26 snapshot.

### 2.2 Bioactivity

| Paste | Locked fact | Revision |
|---|---|---|
| “Database versions, retrieval information, and activity-processing checks are provided in Table S3.” | Table S3 is label/aggregation sensitivity (θ grid + max vs median on ChEMBL 37). | Name ChEMBL 37 in the sentence. Point Table S3 to threshold and aggregation checks only. |
| IC50, Ki, Kd, EC50, or Potency | Also `IC50app` and `Ki app`. | Add both. |
| High-confidence human single-protein subset | Not a typeset eight-pair check; removed from SI. | Delete. |
| Missing measurements not treated as low activity | Correct, but later restated. | Keep once. |

### 2.3 Panel construction

| Paste | Locked fact | Revision |
|---|---|---|
| “Seven target pairs had sufficient single-target-active supply under the strict 6.5/5.5 rule, while EGFR/HER2 was evaluated from the θ=6.0 pool” | EGFR/HER2 failed the 50/50 gate (min class count = 7). PIK3CA/mTOR met that supply gate, but **both** EGFR/HER2 and PIK3CA/mTOR **panels** were sampled from the θ=6.0 pools. | Split supply-gate from panel-pool. Do not imply only EGFR used θ=6.0 sampling. |
| Scaffold cap stated as if it applied after naming all eight pairs | Cap of 5 for EGFR/HER2 and 2 for PIK3CA/mTOR only. The other six primary panels had no extra scaffold cap. Holdout Murcko cap = 3. | Restrict the cap sentence to the two θ=6.0 panels. |
| “Only ligands with valid docking scores at both targets were included in the primary directional analyses.” | True for Table 2 complete-case. Independent GNINA complete-case *n* is smaller and belongs in 2.4.2 / Table S7. | Keep, but do not let it hide the GNINA *n*. |

### 2.4 Docking

Missing locked protocol (needed for a docking paper):

- Binding box = cognate AABB + 5 Å, minimum edge 20 Å (Table S2).
- Shared receptors: JAK1 6N7A for JAK1/JAK2 and JAK1/TYK2; PPARA 6LXA for PPARG/PPARA and PPARA/PPARD.
- `energy_range = 3`.
- Up to nine poses except 4EY7 and 3LXP (eight saved poses).
- 2.0 Å cutoff is lowest-RMSD among **saved** poses (search coverage), not a claim that mode 1 is always near-native.
- Primary Vina knobs: Table S1. Independent GNINA complete-case *n*: Table S7 (EGFR 28/38/32/11; PIK3CA/mTOR 18/13/12/4; JAK1/TYK2 30/32/29/14).

### 2.5 Metrics

| Paste | Problem | Revision |
|---|---|---|
| `summarymin=minAUCD/A(B),AUCD/B(A)` | Garbled. | \( \mathrm{summary}_{\min}=\min\{\mathrm{AUROC}_{D/A}(B),\mathrm{AUROC}_{D/B}(A)\} \) |
| `Smean=SA+SB2` | Chinese comma; missing division. | \( S_{\mathrm{mean}}=(S_A+S_B)/2 \) |
| Same-pocket dual-versus-neither next to \(S_{\mathrm{mean}}\) dual-versus-neither | Two different analyses. | Table S4 / Figure 2A = fixed pocket score. Table 3 / Figure 2C = mean score. Never mix. |
| Bootstrap described as one procedure | Table 2 = pooled non-stratified, B=2000. Table 3 = class-stratified. Dual-versus-neither in Table S4 uses the pooled dual/A-only/B-only resample with neither handled separately. | Split by table. |
| “weaker-arm summary,summarymin” | Concatenated words. | See formula above. |

### 3.2–3.6 numerical / SI mismatches

These numbers **match** the lock and can stay: JAK1/TYK2 Δ 0.444 [0.263, 0.620]; EGFR Δ 0.378 [0.205, 0.547]; directional AUROC 0.345–0.728; `summary_min` 0.345–0.692; PPARG/PPARA 0.649 [0.504, 0.751]; F2/F10 0.345 [0.211, 0.477]; Table 3 range 0.514–0.770; AChE TPSA 0.733 / 0.801; PIK3CA heavy-atom `summary_min` 0.463; pocket Δ EGFR 0.170 [0.060, 0.280], AChE 0.161 [0.037, 0.269]; holdout Δ −0.079 to +0.150; 4JPS/5DXT/4JSX `summary_min`; PPARG RTMScore 0.369 and GNINA CNN 0.500; independent GNINA EGFR dual-versus-neither 0.783 [0.610, 0.922] and dual-versus-B-only 0.220 [0.109, 0.343]; JAK1/TYK2 0.705 and 0.317 [0.183, 0.463]; cluster intervals; BindingDB eligibility failure.

Change or drop:

1. **Delete** “A ChEMBL 37 dump check covering all eight target pairs reproduced the production maximum pChEMBL values…” and “The separate 2026-08-26 ChEMBL API snapshot…”. The article source is ChEMBL 37. Sensitivity is max versus median on those records (EGFR 6/110, AChE 1/95 CHEMBL659, PPARA/PPARD 1/110 CHEMBL121; other five unchanged). That is Table S3, not a second database.
2. Independent GNINA: the paste omits PIK3CA/mTOR while saying “all three target pairs”. Locked Table S7: `summary_min` 0.633 (no bootstrap CI on the min row); weaker arm dual-versus-A-only 0.633 [0.427, 0.825]; dual-versus-neither 0.569 [0.222, 0.889]; n = 18/13/12/4. Add that row, or stop saying all three in prose. Do not cite 0.222 as a PIK3CA `summary_min`; that number is the lower CI of the mean-score dual-versus-neither AUROC.
3. Table 3: range 0.514–0.770 is true but thin. Name Table 3, give EGFR/HER2 0.756 and JAK1/TYK2 0.770, and state PIK3CA/mTOR `n_neither = 4`.
4. Incremental docking ≤ 0.023 is correct; keep Table S5 as the citation.
5. **Drop the 2018 year-split paragraph from Results.** The numbers (JAK1/JAK2 0.653 [0.476, 0.810]; JAK1/TYK2 0.368 [0.200, 0.545]) match `time_split_v1.csv` but that table is archived, not typeset SI. Leaving it in 3.5 implies a SI table that readers will not find.
6. Redocking sentence is correct if 2.0 Å is defined as lowest saved pose.

---

## 3. Repeated content

Keep the **rule** in Methods and the **number** in Results.

| Topic | Now | Keep |
|---|---|---|
| Dual-versus-neither purpose | 2.1 and again 2.5.1 | 2.1: why (comparison population). 2.5.1: how (fixed pocket vs \(S_{\mathrm{mean}}\)). |
| Four states / θ = 6.0 | 2.1, 2.2, 2.3, 3.1 | Definition in 2.2; screening in 2.3; counts in 3.1. |
| Holdouts | 2.6.2, 3.4, 3.5 | Construction in 2.6.2. Results only in 3.5 (Figure 4B). 3.4 should not preview holdout Δ. |
| Cluster bootstrap | 2.5.2 and 3.5 | Method in 2.5.2. Intervals in 3.5 (Figure 6B). |
| “Missing measurements were not treated as low activity” | 2.2 and implied in 2.7 | Once in 2.2. |
| Pair list | 2.3 and 3.1 | Short list in 2.3; 3.1 can say “the eight pairs in Table 1”. |
| No external set | 2.7 and 3.6 | Eligibility rule in 2.7. Failure + Figure 6C,D in 3.6. |

---

## 4. Over-defensive sentences (cut or tighten)

These add caution without adding a docking result. Cut or reduce to one clause.

- 2.1 “dual-target candidate prioritization” — overclaims. This is retrospective experimental-state evaluation of docking scores, not a prospective virtual-screening campaign.
- 2.2 “These labels describe activity states only within the studied target pair at the specified threshold.” — true; one short clause is enough.
- 2.5.1 “This value was used only to summarize the weaker of the two directional results.” — redundant with the formula.
- 2.5.1 “was reported only as a descriptive aggregate ranking … and was not used as the primary directional endpoint.” — keep the second half; drop “only as a descriptive aggregate ranking”.
- 2.6.1 “Because this descriptor was selected on the same panel, comparisons with Vina were treated as descriptive.” — keep (needed). Do not repeat in 3.3.
- 3.4 closing “Additional panel-size, exhaustiveness, and bootstrap sensitivity analyses are provided in the Supporting Information.” — unnamed dump. If kept, name Figure S3 (panel/exhaustiveness) and Table S4 (bootstrap). Do not say “bootstrap” without the table.
- Any “did not change the main directional judgments” — not in this paste, but do not reintroduce it.

---

## 5. English / register

| Paste | Issue | Prefer |
|---|---|---|
| “neither states” | Unnatural | “or neither” / “or inactive at both targets” |
| `A−only` (minus) | Wrong character | A-only |
| Bemis–Murcko vs Bemis-Murcko | Mix | Bemis–Murcko throughout |
| “the two target-aligned directional comparisons” | Heavy | “the two directional comparisons” |
| “PubChem was searched separately as an additional check on paired-data availability.Experimental states” | Missing space; “paired-data” | “paired data”. New sentence. |
| “weaker-arm summary,summarymin” | Concatenation | formula |
| “formed the primary evaluation of docking for dual-target candidate prioritization” | Screening slogan | “were the primary docking endpoints” |
| “to examine the source of the observed discrimination” | Vague | “to test whether discrimination tracked ligand chemistry or the assigned pocket score” |
| “a defined noncovalent small-molecule binding site” | OK | keep |
| “Systems requiring a substantially different docking treatment were excluded.” | Opaque | Name the actual exclusions if they are covalent, peptide, or missing holo sites; otherwise delete. |

---

## 6. Logic (same paper, two stories)

These are the places Methods and Results do not tell one docking story.

1. **Two dual-versus-neither analyses.** 2.1 / 2.5.1 / 3.2 opening use the **same pocket score** (Figure 2A, Table S4). 2.5.1 / 3.2 closing use **\(S_{\mathrm{mean}}\)** (Figure 2C, Table 3). If the reader thinks Figure 2A and Figure 2C are the same comparison, the docking claim collapses. Write two labeled endpoints.
2. **Supply gate versus panel pool.** 2.3 says seven pairs passed 6.5/5.5 and EGFR used θ=6.0. The next paragraph says EGFR **and** PIK3CA/mTOR were sampled from θ=6.0. Both sentences are partly true and together look contradictory. Rewrite as: screening used 6.5/5.5 except EGFR/HER2 (failed 50/50); panel sampling used θ=6.0 for EGFR/HER2 and PIK3CA/mTOR, and 6.5/5.5 for the other six; after sampling, **all** primary labels are θ=6.0.
3. **Primary claim versus ligand-only result.** 2.1 says directional docking is primary, then ligand-only “to examine the source”. 3.3 shows ECFP4 already separates several comparisons and docking adds ≤ 0.023. That is consistent if 2.1 says the chemistry and pocket-swap checks test **whether the directional AUROC is docking-derived**. It is inconsistent if 2.1 sounds like a virtual-screening success narrative.
4. **Matched-pocket Δ versus independent GNINA.** 3.4 reports EGFR/HER2 matched-pocket `summary_min` Δ excluding zero, then independent GNINA dual-versus-B-only = 0.220. Those are different engines. Do not imply GNINA “confirms” the Vina pocket-correspondence result. State: Vina matched-minus-mismatched on the main panels; separate pose generation with GNINA on three pairs.
5. **Holdout in 3.4 and 3.5.** 3.4 already says all seven holdout pocket-Δ intervals include zero. 3.5 then reports holdout `summary_min` levels. Split: 3.4 = pocket assignment (Figure 4A); 3.5 = sample composition including holdout `summary_min` (Figure 4B).
6. **3.2 “0.5” reading.** Intervals excluding 0.5 for PPARG/PPARA (above) and F2/F10 (below) are correct. Do not call the other six “no discrimination”; they are inconclusive at this sample size.
7. **External data.** 2.7 says BindingDB candidate sets were constructed, then 3.6 says no set was formed. Use “searched / filtered” in Methods and “no eligible set remained” in Results. “Constructed” implies a docking set exists.

---

## 7. Suggested section map (keep this docking line)

**2.1** Four experimental states; two directional AUROCs (dual vs A-only on score B; dual vs B-only on score A). Dual-versus-neither is a **population** control (Table S4), not a second primary endpoint. Chemistry and pocket-swap tests ask whether any AUROC is docking-derived.

**2.2** ChEMBL 37; activity types including IC50app/Ki app; max pChEMBL; both-target complete cases; θ = 6.0 states. Table S3 = θ and max/median only.

**2.3** Screening counts → docking-compatible eight pairs in Table 1 order → pool used for each panel → unified θ = 6.0 labels → complete-case scores.

**2.4** Receptors/boxes (S2), Meeko, ETKDGv3/MMFF, Vina 1.2.7, exhaustiveness, poses, mode-1 score, redocking definition, RTMScore/GNINA rescoring vs independent GNINA (S7 *n*).

**2.5** Two directional AUROCs; `summary_min`; sign flip of Vina; Table 2 bootstrap; Table 3 \(S_{\mathrm{mean}}\) and its bootstrap; paired method comparisons.

**2.6** Ligand descriptors + ECFP4 CV (incremental docking); matched vs swapped scores; receptors, holdouts, seeds, exhaustiveness, θ reassignment. No 2018 SI claim.

**2.7** BindingDB/PubChem filters and eligibility rule. Result of the search stays in 3.6.

**3.1** Census (Figure 1C) + eight pairs (Table 1).

**3.2** Fixed-score Δ (Figure 2A, S4) → directional AUROCs (Figure 2B, Table 2) → mean dual-versus-neither (Figure 2C, Table 3).

**3.3** Physicochemical AUROCs and Vina−descriptor Δ (S5; S2 for overlay) → ECFP4 vs ECFP4+docking (Figure 3).

**3.4** Pocket swap (Figure 4A) → GNINA/receptors/rescoring (Figure 5, S7) → redocking (S4/S2). Not holdout `summary_min`.

**3.5** θ grid (Figure 6A, S3) → max/median ChEMBL 37 (S3) → holdout `summary_min` (Figure 4B) → cluster Δ (Figure 6B). No API snapshot. No 2018 paragraph.

**3.6** BindingDB remainder (Figure 6C,D, Table S8). No external docking set.

---

## 8. What not to change

Do not edit `input_snapshot/` hashes or locked CSVs.  
Do not put REST versus dump, 2026-08-26 API, or high-confidence field screens back into Methods, Results, or SI.  
Do not treat Table 3 as a directional endpoint or Table S4 as the two-pocket mean.
