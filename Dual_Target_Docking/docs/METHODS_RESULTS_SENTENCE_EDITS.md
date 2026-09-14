# Sentence-level edits for the pasted Methods 2 + Results 3

Legend: **KEEP** | **CHANGE** | **DELETE** | **ADD**

Replace only the quoted original. Hyphens: write `A-only`, not `A−only`. En dash in Bemis–Murcko.

---

## 2.1 Study design

**CHANGE** `We assigned ligands to dual, A-only, B-only, or neither states based on experimental activity measurements at both targets.`

→ `Ligands were assigned to dual, A-only, B-only, or neither classes from experimental activity at both targets.`

**CHANGE** `Directional evaluation compared dual ligands separately with the A-only and B-only states.`

→ `Directional evaluation compared dual ligands with each single-target selective class, treating dual as the positive class.`

**CHANGE** `Dual and A-only ligands differed experimentally at target B, so target B scores were used for this comparison.`

→ `Dual and A-only ligands both meet the activity threshold at target A. Their experimental difference lies at target B, so pocket B scores were used.`

**KEEP** sense of `Dual and B-only ligands differed at target A, so target A scores were used.`

→ `Dual and B-only ligands differ at target A, so pocket A scores were used.`

**CHANGE** `The four experimental states and the two target-aligned directional comparisons are illustrated in Figure 1A,B.`

→ `The four experimental states and the two directional comparisons are shown in Figure 1A,B.`

**CHANGE** `The two directional comparisons formed the primary evaluation of docking for dual-target candidate prioritization.`

→ `The two directional AUROCs at a fixed target score were the primary docking endpoints.`

**CHANGE** `Dual-versus-neither was analyzed separately to examine how performance changed when the comparison population no longer contained single-target-active compounds.`

→ `A separate fixed-pocket dual-versus-neither comparison asked how AUROC changed when the control class no longer contained single-target-active ligands (Table S4).`

**CHANGE** `Ligand-only and pocket-correspondence analyses were then used to examine the source of the observed discrimination, while additional sensitivity analyses assessed its dependence on data processing, panel composition, receptor choice, and docking implementation.`

→ `Ligand-chemistry baselines and matched versus mismatched pocket scores tested whether that discrimination tracked ligand chemistry or the assigned pocket score. Two-pocket mean dual-versus-neither AUROCs were a conventional descriptive reference only (Table 3). Receptor substitution, independent pose generation, label sensitivity, and unused-pool holdouts assessed sensitivity of the directional results.`

**DELETE** any mention of publication-year subset analyses in 2.1 (not in the typeset SI).

---

## 2.2 Bioactivity data

**CHANGE** `Experimental activity data were obtained from ChEMBL.`

→ `Experimental activities were taken from ChEMBL 37.`

**DELETE** `Database versions, retrieval information, and activity-processing checks are provided in Table S3.`

**CHANGE** `Quantitative records with IC50, Ki, Kd, EC50, or Potency measurements were retained when a pChEMBL value was available.`

→ `Quantitative records with IC50, Ki, Kd, EC50, Potency, IC50app, or Ki app were retained when a pChEMBL value was available.`

**KEEP** `pChEMBL values were used on the standard negative-log molar scale.`

**KEEP** max pChEMBL; salt/mixture; both-target requirement.

**KEEP once** `Missing measurements were not treated as low activity.` Do not repeat later.

**KEEP** the four-state definitions. Fix `A−only` → `A-only`.

**KEEP** `The primary analysis used θ=6.0 for all target pairs.`

**CHANGE** `These labels describe activity states only within the studied target pair at the specified threshold.`

→ `A-only and B-only denote experimental state relative to the current pair and threshold, not proteome-wide selectivity.`

**CHANGE** `Activity-processing sensitivity analyses compared maximum and median aggregation and examined a high-confidence human single-protein subset.`

→ `Activity-processing checks compared maximum and median aggregation on the same ChEMBL 37 records for all eight pairs, with panel membership and docking scores held fixed (Table S3). That check does not replace Table 2.`

**DELETE** the high-confidence human single-protein subset (whole clause). Do not restore API or dump language.

---

## 2.3 Target-pair screening

**KEEP** the 10 dual / 10 A-only / 10 B-only gate and the 6.5/5.5 screening rule.

**CHANGE** `Eight target pairs formed the final evaluation set: EGFR/HER2, PIK3CA/mTOR, AChE/BChE, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD.`

→ `Eight target pairs formed the final evaluation set: EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1).`

**DELETE** `Seven target pairs had sufficient single-target-active supply under the strict 6.5/5.5 rule, while EGFR/HER2 was evaluated from the θ=6.0 pool because it provided adequate directional sample sizes and compatible structural data.`

**ADD in its place**

`EGFR/HER2 did not meet the strict 6.5/5.5 selective-supply criterion (minimum selective-class count = 7). It was retained because it had suitable human holo structures, a cognate-defined docking site, and sufficient dual, A-only, and B-only ligands at θ=6.0.`

**CHANGE** `EGFR/HER2 and PIK3CA/mTOR were sampled from the θ=6.0 pools, whereas the other six pairs were sampled from the strict 6.5/5.5 pools using the class quotas listed in Table 1.`

→ keep this sentence (it is the panel-pool rule, distinct from the supply gate).

**CHANGE** `Within each class, a Bemis–Murcko scaffold was represented at most five times for EGFR/HER2 and twice for PIK3CA/mTOR.`

→ `EGFR/HER2 and PIK3CA/mTOR limited a single Bemis–Murcko scaffold to at most five and two repeats within a class, respectively. The other six primary panels had no additional scaffold cap.`

**KEEP** `After panel construction, experimental states for all primary analyses were assigned uniformly using θ=6.0.`

**KEEP** `Only ligands with valid docking scores at both targets were included in the primary directional analyses.`

**DELETE** `Systems requiring a substantially different docking treatment were excluded.` unless the next sentence names covalent, membrane, or missing-holo cases (see Table S9). Prefer: `Pairs incompatible with noncovalent small-molecule docking in a defined holo site were excluded (Table S9).`

---

## 2.4 Receptor, ligand, docking

**KEEP** Tables 1 and S2; Meeko; cognate-defined sites; PG08-NL exception.

**ADD after site definition**

`Each docking box was the cognate heavy-atom range expanded by 5 Å on each side, with any edge shorter than 20 Å extended to 20 Å (Table S2). JAK1 6N7A was used for both JAK1/JAK2 and JAK1/TYK2. PPARA 6LXA was used for both PPARG/PPARA and PPARA/PPARD.`

**CHANGE** ligand prep: add seed and MMFF cap.

→ `RDKit added hydrogens and generated one three-dimensional conformer with ETKDGv3 (seed in Table S1). Ligands were geometry-optimized with MMFF for at most 200 iterations and converted to PDBQT with Meeko.`

**CHANGE** `Up to nine poses were retained for each ligand–receptor pair.`

→ `At most nine poses were retained for each ligand–receptor pair, with energy_range set to 3 kcal mol−1. AChE 4EY7 and TYK2 3LXP saved eight poses.`

**KEEP** exhaustiveness 16 vs 8; mode-1 Vina score.

**CHANGE** `The RMSD of the top-ranked pose and the lowest RMSD among saved poses were recorded, with 2.0 Å used as the near-native recovery cutoff.`

→ `The RMSD of the top-ranked pose and the lowest RMSD among saved poses were recorded. The 2.0 Å cutoff was applied to the lowest RMSD among saved poses as a search-coverage check, not as a claim that the top-ranked pose was near-native.`

**CHANGE** independent GNINA sentence. After `same receptors, ligand sets, and docking boxes`, **ADD**

`Independent GNINA complete-case counts were 28/38/32/11 for EGFR/HER2, 18/13/12/4 for PIK3CA/mTOR, and 30/32/29/14 for JAK1/TYK2 (dual / A-only / B-only / neither; Table S7). Primary Vina settings are in Table S1.`

**CHANGE** `Detailed scoring settings, docking parameters, and complete-case sample sizes are provided in Table S7.`

→ `Primary Vina settings are in Table S1. Independent GNINA complete-case sample sizes and alternative-receptor results are in Table S7.`

---

## 2.5 Metrics and statistics

**KEEP** dual as positive class; multiply Vina by −1.

**DELETE** the garbled `summarymin=minAUCD/A(B),AUCD/B(A)` and `This value was used only to summarize the weaker of the two directional results.`

**ADD**

`Their lower value was \(\mathrm{summary}_{\min}=\min\{\mathrm{AUROC}_{D/A}(B),\mathrm{AUROC}_{D/B}(A)\}\). It is a descriptive weaker-arm summary, not a separate overall-performance endpoint.`

**CHANGE** the long “same target-specific docking score” paragraph so it cites Table S4 only.

→ `To isolate the effect of changing the control class, the same pocket score was used for dual-versus-single-target-active and dual-versus-neither comparisons (Table S4). Pocket A scores were used for dual-versus-B-only and that dual-versus-neither comparison; pocket B scores were used for dual-versus-A-only and the corresponding dual-versus-neither comparison.`

**DELETE** garbled `Smean=SA+SB2， was reported only as a descriptive aggregate ranking...`

**ADD**

`A two-pocket mean score, \(S_{\mathrm{mean}}=(S_A+S_B)/2\), was used only for the dual-versus-neither AUROCs in Table 3 and Figure 2C. It was not a primary directional endpoint.`

**CHANGE** bootstrap paragraph: split Table 2 vs Table 3.

→ `Point estimates used the full analysis samples. For Table 2, 95% confidence intervals used 2000 ligand-level percentile bootstrap replicates, sampling with replacement from the pooled dual, A-only, and B-only set without class stratification. The same resampled dual ligands were used in both directional comparisons within each replicate. Replicates lacking a required class were omitted. Dual-versus-neither intervals in Table 3 used class-stratified bootstrap. All intervals were pointwise and were not adjusted for multiple comparisons.`

**KEEP** paired resamples for method/pocket comparisons; cluster bootstrap for the two largest Δ (Table S4).

**DELETE** “A class-stratified bootstrap was also performed as a sensitivity analysis…” as a second full procedure in the main text, or reduce to: `A class-stratified bootstrap for Table 2 is archived and did not replace the pooled intervals.`

---

## 2.6 Diagnostics

**KEEP** MW, HAC, cLogP, TPSA, ECFP4 radius 2 / 2048 bits; best-descriptor-on-same-panel is descriptive; incremental AUROC from ECFP4 vs ECFP4+docking.

**KEEP** matched vs swapped scores; Δ on `summary_min`; Table S6.

**KEEP** alternative receptors, holdouts, five seeds, exhaustiveness, θ reassignment.

**DELETE** `A publication-year subset analysis used the earliest source-document year of each ligand and a 2018 cutoff on the existing panels. Results were reported only for target pairs that met the sample requirements for both directional comparisons.`

**DELETE** `Missing scores were not treated as low activity` if already said in 2.2.

---

## 2.7 External data

**CHANGE** `BindingDB was used to construct candidate external evaluation sets for all eight target pairs.`

→ `BindingDB and PubChem were searched for paired experimental data across all eight target pairs.`

**CHANGE** `PubChem was searched separately as an additional check on paired-data availability.Experimental states were assigned using θ=6.0.`

→ `Experimental states were assigned using θ=6.0.`

**KEEP** development-set definition; source/structure/Tanimoto filters; n≥20 and ≥3 sources per class.

**KEEP** `Software versions... Table S1.`

Do not say a docking set was constructed. Eligibility failure belongs in 3.6.

---

## 3.1

**KEEP** census numbers 2,164,618 / 63,790 / 5,253 / 86.

**CHANGE** `Figure 1C summarizes these counts (Figure 1).`

→ `Figure 1C summarizes these counts.`

**KEEP** the 3.1 pair list (already Table 1 order). Optional shorten to `the eight pairs in Table 1`.

---

## 3.2

**CHANGE** `JAK1/TYK2 and EGFR/HER2 showed the largest fixed-score differences. For JAK1/TYK2, replacing B-only with neither increased the JAK1-pocket AUROC by 0.444 [0.263, 0.620]. For EGFR/HER2, the corresponding increase in the EGFR-pocket AUROC was 0.378 [0.205, 0.547].`

→ `EGFR/HER2 and JAK1/TYK2 showed the largest fixed-score differences. Using the EGFR-pocket score, the EGFR/HER2 dual-versus-B-only AUROC rose from 0.430 to 0.808 against neither (difference 0.378 [0.205, 0.547]). Using the JAK1-pocket score, the JAK1/TYK2 difference was 0.444 [0.263, 0.620].`

**KEEP** `For the remaining target pairs, the differences were smaller or their confidence intervals included zero (Figure 2A; Table S4).`

**KEEP** AUROC 0.345–0.728; `summary_min` 0.345–0.692; PPARG 0.649 [0.504, 0.751]; F2 0.345 [0.211, 0.477]; six pairs include 0.5 (Figure 2B; Table 2).

**CHANGE** `the other six target pairs had intervals that included 0.5` — do not add “no discrimination”. Optional: `the other six pairs crossed 0.5`.

**CHANGE** `For the two-pocket mean-score dual-versus-neither comparison, AUROCs ranged from 0.514 to 0.770 across the eight target pairs (Figure 2C; Table 3).`

→ `For the two-pocket mean-score dual-versus-neither comparison, AUROCs ranged from 0.514 to 0.770 (EGFR/HER2 0.756; JAK1/TYK2 0.770). The PIK3CA/mTOR neither sample was n = 4 (Figure 2C; Table 3). This comparison is not the fixed-pocket analysis in Figure 2A.`

---

## 3.3

**KEEP** TPSA 0.733 / 0.801; PIK3CA heavy-atom `summary_min` 0.463.

**CHANGE** `(Figure 3C; Table S5)` — keep Table S5; Figure 3C is the TPSA panel of Figure 3, not a separate figure.

**CHANGE** `The 95% interval for this difference included zero for six of the eight pairs, but not for F2/F10 or JAK1/TYK2 (Figure S2; Table S5).`

→ `The 95% interval for this difference included zero for six of the eight pairs, but not for F2/F10 or JAK1/TYK2 (Table S5). Vina summary_min intervals and best-descriptor points are in Figure S2.`

**KEEP** ECFP4 CV; incremental change ≤ 0.023 (Figure 3B; Table S5).

---

## 3.4

**KEEP** EGFR Δ 0.170 [0.060, 0.280]; AChE 0.161 [0.037, 0.269]; other six include zero.

**CHANGE** `All seven available internal-holdout intervals also included zero, with differences ranging from −0.079 to +0.150 (Figure 4A; Table S6). EGFR/HER2 did not have an unused-pool holdout.`

→ Keep this **only as pocket-Δ**. Do not preview holdout `summary_min` here. Cite Figure 4A, not 4B.

**KEEP** 4JPS / 5DXT / 4JSX numbers (Figure 5B; Table S7).

**KEEP** PPARG RTMScore 0.369 [0.233, 0.475].

**CHANGE** `to 0.500 with GNINA CNN rescoring (Table S7)`

→ `to 0.500 [0.356, 0.623] with GNINA CNN rescoring of the same Vina poses (Table S7).`

**CHANGE** the independent-GNINA paragraph.

Delete the implication that complete results for all three pairs are only EGFR + JAK1.

→ `Independent GNINA pose generation used the same receptors, ligands, and boxes on three pairs (Figure 5A; Table S7). For EGFR/HER2, dual-versus-neither AUROC was 0.783 [0.610, 0.922] (n_neither = 11), whereas dual-versus-B-only was 0.220 [0.109, 0.343]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.517, 0.876] and directional summary_min was 0.317 [0.183, 0.463]. For PIK3CA/mTOR, summary_min was 0.633 (no bootstrap CI on the min row), the weaker arm was dual-versus-A-only 0.633 [0.427, 0.825], and dual-versus-neither was 0.569 [0.222, 0.889] (n = 18/13/12/4). These independent-GNINA runs are a separate pose-generation pipeline, not a confirmation of the Vina matched-minus-mismatched result.`

**CHANGE** `Across five fixed Vina random seeds, the EGFR/HER2 fixed-score task difference remained positive, and the full seed-dependent summarymin distributions are shown in Figure 5C and Table S7.`

→ `Across five fixed Vina random seeds, the EGFR/HER2 fixed-score task difference remained positive. Seed-dependent summary_min distributions are shown in Figure 5C.`

**KEEP** redocking: all 14 receptors had at least one saved pose < 2.0 Å; top-ranked pose did not always (Figure S4; Table S2).

**CHANGE** `Additional panel-size, exhaustiveness, and bootstrap sensitivity analyses are provided in the Supporting Information.`

→ `PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S3.`

**DELETE** unnamed “bootstrap sensitivity in the SI” unless it cites Table S4.

---

## 3.5

**KEEP** θ-grid paragraph (Figure 6A; Table S3).

**DELETE** `A ChEMBL 37 dump check covering all eight target pairs reproduced the production maximum pChEMBL values and the primary directional results.`

**DELETE** `The separate 2026-08-26 ChEMBL API snapshot provided an additional sensitivity check for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR.`

**CHANGE** `Maximum-versus-median aggregation produced similar directional estimates in the corresponding checks (Table S3).`

→ `Relabeling the same ChEMBL 37 records by median rather than maximum pChEMBL flipped 6/110 EGFR/HER2 classes (summary_min 0.430 to 0.424), 1/95 AChE/BChE classes (CHEMBL659; 0.606 to 0.629), and 1/110 PPARA/PPARD classes (CHEMBL121; summary_min remained 0.446). The other five pairs kept class composition and summary_min point estimates (Table S3).`

**KEEP** holdout `summary_min` paragraph here (Figure 4B; Table S6), including PPARG 0.649 [0.504, 0.751] → 0.535 [0.350, 0.717].

**KEEP** cluster intervals (Figure 6B; Table S4).

**DELETE** the entire 2018 publication-year paragraph (`Only JAK1/JAK2 and JAK1/TYK2 met...`).

---

## 3.6

**CHANGE** `BindingDB and PubChem were searched for paired external data across all eight target pairs.` — KEEP.

**KEEP** filters and eligibility failure.

**CHANGE** do not say candidate sets were constructed.

→ `No target pair met the eligibility requirement. No external docking set was therefore formed (Figure 6C,D; Table S8).`

---

## Do not add

REST harvest, dump-versus-API, 2026-08-26 snapshot, high-confidence field screen, 2018 year-split as a typeset SI result.
