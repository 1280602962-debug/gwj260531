# Results (JCIM Articles draft, English)

> Companion to [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) (Chinese authoritative for this rewrite cycle).  
> Numbers from `PRIMARY_METRIC_V2.md`, `PRIMARY_METRIC_CLAIM_GATE.md`, SI Tables S4–S11. No fabricated experiments.  
> **Framing (deliberately not absolute; not a named method product):** not "Docking can/cannot identify dual-target ligands," and not "we developed a novel Dual-target Docking Reliability Assessment Framework (D-DRAF)." Preferred: *evaluating the reliability and limitations of docking-based dual-target recognition* via a systematic benchmarking framework and the DualFourClass-Bench resource. See [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md).

## 3. Results

### 3.1 Construction of a dual-target recognition benchmark: public-data limits on hard-negative supply

Dual-target docking evaluation requires four ligand classes: dual, A-selective, B-selective, and neither. Experimentally defined selective ligands on each arm serve as hard-negative selective ligands for testing whether a score can suppress both single-target arms.

Across 49 audited ChEMBL target pairs under the primary strict rule (dual: both ends ≥ 6.5; selective: active end ≥ 6.5 and opposite ≤ 5.5), only four pairs retained ≥50 hard-negative selective ligands on **both** ends. Despite the large number of ChEMBL target pairs, **balanced dual-target benchmarking was severely constrained by the scarcity of experimentally characterized hard-negative selective ligands**. After excluding metal-dependent HDAC1/HDAC6, three pairs supported reasonably balanced strict panels; EGFR/HER2 entered as a supply-limited case (few strict B-selective ligands), not as a thick panel. The frozen K = 4 set follows this audit, not post-hoc selection of docking-favorable pairs (Methods 2.1–2.3).

### 3.2 Primary analysis uses one unified label rule (θ = 6.0); threshold sensitivity is supporting

To avoid the appearance of "different thresholds for different pairs," Table 2 reports all four pairs under **one common θ = 6.0 rule** (both ends ≥ 6.0 = dual; one end ≥ 6.0 with the other < 6.0 = the corresponding selective class). For EGFR/HER2 and PIK3CA/mTOR, this is the rule already used at construction. For AChE/BChE and PIK3CA/PIK3CB, construction used the stricter 6.5/5.5 rule for supply auditing, but θ = 6.0 gives **identical** ligand classification and AUROC on this data (Table S4) — i.e., labels for these two pairs are insensitive to threshold across the tested grid, a verified fact rather than a selectively reported one.

As a supporting robustness analysis (not a second, competing primary standard), we further relabeled all four pairs at θ ∈ {5.5, 6.5} and the strict 6.5/5.5 rule and recomputed summary_min (Table S4). EGFR/HER2 and PIK3CA/mTOR are more threshold-sensitive: under the strict rule they drop to 0.324 (only 7 B-selective ligands, underpowered) and 0.639 (only 4, underpowered), respectively, versus 0.430 and 0.692 under θ = 6.0 — both declines but no rank change. Ranking trends held across the full grid: PIK3CA/mTOR highest; the other three pairs ≤ 0.61.

### 3.3 Docking shows limited, pair-dependent ability to discriminate true dual-target ligands

A dual-target score must suppress both selective arms. Pooling the two pocket scores (e.g., mean) can let the stronger arm mask failure on the weaker arm. The primary metric is therefore pocket-matched directional AUROC: dual versus A_only scored in pocket B, dual versus B_only scored in pocket A, with summary_min as the smaller arm (Methods 2.6). Scores are \(S=-E_{\mathrm{Vina}}\) (higher better); dual is the positive class.

Two points must be kept distinct. First, **direction-specific discrimination failure**: on EGFR/HER2, pocket-matched dual-versus-B_only AUROC is 0.430, and the weaker arm under a pooled protocol can read near ~0.28; that low value reflects failure of that direction itself, not a pooling arithmetic that "creates" 0.28 from 0.50. Second, **pooling can mask the weak arm**: on the same pair, a pooled summary can approach ~0.50 and appear merely mediocre. Relative to pooling, pocket matching raised point estimates across pairs without changing rank order (Table S6).

**Table 2.** Pocket-matched directional AUROC on the frozen K = 4 set, all four pairs under the single unified θ = 6.0 label rule (Vina). Wrong-pocket, ligand-efficiency, and descriptor baselines are in Supporting Information Table S6; full threshold grid in Table S4.

| Pair | n (dual / A_only / B_only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] |
|---|---:|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] |

Docking signal is generally weak and highly target-pair dependent. Only PIK3CA/mTOR has summary_min point estimate above both 0.5 and the best trivial descriptor baseline (heavy-atom count 0.463), yet its 95% CI lower bound still approaches 0.5. EGFR/HER2 and PIK3CA/PIK3CB do not beat their best descriptor baselines; AChE/BChE (0.606) remains below the TPSA baseline (0.733). RTMScore and GNINA did not change this pattern.

### 3.4 Physicochemical and structural confounding dominates several apparent dual-target signals

The central finding is not that "docking is bad," but that **many apparent dual-target signals are explained by ligand properties**. Relative to the best trivial descriptor, docking–descriptor Δ intervals lie below 0 for EGFR/HER2 and PIK3CA/PIK3CB; AChE/BChE does not clear the gate; PIK3CA/mTOR exceeds heavy-atom count at the point estimate, but the Δ 95% CI still includes 0 (Table S6). Wrong-pocket summary_min values were 0.260, 0.444, 0.349, and 0.602; pocket-matched gains over wrong-pocket exceeded 0.09 on all pairs, underscoring ligand-level confounding.

On AChE/BChE, mean TPSA was ≈ 75 (dual) versus ≈ 51 (hard-negative selective ligands). TPSA alone separated dual from hard negatives with AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA raised dual-versus-B AUROC from 0.606 to 0.807. That increase indicates that **the apparent docking contribution was largely dependent on physicochemical covariates**; the docking odds ratio near 1 (OR ≈ 1.18) should not be read as retained independent directional information. On PIK3CA/mTOR, covariate-adjusted AUROC shifts were smaller (≈ +0.07 to +0.11) with OR ≈ 2.19 and 3.08, suggesting only residual pocket-specific signal, to be read with Δ intervals that still include 0.

Scaffold-grouped ECFP4 baselines further elevate the confounding narrative: fold AUROCs often 0.78–0.91 and frequently above the corresponding docking arms (e.g., EGFR/HER2 dual-versus-B_only: fingerprint 0.85 vs docking 0.43). Labels correlate with chemotype; this alone does not prove pocket-physical specificity of docking scores.

### 3.5 Robustness checks and case-dependent success

Changing PIK3CA/mTOR exhaustiveness from 16 to 8 moved Vina summary_min from 0.692 to 0.660 (Δ ≈ +0.03)—far smaller than between-pair differences. Single-target enrichment against property-matched weak binders (pChEMBL ≤ 5.5) gave AUROC 0.603/0.629 and EF1% 2.04/2.00 on 4L23/4JT6: docking retained limited enrichment capability, not a strong VS engine.

Expanding PIK3CA/mTOR to PM110 is a **stability check**, not independent validation and not an attempt to "rescue" the estimate with a larger panel. PM48 itself is small (18/14/12). On PM110, Vina summary_min was 0.648 [0.51, 0.76] versus 0.692 on PM48 (Δ ≈ −0.04); ranking trend remained consistent.

Across §3.2–3.5, **only PIK3CA/mTOR showed reproducible but modest pocket-related discrimination**; apparent signals on the other three pairs are largely explained by ligand properties or 2D chemotype.

### 3.6 Structural determinants of docking discriminability across pairs (exploratory)

Ligand-level confounding (§3.4) does not by itself explain why discriminability differs so sharply across pairs. We extracted the longest protein chain sequence per target from the frozen receptor PDB files (Biopython `PDBParser`) and computed pairwise **whole-chain sequence identity** within each pair using global BLOSUM62 alignment (`PairwiseAligner`, gap open = −11, extend = −1) as a coarse proxy for overall structural similarity. Script and raw output are in `data/jcim_bench_v0/analysis/structural_context_v1/` (Supporting Information Table S7).

**Table 3.** Within-pair whole-chain sequence identity versus discriminability (summary_min).

| Pair | Identity, % (over alignment length) | Identity, % (over shorter chain) | summary_min (θ = 6.0 / strict) |
|---|---:|---:|---|
| PIK3CA/mTOR | 18.1 | 21.0 | 0.692 / 0.639 |
| PIK3CA/PIK3CB | 40.5 | 43.3 | 0.500 / 0.500 |
| AChE/BChE | 51.9 | 53.1 | 0.606 / 0.606 |
| EGFR/HER2 | 71.4 | 76.6 | 0.430 / 0.324 |

The pattern is counterintuitive: the **most** discriminable pair (PIK3CA/mTOR) has the **lowest** whole-chain identity, while the **least** discriminable pair (EGFR/HER2) has the **highest** identity (ErbB-family kinase domains are highly homologous, consistent with the existence of many pan-ErbB / dual inhibitor chemotypes such as lapatinib and afatinib). This suggests that greater overall similarity between the two targets can make the two pockets more physically indistinguishable to a docking score — the opposite of "more similar targets are easier to recognize as dual."

This proxy has clear limits: it is whole-chain, first-order sequence identity, not a structurally superposed, pocket-residue-level RMSD or interaction-fingerprint (PLIF) similarity. The latter requires validated superposition tools (e.g., TM-align, PyMOL align/super) and pocket-residue correspondence, which were not set up or validated in this round; we do not fabricate a pocket RMSD number here, and leave it to future work. PIK3CA and mTOR have the lowest whole-chain identity yet belong to the PIKK-related superfamily with known local structural homology at the ATP-competitive site — the same structural basis underlying real PI3K/mTOR dual-inhibitor chemotypes (e.g., PI-103, omipalisib). The low whole-chain identity in the table should therefore be read as "divergent overall architecture with a locally accessible, overlapping ATP site," not "dissimilar pockets," consistent with the pose-level clues in §3.7. With only four pairs (n = 4), this section is descriptive; we make no formal correlation or significance claim.

### 3.7 Structural clues for the only reproducible dual-target signal (case-level)

PIK3CA (4L23) and mTOR (4JT6) are ATP-competitive kinase-related pockets; cognate ligands recover near-native poses under protocol checks (Table S3). Existing pose-level failure typology shows hinge/ATP-like duals that can rank well on both ends (type T2) versus cases where rescoring prefers off-hinge poses (type T5). Even on this best pair, ATP-site cross-chemotypes can be misread as dual when both pockets yield geometrically clean, hinge-positive poses. We do not claim a completed shared-residue / PLIF campaign; a dedicated structural-determinants analysis remains future work. The PIK3CA/mTOR advantage should be read as limited directional signal for some chemotypes under a shared ATP recognition framework, not as a validated general dual-target decision rule.

### 3.8 An unresolved robustness gap: panel-composition resampling

The existing bootstrap (Methods 2.6; B = 2000) resamples ligands **within a fixed panel**, answering "how uncertain is this particular docked set," not "would the conclusion still hold with a different, equally sized sample." Answering the latter requires drawing many independent panels (e.g., 1000 draws) from each pair's strict supply pool at the same quota, docking the not-yet-docked pool members, and summarizing the resulting distribution of summary_min. This requires docking beyond the currently frozen score package; §3.9 below reports one such unused-pool draw per pair, which addresses this gap partially but not with a resampled distribution.

### 3.9 Unused-pool holdout (post-panel-freeze validation)

From the strict ChEMBL pool not used in panel construction or protocol tuning, we drew 60 ligands per pair (dual/A_only/B_only = 20/20/20; seed 20260731) for three pairs, docked them under the frozen protocol, and recomputed pocket-matched summary_min (Supporting Information Table S8; `HOLDOUT_VERDICT.md`).

| Pair | Holdout summary_min [95% CI] | Main-panel summary_min | Δ (holdout − main panel) | vs best trivial baseline |
|---|---:|---:|---:|---|
| PIK3CA/mTOR | 0.765 [0.603, 0.891] | 0.692 | +0.073 | beats baseline (Δ = +0.21 vs heavy) |
| AChE/BChE | 0.618 [0.422, 0.759] | 0.606 | +0.012 | narrowly beats baseline (Δ = +0.043 vs cLogP); CI still includes 0.5 |
| PIK3CA/PIK3CB | 0.425 [0.241, 0.618] | 0.500 | −0.075 | does not beat baseline (Δ = −0.266 vs heavy) |

PIK3CA/mTOR's direction is confirmed on unseen ligands, with the bootstrap lower bound above 0.5 and continued superiority over the strongest trivial descriptor. AChE/BChE is close to the main-panel point estimate but its CI still spans 0.5. PIK3CA/PIK3CB shows no usable directional signal. One boron-containing ligand (HOAP_028) failed on both ends (unsupported AutoDock atom type) and was excluded (59/60 ligands analyzed). This holdout shares the same ChEMBL extraction batch as the main panels and is not an independent cross-database validation; it tests whether the scoring rule and protocol generalize to same-rule ligands unseen at construction time.

All three holdout pairs show `wrong_pocket_control_vina` **at or above** `pocket_matched_vina` (PM: 0.788 vs 0.765; AChE/BChE: 0.643 vs 0.618; PIK3CB: 0.520 vs 0.425). To test whether this is a Vina-scoring artifact, we computed a scoring-free geometric proxy directly from the already-docked mode-1 pose coordinates: `contact_count`, the number of ligand heavy atoms within 4.0 Å of any receptor heavy atom (a coarse burial proxy, not a validated PLIF). Using `contact_count` alone to repeat the same own-pocket comparison gave AUROC clearly above chance on all three pairs — dual vs A_only in pocket A: 0.552–0.622; dual vs B_only in pocket B: 0.698–0.714 (`WRONG_POCKET_MECHANISM_VERDICT_V1.md`). Dual ligands were also larger on average than the corresponding hard-negative selective ligands (e.g., AChE/BChE mean heavy-atom count: 34.8 vs 33.8/29.6). This indicates that the holdout's wrong-pocket-not-worse pattern is substantially explained by a **ligand size/burial confound that reproduces at the pose-geometry level**, independent of the Vina energy function — a stronger form of evidence than the 2D descriptor covariates in §3.4. Using the identical control definition, the frozen main panels instead show pocket-matched clearly above wrong-pocket for all four pairs (Table S6); we do not have a resolved explanation for this main-panel-versus-holdout contrast and report it as an open discrepancy rather than smoothing it over.

### 3.10 Structural robustness of the alternate crystal forms (cognate QC + panel re-docking)

Cognate re-docking QC on alternate crystals under the Methods 2.4 protocol: PIK3CA **4JPS** (1LT) best-of-9 = 0.607 Å; **5DXT** (5H5) = 0.624 Å; mTOR **4JSX** (Torin2/17G) = 0.515 Å — **all three pass** the < 2 Å gate (`STRUCTURE_ROBUSTNESS_QC_V1.md`); the chimeric structure 3T8M remains excluded.

Re-docking the frozen PM48 ligands with pocket A replaced by 4JPS or 5DXT (pocket B kept at the original 4JT6 scores) dropped pocket-matched summary_min from the main-panel 0.692 to **0.486** [0.259, 0.692] and **0.505** [0.292, 0.696], respectively; the D-vs-A arm (still scored on unchanged 4JT6) stayed at 0.714, so the decline is concentrated in the D-vs-B arm that depends on the new PIK3CA pocket. Replacing pocket B with 4JSX (pocket A kept at 4L23) gave summary_min = **0.639** [0.418, 0.776] (Δ ≈ −0.05). Under the plan's pre-declared decision rule, the PIK3CA end is recorded as **receptor-dependent**: the favorable signal on 4L23 does not automatically transfer to other cognate-QC-passing PIK3CA crystals; the mTOR end retains a weaker version of the advantage after crystal replacement. Full record in `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`.

### 3.11 A structural mechanism for the receptor dependence: inter-crystal conformational variability (exploratory)

To explain the asymmetry in §3.10 — PIK3CA-end swaps collapse the signal, the mTOR-end swap only mildly reduces it — we performed a rigid-body superposition directly on the already-deposited crystal coordinates (`POCKET_MECHANISM_VERDICT_V1.md`). Pocket residues were defined per structure from its own cognate ligand (heavy-atom distance ≤ 5 Å; 20 residues for PIK3CA, 18 for mTOR, matched by residue number **and** identity with zero mismatches, ruling out a crystallization-mutant explanation). All matched Cα atoms were used for a single Kabsch fit, giving a **global Cα RMSD**; the same transform, restricted to the pocket-residue subset, gave a **local pocket Cα RMSD**.

| Reference (main panel) | Alternate | Matched Cα | Global Cα RMSD | Local pocket Cα RMSD | Cognate-ligand centroid distance |
|---|---|---:|---:|---:|---:|
| 4L23 (PIK3CA) | 4JPS | 982 | **1.486 Å** | 0.867 Å | 2.566 Å |
| 4L23 (PIK3CA) | 5DXT | 862 | **1.441 Å** | 0.343 Å | 2.072 Å |
| 4JT6 (mTOR) | 4JSX | 1054 | **0.454 Å** | 0.467 Å | 2.196 Å |

This quantifies the asymmetry seen at the score level: **these deposited PIK3CA structures differ from one another (global Cα RMSD 1.44–1.49 Å) far more than these mTOR structures differ from one another (0.45 Å)**, tracking the direction and relative magnitude of the AUROC collapse. However, local pocket geometry alone does not fully account for it: on 5DXT, the local pocket Cα RMSD (0.343 Å) is *smaller* than the global RMSD (1.441 Å) — the ATP-site backbone is well conserved — yet summary_min still collapsed to 0.505, essentially matching 4JPS (local pocket RMSD 0.867 Å, summary_min 0.486). **Cα-level pocket conservation is therefore not sufficient to guarantee that pocket-matched discrimination transfers**; side-chain rotamers, protonation states, or docking search-space sensitivity not captured by a Cα-only metric are plausible additional factors, but we have no PLIF- or rotamer-level evidence for this round and do not claim the mechanism is fully resolved. Across all three alternates, the transformed cognate-ligand centroids sit within 2.1–2.6 Å of the main-panel cognate-ligand centroid, indicating that docking still targets the same general ATP-competitive site rather than an unrelated pocket.
