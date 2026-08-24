# Methods (JCIM Articles draft, English)

> Companion to [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md) (Chinese authoritative for this rewrite).  
> Protocol only: counts, cognate RMSDs, AUROCs, and holdout point estimates belong in Results / SI.  
> DualFourClass-Bench is a **four-state curated benchmark with two directional primary tasks**, not a four-class classifier. Call it a curated four-pair panel + evaluation protocol, not a comprehensive suite.

---

## 2. Methods

### 2.1 Data sources and activity curation

Ligand activities used as **experimentally derived activity labels** were retrieved from the public ChEMBL Web API activity endpoint. The target-pair supply audit was frozen on 2026-07-23. pChEMBL converts molar concentration–response measurements (IC50, Ki, Kd, EC50, and related endpoints) to an approximate −log10 activity scale for large-scale integration. Assay types, conditions, and experimental systems are not equivalent; pChEMBL is used here as a curation convenience, not as an absolute affinity measured under one protocol.

When several pChEMBL values existed for the same ligand–target pair, the frozen tables store the **maximum** as the one-to-one representative used for the primary curation. Assay types, conditions, and experimental systems are not equivalent; taking the maximum can inflate a single-assay reading. Activity-aggregation sensitivity was therefore assessed as a **prespecified sensitivity analysis** by re-aggregating all available ChEMBL measurements using the **median** rather than the maximum pChEMBL value. This analysis was performed across all frozen benchmark panels after re-fetching assay-level records from the ChEMBL activity endpoint (Table S29). It was not used to redefine panel membership or docking parameters. Frozen Vina scores were not recomputed. Class assignment was compared at the same θ = 6.0 rule. The frozen `mols_*.json` files still contain only the representative float; median labels live in the A4 tables, not as a rebuilt primary panel. Ligands missing a usable pChEMBL on either target were excluded from analyses that require paired labels.

Salt, solvate, and multicomponent ChEMBL records were split by connected component; the organic fragment with the most heavy atoms was retained.

To test whether the ChEMBL supply gate is an artifact of one database, the frozen target pairs were recounted in BindingDB and PubChem (counts only; no docking and no panel rebuild). BindingDB was queried with REST `getLigandsByUniprots` (cutoff = 1 mM, so weak-end measurements are not truncated); PubChem with PUG REST `protein/accession/…/concise`. Endpoints were restricted to IC50/Ki/Kd/EC50; the representative value was the maximum converted p-activity; class rules matched the strict supply gate in Section 2.2. Identifiers were BindingDB monomerid and PubChem CID, with **no** cross-database structure merge. The primary count used equal-relation measurements (censored `>`/`<` values dropped); treating inequalities as point estimates was a sensitivity only. Counts are reported in Supporting Information Table S12.

### 2.2 Target-pair supply audit and experimental ligand-state definition

Candidate pairs were audited for whether public data can support a strict dual-target benchmark. For each pair A/B, ligands were assigned one of four **experimental states**:

- **dual**: strong activity on both targets;
- **A-only**: strong on A, weak on B;
- **B-only**: strong on B, weak on A;
- **neither**: insufficient activity on both.

A-only and B-only ligands are selectivity hard negatives, not DUD/DUD-E-style assumed decoys.

**Strict supply-audit rule (construction gate, not the sole label for every comparison).** Dual: both pChEMBL values ≥ 6.5. A-only: A ≥ 6.5 and B ≤ 5.5. B-only is symmetric. Neither: both ≤ 5.5. The 5.5–6.5 gray zone was excluded from the strict audit. The rule asks whether both arms have enough selective hard negatives for a reasonably balanced panel. How many pairs pass, and which pairs entered the frozen set, are reported in Results 3.1. Metal-dependent systems (e.g. HDACs) were excluded in advance as unsuitable primary objects for this noncovalent protocol.

**The primary manuscript comparison uses one prespecified θ = 6.0 rule.** Dual: both ends ≥ θ; A-only: A ≥ θ and B < θ; B-only is symmetric; neither: both < θ. Construction may use this single threshold when the strict rule leaves too few selectives to fill quotas. Construction rules were frozen from the supply audit before sampling and are recorded in Table 1. Thresholds were chosen to assemble analyzable quotas, not after inspecting docking scores. As a supporting sensitivity, states were relabeled at θ ∈ {5.5, 6.5} and under the strict 6.5/5.5 rule, and pocket-matched summary_min was recomputed (Table S4). That grid is not a second primary standard. Underpowered cells are flagged in Results.

### 2.3 DualFourClass-Bench panel construction

The resource is a **four-state curated benchmark with two directional primary tasks**. Dual, A-only, B-only, and neither are all retained to describe the experimental space; the prespecified primary endpoint is dual versus A-only and dual versus B-only. Neither does not enter the primary directional AUROCs. This is not a four-class classification benchmark.

Candidate pairs were screened with the strict audit in Section 2.2. The frozen evaluation set comprises PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, and EGFR/HER2. EGFR/HER2 is retained as a **supply-limited case** (`PAIR_ROLES_APPROVED_JCIM.yaml`) and is not treated as equivalent in supply to the other pairs.

For each pair, ligands were drawn from the labeled pool under frozen class quotas and random seed 20260729. Where Bemis–Murcko scaffolds could be computed at draw time, a per-class scaffold cap limited series over-representation: at most two molecules per scaffold in PIK3CA/mTOR (PM48) and at most five in EGFR/HER2. For AChE/BChE and PIK3CA/PIK3CB, SMILES were not yet in the sampling table, so a Murcko cap could not be applied; sampling used class quotas and a deterministic shuffle only. **No additional chemical-diversity constraint was applied.** Post-construction Murcko scaffolds, when later available, are reported with the deposited tables. Final membership, state labels, ChEMBL identifiers, SMILES, and sampling scripts are deposited. Panels were not redrawn after docking scores were seen.

Construction rules were not identical across pairs. AChE/BChE and PIK3CA/PIK3CB were sampled under the strict 6.5/5.5 gate; EGFR/HER2 and PIK3CA/mTOR used θ = 6.0 because the strict gate left too few B-only ligands. Cross-pair AUROCs therefore mix target-pair biology with panel-construction differences (sample size, threshold, series composition, receptor) and are not interpreted as purely intrinsic docking performance.

Quotas and construction labels were as follows. AChE/BChE and PIK3CA/PIK3CB: strict 6.5/5.5, target dual / A_only / B_only / neither = 28 / 28 / 28 / 16 (panel n = 100). EGFR/HER2: existing θ = 6.0 panel (n = 110). PIK3CA/mTOR: θ = 6.0, main comparison panel PM48 (n = 48; constructed 18 / 14 / 12 / 4), on which receptors and the docking protocol were frozen.

Ligand–receptor jobs that failed to yield a score were dropped for that receptor; ligands missing a usable score on either end were omitted from pocket-matched AUROCs that require both scores, so analysis counts can fall below construction quotas (Table 1). AUROC tables are therefore **conditional on compounds the docking engine can process**. Attempted / successful / failed counts, including chemical-coverage failures such as unsupported AutoDock atom type `B`, are reported in Table S27.

An expanded PIK3CA/mTOR panel (historical name PM110) keeps all 48 PM48 ligands and adds molecules under the strict rule, targeting 30 / 30 / 30 / 25. PM110 is a superset of PM48, used to check whether point estimates stay in the same direction after increasing panel size. It is not an independent primary benchmark and not an independent replicate. Cross-pair comparison in the main text uses PM48.

This work does not treat a distribution of non-overlapping re-drawn balanced panels as the robustness readout; leftover hard-negative supply is too limited (quantified in Results). The formal ligand-side check is one unused-pool holdout (Section 2.11). Ligand-level bootstrap (Section 2.8) describes uncertainty inside a fixed panel and is not unused-pool resampling.

**Table 1.** DualFourClass-Bench composition and docking settings (construction rules)

| Pair | Construction labels | PDB (A / B) | Resolution (Å) | Panel n | Analysis n (dual / A_only / B_only) | Vina exhaustiveness |
|------|---------------------|-------------|----------------|-------:|------------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | strict 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.4 Protein structures and binding-site definition

Receptors were PDB entries with experimental structures and a small-molecule cognate ligand: PIK3CA/mTOR, 4L23 / 4JT6 (X6K / PI-103); AChE/BChE, 4EY7 / 4BDS (E20 / THA); PIK3CA/PIK3CB, 4L23 / 2WXF (X6K / 039); EGFR/HER2, 3POZ / 3RCD (03P / TAK-285). Resolutions are in Table 1.

The site was defined from the cognate ligand. An axis-aligned bounding box on cognate heavy atoms was expanded by 5 Å on each axis; any edge shorter than 20 Å was set to at least 20 Å. Box centers and sizes are frozen in JSON and listed in Supporting Information Table S2.

Water and the cognate ligand were removed and Meeko wrote PDBQT. PIK3CA, mTOR, EGFR, and HER2 used hydrogen-containing protein coordinates already in the frozen directories (`mk_prepare_receptor.py --read_pdb`). AChE, BChE, and PIK3CB were extracted from deposited ATOM/TER records (waters and hetero atoms removed) and prepared with `mk_prepare_receptor` (default alternate location A). PDBFixer was not used to rebuild missing atoms, nor Reduce for independent pH-dependent protonation or histidine tautomer enumeration; protonation is part of the frozen preparation. Docking treated noncovalent small-molecule sites only; metals and other cofactors were not extra dockable components in the box.

### 2.5 Cognate redocking quality control

Before production docking, each frozen receptor was redocked with its cognate ligand to test whether the box, receptor preparation, and search settings can generate a near-native pose **within the retained pose ensemble**.

Nine poses were generated and heavy-atom RMSD to the crystal ligand was computed in the docking frame (no protein superposition). PIK3CA/mTOR and EGFR/HER2 used meeko `REMARK SMILES IDX` mapping and the minimum CalcRMS over graph automorphisms; AChE/BChE and PIK3CB used Hungarian matching on heavy atoms. Define

\[
\mathrm{RMSD}_{\mathrm{best9}} = \min_{i=1,\ldots,9} \mathrm{RMSD}_i.
\]

The prespecified pass criterion was \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å.

This QC tests **pose-generation capability**: whether a near-native pose appears among the retained poses. It does **not** require the top-ranked Vina pose (mode 1) to be near-native. Best-of-9 QC and mode-1 scoring are different evaluations.

If the default exhaustiveness failed the gate, search effort was raised to a prespecified fallback without changing the box, receptor, or random seed, and QC was repeated. Production docking therefore used receptor-specific frozen exhaustiveness: 16 for PIK3CA/mTOR and 8 for the other main panels. QC values and mode-1 versus best-of-9 comparisons are in Supporting Information Table S3.

### 2.6 Ligand preparation and molecular docking

Ligands started from frozen ChEMBL SMILES: desalt to the largest organic fragment, add explicit hydrogens in RDKit, embed with ETKDGv3 (seed 20260727), locally optimize with MMFF (at most 200 steps), and convert with default Meeko to PDBQT. Protonation states, tautomers, and conformational ensembles were not systematically enumerated. Schrödinger LigPrep was not used.

Docking used AutoDock Vina 1.2.7 with the default `vina` scoring function. Nine poses were retained per ligand–receptor pair, `energy_range = 3` kcal mol\(^{-1}\), random seed 20260727 (same as ETKDG). Exhaustiveness followed Table 1. Ligand preparation, box rules, and the scoring function were shared across main panels; only receptor coordinates, box numbers, and the prespecified exhaustiveness varied. Full parameters are in Supporting Information Table S1.

### 2.7 Alternative scoring channels

To test dependence on one scoring function, the **same Vina-generated poses** were rescored with RTMScore and GNINA CNN.

RTMScore used public weights `rtmscore_model1` and scored all nine Vina poses; the highest RTMScore in that pocket was kept.

GNINA 1.3.2 ran CNN rescoring on CPU (`--cnn_scoring rescore --minimize`). The final protocol converted each of the nine Vina poses to SDF (Open Babel) and took the highest CNNscore per pocket, matching RTM pose coverage. Mode-1-only GNINA scores are retained as a historical sensitivity control, not the channel readout.

Vina’s primary readout is the mode-1 energy; RTM and GNINA are best-of-9 rescores. The three channels do not aggregate the nine poses in the same way and are **not** a head-to-head docking-engine competition. They are a scoring-channel sensitivity analysis. The primary endpoint is always Vina.

### 2.8 Primary endpoint and statistical analysis

#### 2.8.1 Pocket-matched directional AUROC

Two binary AUROCs were computed per pair. Dual versus A-only used the **pocket B** score:

\[
\mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B),
\]

asking whether docking can use structural information at the nonselective target B to separate dual-actives from A-only ligands that are already potent at A. Dual versus B-only used the pocket A score:

\[
\mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A).
\]

Dual is always the positive class. Neither does not enter these contrasts.

Vina reports \(E_{\mathrm{Vina}}\) (kcal mol\(^{-1}\); more negative is more favorable). Define

\[
S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}
\]

so that larger primary scores mean stronger predicted binding. RTMScore and GNINA CNN scores are already higher-better.

#### 2.8.2 summary_min

The pair summary is the weaker arm:

\[
\mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}).
\]

This is a task-aligned **worst-arm aggregation**, not a new scoring function. The minimum prevents a strong arm from hiding failure on the other; it is not the unique statistically natural aggregator. Arithmetic mean and harmonic mean of the two directional AUROCs are reported as a sensitivity (Table S26). Pair ranking and the EGFR Dual-versus-neither contrast remain in the same direction under all three aggregators. The single primary endpoint is pocket-matched Vina `summary_min` under unified θ = 6.0 (Table 2; PIK3CA/mTOR uses PM48). Prespecified secondary endpoints are the two directional arms, pocket-matched RTMScore, pocket-matched GNINA CNN best-of-9, and the descriptor panel in Section 2.8.3. Robustness / falsification endpoints are the θ grid, PM110, E = 8, unused-pool holdout, receptor replacement, and wrong-pocket (including paired Δ). Exploratory endpoints are ECFP4, contact_count (not a PLIF), and Top-10 hard-negative counts on pooled `vina_mean`. The hierarchy is in Supporting Information Table S16. Pooled `vina_mean` directional AUROC is **not** Table 2.

#### 2.8.3 Physicochemical descriptor controls

A **prespecified** RDKit panel was computed: heavy-atom count, molecular weight, cLogP, and TPSA. Each descriptor was evaluated with the same directional AUROC workflow; **all four are reported** (Table 2; Table S28). The highest AUROC among them is a **best single-descriptor reference** — a descriptive post-hoc maximum, **not** a confirmatory competitor and not a “trivial baseline” hypothesis test. Paired Δ between docking and that reference is not a confirmatory test of “beats the selected best descriptor” (Table S19).

#### 2.8.4 Score-aggregation controls

Pooled means of the two pocket scores, wrong-pocket assignment (Section 2.9.1), and worst-pocket aggregation were computed as auxiliaries, not as the primary endpoint (Table S6).

#### 2.8.5 Bootstrap uncertainty

AUROC and summary_min uncertainty used ligand-level bootstrap: ligands were resampled with replacement, preserving class structure, and both directional AUROCs and summary_min were recomputed. \(B = 2000\), seed 20260729, percentile 95% CI \([P_{2.5}, P_{97.5}]\). Paired contrasts used the **same** resample to form \(\Delta = \mathrm{Metric}_1 - \mathrm{Metric}_2\) (Tables S17, S19). Murcko-scaffold resampling is reported as a control; the text uses ligand-level intervals. Intervals are descriptive. Outside the prespecified primary endpoint, this work does not treat “whether the CI crosses 0.5” as a formal significance test across many pairs and controls.

#### 2.8.6 Benchmark-formulation comparison

As an auxiliary contrast on the **same** frozen Vina scores, a **Dual-versus-neither comparator** (experimental inactives; `vina_mean` and `vina_worst`) and Dual versus all non-duals were computed beside the directional primary endpoint. Dual-versus-neither is a **nonselectivity-controlled comparator** on this panel, not a claim that established dual-target benchmarks use Dual versus neither as their official task. Neither ligands are used here; they still do not enter Table 2. PIK3CA/mTOR neither n = 4 is flagged underpowered. The comparison asks whether omitting selective hard negatives can change the apparent evidence for dual-target recognition; it is not a second primary endpoint and is not a paired significance test (different negative sets; Table 3; Table S22). Single-target analogues—(dual + A-only) versus (B-only + neither) in pocket A, and the symmetric B contrast—are reported only as a Zhou-like backdrop.

### 2.9 Confounder and falsification analyses

#### 2.9.1 Wrong-pocket falsification control

Scores for targets A and B were swapped; ligands, receptors, and all other settings were unchanged. Directional AUROCs and summary_min were recomputed. This is a **falsification control**, not a positive control designed to prove pocket specificity. Matched > wrong on a fixed panel is **not** taken as evidence of pocket-specific signal. Wrong-pocket performance near or above matched-pocket performance counts against a pocket-specific reading. Holdout reversal further means wrong-pocket is **not a reliable universal negative control under panel shift**.

#### 2.9.2 Ligand-efficiency normalization

Each pocket score was divided by heavy-atom count, \(S_{\mathrm{LE}} = S_{\mathrm{dock}} / N_{\mathrm{heavy}}\), and directional AUROCs and summary_min were recomputed.

#### 2.9.3 Potency- and size-matched subsets

Subsets with \(|\Delta\mathrm{pChEMBL}| \leq 0.5\) or \(|\Delta N_{\mathrm{heavy}}| \leq 2\) were formed and directional AUROCs recomputed. Matching reduces n; the analysis asks whether the direction changes, not whether a small subset is independent strong evidence (Table S5).

#### 2.9.4 Covariate-adjusted analysis

Logistic regression compared

\[
\mathrm{Model}_1:\ Y \sim S_{\mathrm{dock}}, \qquad
\mathrm{Model}_2:\ Y \sim S_{\mathrm{dock}} + N_{\mathrm{heavy}} + \mathrm{TPSA},
\]

where \(Y\) is the dual versus selective-hard-negative label. scikit-learn `LogisticRegression` was used (\(C = 1.0\), `max_iter = 2000`). Model AUROC, the docking coefficient, and its odds ratio are reported. The question is residual discrimination after size and polarity, not a new primary predictor.

#### 2.9.5 Two-dimensional chemical baseline

Morgan/ECFP4 fingerprints (radius 2, 2048 bits) with the same logistic settings provided a ligand-only chemical baseline. Evaluation used Bemis–Murcko scaffold `GroupKFold` with \(K = \min(5, N_{+}, N_{-}, N_{\mathrm{scaffold}})\) and at least two folds, so the same scaffold does not span train and test. High CV AUROC therefore means discrimination remains when molecules from the same Murcko scaffold are not shared between folds; it is **not** target-external generalization. On PIK3CA/mTOR, \(n_{\mathrm{scaffolds}} \approx n\), so the split is nearly leave-one-scaffold. Random `StratifiedKFold` is a leakage check only (Table S20). Incremental models (physchem, ECFP4, docking, and combinations) use the same split; logistic docking AUROC is not the rank AUROC in Table 2 (Table S24). Nearest-neighbor ECFP4 Tanimoto matching of A-only/B-only ligands to duals is reported at T ≥ 0.3 / 0.4 / 0.5 because T ≥ 0.7 matching is empty on these panels (Table S23). T ≥ 0.3 is a **similarity-constrained subset**, not a chemically matched analogue set.

#### 2.9.6 Scoring-independent contact count

On frozen Vina **mode-1** poses, a scoring-free geometric descriptor was the number of ligand heavy atoms within 4.0 Å of any receptor heavy atom:

\[
N_{\mathrm{contact}} = \#\{i:\ \min_j d_{ij} \le 4.0\,\text{Å}\}.
\]

This descriptor does not use the docking energy function. \(N_{\mathrm{contact}}\) was used for dual versus A-only in pocket A and dual versus B-only in pocket B, isomorphic to the wrong-pocket pocket-wise comparison, as a geometric confounder control (size/burial). The 4.0 Å cutoff is coarse and **not** a PLIF. Magnitude agreement with Vina wrong-pocket is not assumed (Table S11).

#### 2.9.7 Cross-pair sequence identity (exploratory)

The longest protein chain was read from each frozen `*_protein.pdb` (standard amino-acid ATOM records). Pairwise global identity used Biopython `PairwiseAligner` (BLOSUM62, gap open = −11, extend = −1), normalized by alignment length or by the shorter chain (Table S7). The quantity is a coarse whole-chain similarity proxy, not pocket residue correspondence, pocket RMSD, or a PLIF.

### 2.10 Single-target enrichment reference

On PIK3CA 4L23 and mTOR 4JT6, single-target active versus weak-active sets were built. Actives: pChEMBL ≥ 6.5. Weak actives: measured on the same target with pChEMBL ≤ 5.5 and property-matched to actives within ±50 Da (MW), ±1.5 (cLogP), and ±25 Å² (TPSA). MW and logP windows follow common property-matched decoy practice (Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594); TPSA is an added polarity match. Target size was about 50 actives and 150 weak actives. Preparation, receptor, box, and Vina settings matched the PIK3CA/mTOR main panel (exhaustiveness = 16). AUROC, EF1%, and EF5% are reported as a single-target enrichment backdrop, not as a substitute for dual-target summary_min.

### 2.11 Unused-pool holdout

To test dependence on the exact frozen panel members, all ChEMBL entries used in the main panels and in PM110 were excluded from the strict-label pool, and an **unused-pool, panel-external holdout** was drawn from what remained. This is not cross-database or cross-assay external validation: ligands still come from the same ChEMBL harvest, the same target pairs, and the same label rules.

Holdout was built only for pairs whose unused pool could supply 20 dual, 20 A-only, and 20 B-only ligands: PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB. EGFR/HER2 was not eligible for an equivalent draw and was not patched with a non-equivalent sample. PIK3CA/mTOR exclusion used the PM110 superset, which covers PM48. The draw used `HOLDOUT_SEED = 20260731` (distinct from the construction seed) and a Murcko cap of three members per state class. The list was frozen before docking.

The holdout did not enter panel construction, protocol tuning, or primary-endpoint choice. Receptor, box, ligand preparation, exhaustiveness, scoring, and statistics matched the main benchmark, including `summary_min` and ligand-level bootstrap. Jobs that produced no Vina score were dropped as in Section 2.3. Descriptor, wrong-pocket, matched-subset, and contact-count analyses were repeated on the holdout (Tables S8, S13). Matched-subset diagnostics do not rewrite Table S8.

### 2.12 Receptor-structure sensitivity analysis

To test sensitivity of the benchmark conclusion to receptor choice, alternate crystals were required, **before scores were seen**, to (i) have a polymer entity that is the true target protein (no chimeras or off-target scaffolds), (ii) contain a small-molecule cognate in the ATP or target site, (iii) have acceptable resolution, and (iv) pass the same cognate redocking QC as Section 2.5. Structures actually docked were PIK3CA 4JPS and 5DXT and mTOR 4JSX. This is a **receptor-structure sensitivity analysis** (a receptor-realization effect), not a robustness check and not a certification of a “more correct” crystal. The purpose was to quantify sensitivity of the dual-target discrimination endpoint to receptor realization rather than to identify a superior receptor structure. PIK3CA/mTOR is not treated as a prespecified structure-invariant positive case.

Replacement was **one pocket at a time**. On PIK3CA/mTOR (PM48), 4JPS/5DXT replaced pocket A while pocket B kept frozen 4JT6 scores; 4JSX replaced pocket B while pocket A kept frozen 4L23 scores. Exhaustiveness was 16, matching the PM48 main panel. On PIK3CA/PIK3CB, the same prepared 4JPS and 5DXT receptors replaced pocket A while pocket B kept frozen 2WXF scores; exhaustiveness was 8, matching that main panel. New boxes used that crystal’s own cognate ligand and the Section 2.4 AABB rule. Ligand preparation, seed (20260727), scoring function, and the primary endpoint matched the corresponding main analysis. Jobs that produced no Vina score were dropped as in Section 2.3; attempted / successful / failed counts are reported with the swap tables.

As an exploratory, docking-free geometric control, rigid superposition used Biopython `Superimposer` on residue-matched Cα atoms of the longest chain (Kabsch). Pocket residues were those within 5 Å of the reference cognate heavy atoms; local pocket RMSD used the **same** transform, with no second local fit. Alternate cognate centroids were projected in that frame. Matched Cα counts can differ, so global RMSDs are not equal-coverage comparisons. The set of alternates is small; Cα RMSD is not assumed to explain AUROC changes quantitatively (Table S10).

### 2.13 Software and data availability

Analyses ran under Python 3 with RDKit 2026.3.1, meeko 0.7.1, AutoDock Vina 1.2.7, GNINA 1.3.2, and RTMScore (`rtmscore_model1`); Open Babel converted Vina poses to SDF. Superposition and chain alignment used Biopython. AUROC, logistic regression, and cross-validation used NumPy, SciPy, scikit-learn, and pandas. Panels, scores, scripts, and parameter tables will be released with the public pack (Data and Software Availability).
