# Figure captions (JCAMD)

Lettering is 8 pt sans-serif (Liberation Sans; Arial/Helvetica metric-compatible substitute — Arial is not licensed on the builder). All numbers below are copied from archived files under `data/`; they are the same values asserted by `scripts/plot_jcamd_publication_figures.py`. Panel letters sit outside the data area. Captions are **not** drawn on the artwork.

## Figure 1. Locked URAT1 docking readout is a weak activity retriever and is not pose-accurate.

**a** Enrichment factor at the top 1% (EF@1%) for protocols P0–P5 on the pre-registered TrueDecoy and RandomDecoy benchmarks (9DKB). Bars are point estimates; whiskers are molecule-resampled bootstrap 95% percentile intervals (1,000 draws) from `data/si/protocol_enrichment_ci/protocol_ef_ci.csv`. The dashed line is chance (EF = 1). P2 (gnina CNNaffinity) is the production readout locked by the pre-registered rule (green band): TrueDecoy EF@1% = 2.59 (12/51, hypergeometric *p* = 0.0016), RandomDecoy EF@1% = 0.22 (1/51, 95% CI 0.00–1.04). P0 is the pre-registered negative-control readout (CNNscore), not a salvage protocol. **b** lesinurad self-docking on 9DKB at exhaustiveness = 32, `num_modes` = 9 (`data/redock_smoke/redock_results_lesinurad_9DKB.csv`). Top-1 heavy-atom RMSD versus the best RMSD in the nine-pose ensemble. The dashed line is the 2 Å Top-1 gate. P2 Top-1 RMSD = 4.16 Å (fails the gate); every protocol’s ensemble minimum is ≤ 1.00 Å.

## Figure 2. Transfer of the locked P2 readout onto the clinical library is an audit, not a dual-node nomination.

**a** Compound counts at each frozen funnel step: ChEMBL clinical library *n* = 8,319; NLRP3 ML shrink q_N ≥ 0.5, *n* = 1,588; P2 dual-dock complete-case table *n* = 1,580 (1,579 rows with valid dual scores; one empty-pose row, tauroselcholic acid, retained with percentile 0); dual-structure gate S_U ≥ 90 and S_N,dock ≥ 90, *n* = 51; chemistry-filtered list (Veber + Ro5 HBD/HBA/logP, MW 200–550 Da, macrolide demotion), *n* = 7. **b** S_U versus S_N,dock for the 1,579 valid dual-score rows. The dashed box is the τ = 90 gate. Red diamonds: raw docking Pareto front (*n* = 4; Idremcinal, Alemcinal, Cethromycin, Zamzetoclax — macrolide/erythromycin audit, not a follow-up list). Blue circles: chemistry-filtered list (*n* = 7). Black squares: gout-related drugs that are in the complete-case table (lesinurad, verinurad, colchicine). Names are omitted from the scatter so they do not cover points; identities are in Fig. 3, Fig. S4, and `figures/generated/tables/`.

## Figure 3. Production P2 poses sit in the crystal cavity but lesinurad has lost the Arg477 salt bridge.

**a** Shortest carboxylate-oxygen to Arg477 distance in the production first pose, reported only for molecules with an archived `acid_arg477` value (`data/si/pose_qc/pose_qc_table.csv`): lesinurad 14.20 Å, verinurad 2.86 Å, GSK-3008348 3.19 Å. The other six chemistry-filtered molecules are not carboxylic acids and have no value in that column. **b** URAT1 centre-of-mass displacement relative to co-crystal lesinurad for the seven chemistry-filtered molecules plus the two URAT1 controls. The dashed line is the 6 Å COM in-pocket cutoff used in pose QC. **c** NLRP3 COM displacement relative to co-crystal NP3-146/RM5 for the same seven molecules. lesinurad and verinurad were not archived on 7ALV. All seven have `both_in_pocket = True` and zero 2.2 Å clashes (`pose_qc_dual.csv`). These geometries do not constitute binding-mode or affinity evidence.

## Figure 4. The TrueDecoy active set is analog-biased and excludes textbook URAT1 drugs; RandomDecoy is not a near-neighbour leak.

**a** Membership of named URAT1-related compounds in the 469-compound p≥6 TrueDecoy active set, by ChEMBL identifier (library IDs plus the identifiers listed in `docs/DATA_FACT_CHECK.md`; isobavachin from `data/raw/URAT1_CHEMBL_cf12.csv`). lesinurad, benzbromarone, dotinurad and probenecid are absent; verinurad, puliginurad, SHR-4640 and isobavachin are present. **b** The five most frequent Murcko scaffolds among the 469 actives (118 unique scaffolds). Scaffold 1 accounts for 127/469 (27.1%); SMILES are in `tables/table_s_top_scaffolds.csv`. **c** Maximum Tanimoto similarity of each of 4,690 RandomDecoy molecules to any TrueDecoy active. None exceed 0.5; scaffold overlap with actives is 0 (`data/si/decoy_leakage_audit/`). **d** The same metric for 80 experimental weak actives (designed hard negatives): 14 have TC > 0.5 and none have TC > 0.85.

## Figure S1. Gate-threshold sensitivity (does not replace τ = 90).

Number of compounds passing the dual-structure gate and the subsequent chemistry filter as a function of percentile threshold τ (`data/si/nomination_sensitivity/gate_counts.csv`). The production analysis is locked at τ = 90 (51 → 7). Wider gates are monotonic sensitivity, not a second production shortlist.

## Figure S2. ROC-AUC for P0–P5 on TrueDecoy and RandomDecoy.

Point estimates and bootstrap 95% percentile intervals from the same file as Fig. 1a. The dashed line is AUC = 0.5. P2 TrueDecoy AUC = 0.580. Absolute discrimination remains modest; AUC is not the locking metric.

## Figure S3. NLRP3 assay-conditioned classifier, scaffold-grouped CV.

Per-fold out-of-fold AUROC and AUPRC from `docs/MODEL_TRAINING_SUMMARY.json` (five folds). Dashed lines: pooled AUROC = 0.893, AUPRC = 0.914. EF@10% = 1.57 at prevalence ≈ 60%. No ROC/PR curve is drawn because the OOF prediction table is not in this archive.

## Figure S4. Known gout-related drugs that reached the dual-dock table.

P2 percentiles inside the archived 1,580-row denominator. None meet S_U ≥ 90 and S_N,dock ≥ 90. benzbromarone, dotinurad, allopurinol, febuxostat and probenecid have q_N = 0 and never entered the 1,588 docking pool (not plotted).

## Figure S5. NLRP3 q_N on the 8,319-compound clinical library.

Histogram of predicted P(active). The dashed line is the pre-registered shrink threshold 0.5 (*n* = 1,588).
