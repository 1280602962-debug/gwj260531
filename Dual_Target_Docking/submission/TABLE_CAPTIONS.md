## Table 1. Target pairs and experimental activity classes

靶点对与实验活性类别

Counts are target-pair observations. A molecule may occur in more than one pair. Final activity labels use pChEMBL >= 6.0. Two ineligible observations and the AB_040 alias are excluded.

## Table 2. Dual-ligand recovery in matched top-10% candidate lists

共同成员上前 10% 候选名单的双靶配体检出

Every comparison uses its own common population with valid scores on both targets. k = ceil(0.10 N). S4 versus S1 is a supporting comparison. Class composition and directional changes are given in Table S10.

## Table S1. Evaluation-set construction and activity provenance

评价集构建与实验活性来源

All final labels use pChEMBL >= 6.0. Provenance flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied maxima at both endpoints. Across all pairs, cross-document aggregation is present in 81 observations and endpoint heterogeneity in 123 (union 168). The full source crosswalk, raw activity rows, selection evidence and exceptions are included in Supplementary_Data.zip. CAP5 is a post hoc sensitivity analysis and does not redefine initial selection.

## Table S2. Receptor structures and docking boxes

受体结构与对接盒参数

Coordinates and dimensions are in angstrom. Waters and cocrystal ligands are removed from the receptor. Box centers are bounding-box midpoints; padding is 5 A with minimum edge 20 A. Chain identifiers are author chain identifiers. Formal prepared files and preparation exceptions, including observed-fragment templates, are documented in the accompanying source files. Alternative structures cover only the six conditions used in the current analysis; no Q110 alternative-receptor result is included.

## Table S3. Cocrystal-ligand redocking across five random seeds

五随机种子的共晶配体重对接

Seeds: 17, 29, 42, 71 and 101. Symmetry-aware heavy-atom RMSD is calculated in the receptor coordinate frame. Top-ranked-pose success and success among saved poses are separate outcomes. The summary uses recorded RMSDs without recomputing coordinates; every seed-level record is in Supplementary_Data.zip.

## Table S4. Analysis populations and score availability

分析人群与评分可用性

Panel A counts ligand-direction records. Panel B aggregates the four activity classes; class-specific missingness, exclusions and ligand identifiers remain in the complete source files. D_vs_A means dual versus A-only using target B; D_vs_B uses target A. Missing fold assignments are excluded from chemistry evaluation. Paired method comparisons use their own complete-case populations, reported in Table S5B.

## Table S5. Directional discrimination and paired scoring-workflow comparisons

方向性判别及评分工作流配对比较

AUROC B compares dual with A-only; AUROC A compares dual with B-only. S1 Vina; S2 fixed-pose GNINA; S3 fixed-pose RTMScore; S4 maximum saved-pose CNNscore; S5 independent GNINA. All scores have higher-better direction. Each paired contrast is recalculated on its own shared members; it need not equal subtraction of the method-specific estimates in panel A. CI: 95% class-stratified percentile bootstrap, 10,000 resamples, seed 271828. S4-S2 isolates saved-pose reranking; the supporting S4-S1 link in Table S10 has no added confidence interval.

## Table S6. Negative-class and pocket-correspondence controls

负类组成与口袋对应对照

Panel A difference = neither minus single-target negative class. Panel B difference = corresponding minus noncorresponding pocket. All estimates use S1. Individual AUROC confidence intervals and ligand identifiers are retained in the complete source CSVs.

## Table S7. Out-of-fold performance of ligand-chemistry models

配体化学模型的折外判别表现

All models and S1 are evaluated on identical members within each direction. The increment is AUROC(ECFP4 + S1) minus AUROC(ECFP4). Existing frozen five-fold assignments are retained; Q110 uses sorted canonical Murcko scaffolds assigned cyclically to folds 0-4. Standardization is fitted in training folds. Full descriptors, scaffolds, similarities, folds and OOF predictions are supplied as data files; no new model is fitted during table preparation.

## Table S8. Random-seed and receptor-structure sensitivity

随机种子与受体结构敏感性

Five Vina seeds: 17, 29, 42, 71 and 101. Each seed uses its own valid population. Receptor comparisons use paired common members and replace only one target. Full direction-specific estimates, confidence intervals and members are in the source data. PIK3CA/mTOR Q110 has no alternative-receptor calculation.

## Table S9. Chemical-diversity and activity-threshold sensitivity

化学多样性与活性阈值敏感性

CAP5 is post hoc supporting sensitivity: at most five ligands per nonempty Murcko scaffold within pair and activity class; acyclic ligands are exempt. Six ligands are excluded without replacement. Complete CAP5 method comparisons, chemical models, controls and ranking outputs are retained in Supplementary_Data.zip. The PPARA/PPARD S5 dual-versus-A-only CI changes from [0.4921875, 0.767578125] to [0.5030241935, 0.7772177419]; this boundary change is retained, without a general significance claim. The strict Q110 subset has 56 members (25 dual, 12 A-only, 14 B-only, 5 neither). Main and subset populations differ; no new paired CI for their difference is calculated.

## Table S10. Joint-ranking composition and directional concordance

联合排序候选组成与方向性变化关系

Panel A uses each method's own members and is descriptive. Table 2 and panel B use separate matched populations for every comparison. Joint rank = max(rank A, rank B), with average single-target ranks; ties are resolved by rank sum then global molecular ID. Both differences are new scheme minus S1. Opposite signs define concordant improvement or worsening; equal nonzero signs are discordant; either exact zero defines no change. S4-S1 is supporting, using verified identical directional members and no additional CI. Complete paired composition, per-ligand ranks and all CAP5 ranking results are supplied as data files.
