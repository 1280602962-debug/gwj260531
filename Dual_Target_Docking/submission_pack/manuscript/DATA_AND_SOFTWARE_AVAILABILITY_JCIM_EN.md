# Data and Software Availability

## Data and Software Availability

Benchmark membership, experimental-state labels, receptor and docking-box definitions, per-ligand docking scores, analysis tables, and all scripts used to regenerate the reported statistics and figures are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531.

Typeset Supporting Information contains **Tables S1–S13**. The tables cover software and seeds; docking boxes and cognate RMSD; threshold and pChEMBL-aggregation sensitivity; fixed-score control-class contrast; ECFP4 increment; matched versus mismatched pockets; unused-pool holdout; PIK3CA/mTOR crystal substitution; independent GNINA and five-seed Vina; document-cluster uncertainty and detectable-effect simulation; BindingDB external eligibility; 2018 year-split counts; and EGFR/HER2 operating points. The map from the former working Tables S1–S54, including the old multi-seed Table S54, is `Dual_Target_Docking/data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv`.

Per-ligand scores, holdout membership, multi-seed long tables, property-caliper matches, chemotype-matched selectives, aggregation means, complete-case and assay-context ledgers, the J0 candidate-pair census, historical BindingDB REST counts, leave-cognate-out, occupancy snapshots, contact counts, whole-chain sequence identity, and the MCL1/Bcl-xL applicability archive (not Table 2) are listed in **Note S14**. These files will ship with code and SHA-256 checksums in a GitHub Release. A Zenodo DOI will be issued from a tagged snapshot, not from this moving branch.

`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` indexes the principal numerical results and their source CSVs. SHA-256 checksums of manuscript-facing tables are in `REVISION_CHECKSUM_MANIFEST_v1.csv`. The native-slice contract is `protocol/external_slice_contract.yaml`; the evaluation contract is `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`. The ChEMBL supply audit was frozen on 2026-07-23; the high-confidence activity view was fetched on 2026-08-26; the BindingDB-native archive lock is release 202609. The analysis environment and zero-docking reproduction commands are documented in the repository README. The BindingDB TSV archives themselves are not redistributed; CI checks committed CSVs only.

