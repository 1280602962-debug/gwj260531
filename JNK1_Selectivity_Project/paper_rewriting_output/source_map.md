# Source Map — Materials → Manuscript Units

| Source ID | Local path | Content summary | Maps to manuscript unit | Claim it can support | Must not over-claim |
|-----------|------------|-----------------|-------------------------|----------------------|---------------------|
| REF007 | `docs/JNK1_PROJECT_REPORT.md` | Full project narrative, funnel, pivot, MD, purchase | Intro / Methods / Results / Discussion backbone | Pipeline design; selectivity-method failures; purchase rationale | “Selective hit discovered” |
| REF006 | `docs/REFERENCES.md` | Curated literature + PDB table | Intro / Related structural background | Field mechanisms (Ile106/Leu; Gly87 trial) | Invented PDBs |
| REF034 | `results/model_comparison/MODEL_COMPARISON_REPORT.md` | XGBoost vs Chemprop | Results — ML | Holdout/CV R²; model choice | Isoform selectivity from ML |
| REF067 | `results/ml_external_validation/*` | Decoy FPR / EF1% | Results — ML specificity | F1@6.0 high recall, low specificity | “ML is selective” |
| REF069 | `results/docking_validation/validation_report.md` | Redock + Δsel calibration | Results — docking | Pose grids OK; Δsel ~43% | Docking selects isoform |
| REF068 | `results/docking_validation/benchmark_mmgbsa_calibration.md` | MM-GBSA thresholds | Methods / Results | Potency gate calibration | Δsel_MMGBSA selectivity |
| Table27 | `docs/popular_science/data_tables/27_MD16_选择性排序与报价.csv` | MD bias scores + HIT IDs + prices | Results — shortlist / purchase | 690/2157 IDs; MD ranks | MD = selectivity |
| Purchase | `data/purchase/purchase_after_md.csv` | 10-compound purchase design | Methods — experimental design | Assay panel plan | Actual IC50 (none yet) |
| Shortlist | `data/shortlist/md_shortlist_final.csv` | 16 MD compounds + QikProp | Methods / Results — ADMET/MD input | HIT103871685=690; HIT101201113=2157 | |
| C1 | `results/chemotype_novelty/` | ECFP4/Murcko novelty for 690/2157 | Results — chemotype | maxTc vs refs/ChEMBL | Fingerprint ≠ pharmacophore proof |
| C4 | `results/assay_analysis/` + `scripts/c4_preregistered_ic50_analysis.py` | Locked IC50/SI rules | Methods / Results — wet-lab | RQ-A/RQ-B endpoints | Post-hoc SI redefinition |
| C5 | `results/selectivity_autopsy/` | Δsel/Gly87/ML failure table | Results — method autopsy | Purchase decoupling | Claiming predictors worked |
| C2/C3 proto | `docs/protocols/C2_C3_pose_md_replica_protocol.md` | Redock/MD replica plan | Methods | Pose credibility | MD = selectivity |
| C2 Vina | `results/pose_consensus/` | Multi-seed Vina consensus | Results — pose QC | 690 pass-all; 2157 JNK2 weak | Vina score = Glide rank |
| C7 | `results/purchase_risk/` | PAINS/physchem | Methods / SI | No PAINS on 690/2157 | Alert ≠ artifact proof |
| C11 | `results/c11_2231_comparison/` | Unbought 2231 table | Discussion | Opportunity cost | 2231 inactive claim |
| Blueprints | `paper_rewriting_output/section_blueprints.md` | Section plan | Writing | Option A structure | |
| Draft | `paper_rewriting_output/draft_intro_methods_rqc_en.md` | EN draft Intro/RQ-C | Writing | Paste-ready prose | IC50 overclaim |
| MD2231 | `results/md_2231_200ns/` | Extended MD for 2231 | Results — optional / SI | Directional JNK1 RMSD trend | Confirm selectivity; note **not purchased** |
| Figs | `docs/popular_science/figures/` | Funnel, redock, direction, hinge | Main figures | Visual evidence of funnel + failures | |
| Workflow | `docs/JNK1_selectivity_screening_workflow.md` | Original design intent | Intro historical design | Why JNK1 was chosen | Current computational success |

## Purchase-critical ID check (locked)

| Library ID | HIT ID | Role in paper |
|------------|--------|---------------|
| 690 | HIT103871685 | New candidate — activity/pan anchor |
| 2157 | HIT101201113 | New candidate — weak JNK1-bias hypothesis / G1 |
| E1 | literature | Positive control (in hand) |
| CC-90001 | literature | Positive / near-pan control (in hand) |
