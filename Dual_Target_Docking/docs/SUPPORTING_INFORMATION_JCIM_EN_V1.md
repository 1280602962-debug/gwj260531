# Supporting Information (English, compressed submission draft)

**Article:** A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition  
**Hierarchy:** Primary / Prespecified sensitivity / Post-hoc exploratory — see `docs/ANALYSIS_HIERARCHY_V1.md`.  
**Archive policy:** Run logs, per-ligand long tables, historical receptor candidates, and full filter streams remain in the GitHub repository; this SI retains tables that answer a distinct reviewer question.

---

## S1. Data sources and label rules

| SI table | Content | Role |
|---------|---------|------|
| S1 | Software versions, seeds, Vina/GNINA/RTM settings | Reproducibility |
| S4 | θ sensitivity (5.5 / 6.0 / 6.5 / strict) | Prespecified sensitivity |
| S12 | BindingDB/PubChem count-only supply check | Prespecified sensitivity |
| S29, S35, S36 | Max vs median; measurement frequency; high-confidence view | Sensitivity / post-hoc label audit |
| S37 | Complete-case coverage and document concentration | Primary limitation |

## S2. Frozen panel composition

| SI table | Content | Role |
|---------|---------|------|
| Table 1 (main) | K=4 composition and exhaustiveness | Primary |
| S27 | Docking success/fail by class and properties | Prespecified sensitivity |
| S38 | Class-wise chemistry summary | Post-hoc exploratory |
| S44 | θ=6.0 pair census (not docked) | Post-hoc exploratory |

## S3. Receptors, boxes, and cognate QC

| SI table | Content | Role |
|---------|---------|------|
| S2 | Box definitions | Primary protocol |
| S3 | Cognate ranked RMSD re-audit (topology-checked where available) | Prespecified / reconstructed QC |
| S9–S10, S30 | Alternate receptors and pocket RMSD | Prespecified sensitivity |
| S33 | Geometric occupancy snapshot | Post-hoc exploratory |

## S4. Primary AUROCs and intervals

| SI table | Content | Role |
|---------|---------|------|
| Tables 2–3 (main) | Directional AUROC; Dual-vs-neither contrast | Primary / formulation contrast |
| S6, S17, S19 | Wrong-pocket, paired deltas, descriptor contrasts | Sensitivity |
| S22, S25, S26, S28, S34 | Formulation, Top-10, aggregation, descriptors, fixed-score neither | Mix of primary support and sensitivity |
| S31 | Detectable-effect simulation | Prespecified sensitivity |
| S32 | Independent GNINA pose generation | Prespecified sensitivity |
| S39–S40 | Document-blocked CV | Prespecified sensitivity |
| `document_cluster_bootstrap_v1.csv` | Document-cluster CIs | Prespecified sensitivity |
| `scaffold_cluster_bootstrap_v1.csv` | Scaffold-cluster CIs | Prespecified sensitivity |

**EGFR/HER2 weak arm (Dual vs B-only = 0.430):** ligand CI [0.282, 0.578]; document-cluster [0.321, 0.617]; scaffold-cluster [0.278, 0.595]. All span 0.5.

## S5. Chemistry and source controls

| SI table | Content | Role |
|---------|---------|------|
| S5, S20, S23–S24 | Matched subsets; scaffold vs random; chemotype proximity; ECFP+docking | Sensitivity / exploratory |
| S45–S47 | Property calipers; AND filter; full-map ECFP4 | Post-hoc exploratory |

## S6. Docking-failure sensitivity

See Table S27 and Methods missing-data bounds. Failures concentrate among large/flexible ligands; pair-level conclusions were not reversed under rank-extreme bounds.

## S7. Receptor sensitivity

Tables S9, S30 and Figure 5: same PIK3CA alternate crystals lower PIK3CA/mTOR `summary_min` (0.692 → 0.486/0.505) and raise PIK3CA/PIK3CB (0.500 → 0.691/0.685).

## S8. Document / time split

| SI table | Content | Role |
|---------|---------|------|
| S39–S41 | Document-blocked CV; time split | Prespecified sensitivity |
| S42 | Assay-context metadata review | Post-hoc; no label flips |

PIK3CA/mTOR Dual vs B-only under document blocking: **not stably estimable**.

## S9. BindingDB supply freeze

| SI table | Content | Role |
|---------|---------|------|
| S43 | Historical REST independence counts | Historical supply |
| S48–S49 | Native archive flow and gate summary | Prespecified external audit (**zero pairs**; not docked) |

## S10. Reproducibility

| Item | Location |
|------|----------|
| Evaluation contract | `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json` |
| Checksum manifest | `data/jcim_novelty_v0/tables/REVISION_CHECKSUM_MANIFEST_v1.csv` |
| Validator | `data/jcim_novelty_v0/scripts/validate_revision_v1.py` |
| Manuscript assembly | `docs/assemble_manuscript_en.py` |
| Analysis hierarchy | `docs/ANALYSIS_HIERARCHY_V1.md` |
| MCL1 formal demotion | `data/mcl1_bclxl_panel_v0/analysis/MCL1_BCLXL_FORMAL_DEMOTION_V1.md` |

**MCL1/Bcl-xL:** exploratory archive only after Option B demotion; not a main SI performance table.

## Figure list (SI)

| Figure | Content |
|--------|---------|
| S1 | Threshold / scorer / panel-size sensitivities |
| S3 | Paired delta / bootstrap panels |
| S4 | Pocket-matched forest plot |
| S5 | Unused-pool holdout |
| S6 | Detectable-effect heatmap |
| S7 | Census / AND / full-map (exploratory) |
| S8 | BindingDB-native slice (gate fail) |

## Note on Chinese working draft

`SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` remains a historical working archive. This English compressed SI is the submission-facing outline; detailed numeric cells continue to be drawn from the cited CSV paths until the SI PDF is typeset.
