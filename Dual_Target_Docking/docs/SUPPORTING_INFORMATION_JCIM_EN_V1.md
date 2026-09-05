# Supporting Information (English, compressed submission draft)

**Article:** A Three-Pair Formulation Audit of Docking-Based Dual-Target Recognition  
**Hierarchy:** Primary / Prespecified sensitivity / Post-hoc exploratory — see `docs/ANALYSIS_HIERARCHY_V1.md`.  
**Archive policy:** Run logs, per-ligand long tables, historical receptor candidates, and full filter streams remain in the GitHub repository; this SI retains tables that answer a distinct reviewer question.

## S0. Receptor-identity failure: PIK3CA/PIK3CB (withdrawn from the main text)

A fourth candidate pair, PIK3CA/PIK3CB, passed the same ChEMBL supply screen as the three primary pairs and was docked (4L23/2WXF, exhaustiveness 8, n_panel = 100, n_scored 28/27/28) before a post-hoc receptor-identity audit (`data/jcim_chembl_universe_v0/analysis/RECEPTOR_IDENTITY_AUDIT_V1.md`) found that PDB entry 2WXF resolves to **murine PIK3CD** (SIFTS UniProt O35904, gene *Pik3cd*, *Mus musculus*, entity coverage 0.998; RCSB title: "the murine class IA PI 3-kinase p110delta in complex with PIK-39"), not human PIK3CB (P42338), which has no deposited PDB structure. The A pocket (PIK3CA 4L23) is correct; only the B pocket is wrong. Cognate redocking QC could not detect this (best-of-9 RMSD 0.405 Å) because the cognate ligand 039 (PIK-39) belongs to the same mouse protein, so the QC pass is a true statement about pose recovery on the docked receptor and not evidence of correct protein identity.

Numbers computed on this pair before withdrawal, retained here only as a documented failure case and **not** as primary or sensitivity evidence for dual-target recognition:

- Primary Vina, unified θ = 6.0: pocket-matched `summary_min` 0.500 [0.350, 0.650] (dual vs A_only 0.691; dual vs B_only 0.500); Dual vs neither 0.559 [0.373, 0.746] (n_neither = 16); Dual vs all non-duals 0.556 [0.437, 0.672].
- Receptor-structure sensitivity (same PIK3CA 4JPS/5DXT crystals used for PIK3CA/mTOR in Table S30, B pocket frozen at 2WXF): `summary_min` rose from 0.500 to 0.691 and 0.685 — the opposite direction from PIK3CA/mTOR's drop under the identical PIK3CA substitution. Because the B end is a nonhuman off-target receptor, this contrast is reported only as a caution that receptor substitution results are not interpretable once receptor identity is wrong, not as evidence that receptor effects are pair-direction-dependent in general.
- Five-seed Vina (Table S54 protocol): `summary_min` median 0.478 (range 0.468–0.502) across five frozen seeds.
- GNINA best-of-9 rescore of the same Vina poses: 0.554 (mode01) / 0.533 (best9), both inside the Vina bootstrap CI.
- Docking-failure rank-extreme check: one failed A-only ligand, using the available pocket score left `summary_min` at 0.500.
- Receptor-swap holdout: this pair was included in the unused-pool holdout resample before withdrawal; its holdout numbers are likewise not carried into the primary Figure 5 comparison.

**Lesson for the protocol:** identity (UniProt accession + source organism) must be verified against the entity-level SIFTS mapping *before* cognate-ligand redocking QC is run, not only after. Cognate RMSD tests pose recovery on the receptor as given; it cannot substitute for a protein-identity check. This is now folded into the Layer-2 site-verification checklist (`SITE_VERIFICATION_CHECKLIST_V1.md`) used for all pairs added after this finding.

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
| Table 1 (main) | K=3 composition and exhaustiveness (PIK3CA/PIK3CB withdrawn; S0) | Primary |
| S27 | Docking success/fail by class and properties | Prespecified sensitivity |
| S38 | Class-wise chemistry summary | Post-hoc exploratory |
| S44 | θ=6.0 pair census (not docked) | Post-hoc exploratory |

## S3. Receptors, boxes, and cognate QC

| SI table | Content | Role |
|---------|---------|------|
| S2 | Box definitions | Primary protocol |
| S3 | Cognate ranked RMSD re-audit (topology-checked where available) | Prespecified / reconstructed QC |
| S9–S10, S30 | Alternate receptors and pocket RMSD (PIK3CA/mTOR only; the parallel PIK3CA/PIK3CB arm is S0, not sensitivity evidence) | Prespecified sensitivity |
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
| S54 | Five-seed Vina sensitivity (`AUC(vina_mean)` Dual vs neither) | Prespecified sensitivity |
| `document_cluster_bootstrap_v1.csv` | Document-cluster CIs | Prespecified sensitivity |
| `scaffold_cluster_bootstrap_v1.csv` | Scaffold-cluster CIs | Prespecified sensitivity |

**EGFR/HER2 weak arm (Dual vs B-only = 0.430):** ligand CI [0.282, 0.578]; document-cluster [0.321, 0.617]; scaffold-cluster [0.278, 0.595]. All span 0.5.

Table S54 reports four additional frozen Vina seeds (20260811–20260814) with production seed 20260727. Ligands, receptors, boxes, exhaustiveness, retained modes, and analysis rules were unchanged. Dual versus neither uses per-ligand `vina_mean`, matching Table 3; the primary seed recovered 0.756 / 0.649 / 0.559 / 0.514. Directional `summary_min` medians (ranges) were EGFR/HER2 0.373 (0.321–0.430), AChE/BChE 0.599 (0.553–0.606), and PIK3CA/mTOR 0.704 (0.676–0.726). The EGFR/HER2 formulation gap was positive at all five seeds. (PIK3CA/PIK3CB was also seeded five ways before its withdrawal, median 0.478 (0.468–0.502); see S0 — not primary.) PIK3CA/mTOR Dual versus neither remains underpowered (n_neither = 4) and is not interpreted as a reversal. Source: `data/jcim_multiseed_v0/tables/multiseed_auroc_aggregate_v2.csv`. The v1 Dual-versus-neither column used `mean(AUC_A, AUC_B)` and is not cited.

## S5. Chemistry and source controls

| SI table | Content | Role |
|---------|---------|------|
| S5, S20, S23–S24 | Matched subsets; scaffold vs random; chemotype proximity; ECFP+docking | Sensitivity / exploratory |
| S45–S47 | Property calipers; AND filter; full-map ECFP4 | Post-hoc exploratory |

## S6. Docking-failure sensitivity

See Table S27 and Methods missing-data bounds. Failures concentrate among large/flexible ligands; pair-level conclusions were not reversed under rank-extreme bounds.

## S7. Receptor sensitivity

Tables S9, S30 and Figure 4B: same PIK3CA alternate crystals lower PIK3CA/mTOR `summary_min` (0.692 → 0.486/0.505). A parallel swap on the withdrawn PIK3CA/PIK3CB pair raised its point estimate (0.500 → 0.691/0.685), but that pair's B pocket is a receptor-identity failure (S0), so the contrast is not used as evidence of pair-direction-dependence.

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
| Statistical lock | `docs/STATISTICAL_LOCK_V1.md` |
| MCL1 formal demotion | `data/mcl1_bclxl_panel_v0/analysis/MCL1_BCLXL_FORMAL_DEMOTION_V1.md` |

**MCL1/Bcl-xL (formally demoted; exploratory repository archive).** Prespecified LC6 topology-aware pose-gold was not established for 3WIY/3WIZ, so this pair is not a fifth main pair and is not used as confirmatory evidence. A frozen 24/24/24/24 panel was docked as an applicability archive: 93/96 ligands scored on both pockets. Vina Dual versus neither was 0.628 [0.462, 0.786], Dual versus A-only 0.793 [0.655, 0.915], Dual versus B-only 0.609 [0.439, 0.776], worst-arm AUROC 0.609 (Tables S50–S51, S53).[17] Comparator notes versus Zhou 2013, DUD-E, LIT-PCBA, CASF-2016, and DOCKSTRING are in Table S52.

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
