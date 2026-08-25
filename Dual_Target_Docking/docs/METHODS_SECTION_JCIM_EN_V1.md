# Methods (JCIM Articles draft, English)

## 2. Methods

### 2.1 Data and experimental-state definition

Ligand activities used as experimentally derived labels were retrieved from the public ChEMBL Web API activity endpoint. The target-pair supply audit was frozen on 2026-07-23. pChEMBL converts standardized quantitative potency or affinity measurements (e.g. IC50, EC50, Ki, Kd, and Potency) to an approximate −log10 activity scale. Assay types, conditions, and experimental systems are not equivalent; pChEMBL is used here as a unified curation scale across heterogeneous assays.

When several pChEMBL values existed for the same ligand–target pair, the **maximum** was used as the one-to-one representative for primary curation. Activity-aggregation sensitivity was assessed by re-fetching assay-level records and replacing the maximum with the **median** of repeated measurements under the same θ = 6.0 rule, without changing panel membership, docking parameters, or Vina scores. API-refetched max-versus-median estimates are reported as a label-aggregation sensitivity alongside Table 2 (Table S29). Ligands missing a usable pChEMBL value on either target were excluded from analyses requiring paired labels. Salt, solvate, and multicomponent records were split by connected component, retaining the organic fragment with the most heavy atoms.

For each pair A/B, ligands were assigned one of four experimental states: **dual** (strong on both), **A-only** (strong on A, weak on B), **B-only** (strong on B, weak on A), and **neither** (insufficient on both). A-only and B-only ligands are selectivity hard negatives.

The **strict supply-audit rule** (construction gate) was: dual, both pChEMBL ≥ 6.5; A-only, A ≥ 6.5 and B ≤ 5.5; B-only symmetric; neither, both ≤ 5.5. The 5.5–6.5 gray zone was excluded from the strict audit. Metal-dependent systems (e.g. HDACs) were excluded in advance. **Primary AUROCs use one prespecified θ = 6.0 rule:** dual, both ends ≥ θ; A-only, A ≥ θ and B < θ; B-only symmetric; neither, both < θ. Construction rules were frozen from the supply audit before sampling (Table 1). Thresholds were chosen to assemble analyzable panels, not after inspecting docking scores. Relabeling at θ ∈ {5.5, 6.5} and under the strict 6.5/5.5 rule is reported as sensitivity (Table S4). Underpowered cells are flagged in Results.

As a count-only check of the ChEMBL supply gate, the frozen pairs were recounted in BindingDB and PubChem without docking or panel rebuild (Table S12). Endpoints were restricted to IC50/Ki/Kd/EC50; identifiers were BindingDB monomerid and PubChem CID, with no cross-database structure merge. The primary count used equal-relation measurements.

### 2.2 Benchmark construction

DualFourClass-Bench retains four experimental states and two directional primary tasks: dual versus A-only and dual versus B-only. Neither is curated to describe the experimental space but does not enter the primary directional AUROCs.

Candidate pairs were screened with the strict audit in Section 2.1. The frozen evaluation set comprises PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, and EGFR/HER2. EGFR/HER2 is retained as a supply-limited case. Ligands were drawn under frozen class quotas and random seed 20260729. Where structures were available at sampling, a per-class Bemis–Murcko scaffold cap limited series over-representation: at most two molecules per scaffold in PIK3CA/mTOR (PM48) and at most five in EGFR/HER2. AChE/BChE and PIK3CA/PIK3CB used class quotas and a deterministic shuffle. Panels were not redrawn after docking scores were seen.

AChE/BChE and PIK3CA/PIK3CB were sampled under the strict 6.5/5.5 gate (target 28 / 28 / 28 / 16; panel n = 100). EGFR/HER2 (n = 110) and PIK3CA/mTOR PM48 (n = 48; constructed 18 / 14 / 12 / 4) used θ = 6.0 because the strict gate left too few B-only ligands. Cross-pair AUROCs therefore mix target-pair biology with panel-construction differences. Ligand–receptor jobs that failed to yield a score were dropped; analysis counts can fall below construction quotas (Table 1; Table S27). An expanded PIK3CA/mTOR panel (PM110) keeps all 48 PM48 ligands as a nested size check.

**Table 1.** DualFourClass-Bench composition and docking settings. Construction labels record the supply/panel-building rule for each pair. All primary AUROCs in Tables 2–3 use unified θ = 6.0 experimental-state labels.

| Pair | Construction labels | PDB (A / B) | Resolution (Å) | Panel n | Analysis n (dual / A_only / B_only) | Vina exhaustiveness |
|------|---------------------|-------------|----------------|-------:|------------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | strict 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.3 Receptor preparation and docking protocol

Receptors were PDB entries with a small-molecule cognate ligand: PIK3CA/mTOR, 4L23 / 4JT6 (X6K / PI-103); AChE/BChE, 4EY7 / 4BDS (E20 / THA); PIK3CA/PIK3CB, 4L23 / 2WXF (X6K / 039); EGFR/HER2, 3POZ / 3RCD (03P / TAK-285). The site was defined from the cognate ligand. An axis-aligned bounding box on cognate heavy atoms was expanded by 5 Å on each axis; any edge shorter than 20 Å was set to at least 20 Å (Table S2). Water and the cognate ligand were removed and Meeko wrote PDBQT. PIK3CA, mTOR, EGFR, and HER2 used hydrogen-containing protein coordinates already in the frozen directories (`mk_prepare_receptor.py --read_pdb`). AChE, BChE, and PIK3CB were extracted from deposited ATOM/TER records and prepared with `mk_prepare_receptor` (default alternate location A). PDBFixer and Reduce were not used for independent rebuilding or tautomer enumeration. Docking treated noncovalent small-molecule sites only.

Before production docking, each frozen receptor was redocked with its cognate ligand. Nine poses were generated and heavy-atom RMSD to the crystal ligand was computed in the docking frame. The prespecified pass criterion was \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å, testing pose-generation capability among retained poses rather than requiring the top-ranked Vina pose to be near-native. If default exhaustiveness failed the gate, search effort was raised to a prespecified fallback. Production docking therefore used exhaustiveness 16 for PIK3CA/mTOR and 8 for the other main panels (Table S3).

Ligands started from frozen ChEMBL SMILES: desalt to the largest organic fragment, add explicit hydrogens in RDKit, embed with ETKDGv3 (seed 20260727), locally optimize with MMFF (at most 200 steps), and convert with default Meeko to PDBQT. Protonation states, tautomers, and conformational ensembles were not systematically enumerated. Docking used AutoDock Vina 1.2.7 with the default `vina` scoring function, nine retained poses, `energy_range = 3` kcal mol\(^{-1}\), and random seed 20260727 (Table S1). To test scoring-function dependence, the same Vina poses were rescored with RTMScore (`rtmscore_model1`, best of nine) and GNINA 1.3.2 CNN (`--cnn_scoring rescore --minimize`, best of nine after Open Babel SDF conversion). Vina’s primary readout is the mode-1 energy; RTM and GNINA are best-of-9 rescores. The primary endpoint remains Vina.

### 2.4 Primary endpoint and statistics

Throughout this Article, “dual-target recognition” names this computational discrimination task. Two binary AUROCs were computed per pair. Dual versus A-only used the pocket B score, \( \mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B) \). Dual versus B-only used the pocket A score, \( \mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A) \). Dual is always the positive class. Vina reports \(E_{\mathrm{Vina}}\) (kcal mol\(^{-1}\); more negative is more favorable); \(S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}\).

The worst-direction discrimination summary is \( \mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}) \). It is a conservative summary of two AUROCs. Arithmetic, geometric, and harmonic means are aggregation sensitivities (Table S26). The single primary endpoint is pocket-matched Vina `summary_min` under unified θ = 6.0 (Table 2; PIK3CA/mTOR uses PM48). A prespecified RDKit panel (heavy-atom count, molecular weight, cLogP, and TPSA) was evaluated with the same directional workflow; the highest AUROC among them is a best single-descriptor reference (Tables 2, S28, S19). Dual versus neither (experimental inactives; `vina_mean`) and Dual versus all non-duals were computed as auxiliary formulation contrasts on the same frozen scores (Table 3; Table S22). PIK3CA/mTOR neither n = 4 is flagged underpowered.

AUROC and summary_min uncertainty used ligand-level bootstrap: ligands were resampled with replacement, preserving class structure (\(B = 2000\), seed 20260729, percentile 95% CI). Paired contrasts used the same resample (Tables S17, S19). Intervals are descriptive.

### 2.5 Confounder, holdout, and receptor-sensitivity analyses

Scores for targets A and B were swapped as a falsification control, leaving ligands, receptors, and all other settings unchanged. Directional AUROCs were also recomputed after ligand-efficiency normalization (\(S_{\mathrm{dock}}/N_{\mathrm{heavy}}\)) and on potency- or size-constrained subsets (\(|\Delta\mathrm{pChEMBL}| \leq 0.5\); \(|\Delta N_{\mathrm{heavy}}| \leq 2\)). Logistic models compared docking alone with docking plus heavy-atom count and TPSA. Morgan/ECFP4 fingerprints (radius 2, 2048 bits) with logistic regression provided a ligand-only chemical baseline under Bemis–Murcko scaffold `GroupKFold` (Tables S5, S20, S23, S24). Nearest-neighbor Tanimoto subsets were diagnostic. A coarse contact count and whole-chain sequence identity were exploratory only (Tables S7, S11).

To test dependence on exact panel membership, an unused-pool holdout was drawn after excluding all ChEMBL entries used in the main panels and PM110. Ligands still come from the same ChEMBL harvest, target pairs, and label rules. Holdout was built for PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB (20 dual / 20 A-only / 20 B-only; `HOLDOUT_SEED = 20260731`); EGFR/HER2 was not eligible. Receptor, box, ligand preparation, exhaustiveness, scoring, and statistics matched the main benchmark (Tables S8, S13).

Receptor-structure sensitivity used alternate crystals that, before scores were seen, (i) matched the true target protein, (ii) contained a small-molecule cognate in the ATP or target site, (iii) had acceptable resolution, and (iv) passed the same cognate redocking QC. Structures docked were PIK3CA 4JPS and 5DXT and mTOR 4JSX. Replacement was one pocket at a time: on PIK3CA/mTOR (PM48), 4JPS/5DXT replaced pocket A while pocket B kept frozen 4JT6 scores, and 4JSX replaced pocket B while pocket A kept frozen 4L23 scores (exhaustiveness 16). On PIK3CA/PIK3CB, the same 4JPS and 5DXT receptors replaced pocket A while pocket B kept frozen 2WXF scores (exhaustiveness 8). Rigid Cα superposition was an exploratory geometric control (Table S10).

Analyses ran under Python 3 with RDKit 2026.3.1, meeko 0.7.1, AutoDock Vina 1.2.7, GNINA 1.3.2, and RTMScore. Panels, scores, scripts, and parameter tables are available in the public repository (Data and Software Availability).
