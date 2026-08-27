# DualFourClass analysis hierarchy (frozen before further manuscript claims)

Updated: 2026-08-27  
Scope: JCIM Article on DualFourClass-Bench. Every Results sentence, table title, and figure caption must use one of the three labels below. Post-hoc analyses must not be phrased as confirmatory primary evidence.

## Primary (confirmatory for the article’s core claim)

1. **Four frozen pairs:** PIK3CA/mTOR (PM48), AChE/BChE, PIK3CA/PIK3CB, EGFR/HER2.
2. **Primary label:** θ = 6.0 four-state assignment (dual / A-only / B-only / neither).
3. **Primary scores:** AutoDock Vina mode-1 pocket scores on the frozen receptors/boxes.
4. **Primary estimands:** two pocket-matched directional AUROCs  
   - Dual vs A-only in pocket B  
   - Dual vs B-only in pocket A  
   Both arms are always reported.
5. **`summary_min`:** prespecified *descriptive* worst-arm summary only; not a calibrated score and not the sole inferential quantity.
6. **Primary uncertainty (minimum):** ligand-level bootstrap CIs on the two directional arms (and on `summary_min` where reported), with underpowered cells flagged.

Main-text Tables 2–3 and Figures 3–4 are primary unless explicitly labeled otherwise.

## Prespecified sensitivity (planned; not claim-expanding)

Report as sensitivity; do not elevate to “robustness proof.”

- θ ∈ {5.5, 6.5} and strict 6.5/5.5 relabeling on the same scores  
- Docking-failure arm-available AUROC and rank-extreme missing-data bounds  
- Receptor-realization swaps (4JPS/5DXT/4JSX) chosen before AUROC  
- Scaffold-grouped and document-blocked grouping of the same primary scores  
- Ligand / document-cluster / scaffold-cluster bootstrap intervals on primary arms  
- Dual vs neither on the same frozen scores as a *formulation contrast* (changes negative class; not a paired test of one estimand)  
- Independent GNINA pose-generation on EGFR/HER2 and PIK3CA/mTOR (formulation persistence check, not engine bake-off)  
- Unused-pool holdout (internal panel-membership only)  
- Pre-frozen 2018 literature-year split (supply failure; not external validation)  
- BindingDB-native archive rebuild under the pre-frozen external gate (zero pairs; not docked)

## Post-hoc exploratory (diagnostic; never confirmatory)

Must be labeled exploratory in Methods, Results, and captions.

- AND-like dual pocket filters and operating-point grids  
- θ = 6.0 17-pair / 49-pair label census (supply diagnosis; no extra docking)  
- Full-map ligand-only ECFP4 models beyond the docking panels  
- MCL1/Bcl-xL panel docking and AUROCs after **formal demotion** (Option B; not pose-gold validated)  
- Geometric PLIF / occupancy snapshots and contact-count baselines  
- Property caliper matching, measurement-frequency audits, and related chemistry dumps that were not the primary endpoint  
- Aggregation alternatives beyond the designated `summary_min` (arithmetic / geometric / harmonic) when used only to show rank stability

## External validation

**Not achieved.** BindingDB supply under the pre-frozen gate failed. Do not claim database-external or prospective utility. Acceptable manuscript stance: data-constrained formulation audit with an explicit failed external-supply audit.

## One-sentence claim ceiling

Dual-versus-neither docking evaluation can overlook single-target selectives; on the present four pairs that formulation effect is pair-, chemistry-, and receptor-dependent and does not establish target-general docking performance.
