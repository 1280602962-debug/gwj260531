# Zero-docking JCIM formulation upgrades (v1)

Not a docking scale-up. Not wet-lab. Not BindingDB docking. Not external validation.

- Audited pairs with both cached maps (aliases dropped): **49**.
- θ = 6.0 directional gate (dual/A/B each n≥10): **17** pairs.
- θ = 6.0 formulation gate (also neither n≥10): **17** pairs.
- Of directional pairs, metal-enzyme-risk: **1**.
- Docked evaluation remains K = 4.

Property-caliper matching uses z-scored MW/cLogP/TPSA/heavy on the frozen scored panels.
AND-filter tables score Dual+A-only+B-only libraries at Dual-percentile cuts of `vina_worst`/`vina_mean`.
Ligand-only ECFP uses the full ChEMBL maps of the four frozen pairs, capped at 250/class, scaffold GroupKFold.

Do not write that docking was evaluated on the census pairs.
Do not write that ligand-only full-map AUROC replaces Table 2.
