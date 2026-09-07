# Five-pair local channels (same article, not Table 2)

Scores were already computed locally. This file is **statistics only**: Table-2-comparable ligand-level non-stratified `summary_min` CIs (B=2000, seed 20260729 + SHA offset), wrong-pocket paired Δ, and property caliper (1 SD). CNN joined by **ligand ID**, not `pair_hint`. Primary CNN readout = `cnn_affinity`. Holdout is unused-pool, **not** external validation.

Production Vina (`vina_20260727`) **reproduced** `tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv` (same estimand, same bootstrap seeds).

## Primary Vina (seed 20260727) vs local channels (`summary_min` [95% CI])

| pair | Vina | five-seed per-seed median (range) | per-ligand median | RTM | CNN affinity | Holdout | indep. GNINA |
|---|---|---|---|---|---|---|---|
| F2/F10 | 0.3448 [0.2109, 0.4773] | 0.3659 (0.3448–0.3851) | 0.3619 [0.2211, 0.4928] | 0.5907 [0.4295, 0.6874] | 0.5756 [0.4182, 0.6887] | 0.3921 [0.2137, 0.5733] | not run |
| JAK1/TYK2 | 0.3649 [0.2306, 0.503] | 0.3955 (0.3649–0.3994) | 0.3715 [0.2375, 0.5121] | 0.6099 [0.449, 0.7096] | 0.499 [0.3463, 0.6471] | 0.475 [0.2816, 0.66] | 0.3172 [0.1826, 0.4632] |
| JAK1/JAK2 | 0.5884 [0.4444, 0.7246] | 0.5879 (0.5742–0.5918) | 0.5791 [0.4325, 0.7136] | 0.6533 [0.5109, 0.7857] | 0.5801 [0.4393, 0.7049] | 0.6194 [0.4196, 0.7493] | not run |
| PPARG/PPARA | 0.6492 [0.5045, 0.7508] | 0.6514 (0.6489–0.691) | 0.6719 [0.5202, 0.7546] | 0.3691 [0.233, 0.4753] | 0.5 [0.3559, 0.6229] | 0.535 [0.3499, 0.7167] | not run |
| PPARA/PPARD | 0.4463 [0.2958, 0.5841] | 0.4536 (0.4463–0.4688) | 0.4492 [0.2991, 0.5867] | 0.5625 [0.4113, 0.6847] | 0.5947 [0.4507, 0.7024] | 0.445 [0.2412, 0.5585] | not run |

## Completeness that matters

- Production Vina / RTM / CNN / per-ligand five-seed median use the same both-end set as the zero-dock stack (F2 107, JAK1/TYK2 109, JAK1/JAK2 110, PPARG 109, PPARA/PPARD 110).
- Later Vina seeds recovered a few production misses (`F2F10_105/106`, `J1TYK2_092`, `PGPA_030`). Per-seed AUROCs are therefore not on an identical ligand set. The per-ligand median requires all five seeds and stays on the production intersection. Do not prefer a later seed because it scored more ligands.
- Independent GNINA JAK1/TYK2: 210/220 jobs ok; four CG0-fail ligands (`J1TYK2_033/034/039/065`) plus `J1TYK2_092` timeout are dropped (not imputed).
- Holdout: `HOF2F10_045` both-end timeout (dual; 19/20/20); `HOPGPA_001` missing 9V8H (A-only; 20/19/20). JAK1/JAK2 remains 20/20/18.
- CNN `pair_hint` concatenates shared receptors; join is by ligand ID. Duplicate ligand+target rows: 0.

## What changes relative to production Vina

- Five-seed per-seed `summary_min` ranges do not cross 0.5 on any pair. Primary seed 20260727 is kept.
- Independent GNINA on JAK1/TYK2: `summary_min` 0.3172 [0.1826, 0.4632], Dual-vs-neither 0.7048 (Vina 0.3649 [0.2306, 0.503], Dual-vs-neither 0.7696). The formulation gap (weak min, stronger Dual-vs-neither) remains.
- RTM on PPARG/PPARA: `summary_min` 0.3691 [0.233, 0.4753] vs Vina 0.6492 [0.5045, 0.7508]; Dual-vs-neither 0.817 vs 0.6853. A high Dual-vs-neither with a collapsed min is the same negative-class failure mode, not a PPAR success.
- RTM vs Vina `summary_min`: lifts F2/F10 (0.5907 vs 0.3448), JAK1/TYK2 (0.6099 vs 0.3649), JAK1/JAK2 (0.6533 vs 0.5884), PPARA/PPARD (0.5625 vs 0.4463); drops PPARG/PPARA (0.3691 vs 0.6492).
- CNN affinity (primary CNN readout): JAK1/TYK2 `0.499` [0.3463, 0.6471]; PPARG/PPARA `0.5` [0.3559, 0.6229]. `cnn_score` is sensitivity only and is not promoted.
- Holdout PPARG/PPARA `summary_min` 0.535 [0.3499, 0.7167] vs main-panel 0.6492 [0.5045, 0.7508]. The main-panel PPARG min does **not** replicate on the unused-pool draw (membership-sensitive). Holdout is not external validation.
- Holdout JAK1/JAK2 stays same-direction (0.6194 [0.4196, 0.7493]; main 0.5884 [0.4444, 0.7246]). Drawn 20/20/18; do not relax the Murcko cap.
- Holdout stays weak (<0.5) on: F2/F10 0.3921, JAK1/TYK2 0.475, PPARA/PPARD 0.445.
- Wrong-pocket Δ CI excludes 0 on: rtm_best9 JAK1/TYK2, gnina_cnn_affinity F2/F10, gnina_cnn_score F2/F10. Those are **same-pose rescores**, not independent docking. Production Vina Δ CIs all include 0. Do not read rescore Δ as pocket-geometry proof. `cnn_score` on F2/F10 is sensitivity only.

## Still true

Internal four-state formulation audit. Not database-external. Does not restock Table 2 or change the title. Does not hard-dock BindingDB. JAK1/JAK2 holdout remains 20/20/18. Do not promote `cnn_score`. Do not average channels into one master AUROC.
