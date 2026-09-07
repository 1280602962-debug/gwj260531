# Five-pair results merged (zero-dock + dump-gated + local channels)

Same article, 8-row destination after withdrawing PIK3CA/PIK3CB (`PROJECT_IDENTITY_LOCK_V1.md`). Five pairs added after the ChEMBL 37 census; not a 2026-07-23 freeze. **Does not restock Table 2 or retitle.**

This note merges three already-run layers. It does not re-dock.

## 1. Primary numbers (production Vina, θ = 6.0)

Source: `FIVE_PAIR_STACK_V1.md`. Estimand = ligand-level non-stratified `summary_min` of Dual-vs-A (pocket B) and Dual-vs-B (pocket A). Reproduced by this channel run.

| pair | n D/A/B | summary_min [95% CI] | Dual-vs-neither | best descriptor |
|---|---:|---|---:|---|
| F2/F10 | 31/32/32 | 0.3448 [0.2109, 0.4773] | 0.5188 | clogp 0.5151 |
| JAK1/TYK2 | 31/32/32 | 0.3649 [0.2306, 0.503] | 0.7696 | clogp 0.5796 |
| JAK1/JAK2 | 32/32/32 | 0.5884 [0.4444, 0.7246] | 0.7299 | heavy 0.5781 |
| PPARG/PPARA | 32/31/32 | 0.6492 [0.5045, 0.7508] | 0.6853 | tpsa 0.6274 |
| PPARA/PPARD | 32/32/32 | 0.4463 [0.2958, 0.5841] | 0.5647 | clogp 0.5635 |

Only PPARG/PPARA has a `summary_min` CI entirely above 0.5. JAK1/JAK2 point > 0.5 but CI includes 0.5. The other three are anti-directional or straddle 0.5. Dual-vs-neither is **not** the primary claim (it can look good while a weak arm fails).

## 2. Dump-gated and leftover holdout IDs

Source: `FIVE_PAIR_DUMP_GATED_V1.md`. sqlite SHA matches the freeze.

- Dump max vs panel pChEMBL: **0/545** mismatches; **0** missing dump ends.
- Leftover vs frozen panel summary: exact MATCH. Holdout IDs frozen (seed 20260731, Murcko cap 3). JAK1/JAK2 **20/20/18**.
- 1 max→median class flip: PPARA/PPARD CHEMBL121. Do not switch labels.
- 2018 time split (earliest document year): reportable only on JAK1/TYK2 (min 0.3681) and JAK1/JAK2 (min 0.6531). Shared JAK1; **not** external validation. 2020 fails n≥10 everywhere.

## 3. BindingDB / PubChem count-only

Source: `FIVE_PAIR_CROSSDB_V1.md`. `equal_only` does **not** flip the ≥50 hard-neg gate. Do not hard-dock BindingDB.

## 4. Local score channels (this run)

Source: `FIVE_PAIR_LOCAL_CHANNELS_V1.md` and `tables/five_pair_local_channels_v1/`.

## Primary Vina (seed 20260727) vs local channels (`summary_min` [95% CI])

| pair | Vina | five-seed per-seed median (range) | per-ligand median | RTM | CNN affinity | Holdout | indep. GNINA |
|---|---|---|---|---|---|---|---|
| F2/F10 | 0.3448 [0.2109, 0.4773] | 0.3659 (0.3448–0.3851) | 0.3619 [0.2211, 0.4928] | 0.5907 [0.4295, 0.6874] | 0.5756 [0.4182, 0.6887] | 0.3921 [0.2137, 0.5733] | not run |
| JAK1/TYK2 | 0.3649 [0.2306, 0.503] | 0.3955 (0.3649–0.3994) | 0.3715 [0.2375, 0.5121] | 0.6099 [0.449, 0.7096] | 0.499 [0.3463, 0.6471] | 0.475 [0.2816, 0.66] | 0.3172 [0.1826, 0.4632] |
| JAK1/JAK2 | 0.5884 [0.4444, 0.7246] | 0.5879 (0.5742–0.5918) | 0.5791 [0.4325, 0.7136] | 0.6533 [0.5109, 0.7857] | 0.5801 [0.4393, 0.7049] | 0.6194 [0.4196, 0.7493] | not run |
| PPARG/PPARA | 0.6492 [0.5045, 0.7508] | 0.6514 (0.6489–0.691) | 0.6719 [0.5202, 0.7546] | 0.3691 [0.233, 0.4753] | 0.5 [0.3559, 0.6229] | 0.535 [0.3499, 0.7167] | not run |
| PPARA/PPARD | 0.4463 [0.2958, 0.5841] | 0.4536 (0.4463–0.4688) | 0.4492 [0.2991, 0.5867] | 0.5625 [0.4113, 0.6847] | 0.5947 [0.4507, 0.7024] | 0.445 [0.2412, 0.5585] | not run |

## 5. What the channels do to the story

- Five-seed per-seed `summary_min` ranges do not cross 0.5 on any pair. Primary seed 20260727 is kept.
- Independent GNINA on JAK1/TYK2: `summary_min` 0.3172 [0.1826, 0.4632], Dual-vs-neither 0.7048 (Vina 0.3649 [0.2306, 0.503], Dual-vs-neither 0.7696). The formulation gap (weak min, stronger Dual-vs-neither) remains.
- RTM on PPARG/PPARA: `summary_min` 0.3691 [0.233, 0.4753] vs Vina 0.6492 [0.5045, 0.7508]; Dual-vs-neither 0.817 vs 0.6853. A high Dual-vs-neither with a collapsed min is the same negative-class failure mode, not a PPAR success.
- RTM vs Vina `summary_min`: lifts F2/F10 (0.5907 vs 0.3448), JAK1/TYK2 (0.6099 vs 0.3649), JAK1/JAK2 (0.6533 vs 0.5884), PPARA/PPARD (0.5625 vs 0.4463); drops PPARG/PPARA (0.3691 vs 0.6492).
- CNN affinity (primary CNN readout): JAK1/TYK2 `0.499` [0.3463, 0.6471]; PPARG/PPARA `0.5` [0.3559, 0.6229]. `cnn_score` is sensitivity only and is not promoted.
- Holdout PPARG/PPARA `summary_min` 0.535 [0.3499, 0.7167] vs main-panel 0.6492 [0.5045, 0.7508]. The main-panel PPARG min does **not** replicate on the unused-pool draw (membership-sensitive). Holdout is not external validation.
- Holdout JAK1/JAK2 stays same-direction (0.6194 [0.4196, 0.7493]; main 0.5884 [0.4444, 0.7246]). Drawn 20/20/18; do not relax the Murcko cap.
- Holdout stays weak (<0.5) on: F2/F10 0.3921, JAK1/TYK2 0.475, PPARA/PPARD 0.445.
- Wrong-pocket Δ CI excludes 0 on: rtm_best9 JAK1/TYK2, gnina_cnn_affinity F2/F10, gnina_cnn_score F2/F10. Those are **same-pose rescores**, not independent docking. Production Vina Δ CIs all include 0. Do not read rescore Δ as pocket-geometry proof. `cnn_score` on F2/F10 is sensitivity only.

## 6. Honest combined reading

- The five new pairs are **not** a second PPAR-success paper. PPARG/PPARA is the only main-panel Vina min whose CI clears 0.5, and that number is membership-sensitive (holdout) and formulation-sensitive (RTM).
- JAK1/JAK2 is the most internally consistent Vina pair (main + holdout same direction; five-seed stable), but the CI still includes 0.5 and it is not an independent confirmation of JAK1/TYK2.
- F2/F10, JAK1/TYK2, and PPARA/PPARD stay weak-armed under Vina. RTM/CNN can move a point estimate without turning the stack into a pocket-matched dual classifier.
- Independent GNINA on the gap pair does not close the EGFR-like gap.
- Unused-pool holdout and 2018 year-split are **internal** confirmations. The frozen BindingDB external contract still has 0 pairs. Do not relax it.
- Stack completeness is now a **precondition** for a later 8-row Table 2 draft. That retitle/restock is a later manuscript step, not this commit.

## 7. Still forbidden

Restock Table 2 today; retitle; average eight pairs; hard-dock BindingDB; promote `cnn_score`; re-draw holdout IDs; relax Murcko cap; run independent GNINA on the other four new pairs; LigPrep / seed 42; CTSK ordinary Vina; PIK3CB 2Y3A swap.
