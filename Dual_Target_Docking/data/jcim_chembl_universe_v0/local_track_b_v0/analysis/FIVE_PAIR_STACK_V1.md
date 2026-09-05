# Five-pair zero-dock stack (same article, not Track B)

Destination: 8-row main table after withdrawing PIK3CA/PIK3CB (`PROJECT_IDENTITY_LOCK_V1.md`). These five pairs were added after the ChEMBL 37 census; they were not frozen on 2026-07-23. This run does **not** restock Table 2 or change the title.

Bootstrap for `summary_min` is the Table 2 estimand: ligand-level non-stratified resample of the dual+A-only+B-only pool, B=2000, seed 20260729 + SHA-256 pair offset. Existing `track_b_directional_auroc_v1.csv` CIs remain class-preserving and are not Table 2 intervals.

## Table-2-comparable Vina (θ = 6.0, both-end scores)

| pair | n D/A/B | D/A (B) | D/B (A) | summary_min [95% CI] | Dual-vs-neither | best descriptor |
|---|---:|---:|---:|---|---:|---|
| F2/F10 | 31/32/32 | 0.4133 | 0.3448 | 0.3448 [0.2109, 0.4773] | 0.5188 | clogp 0.5151 |
| JAK1/TYK2 | 31/32/32 | 0.5751 | 0.3649 | 0.3649 [0.2306, 0.503] | 0.7696 | clogp 0.5796 |
| JAK1/JAK2 | 32/32/32 | 0.5884 | 0.7275 | 0.5884 [0.4444, 0.7246] | 0.7299 | heavy 0.5781 |
| PPARG/PPARA | 32/31/32 | 0.6492 | 0.7061 | 0.6492 [0.5045, 0.7508] | 0.6853 | tpsa 0.6274 |
| PPARA/PPARD | 32/32/32 | 0.6465 | 0.4463 | 0.4463 [0.2958, 0.5841] | 0.5647 | clogp 0.5635 |

## Holdout leftover (counts only; IDs need sqlite)

| pair | leftover D/A/B | 20/20/20 | thin |
|---|---:|---:|---:|
| F2/F10 | 312/76/245 | 1 | 0 |
| JAK1/TYK2 | 1874/59/80 | 1 | 0 |
| JAK1/JAK2 | 5953/76/21 | 1 | 1 |
| PPARG/PPARA | 408/50/59 | 1 | 0 |
| PPARA/PPARD | 187/50/68 | 1 | 0 |

## Still blocked or local

- **max_vs_median_pchembl** — `blocked_no_sqlite`: panel CSV stores one pChEMBL per end; repeat-record graph needs ChEMBL 37 sqlite
- **document_year_split** — `blocked_no_sqlite`: no document.year on panels; report counts only after dump join; AUROC only if dual/A/B each n≥10 after 2018
- **document_cluster_bootstrap** — `blocked_no_sqlite`: no document_id; scaffold-cluster bootstrap was computed instead
- **document_blocked_cv** — `blocked_no_sqlite`: same missing document_id; ECFP4 used Bemis–Murcko GroupKFold
- **bindingdb_pubchem_count_only** — `blocked_no_cache`: jcim_supply_crossdb_v0 caches only the original K=4 UniProts; five new pairs need a new count-only fetch (no Docker)
- **holdout_panel_ids** — `blocked_no_sqlite`: leftover counts are known; 20/20/20 member lists need the dump + HOLDOUT_SEED=20260731
- **five_seed_vina** — `local_recompute`: user will submit locally; seeds 20260727 + 20260811–20260814; see LOCAL_RECOMPUTE_PACK_V1.md
- **rtm_best_of_9** — `local_recompute`: poses gitignored; regenerate 9 modes then rtmscore_model1
- **gnina_cnn_rescore** — `local_recompute`: same poses; --cnn_scoring rescore --minimize
- **independent_gnina_search** — `local_rule_subset`: JAK1/TYK2 only (EGFR-like formulation gap)

## ECFP4 scaffold GroupKFold (docking increment)

| pair | contrast | ECFP4 | ECFP4+docking | Δ | docking rank |
|---|---|---:|---:|---:|---:|
| F2/F10 | D_vs_A | 0.9365 | 0.9294 | -0.0071 | 0.4133 |
| F2/F10 | D_vs_B | 0.6653 | 0.6865 | 0.0212 | 0.3448 |
| JAK1/TYK2 | D_vs_A | 0.8458 | 0.8397 | -0.0061 | 0.5751 |
| JAK1/TYK2 | D_vs_B | 0.9022 | 0.9032 | 0.001 | 0.3649 |
| JAK1/JAK2 | D_vs_A | 0.9146 | 0.916 | 0.0014 | 0.5884 |
| JAK1/JAK2 | D_vs_B | 0.9678 | 0.9688 | 0.001 | 0.7275 |
| PPARG/PPARA | D_vs_A | 0.8327 | 0.8196 | -0.0131 | 0.6492 |
| PPARG/PPARA | D_vs_B | 0.668 | 0.6758 | 0.0078 | 0.7061 |
| PPARA/PPARD | D_vs_A | 0.9316 | 0.9277 | -0.0039 | 0.6465 |
| PPARA/PPARD | D_vs_B | 0.8584 | 0.835 | -0.0234 | 0.4463 |

Independent GNINA search is **not** in this run. JAK1/TYK2 is the only new pair that qualifies under the original formulation-gap rule.
