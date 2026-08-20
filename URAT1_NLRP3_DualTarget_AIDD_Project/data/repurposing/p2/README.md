# Archived gnina P2 funnel tables

Slim copies from `docking_export_20260820/` (poses stay in the export).  
Engine: **gnina P2** (CNNaffinity, exhaustiveness 32, `num_modes=1`, CNN rescore).  
Leftover `glide_score_xp` aliases were stripped; they were identical to `dock_score`.

| File | n | Role |
|------|---|------|
| `docking_9dkb_gnina.csv` / `docking_7alv_gnina.csv` | 1583 rows; 1582 docked | Per-target P2 scores |
| `pareto_merged_scores.csv` | **1580** dual-success complete case | Percentiles and Pareto audit |
| `pareto_shortlist.csv` | 4 | Raw non-dominated front (macrolides; audit only) |
| `filters_pool.csv` / `admet_pool.csv` | 1580 | PAINS/Brenk and rule-based drug-likeness |
| `nominated_candidates.csv` | 51 | Dual-dock gate \(S_U\ge90\) and \(S_{N,\mathrm{dock}}\ge90\) |
| `nominated_shortlist_diverse.csv` | 7 preferred | Chemistry nomination; follow-up is GSK-3008348 + Vecabrutinib |
| `novelty_pool.csv` | 1580 | NN Tanimoto vs known actives |
| `pareto_robustness/` | — | Top-k% / tau gates / bootstrap front membership |
| `funnel_snapshot.json` | — | Counts for Methods |

Reproduce nomination: `python3 scripts/14_candidate_nomination.py --tau 90`.
