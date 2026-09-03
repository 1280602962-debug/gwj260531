# Leave-cognate-out sensitivity

Source: `data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv`; production seed 20260727.
The analysis removes the exact co-crystallized ligand when it is a member of the main panel, then recomputes the two Table 2 directional AUROCs, their descriptive `summary_min`, and the Table 3 Dual-versus-neither `vina_mean` AUROC.
It tests sensitivity to the single exact cognate molecule; it does not remove cognate-like chemotypes and is not described as a train/test leakage analysis.

| pair | removed cognate | n dual before/after | D/A before→after | D/B before→after | summary_min before→after | D/neither before→after |
|---|---|---:|---:|---:|---:|---:|
| EGFR/HER2 | TAK-285 (EH40_01; CHEMBL1614725) | 28/27 | 0.6664→0.6637 | 0.4297→0.4167 | 0.4297→0.4167 | 0.7560→0.7531 |
| PIK3CA/mTOR | PI-103 (PM48_01; CHEMBL573339) | 18/17 | 0.7143→0.7017 | 0.6921→0.6740 | 0.6921→0.6740 | 0.5139→0.5000 |

Removing either exact cognate ligand produced only small point-estimate changes and did not alter the qualitative interpretation. Residual receptor–chemotype favorability is not excluded by this one-ligand sensitivity.
