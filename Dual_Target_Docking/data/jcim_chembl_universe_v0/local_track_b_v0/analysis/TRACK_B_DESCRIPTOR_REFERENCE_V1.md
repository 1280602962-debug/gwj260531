# Track B four-descriptor reference (same scored ligands as Vina)

Prespecified panel: heavy-atom count, MW, cLogP, TPSA. Largest organic fragment, RDKit.
Estimand matches Table 2 companion columns. Best single descriptor = highest `summary_min`.
Descriptive reference only. **Does not replace Table 2.** Zero docking.

| pair | Vina summary_min | best descriptor | descriptor summary_min | Δ (Vina − desc) |
|------|-----------------:|-----------------|-----------------------:|----------------:|
| F2/F10 | 0.3448 | clogp | 0.5151 | -0.1703 |
| JAK1/TYK2 | 0.3649 | clogp | 0.5796 | -0.2147 |
| JAK1/JAK2 | 0.5884 | heavy | 0.5781 | 0.0103 |
| PPARG/PPARA | 0.6492 | tpsa | 0.6274 | 0.0218 |
| PPARA/PPARD | 0.4463 | clogp | 0.5635 | -0.1172 |

Vina does not beat the best single descriptor by more than +0.022. PPARG/PPARA TPSA `summary_min` is 0.627 (Vina 0.649). JAK1/TYK2 Dual-versus-neither is 0.770 on Vina mean and 0.810 on heavy-atom count.

Artifacts: `tables/track_b_descriptor_summary_v1.csv`, `tables/track_b_descriptor_auroc_v1.csv`, `tables/track_b_descriptor_by_class_v1.csv`.
