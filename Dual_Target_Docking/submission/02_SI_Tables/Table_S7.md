# Table S7. Out-of-fold performance of ligand-chemistry models

配体化学模型的折外判别表现

| Target pair | Direction | Common OOF n | Physicochemical AUROC | ECFP4 AUROC | ECFP4 + S1 AUROC | Matched S1 AUROC | ECFP4 + S1 increment |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | D_vs_A | 65 | 0.6757 | 0.8002 | 0.7838 | 0.6274 | -0.0164 |
| EGFR/HER2 | D_vs_B | 59 | 0.7592 | 0.8790 | 0.8767 | 0.2408 | -0.0023 |
| JAK1/JAK2 | D_vs_A | 64 | 0.6030 | 0.9185 | 0.9199 | 0.5908 | 0.0015 |
| JAK1/JAK2 | D_vs_B | 64 | 0.9014 | 0.9678 | 0.9707 | 0.7485 | 0.0029 |
| JAK1/TYK2 | D_vs_A | 63 | 0.5655 | 0.8488 | 0.8438 | 0.6371 | -0.0050 |
| JAK1/TYK2 | D_vs_B | 63 | 0.8407 | 0.9042 | 0.9032 | 0.4032 | -0.0010 |
| PIK3CA/mTOR | D_vs_A | 64 | 0.6094 | 0.8213 | 0.8428 | 0.7363 | 0.0215 |
| PIK3CA/mTOR | D_vs_B | 64 | 0.2939 | 0.7285 | 0.7256 | 0.5801 | -0.0029 |
| AChE/BChE | D_vs_A | 52 | 0.6815 | 0.8993 | 0.8948 | 0.6341 | -0.0044 |
| AChE/BChE | D_vs_B | 55 | 0.7884 | 0.8201 | 0.8148 | 0.6442 | -0.0053 |
| F2/F10 | D_vs_A | 63 | 0.6139 | 0.9405 | 0.9335 | 0.4204 | -0.0071 |
| F2/F10 | D_vs_B | 63 | 0.4889 | 0.6633 | 0.6996 | 0.3740 | 0.0363 |
| PPARG/PPARA | D_vs_A | 63 | 0.7198 | 0.8206 | 0.8065 | 0.6764 | -0.0141 |
| PPARG/PPARA | D_vs_B | 64 | 0.7158 | 0.6670 | 0.6611 | 0.6021 | -0.0059 |
| PPARA/PPARD | D_vs_A | 64 | 0.4375 | 0.9316 | 0.9287 | 0.6787 | -0.0029 |
| PPARA/PPARD | D_vs_B | 64 | 0.3115 | 0.8740 | 0.8789 | 0.5107 | 0.0049 |


All models and S1 are evaluated on identical members within each direction. The increment is AUROC(ECFP4 + S1) minus AUROC(ECFP4). Existing frozen five-fold assignments are retained; Q110 uses sorted canonical Murcko scaffolds assigned cyclically to folds 0-4. Standardization is fitted in training folds. Full descriptors, scaffolds, similarities, folds and OOF predictions are supplied as data files; no new model is fitted during table preparation.
