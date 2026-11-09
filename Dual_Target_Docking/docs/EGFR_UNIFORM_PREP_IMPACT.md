# EGFR uniform RDKit/Meeko rebuild — impact

Date: 2026-09-19  
Decision: CASE 2 promoted. The uniform table is the only current EGFR production truth.

## Docking campaign

- Ligand prep: RDKit ETKDGv3 + MMFF≤200 + Meeko 0.7.1; 110/110 PDBQT.
- Vina 1.2.7 five-seed, corrected boxes, timeout 600 s skip: 1090 ok / 10 `timeout_skipped` (EH40_31, all seeds/pockets).
- Independent GNINA: 185 ok / 33 `timeout_skipped` / 2 fail (`EH120_109`, invalid PDBQT atom type `CG0`).
- Both-pocket Vina production seed 20260727: 109/110. Primary complete-case after activity gate: 28 / 37 / 31 / 12 (n=108).

## Numeric impact (new numbers win)

| Channel | Historical deposited | Uniform current |
|---------|----------------------|-----------------|
| Table 2 `summary_min` | 0.3237 [0.1953, 0.4744]; n=28/37/32 | **0.3341 [0.1970, 0.4713]; n=28/37/31** |
| Five-seed range | 0.3237–0.3504 | **0.3168–0.3491** |
| Independent GNINA `summary_min` | 0.265; n≈108 | **0.2269 [0.1038, 0.3731]; n=89 (20/33/26/10)** |
| Pocket-A fixed-score Δ | 0.4621 | **0.4457 [0.2594, 0.6319]** |
| Matched−mismatched | 0.056 [−0.037, 0.155], includes 0 | **0.1071 [0.0058, 0.2201], excludes 0** |

Do not keep 0.3237 by reverting inputs. Historical ablation and corrected-box five-seed tables remain on disk as historical / unavailable current truth.
