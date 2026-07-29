# GNINA rescore — STATUS: DONE (mode_01 CNN minimize/rescore)

**Updated:** auto after JCIM gnina orchestrator

## Binary
`/mnt/d/CADD paper exercise/gnina/bin/gnina` (v1.3.2), CPU `--no_gpu`

## Protocol
- Input: Vina `mode_01.pdbqt` → SDF (obabel)
- `gnina --cnn_scoring rescore --minimize --seed 20260727 --cpu 1`
- Packs: AChE/BChE, PIK3CA/PIK3CB, PM48 RDKit, EGFR if available
- Tables: `tables/scores_gnina_long.csv`, `tables/scores_gnina_best.csv`
