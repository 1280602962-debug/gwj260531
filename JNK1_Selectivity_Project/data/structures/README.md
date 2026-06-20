# JNK Ensemble Docking Structures

Curated **6-receptor panel** for F3 ensemble docking (Type I / DFG-in).

| File | Description |
|------|-------------|
| [`docking_ensemble_pdb.csv`](docking_ensemble_pdb.csv) | Master table (KLIFS + RCSB metadata) |
| [`../../config/docking_ensemble.yaml`](../../config/docking_ensemble.yaml) | Machine-readable config for docking scripts |

## Selected panel

| Isoform | Primary | Secondary |
|---------|---------|-----------|
| JNK1 | 3ELJ | 4L7F |
| JNK2 | 3E7O | 7N8T |
| JNK3 | 3TTI | 4WHZ |

## Scoring

```text
Score_JNK1 = mean(dock_3ELJ, dock_4L7F)
Score_JNK2 = mean(dock_3E7O, dock_7N8T)
Score_JNK3 = mean(dock_3TTI, dock_4WHZ)
Selectivity = Score_JNK1 - max(Score_JNK2, Score_JNK3)
```

Run co-crystal redocking and benchmark direction checks before screening `top500_diverse.csv`.
