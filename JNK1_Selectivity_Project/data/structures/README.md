# JNK Ensemble Docking Structures

Curated **5-receptor panel** for F3 ensemble docking (Type I / DFG-in).

| File | Description |
|------|-------------|
| [`docking_ensemble_pdb.csv`](docking_ensemble_pdb.csv) | Master table (KLIFS + RCSB metadata) |
| [`../../config/docking_ensemble.yaml`](../../config/docking_ensemble.yaml) | Machine-readable config for docking scripts |

## Selected panel

| Isoform | Primary | Secondary | Note |
|---------|---------|-----------|------|
| JNK1 | 3ELJ | 4L7F | 2-structure ensemble |
| JNK2 | 3E7O | — | **sole** drug co-crystal (chain A) |
| JNK3 | 3TTI | 4WHZ | 2-structure ensemble |

## Scoring

```text
Score_JNK1 = mean(dock_3ELJ, dock_4L7F)
Score_JNK2 = dock_3E7O
Score_JNK3 = mean(dock_3TTI, dock_4WHZ)
Selectivity = Score_JNK1 - max(Score_JNK2, Score_JNK3)
```

Run co-crystal redocking and benchmark direction checks before screening `top500_diverse.csv`.
