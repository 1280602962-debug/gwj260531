# URAT1 TrueDecoy Benchmark

Property-matched enrichment set for URAT1 docking / scoring validation, framed after
**Gu et al., *Nat. Mach. Intell.* 2025** (TrueDecoy vs RandomDecoy).

## Design

| Item | Choice |
|------|--------|
| Actives | Curated URAT1, `pActivity >= 6.0` (n = 469) |
| Decoy pool | Distill subset D (unlabeled diversity negatives) |
| Matching | Round-robin MW / logP / TPSA / HBD / HBA / rotatable bonds (DUD-E-inspired windows); 1.5× relaxed top-up |
| Near-analog filter | Max Morgan TC to any active ≤ 0.5 |
| Target ratio | 1 : 30 (active : decoy) |
| Achieved TrueDecoy ratio | 1 : 12.95 (n_decoy = 6073) |
| Actives with ≥1 decoy | 460 / 469 |
| Decoys / covered active | min 1, median 10.0, mean 13.2, max 30 |
| RandomDecoy companion | Same size random sample from pool (n = 6073) |
| Seed | 42 |

### Property windows (strict)

| Descriptor | Window |
|------------|--------|
| MolWt | ±40 Da |
| MolLogP | ±1.0 |
| TPSA | ±25 Å² |
| NumHDonors | ±1 |
| NumHAcceptors | ±2 |
| NumRotatableBonds | ±2 |

## Files

| File | Content |
|------|---------|
| `actives.csv` | Potent URAT1 actives + descriptors |
| `true_decoys.csv` | Property-matched decoys |
| `true_decoy_benchmark.csv` | Combined set with `label` (1=active, 0=decoy) and `decoy_class=true` |
| `random_decoys.csv` | Random unmatched decoys (same n when possible) |
| `random_decoy_benchmark.csv` | Combined RandomDecoy set |
| `matching_assignments.csv` | Active→decoy pairs with property distance / window pass |
| `summary.json` | Counts, windows, property stats |

## Usage notes

- **TrueDecoy** is the harder test: physics-based docking/scoring should be judged here.
- **RandomDecoy** mirrors easier / library-like enrichment (subset A vs raw D is closer to this).
- Do **not** train ML models on these decoy labels; they are putative inactives for enrichment only.
- Pool size (~8k subset D) limits the achievable ratio below the 1:30–50 ideal; rebuild with a larger `--pool` if needed.

## Rebuild

```bash
python3 scripts/build_urat1_true_decoy.py \
  --pactivity-min 6.0 \
  --ratio 30 \
  --max-sim-to-active 0.5 \
  --seed 42
```
