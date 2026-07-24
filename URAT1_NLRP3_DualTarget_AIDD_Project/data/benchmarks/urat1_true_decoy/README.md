# URAT1 TrueDecoy / RandomDecoy Benchmark

Framed after **Gu et al., *Nat. Mach. Intell.* 2025** (hard negatives for protocol
selection + library-random companion). This is a **single-target URAT1** set — not a
replica of Gu's multi-target BindingDB TrueDecoy.

## Local rebuild note (taosu library)

- **Actives / experimental inactives**: GitHub `urat1_curated.csv` (same thresholds as official protocol).
- **Decoy pool**: `taosu_20210823_100w_asteroid_murcko_protonized.csv` (~1M), reservoir-sampled 200k (seed=42), property-envelope filtered, capped at 60k for matching.
- Built with `scripts/build_urat1_true_decoy.py` (Gu *Nat Mach Intell* 2025 framing).


## Design

| Item | Choice |
|------|--------|
| Actives | Curated URAT1, `pActivity >= 6.0` (n = 469) |
| TrueDecoy negatives | (1) experimental weak/inactives `pActivity < 5.0` (n = 80); (2) property-matched unlabeled from distill subset D (n = 4610) |
| Matching | Round-robin MW / logP / TPSA / HBD / HBA / rotatable bonds; 1.5× relaxed top-up |
| Near-analog filter (matched only) | Max Morgan TC to any active ≤ 0.5 |
| Target ratio | 1 : 10 (active : TrueDecoy negative) |
| Achieved TrueDecoy ratio | 1 : 10.00 (n_decoy = 4690) |
| Actives with ≥1 matched decoy | 469 / 469 |
| Matched decoys / covered active | min 1, median 10.0, mean 9.8, max 10 |
| RandomDecoy | Gu-style random draw from **remaining** subset D only; target n = TrueDecoy negatives; **zero SMILES overlap** with TrueDecoy |
| Achieved RandomDecoy ratio | 1 : 10.00 (n_decoy = 4690); True∩Random decoy overlap = 0 |
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
| `experimental_inactives.csv` | Curated weak/inactives used as TrueDecoy negatives |
| `true_decoys.csv` | All TrueDecoy negatives (`decoy_source` = experimental_inactive \| property_matched) |
| `true_decoy_benchmark.csv` | Actives + TrueDecoy negatives (`label` 1/0) |
| `random_decoys.csv` | Non-overlapping random library decoys |
| `random_decoy_benchmark.csv` | Actives + RandomDecoy negatives |
| `matching_assignments.csv` | Active→property-matched decoy pairs |
| `summary.json` | Counts, windows, overlap check |
| `unique_docking_pool.csv` | Unique SMILES across both benchmarks (dock once) |

## Usage notes

- **TrueDecoy** is the harder test (experimental weak + property-matched).
- **RandomDecoy** is the easier / VS-like control; must not overlap TrueDecoy negatives.
- Dock **`unique_docking_pool.csv` once**, then join scores into each benchmark for EF/AUC.
- Do **not** train ML models on these decoy labels.

## Rebuild

```bash
python3 scripts/build_urat1_true_decoy.py \
  --pactivity-min 6.0 \
  --inactive-pactivity-max 5.0 \
  --ratio 10 \
  --max-sim-to-active 0.5 \
  --seed 42
```
