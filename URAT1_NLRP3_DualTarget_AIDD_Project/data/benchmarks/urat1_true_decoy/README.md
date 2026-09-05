# URAT1 TrueDecoy / RandomDecoy Benchmark

Framed after **Gu et al., *Nat. Mach. Intell.* 2025** (hard negatives for protocol
selection + library-random companion). This is a **single-target URAT1** set — not a
replica of Gu's multi-target BindingDB TrueDecoy.

> **Current official set (2026-07-24):** taosu commercial-library pool, ratio **1:10**,
> True∩Random decoy overlap **= 0**.  
> **Superseded:** older builds from distill subset D (~8k) with overlapping RandomDecoys
> and incomplete 1:10 fill — do **not** use those numbers or files for protocol selection.

## Pool construction (taosu)

| Step | Detail |
|------|--------|
| Source library | `taosu_20210823_100w_asteroid_murcko_protonized.csv` (~1M) |
| Sample | Reservoir sample 200k (seed=42) |
| Envelope filter | Property envelope around actives (1.5× windows), cap **60k** |
| Meta | `taosu_pool_prefilter_meta.json` |

Actives / experimental inactives come from GitHub `data/processed/urat1_curated.csv`.

## Design

| Item | Choice |
|------|--------|
| Actives | Curated URAT1, `pActivity >= 6.0` (n = **469**) |
| TrueDecoy negatives | (1) experimental weak/inactives `pActivity < 5.0` (n = **80**); (2) property-matched from taosu prefiltered pool (n = **4610**) |
| Matching | Round-robin MW / logP / TPSA / HBD / HBA / rotatable bonds; **all 4610 strict-window** |
| Near-analog filter (matched only) | Max Morgan TC to any active ≤ 0.5 |
| Target / achieved TrueDecoy ratio | **1 : 10** (n_decoy = **4690**) |
| Actives with ≥1 matched decoy | 469 / 469 |
| Matched decoys / covered active | min 1, median 10.0, mean 9.8, max 10 |
| RandomDecoy | Gu-style random from **remaining** taosu pool only; n = **4690**; **zero SMILES overlap** |
| Achieved RandomDecoy ratio | **1 : 10**; True∩Random decoy overlap = **0** |
| Unique docking pool | **9849** SMILES (`unique_docking_pool.csv`) — dock once |
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
| `taosu_pool_prefilter_meta.json` | Pool sampling / envelope provenance |

## Usage notes

- **TrueDecoy** is the harder test (experimental weak + property-matched).
- **RandomDecoy** is the easier / VS-like control; must not overlap TrueDecoy negatives.
- Dock **`unique_docking_pool.csv` once**, then join scores into each benchmark for EF/AUC.
- Do **not** train ML models on these decoy labels.
- Distill **subset D** is a separate retrospective resource; it is **not** the current protocol-selection decoy pool.

## Rebuild

```bash
# Requires a large prefiltered pool CSV (taosu-derived), not distill_subset_d.csv
python3 scripts/build_urat1_true_decoy.py \
  --actives data/processed/urat1_curated.csv \
  --pool /path/to/taosu_pool_prefiltered.csv \
  --pactivity-min 6.0 \
  --inactive-pactivity-max 5.0 \
  --ratio 10 \
  --max-sim-to-active 0.5 \
  --seed 42
```
