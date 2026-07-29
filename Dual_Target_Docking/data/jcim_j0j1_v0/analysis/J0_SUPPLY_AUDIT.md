# J0 — Expanded ChEMBL strict hard-negative supply audit

**Date:** 2026-07-29  
**Script:** `scripts/run_j0_supply_audit.py`  
**Tables:** `tables/j0_strict_label_supply.csv`, `j0_candidate_pairs.csv`, `j0_fetch_queue.csv`  
**Rule:** dual = both ≥6.5; A_only = A≥6.5 & B≤5.5; B_only = converse; gray = both measured else.  
**Panel flags:** `supports_strict_panel` = both hard-neg classes ≥50; `supports_thin_panel` = ≥20.

## Scope

| Item | Count |
|------|------:|
| Candidate pairs listed | 53 (incl. private holdout + aliases) |
| **Audited now** (both `mols_*.json` present) | **49** |
| Fetch queue (new targets; ChEMBL API was 500) | 22 |
| Explicit exclude | NLRP3/JNK1 (private holdout) |

> ChEMBL Web API returned HTTP 500 during this run. Audited pairs use **already-cached** dictionaries only. Additional literature pairs that need new targets are queued — **not** fabricated.

## One-line result

**Among 49 auditable pairs, only 4 meet strict hard-neg ≥50 on both sides; after dropping metal/isozyme HDAC1/HDAC6, only 3 remain usable for a conventional docking Tier-S panel — confirming public-data supply is the bottleneck, not docking budget.**

## Strict-pass pairs (Y)

| pair | both | strict A/B | min HN | caveat |
|------|-----:|------------|-------:|--------|
| HDAC1/HDAC6 | 3987 | 93/494 | 93 | **Zn metal + isoform** — not Tier-S docking primary |
| **PIK3CA/MTOR** | 2713 | 80/81 | 80 | best development pair (already docked) |
| **ACHE/BCHE** | 2537 | 189/78 | 78 | best *new* conventional pair |
| **PIK3CA/PIK3CB** | 1990 | 56/67 | 56 | isoform — narrative “too close”; useful control seat |

## Thin panel only (T: 20≤minHN<50)

| pair | min HN | note |
|------|-------:|------|
| PIK3CB/MTOR | 25 | pathway; thin |

## Named failures (selected)

| pair | strict A/B | min HN | meaning |
|------|------------|-------:|---------|
| EGFR/HER2 | 39/**7** | 7 | **cannot** support thick strict four-class; case/supply-limited |
| MCL1/BCL2L1 | 41/12 | 12 | pose-gold PPI; Tier T |
| EGFR/MET | 17/10 | 10 | below thin |
| BRD4/HDAC* / JAK2/HDAC1 | ≤3 | 0–3 | metal + no supply |
| PARP1/MET, CDK6/BRD4 | ≤1 | 0 | essentially empty paired set |

## Forest (top 15 by min_strict_hardneg)

See CSV for full 49. Ranking is dominated by homologous hydrolases / PI3K pathway / HDAC isoforms; most literature “hot” epigenetic hybrids fail the strict hard-neg gate.

## Implication for C2 (paper)

Public DualFourClass evaluation is **supply-limited**: expanding docking cannot create EGFR strict B_only molecules that ChEMBL does not contain. JCIM-scale K=4 must be chosen from the tiny Y/T set (plus intentional Tier-T cases), not from literature popularity alone.

## Reproduce

```bash
python3 Dual_Target_Docking/data/jcim_j0j1_v0/scripts/run_j0_supply_audit.py
```

When ChEMBL recovers, fetch `j0_fetch_queue.csv` targets into `data/public_pair_selection/mols_*.json` and re-run.
