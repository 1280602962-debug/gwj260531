# Forced-recovery URAT1 A1 gate table (seed 42, free dock)

- Protocol: gnina exhaustiveness=32, num_modes=9, CNNscore-selected pose, Arg threshold **7.7027 Å** (A1), pocket COM ≤ 6.0 Å vs lesinurad crystal.
- **pass_A1_carboxylate_gate** = carboxylate oxygens present AND Arg≤7.7027 AND COM≤6.
- Phenol-class ligands (benzbromarone, dotinurad): carboxylate gate reported separately; see `acid_or_phenolate_arg477_min_A`.
- This table is an **audit**, not a new scientific lock.

| ligand | prep CO2 | Arg CO2 (Å) | pass Arg A1 | COM (Å) | pass pocket | **pass A1 gate** | phenolate/acid Arg (Å) | best CO2 Arg any mode (Å) |
|---|:---:|---:|:---:|---:|:---:|:---:|---:|---:|
| lesinurad | yes | 14.11 | no | 1.63 | yes | **no** | 14.11 | 2.99 |
| benzbromarone | no | — | no | 1.17 | yes | **no** | — | — |
| dotinurad | no | — | no | 0.95 | yes | **no** | — | — |
| verinurad | yes | 12.30 | no | 0.75 | yes | **no** | 12.30 | 2.88 |
| probenecid | yes | 14.24 | no | 2.69 | yes | **no** | 14.24 | 12.88 |
| puliginurad | yes | 14.55 | no | 1.62 | yes | **no** | 14.55 | 2.86 |
| SHR-4640 | yes | 14.82 | no | 3.38 | yes | **no** | 14.82 | 3.04 |
| GSK-3008348 | yes | 2.98 | yes | 0.65 | yes | **yes** | 2.98 | 2.95 |

**Summary:** 1 / 8 pass A1 carboxylate gate; 6 / 8 prepared with carboxylate.

CSV: `data/campaigns/c1/03_forced_recovery/forced_recovery_A1_gate_table_seed42.csv`
