# Layer-3 cognate QC — Track B eight receptors

Protocol: Vina 1.2.7, seed 20260727, num_modes 9, energy_range 3; box = cognate AABB + 5 Å/axis (min 20 Å).
Gate: best-of-9 heavy-atom RMSD < 2.0 Å at E=8 (fallback E=16 not needed).
9V8H keeps PG08-NL peptide (chain B). **Production docking not started.**

**Result: 8/8 PASS**

| protein | PDB | cognate | E_pass | top-1 RMSD | best-of-9 RMSD | best mode | status |
|---------|-----|---------|-------:|-----------:|---------------:|----------:|--------|
| F2 | 4UDW | N6L | 8 | 0.382 | 0.382 | 1 | PASS |
| F10 | 2JKH | BI7 | 8 | 0.658 | 0.658 | 1 | PASS |
| JAK1 | 6N7A | KEV | 8 | 0.459 | 0.459 | 1 | PASS |
| TYK2 | 3LXP | IZA | 8 | 0.197 | 0.197 | 1 | PASS |
| JAK2 | 8BXH | C87 | 8 | 4.064 | 0.807 | 2 | PASS |
| PPARG | 9V8H | BRL | 8 | 6.493 | 1.459 | 3 | PASS |
| PPARA | 6LXA | EPA | 8 | 7.508 | 1.098 | 5 | PASS |
| PPARD | 5U3Q | 7UJ | 8 | 1.452 | 1.452 | 1 | PASS |
