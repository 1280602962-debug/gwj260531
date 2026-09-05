# A2 reference validation (seed 42)

- Arg threshold: **7.7027 Å**
- A1 carboxylate refs pass: **1/6**
- A2 carboxylate refs pass: **5/6**

| ligand | rule | keep | Arg477 (Å) | best-acid among modes (Å) | status |
|--------|------|-----:|-----------:|--------------------------:|--------|
| lesinurad | a1 | False | 14.11 | 2.99 | a1_cnnscore_top1 |
| lesinurad | a2 | True | 3.93 | 2.99 | a2_geometry_then_cnnscore |
| verinurad | a1 | False | 12.30 | 2.88 | a1_cnnscore_top1 |
| verinurad | a2 | True | 3.14 | 2.88 | a2_geometry_then_cnnscore |
| probenecid | a1 | False | 14.24 | 12.88 | a1_cnnscore_top1 |
| probenecid | a2 | False | nan | 12.88 | a2_no_geometry_pass |
| puliginurad | a1 | False | 14.55 | 2.86 | a1_cnnscore_top1 |
| puliginurad | a2 | True | 2.92 | 2.86 | a2_geometry_then_cnnscore |
| SHR-4640 | a1 | False | 14.82 | 3.04 | a1_cnnscore_top1 |
| SHR-4640 | a2 | True | 3.14 | 3.04 | a2_geometry_then_cnnscore |
| GSK-3008348 | a1 | True | 2.98 | 2.95 | a1_cnnscore_top1 |
| GSK-3008348 | a2 | True | 2.98 | 2.95 | a2_geometry_then_cnnscore |
