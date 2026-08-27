# Acid track operations (Amendment A2 + A2b)

## Status (2026-08-27)

| Step | Status |
|------|--------|
| A1 freeze (24 exploratory) | DONE — `acid_dual_a1_frozen/` |
| Amendment A2 (URAT1 geometry-first) | DONE |
| A2 reference validation | DONE — 5/6 carboxylate refs |
| URAT1 acid-gate benchmark | DONE — A1 OR≈3.2; A2 OR≈1.0 |
| A2 clinical seed42 | DONE — 59 dual loose keep |
| **A2b NLRP3 IFP / overlap / key residues** | **DONE (seed42)** — loose 78 → structural 74; dual structural 56 |
| NLRP3 known-ligand panel | IN PROGRESS — `05_metrics/nlrp3_structural_panel/` |
| Seeds 43/44 dual dock | IN PROGRESS |
| Competition shortlist | DONE (seed42 + IFP annotation) |
| MD (L7) | CLOSED |

## NLRP3 gates

| Gate | Definition | Use |
|------|------------|-----|
| Loose `keep_nlrp3_pose` | COM≤6 Å vs NP3-146 + CNNscore≥0.5 | sensitivity / funnel annotation |
| Structural `keep_nlrp3_structural` | + overlap≥0.50 + IFP Jaccard≥0.50 + key contacts≥5/7 + no clash | nomination preference |

Self-dock NP3-146 seeds 42/43/44: RMSD ≈ 0.82 / 0.67 / 0.68 Å.

## PF-04620110 note

Passes loose NLRP3 pose; **fails** structural IFP (overlap 0.48, IFP 0.45, key 4/7).
Keep as pathway-anchored primary with explicit claim:
`experimental NLRP3-pathway evidence + pocket-accessible pose ≠ NACHT co-crystal mode match`.

## Do not

- Overwrite `acid_dual_a1_frozen/`
- Claim docking affinity or validated NLRP3 binding
- Open MD before multiseed stability + authorization
