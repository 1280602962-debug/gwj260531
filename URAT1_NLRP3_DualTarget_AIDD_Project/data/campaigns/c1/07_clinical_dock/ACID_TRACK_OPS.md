# Acid track operations (Amendment A2 + A2b)

## Status (2026-08-28)

| Step | Status |
|------|--------|
| A1 freeze (24 exploratory) | DONE — `acid_dual_a1_frozen/` |
| Amendment A2 (URAT1 geometry-first) | DONE |
| A2 reference validation | DONE — 5/6 carboxylate refs |
| URAT1 acid-gate benchmark | DONE — A1 OR≈3.2; A2 OR≈1.0 |
| A2 clinical seed42 | DONE — 59 dual loose keep |
| **A2b NLRP3 IFP / overlap / key residues** | **DONE seeds 42/43/44** — dual structural 56 / 57 / 53 |
| NLRP3 known-ligand panel | DONE — positives 10/10 vs bg 11/20; Fisher p≈0.013 |
| Seeds 43/44 dual dock | DONE — 311/312 each; `REP_07837` NLRP3 timeout |
| Multi-seed ≥2/3 stability | DONE — eligible 40 after audit |
| Competition shortlist | DONE — `acid_shortlist_a2_competition.csv` |
| MD (L7) | CLOSED (`md_authorized=false`) |

## NLRP3 gates

| Gate | Definition | Use |
|------|------------|-----|
| Loose `keep_nlrp3_pose` | COM≤6 Å vs NP3-146 + CNNscore≥0.5 | sensitivity / funnel annotation |
| Structural `keep_nlrp3_structural` | + overlap≥0.50 + IFP Jaccard≥0.50 + key contacts≥5/7 + no clash | nomination preference |

Self-dock NP3-146 seeds 42/43/44: RMSD ≈ 0.82 / 0.67 / 0.68 Å.

## Multi-seed dual keep

| seed | URAT1 | NLRP3 loose | dual loose | NLRP3 struct | dual struct |
|------|------:|------------:|-----------:|-------------:|------------:|
| 42 | 121 | 78 | 59 | 74 | 56 |
| 43 | 120 | 77 | 59 | 74 | 57 |
| 44 | 121 | 83 | 61 | 75 | 53 |

≥2/3 dual loose = 59; ≥2/3 dual structural ≈ 54.

## PF-04620110 note

Dual-pass 2/3 seeds. **Fails structural IFP on all three seeds** (seed42: overlap 0.48 / IFP 0.45 / keys 4).
Keep as pathway-anchored primary with explicit claim:
`experimental NLRP3-pathway evidence + pocket-accessible pose ≠ NACHT co-crystal mode match`.

## Do not

- Overwrite `acid_dual_a1_frozen/`
- Claim docking affinity or validated NLRP3 binding
- Open MD before explicit shortlist authorization
