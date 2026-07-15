# Research Dossier — Target Scene (Journal / CADD)

## Venue Requirements (candidate cluster)

| Venue | Fit for this project | Key expectations |
|-------|----------------------|------------------|
| **JCIM** (ACS) | Best fit if **method + negative benchmark** is the headline | Rigorous validation, data/code availability, clear claim boundaries; computational methods papers welcome honest negatives |
| **Molecules / IJMS** (MDPI) | Fit for **pipeline + experimental IC50** narrative | Faster review; require clear novelty statement; kinome caveat must be explicit |
| **ChemMedChem** | Fit if **new binders with chemotype novelty + IC50** | Stronger chemistry/biology expectation; n=2 is thin unless chemotype is clearly new |
| **Front. Pharmacol. / RSC Adv.** | Fallback if novelty is mainly workflow + honesty | Accept applied CADD with wet-lab closure |

**Decision rule (PaperSpine):** lock venue **after** `confirmed_contribution.md`. If Core contribution = “selectivity predictors fail on 98% identical ATP sites,” prefer **JCIM**. If Core = “commercial-library JNK-family binders with IC50,” prefer **Molecules/ChemMedChem**.

## Review Criteria (what reviewers will score)

1. **Is the contribution falsifiable and bounded?** Over-claiming “JNK1-selective discovery” will reject.
2. **Are negative results quantified?** Δsel 43%, Gly87 non-discriminative, ML F1=0 must be front-stage, not buried.
3. **Is wet-lab closure adequate?** n=2 new + 2 positives is minimal; kinome panel absence must be Future Work, not hidden.
4. **Reproducibility:** SMILES, PDB IDs, thresholds, scripts, purchase IDs (HIT103871685 / HIT101201113).
5. **Separate selectivity types:** isoform vs kinome — confuse them → desk rejection risk.

## Accepted Paper Patterns (CADD journals)

| Pattern | Example use here |
|---------|------------------|
| Pipeline paper with experimental validation | ML→dock→ADMET→MD→IC50 |
| Method failure / cautionary benchmark | Isoform score-difference pitfalls on homologous kinases |
| Prospective library campaign | Enamine/Chembridge shortlist + purchase |
| Calibration with literature controls | E1, CC-90001, SP600125 |

## Constraints for This Paper (hard, from materials)

- Project **already pivoted**: computational JNK1 selectivity is **not** the decision driver.
- Purchase: **690 + 2157** only (new); positives **E1 + CC-90001** in hand.
- Strongest MD JNK1-bias hypothesis molecule **2231 not purchased** → cannot claim prospective test of best bias hypothesis.
- ATP-site identity ~98%; literature structural rules favor **JNK2/3**, not JNK1.
- No kinome panel yet.
- Single-replica MD; 2231 extended MD used ligand restraints.

## Materials Authority

Primary local evidence: `docs/JNK1_PROJECT_REPORT.md`, `results/model_comparison/`, `results/docking_validation/`, `results/ml_external_validation/`, `data/purchase/`, `docs/popular_science/data_tables/27_*.csv`.
