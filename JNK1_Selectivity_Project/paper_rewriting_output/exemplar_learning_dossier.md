# Exemplar Learning Dossier (flash tier)

> Exemplars teach **structure and rhetoric**, not data. Project numbers come only from local materials.

## Exemplar Inventory

| ID | Exemplar / class | Why learn from it | Transferable pattern |
|----|------------------|-------------------|----------------------|
| E1 | Bennett et al., *J Med Chem* 2021 (CC-90001) | Field precedent for JNK1-biased clinical candidate | Biology motivation → chemotype SAR → assay panel; **no public co-crystal** → honest structural limits |
| E2 | Pan et al., *J Med Chem* 2024 (E1) | Recent JNK1 enzyme-biased chemotype | Report IC50 matrix across isoforms; SBDD without over-claiming structure |
| E3 | Typical JCIM “prospective VS + experiment” papers | Venue pattern | Funnel figure → validation tables → purchased hits → IC50 |
| E4 | Cautionary / negative computational studies in JCIM/JCIM-like venues | Method honesty | Benchmark table showing when docking/ML **fails**; claim “do not use X for Y” |
| E5 | Kinase selectivity kinome-scan papers | Anticipated reviewer ask | Even if not done, learn how they **bound** isoform vs kinome claims |
| E6 | Local project report `JNK1_PROJECT_REPORT.md` | Author’s own structured narrative | Single-funnel mainline + failed selectivity exploration as separate section |

## Structural Patterns to Imitate

1. **Introduction ladder:** disease/biology → pan-JNK toxicity concern → computational hope → **homology barrier** → RQ becomes “can we enrich family binders *and* can common selectivity predictors work?”
2. **Methods:** justify each gate (F1 recall vs specificity; Glide redock; MD QC thresholds; ADMET).
3. **Results order:** (A) ML potency models; (B) docking pose validation; (C) **selectivity-method failure benchmark**; (D) shortlist/MD; (E) purchase + planned/actual IC50.
4. **Discussion:** return to contribution; separate isoform vs kinome; state what IC50 can/cannot prove with n=2.

## Rhetorical Patterns

- Prefer **“enrichment / binder / pose-credible candidate”** over **“selective inhibitor”**.
- Put negative results in Results, not only Limitations.
- Use G3 controls (E1, CC-90001) as **assay calibration**, not as proof that MD predicts selectivity (MD can mis-rank E1).

## Language Patterns (English manuscript)

- Allowed strong: “Δsel direction accuracy was 43% on the literature benchmark.”
- Soften: “MD hinge asymmetry suggests a testable JNK1-bias hypothesis.”
- Forbid without IC50: “We discovered JNK1-selective inhibitors.”
