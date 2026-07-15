# SOTA Gap Map

| Candidate Contribution | What SOTA Already Does | User Evidence (local) | Real Gap | Claim Strength | Risk |
|------------------------|------------------------|-----------------------|----------|----------------|------|
| **C1. Critical benchmark: common computational isoform-selectivity predictors fail on JNK1/2/3** | Many VS papers use docking score differences / IFP / ML Δ as selectivity filters without hard negative benchmarks | Δsel direction accuracy **43%**; Gly87 occupancy **non-discriminative**; ML selectivity **F1=0** (n_pos=8); ATP ~98% identity | Quantified, multi-method failure on literature-calibrated JNK panel — still under-reported as a **primary** result | **High** (if framed as analysis/benchmark) | Reviewer: “known that homologous kinases are hard” — answer with **quantified** multi-method table + purchase decoupling |
| **C2. End-to-end commercial-library pipeline to pose-credible JNK-family shortlist with wet-lab design** | Many VS pipelines exist; fewer with full ML→dock→ADMET→MD→purchase + literature controls | Funnel 4979→157→25→16; MD QC G1 3/4; purchase HIT690/2157 + E1/CC-90001 | Prospective, documented purchase for **family-binder enrichment** (not selectivity) | **Medium** | n=2 thin; need IC50; chemotype novelty uncertain |
| **C3. Discovery of JNK1-selective inhibitors** | Hard; literature structural rules favor JNK2/3 (Ile106/Leu); few co-crystal explanations for JNK1 bias | No IC50 yet; best MD-bias molecule 2231 **not purchased**; 2157 Δsel_dock **negative** | **Not supported by current design** | **Low / invalid as Core** | Over-claim → rejection |
| **C4. MD can rank isoform selectivity** | Some MD papers claim occupancy/RMSD asymmetry as selectivity | G3: E1 hinge **mis-ranks**; SP600125 active with low hinge | Project evidence **refutes** MD-as-selectivity | Use only as **negative / caution** | If claimed positive → false |
| **C5. Kinome-selective JNK inhibitors** | Requires panel / counter-screens | **No data** | Future work only | None now | Common hinge chemotypes → prior risk of polypharmacology |

## Gap Summary

1. **Real, publishable gap:** honest, multi-method demonstration that **isoform selectivity should not be purchased from Δsel/Gly87/ML labels** on JNK-like sites — with a still-useful **family-binder enrichment pipeline**.
2. **Not a gap you currently close:** JNK1-selective discovery; kinome selectivity; MD-based isoform adjudication.
3. **Evidence missing for strong chemistry claim:** IC50 for 690/2157; optional p38/ERK counter-screen; chemotype-vs-known-JNK similarity quantification beyond Tc~0.22 to E1/Q63.
