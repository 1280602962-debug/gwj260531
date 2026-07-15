# Draft prose (Option A) — Intro / Methods notes / RQ-C Results

*Working draft for editing. IC50 Results intentionally omitted until assay data.*

---

## Introduction (draft)

c-Jun N-terminal kinases (JNK1/2/3) are stress-activated MAPKs implicated in fibrosis and inflammation. Clinical and recent medicinal chemistry efforts, including CC-90001 and compound E1, have motivated JNK1-biased inhibition for idiopathic pulmonary fibrosis while avoiding pan-JNK liabilities. However, the ATP-binding sites of JNK isoforms are nearly identical, so isoform preference remains difficult both chemically and computationally.

Structure-based virtual screening often uses docking score differences, residue-occupancy heuristics, or machine-learning selectivity labels as proxies for isoform preference. Related kinases have seen successes with interaction fingerprints or tailored selectivity benchmarks, and free-energy methods can sometimes exploit error cancellation between homologs. Whether inexpensive docking- and ML-based filters are reliable **purchase criteria** for JNK1 versus JNK2/3 has been less systematically stress-tested on literature-calibrated panels.

Here we report an end-to-end commercial-library pipeline—machine-learning family-activity gating, Glide docking, ADMET filtering, and MD pose quality control—designed to enrich **JNK-family binders** for enzymatic testing. In parallel, we evaluate common computational isoform-selectivity predictors on a literature JNK benchmark and show that they fail usability thresholds for directional isoform ranking. Purchase of two shortlisted candidates (library IDs 690 and 2157) was therefore decoupled from those selectivity filters and reserved for prospective JNK1/2/3 IC50 measurement alongside E1 and CC-90001 controls under a pre-registered analysis plan.

---

## Methods notes (software / docking)

**License-compliant option A (preferred if authorized):**  
Molecular docking was performed with Glide XP (Schrödinger Suite, version as used institutionally). Protein preparation and grid generation followed the project ensemble (JNK1: 3ELJ; JNK2: 3E7O; JNK3: 3TTI; secondary structures for benchmark Δsel as listed in SI).

**Open pose-consensus option B (executed in this archive):**  
To assess multi-seed pose stability independent of single Glide runs, purchased ligands were re-docked with AutoDock Vina 1.2.5 (exhaustiveness 16; seeds 1–3) into cognate-ligand-centered boxes on meeko-prepared chain-A receptors. Pairwise heavy-atom RMSD among top poses was used only as a **pose consensus** metric and not to re-rank the historical shortlist.

Selectivity metrics tested on the literature benchmark (not used as purchase hard gates):  
(i) Δsel_dock = min(score_JNK2, score_JNK3) − score_JNK1;  
(ii) Gly87 (KLIFS b.l.37) occupancy heuristic;  
(iii) ML selective-class labels trained on sparse ChEMBL paired actives.  
The ML family gate used p_family = max(pred_JNK1, pred_JNK2, pred_JNK3) ≥ 6.0 as a high-recall filter only.

Pre-registered wet-lab endpoints (C4): primary success if either purchased compound shows IC50 ≤ 10 µM on any JNK isoform; JNK1 preference claimed only if SI_J2 ≥ 3 and SI_J3 ≥ 3.

---

## Results — Computational selectivity autopsy (RQ-C)

On the literature benchmark panel, Glide-derived Δsel_dock did not meet a 55% direction-accuracy usability threshold (archived VSW single-PDB accuracy 43% [3/7]; ensemble recomputation from the project docking table 25% [2/8]). Among key controls (SP600125, TCS JNK 6O, CC-930, E1), directional agreement was 50% (2/4). A Gly87 occupancy heuristic returned occ_JNK1 = True for all five tested benchmarks (ligand–Gly87 distances 0.59–1.18 Å) and did not separate JNK1-preferring profiles from pan or reverse isoform profiles. An ML selective-class model trained on sparse ChEMBL positives yielded test F1 = 0. By contrast, the family-activity ML gate recovered 9/9 literature benchmarks at p_family ≥ 6.0 but passed 95.3% of 10,000 external decoys, confirming high recall with low specificity and motivating docking/MD filtering rather than ML selectivity claims.

Accordingly, Δsel, Gly87, and ML selective labels were retained as negative controls and were **not** used as hard gates for purchasing 690 and 2157. Shortlisting prioritized pose-credible family-binder enrichment (MD overall pass, grade A) over MD-predicted JNK1 bias. The strongest MD-bias molecule (2231; grade C; overall MD fail) was documented in silico as an opportunity-cost comparator but not purchased.

Chemotype novelty audit (ECFP4) placed both purchased molecules far from E1, CC-90001, and SP600125 (maxTc ≈ 0.23 versus literature references; ≈ 0.27 versus the curated ChEMBL JNK pool). PAINS filters were negative for 690 and 2157. Multi-seed Vina redocking showed pose consensus (pairwise RMSD ≤ 2 Å for ≥66% of seed pairs) for 690 on all three isoforms and for 2157 on JNK1 and JNK3, with weaker consensus on JNK2 for 2157.

*(Insert IC50 Results here after C4 unblinding.)*

---

## Contribution check

| Item | OK? |
|------|-----|
| Core claim = predictor failure + family pipeline | Yes |
| No selective-hit discovery claim | Yes |
| Purchase decoupling stated | Yes |
| IC50 not over-claimed | N/A pending |
