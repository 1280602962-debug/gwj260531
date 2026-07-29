# AGENT COMMAND — JCIM Strengthening (Wave 1 first)

Authorization follows [`JCIM_STRENGTHENING_PLAN_V1.md`](JCIM_STRENGTHENING_PLAN_V1.md).  
Claim ceiling: [`../data/jcim_bench_v0/CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md).

## Do now (Wave 1, cloud, zero docking)

Create pack `data/jcim_strengthen_w1_v0/` with scripts + tables:

1. **W1.1** — Resume/complete ChEMBL fetch for 22 missing targets from `jcim_j0j1_v0`; refresh strict supply audit.  
2. **W1.2** — Aggregation sensitivity (pooled / mean / harmonic / summary_min) on assembled K=4.  
3. **W1.3** — Murcko/cluster bootstrap CIs for summary_min (vina + best baseline).  
4. **W1.4** — Scaffold-stratified metrics + AChE scaffold overlap / TPSA-within-scaffold.  
5. **W1.5** — Covariate-controlled comparison (bin by heavy or TPSA; residualize).  
6. **W1.6** — Label noise ceiling for all K=4 pairs.  
7. **W1.7** — Dock-fail MNAR table.  
8. **W1.8** — Failure taxonomy table template filled for 4 pairs (cases from existing ligands).  
9. **W1.9** — ECFP4 + logistic/RF ligand baseline with nested CV on exploration panels only (no holdout peeking).  
10. **W1.10** — Write `data/protocols/PROTOCOL_BENCH_V1.md`.

Also write `data/jcim_strengthen_w1_v0/analysis/WAVE1_VERDICT.md` mapping each reviewer critique → table/figure.

## Do NOT

- Expand EGFR docking  
- Retune primary score / Track B competition  
- Mix LigPrep into primary tables  
- Start Wave 2 docking from cloud unless user explicitly authorizes local docking ops

## After Wave 1

Stop for user decision on Wave 2 docking budget (PM expand ± structure ± single-target enrichment).
