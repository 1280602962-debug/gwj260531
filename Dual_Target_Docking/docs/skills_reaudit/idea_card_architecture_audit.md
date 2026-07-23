# Architecture-Stratified Dual Docking Audit

## Motivation

Dual-target docking almost always docks one ligand twice independently, then aggregates scores. This treats a single connected molecule as two unrelated ligands. That assumption may hold for compact fused/merged chemotypes but fails for linked hybrids whose viable poses are coupled by the linker. Score fusion cannot repair missing coupled-pose information. The field lacks an architecture-stratified diagnostic that isolates this failure, and a protocol-level consistency gate before dual ranking.

## Method

1. Curate architecture labels (fused/merged/linked) for paired dual-activity molecules and Tier-A both-end co-crystals.
2. Freeze an independent dual-docking protocol (GNINA/Vina + PoseBusters).
3. Report pose (where available), enrichment, and dual-vs-single metrics **stratified by architecture**; include A-only/B-only hard negatives (measured weak, never untested-as-inactive).
4. Define connectivity/ensemble consistency checks (linker strain, pocket clash, shared-graph pose plausibility) as **gates**, not as a new generative sampler.
5. Primary evidence: independent docking + fusion looks acceptable on fused/merged but collapses on linked; the consistency gate restores dual-vs-single ranking; kill-switch = gate off returns to failure.
6. Optional holdout: laboratory fused vs linked series, ranking concordance with dual biochemical labels (cell only as secondary holdout, not binding gold).
