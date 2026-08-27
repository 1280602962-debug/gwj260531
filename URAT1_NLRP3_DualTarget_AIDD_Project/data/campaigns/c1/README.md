# Campaign C1 output root

**Active path (2026-08-26 Amendment A1):** Acid track OPEN · Rank track CLOSED · L3 full-decoy forbidden.

Locks: `config/campaign_c1.yaml`  
Amendment: `00_preregistration/AMENDMENT_A1_Arg477_crystal_relative.yaml`  
Gate file: `05_metrics/pass_fail.json`

Do not overwrite `data/repurposing/p2/`.
Claim language: **acid-pose dual-node hypotheses** (not docking-rank activity retrieval).

| Dir | Status |
|-----|--------|
| `00_preregistration/` | A1 accepted |
| `01_ligand_prep/` | forced recovery + acid clinical prep |
| `02_selfdock/` | L2 free + L2b constrained done |
| `05_metrics/pass_fail.json` | Rank FAIL / Acid OPEN |
| `07_clinical_dock/acid_pool/` | acid clinical lists |
| `07_clinical_dock/acid_dual_a1_frozen/` | A1 exploratory 24 dual keep (frozen) |
| `07_clinical_dock/acid_dual_a2/` | A2 geometry-first rescoring / multiseed |
| `05_metrics/acid_gate_retrospective_benchmark/` | A1 vs A2 OR on acid actives vs decoys |
| `08_nomination/` | A2 competition shortlist; MD not authorized |
