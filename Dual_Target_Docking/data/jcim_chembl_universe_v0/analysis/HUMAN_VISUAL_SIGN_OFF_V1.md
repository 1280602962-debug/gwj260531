# Layer 2 visual / literature sign-off — complete for the ten new PDBs

Signed 2026-09-04 against RCSB, UniProt, and primary structure papers (`LAYER2_LITERATURE_SIGN_OFF_V1.md`). Log: `tables/site_verification_log_v1.csv` — every new end is `PASS`.

This page is no longer a to-do list. It records **what was confirmed** and **what still cannot start**.

## Outcome

**All ten PDBs remain frozen (`RECEPTOR_FREEZE_V1.md`).** None will be swapped. CTSK/CTSS crystals are **structure records only** — they are out of the unified Vina campaign (`DOCKING_PLAN_V1.md`). The covalent ligand-prep note (`COVALENT_LIGAND_PREP_V1.md`) stays on file in case a later written amendment opens a separate covalent arm; it is not a Track B job. Numbering must keep construct / resolved / UniProt columns separate (`tables/receptor_span_registry_v1.csv`).

| # | Protein (UniProt) | PDB / CCD | Layer-2 | Receptor freeze | Ligand for Vina |
|---|-------------------|-----------|---------|-----------------|-----------------|
| 1 | Thrombin F2 **P00734** | **4UDW / N6L** | PASS | keep | N6L as deposited (S1, chain H, PDB 364–621) |
| 2 | Factor Xa F10 **P00742** | **2JKH / BI7** | PASS | keep | BI7 as deposited (S1/S4, heavy 235–475) |
| 3 | JAK1 **P23458** | **6N7A / KEV** | PASS | keep | KEV in JH1 ATP; construct ~854–1154 ≠ UniProt JH1 875–1153 |
| 4 | TYK2 **P29597** | **3LXP / IZA** | PASS | keep | IZA=CMP-6 in JH1 ATP; construct 888–1182 / resolved ~1178 / JH1 897–1176 |
| 5 | JAK2 **O60674** | **8BXH / C87** | PASS | keep | C87=momelotinib in JH1 ATP |
| 6 | PPARG **P37231** | **9V8H / BRL** | PASS | keep | BRL=rosiglitazone in LBD; **note PG08-NL peptide** (ternary, not binary) |
| 7 | PPARA **Q07869** | **6LXA / EPA** | PASS | keep | EPA in LBD 200–468 |
| 8 | PPARD **Q03181** | **5U3Q / 7UJ** | PASS | keep | 7UJ=specific agonist 1 in LBD; not PEG |
| 9 | Cathepsin K **P43235** | **4X6H / I37** (not 3XT) | PASS | keep **4X6H** | **I37 pre-reaction nitrile**; do not cut 3XT |
| 10 | Cathepsin S **P25774** | **9GJ2 / KH0** | PASS | keep **9GJ2** | **reconstruct pre-reaction α-ketoamide 13b**; do not cut KH0; no journal paper yet |

Frozen pairs (already docked): 4L23/X6K, 4JT6/X6K, 4EY7/E20, 4BDS/THA remain as before.

## Still blocked after this PASS

1. **Layer 3** cognate best-of-9 RMSD on the **eight ordinary new receptors** — needs local Vina.
2. **Production docking** of the **five** ordinary pairs in `DOCKING_PLAN_V1.md`.
3. Do **not** Vina-dock CTSK/CTSS (covalent). Do **not** re-dock PIK3CA/PIK3CB. Do **not** dock CREBBP/BRD4, GPCRs, or SLC6.

Next compute that can run in-cloud without Vina: extract four-state ligand panels for the five ordinary pairs from the local ChEMBL 37 SQLite and RDKit/meeko-prep. Vina itself stays local.
