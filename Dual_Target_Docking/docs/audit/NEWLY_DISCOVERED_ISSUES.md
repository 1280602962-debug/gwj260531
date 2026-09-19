# NEWLY DISCOVERED ISSUES

Date: 2026-09-19  
These items were not treated as pre-known just because earlier freeze memos exist. Each was re-found from files or independent replay.

## NEWLY_DISCOVERED_ISSUES (non-empty)

1. **SI Table S6 prose contradicts its own table (P1)**  
   EN SI (and submission_pack copy) still says: “On the main panels only AChE/BChE has a 95% CI excluding 0; EGFR/HER2 includes 0.”  
   The same section’s table: EGFR/HER2 main 0.107 [0.006, 0.220], CI excludes 0 = **yes**.  
   Main-text Figure 4 already states both AChE and EGFR exclude 0.  
   Evidence: `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md` lines 207 vs 211; `results/canonical/matched_minus_mismatched.csv`.

2. **Manuscript cluster CIs are stale (P1)**  
   EN/ZH Results still typeset EGFR cluster intervals [0.235, 0.665] (scaffold) and [0.125, 0.644] (document).  
   Canonical: [0.2316, 0.6374] and [0.1172, 0.6356] → typeset [0.232, 0.637] and [0.117, 0.636].  
   SI Table S4 is already correct. Main text was not updated with the same substitution.  
   Qualitative “excludes 0” is unchanged.

3. **SI numeric-rule footnote still cites EGFR TPSA 0.4275 (P1)**  
   Current `descriptor_baselines.csv` EGFR TPSA D-vs-B / summary_min = 0.4389.  
   0.4275 is a leftover typeset note, not Table 2 (Table 2 is docking AUROC).

4. **Figure 5C lock still uses fuzzy `five_seed_comparable` (P1)**  
   Canonical five-seed CSV now has `same_protocol_as_primary` / `same_membership_as_primary`.  
   AChE: protocol=1, membership=0, but `comparable_to_current_primary=1` and plotted `five_seed_comparable=1`.  
   That is the exact ambiguity the writing-freeze asked to stop using.

5. **`GNINA_SOURCES` default still points at historical EGFR GNINA (P1)**  
   Runtime override exists. A direct import of the dict would load the wrong file.

6. **Root freeze-check docs still state 0.3237 as current Table 2 (P1 hygiene)**  
   superseded historical snapshots now under `docs/archive/`: `PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md`, `PR39_SCIENTIFIC_DATA_AUDIT.md`, `SUBMISSION_ONLY_CLEANUP_REPORT.md`.  
   Not manuscript body, but they sit in `docs/` without an archive/ prefix and can be mistaken for current instructions.

7. **Legacy EGFR box JSON files remain (P2 / P1 confusion risk)**  
   `3POZ_box.json` / `3RCD_box.json` were current-facing at audit time; they are now superseded historical snapshots in `data/egfr_her2_panel120_v0/boxes/archive/`. Current scripts use `*_box_corrected.json`.

8. **Duplicate identical score master (P2)**  
   `data/processed/current_score_master.csv` was an identical duplicate of `results/canonical/current_score_master.csv` at audit time (superseded). Current pointer: `data/processed/CURRENT_SCORE_MASTER.md`.

## Scans that did **not** find a new scientific defect

- Intra-pair duplicate ligand_id / SMILES / ChEMBL: 0  
- Primary θ=6.0 class mismatch: 0  
- Unresolved_missing_arm in primary: 0  
- Eligibility CSV performance tokens: 0  
- Negative `score_S` in primary: 0  
- Extreme scores outside 2–20: 0  
- Score-extraction sample mismatches: 0  
- Independent statistic replay discrepancies >1e-4: 0 (after audit MM pairing was corrected; first audit CI error was in the audit script, not in production)  
- `except: pass` writing primary AUROC: 0  
- Hard-coded 0.3237 in analysis scripts: 0  

## First-pass false P0 (documented so it is not reused)

An audit-only matched/mismatched CI that resampled matched and mismatched **separately** produced 16 CI mismatches and would have moved AChE’s interval across 0.  
Production uses one paired ligand draw per replicate. After the audit implementation used the same draw, all 16 CIs MATCH.  
That was an audit bug, not a project P0.
