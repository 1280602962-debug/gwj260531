# Issue register — writing freeze (2026-09-17)

Authoritative numbers live in `results/canonical/`. This register records confirmed data/calculation issues that were fixed, unresolved items that remain in the published claim set, and items that were not treated as inactive or fabricated.

Do not restore pre-adjudication AUROCs to keep old Table 2 values.

## Closed (fixed in this freeze)

| ID | Item | Decision | Effect on primary analysis |
|---|---|---|---|
| EV1 | Pictilisib PM48_05 (CHEMBL521851) PIK3CA 10.00 / mTOR 8.52 | Excluded the source-adjudicated rows; remaining qualified max 9.12 / 7.32 | Still dual; pair n_scored unchanged |
| EV2 | Alpelisib PM48_22 (CHEMBL2396661) mTOR 5.83 cellular | Excluded the cellular A-type row; biochemical 5.52 retained | Still A-only |
| EV3 | EH120_059 (CHEMBL305194) HER2 4.31 cellular | Excluded; HER2 arm unresolved. Missing opposite arm is not inactive | Dropped from activity-eligible set; EGFR A-only 38→37; D vs A 0.6607→0.6564 |
| EV4 | EH40_05 HER2 7.5; EH120_066 HER2 5.58 | Applied committed `machine_human_exclude` rows | No class change on remaining qualified records |
| EV5 | AB_087 BChE IC50/Ki | Applied `machine_human_exclude`; AChE arm remains; BChE unresolved | Unresolved missing arm, not neither-as-inactive; eligible neither 15→14; primary D/A/B 27/26/28 |
| EV6 | PM48_04 PIK3CA/mTOR source rows | Applied source adjudication | Still dual |
| EV7 | AChE max/median denominator mixed dump-missing AB_056 into 1/96 | Max/median uses the same qualified intersection (94 dump-gated ligands); AB_056 stays in primary as `panel_pchembl_no_audit_rows` | Flips 1/94 (CHEMBL659), not 1/96 |
| EV8 | EGFR max/median 6/110 | Same qualified intersection after EH120_059 exclusion | Flips 5/109 |
| EV9 | JAK1/TYK2 document-cluster CI copied from an old JAK mapping | Not recomputed; map and ChEMBL 37 sqlite unavailable | Table S4 / Fig S4 omit that interval; do not copy [−0.034, 0.682] |
| EV10 | Table S2 described as chain/altLoc/residue-template | Table S2 is PDB / cognate / resolution / box | Methods and SI wording updated |
| EV11 | 4JT6 treated as if it met the ≤3.5 Å census gate | Census gate remains ≤3.5 Å; selected mTOR 4JT6 is 3.60 Å | Stated as supply vs selected receptor |
| EV12 | Figure-lock EGFR D vs A = 0.661 | Regenerated from eligible master | 0.656 [0.520, 0.791]; n = 28/37/32 |

## Layering (do not collapse)

- Raw: ChEMBL 37 / panel pChEMBL as deposited.
- Adjudicated: `data/processed/activity_adjudication/` after `adjudicate_activity_records.py`.
- Master: `results/canonical/current_score_master.csv` with `activity_eligible` and `complete_case`.
- Primary Table 2: `activity_eligible == 1` and complete both-end scores.

## Unresolved (limitations, not silent fixes)

| ID | Item | Why it remains | Publication handling |
|---|---|---|---|
| U1 | Not every max-determining ChEMBL record has a paper-level full-text audit | That would invent a new review, not freeze the committed high-impact exclusions | Stated as a limitation; remaining pairs use dump-gated max/median |
| U2 | JAK1/TYK2 document-cluster grouping map not deposited | Cannot rebuild groups from current scores | “not recomputed”; previous interval not used |
| U3 | EGFR/HER2 has no unused-pool holdout | Supply-limited pair | Figure 5 / Table S6 dagger; not filled |
| U4 | No pair passed the independent external docking gate | Eligibility screen, not external validation | Figure S5; 0/8 |
| U5 | Independent GNINA missing one EGFR neither (n=11 vs Vina 12) | Engine did not return both-end scores | Table S7 reports 28/37/32/11 |
| U6 | AB_056 has no dump assay rows | Kept in primary; excluded from dump-gated max/median | Denominator 94 vs eligible 95 |

## Explicit non-actions

- Did not recode missing opposite-arm ligands as inactive.
- Did not delete all EC50/functional assays.
- Did not change the research question, endpoints, or eight-pair set.
- Did not run new MD/FEP or production docking.
- Did not preserve old AUROCs by skipping the eligible filter.
