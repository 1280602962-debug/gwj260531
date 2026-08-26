# Campaign C1 output root

Do not write here until a C1 job is running locally.
Do not copy or overwrite `data/repurposing/p2/`.

Plan: [`docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md`](../../../docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md)  
Locks: [`config/campaign_c1.yaml`](../../../config/campaign_c1.yaml)  
GNINA engine: [`config/docking_c1.yaml`](../../../config/docking_c1.yaml)

Layout (created by L0–L8, not in advance):

| Dir | When |
|-----|------|
| `00_preregistration/` | copy of yaml + timestamp at L0 |
| `01_ligand_prep/` | carboxylate PDBQT |
| `02_selfdock/` | lesinurad@9DKB, NP3-146@7ALV |
| `03_forced_recovery/` | textbook URAT1 drugs |
| `04_decoy_dock_9dkb/` | TrueDecoy + RandomDecoy |
| `05_metrics/pass_fail.json` | Rank-track gate; required before clinical docking |
| `06_nlrp3_benchmark/` | site-related set |
| `07_clinical_dock/` | only after `pass_fail.json` |
| `08_nomination/` | Acid / Rank shortlists |
| `09_md/` | after shortlist freeze |
