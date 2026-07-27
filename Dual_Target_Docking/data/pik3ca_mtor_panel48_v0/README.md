# pik3ca_mtor_panel48_v0 freeze

Local dual-target docking + RTMScore ablation pack for PIK3CA (4L23) / mTOR (4JT6) panel48.

- Protocol: exhaustiveness=16, seed=20260727, n_modes=9
- Cognate QC (PI-103 / X6K): Go @ E=16 — see `analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md`
- Vina: 96/96 jobs (4JT6/PM48_34 has 8 valid modes)
- RTMScore model1 best-of-9 + ablation tables under `tables/`
- Rebuild: `scripts/run_rtm_and_ablation.py`
- Primary ablation arm: `rtm_min_z` (AUROC Dual vs rest ≈ 0.685)

## Class roster (named anchors)

**Dual (18)** — must dock both ends as positives  
**PI-103** (`PM48_01`), Omipalisib, Gedatolisib, Dactolisib, Pictilisib, Apitolisib, Torin2, PF-04691502, Bimiralisib, Torin1, Sapanisertib, Buparlisib, Vistusertib, VS-5584, Paxalisib, Voxtalisib, Samotolisib, ZSTK-474.

**A_only (14)** — PIK3CA-strong / mTOR-weak hard negatives  
Includes **Taselisib**, **Alpelisib**, AZD-6482, Sonolisib + ChEMBL fillers.

**B_only (12)** — mTOR-strong / PIK3CA-weak hard negatives  
Includes **AZD-8055**, **Ku-0063794**, **WYE-132**, **OSI-027** + fingerprint-diverse fillers (WYE-analog pileup avoided).

**Neither (4)** — both ends measured weak.

## Selection rules used

- Paired ChEMBL only (`mols_PIK3CA.json` ∩ `mols_MTOR.json`)
- MW 180–750; heavy atoms ≤ 55; small molecule only
- Exclude rapalog / PROTAC name hits
- ≤2 molecules per Murcko scaffold per class
- Non-seed fillers: Morgan Tanimoto &lt; 0.55 vs already-selected in-class
- PROTACs / rapalogs not in panel

## Docs

完整操作单见：[`../../docs/PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md`](../../docs/PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md)

## Explicitly not frozen here

- Full 2713 paired set  
- Architecture Morphy labels (leave `unknown` unless literature-backed)
