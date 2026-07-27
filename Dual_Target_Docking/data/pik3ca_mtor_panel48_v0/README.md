# PIK3CA / mTOR docking panel freeze — `pik3ca_mtor_panel48_v0`

**Status:** molecule list frozen for protocol-transfer docking (not full ChEMBL scale).

| Item | Value |
|------|-------|
| Freeze ID | `pik3ca_mtor_panel48_v0` |
| N | **48** |
| Classes | dual 18 / A_only 14 / B_only 12 / neither 4 |
| Targets | PIK3CA (A) + mTOR (B) |
| Receptors (next step) | **4L23** + **4JT6** |
| Pose gold | **PI-103** = `PM48_01` / CHEMBL573339 / PDB X6K |
| Label rule | pChEMBL ≥ 6 active; measured &lt;6 weak; untested ≠ inactive |
| Panel CSV | [`tables/panel_v0_48.csv`](tables/panel_v0_48.csv) |
| Builder | [`../../scripts/build_pik3ca_mtor_panel48.py`](../../scripts/build_pik3ca_mtor_panel48.py) |

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

## What to do next (local docking)

1. LigPrep the 48 SMILES (same settings as EGFR/HER2 panel40) → 1 PDBQT each  
2. Prepare receptors 4L23 / 4JT6; boxes from PI-103 (X6K)  
3. Cognate QC: PI-103 dual-end RMSD &lt; 2 Å  
4. Dock all 48 × 2 ends with frozen Vina seed `20260727`  
5. RTMScore best-of-9 → `rtm_min` / `rtm_min_z` ablation (same arms as panel40)

## Explicitly not frozen here

- Full 2713 paired set  
- Exhaustiveness final value (follow EGFR/HER2 v0.1 choice)  
- Architecture Morphy labels (leave `unknown` unless literature-backed)
