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

完整操作单见：[`../../docs/PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md`](../../docs/PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md)

**顺序硬约束：先 PI-103（PM48_01）双端 cognate QC，Go 后再 48×2。**

给定本机输入：
- LigPrep：`...\Maestro doc\pik3ca_mtor_panel48_v0_ligprep\pik3ca_mtor_panel48_v0_ligprep-out.maegz`
- 蛋白：`4L23_PIK3CA_prepared.pdb` · `4JT6_mTOR_prepared.pdb`
- 协议：Vina 1.2.7，`seed=20260727`，`exhaustiveness=8`，`n_modes=9`
- 盒子：共晶 **X6K** AABB + 5 Å（min edge 20 Å）
- QC：两端重原子 `best_of_9 RMSD < 2 Å`

## Explicitly not frozen here

- Full 2713 paired set  
- Exhaustiveness final value (follow EGFR/HER2 v0.1 choice)  
- Architecture Morphy labels (leave `unknown` unless literature-backed)
