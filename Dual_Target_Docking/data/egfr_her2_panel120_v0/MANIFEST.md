# MANIFEST — egfr_her2_panel120_v0 (Stage-1 expand)

- theta_pchembl: **6.0**
- N: **110**
- class counts: {'dual': 28, 'A_only': 38, 'B_only': 32, 'neither': 12}
- Murcko max per (class,scaffold): **5** (observed max=5)
- retained panel40: 40
- new ligands: 70
- docking protocol: 3POZ/3RCD, E=8, seed=20260727, n_modes=9, RTM best-of-9
- ligand prep: panel40 LigPrep reused; new ligands RDKit ETKDG + meeko (documented)
- architecture: not used as selection filter
- warning flags: diagnostic only (not gated into score)

## Stage-1 gate (S1)

- AUROC vina_mean = 0.551; rtm_min_z = 0.512
- ΔAUROC(rtm_min_z − vina_mean) = −0.039; 95% CI [−0.166, +0.085] (B=2000)
- **Verdict: No-Go (S1)** — CI includes 0; see `analysis/STAGE1_VERDICT.md`
