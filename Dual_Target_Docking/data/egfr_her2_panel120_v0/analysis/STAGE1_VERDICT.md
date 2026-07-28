# STAGE1 VERDICT — EGFR/HER2 panel120 expand (S1 gate)

- N = **110** (dual=28)
- Protocol: 3POZ/3RCD, E=8, seed=20260727, n_modes=9, RTM best-of-9
- AUROC vina_mean = **0.551**; rtm_min_z = **0.512**
- ΔAUROC(rtm_min_z − vina_mean) = **-0.039** 95% CI [-0.166, +0.085] (B=2000)
- Top10 hardneg: vina=9 (A=5, B=4); rtm=7 (A=5, B=2)

## Verdict: **No-Go (S1)**

ΔAUROC CI includes 0 at N≈110. Do **not** stack complex methods yet. Recommended downgrade: keep DualFourClass-Bench + failure typology as the publishable claim; optional further EGFR expansion only if effect estimate remains large; otherwise accept diagnosis-paper route (JCIM gap scenario 甲).

### Notes
- Panel40 poses reused; new ligands RDKit+meeko (documented).
- Architecture not used as selection filter; flags not gated into score.
- No clash/shortfall retune for this gate.
