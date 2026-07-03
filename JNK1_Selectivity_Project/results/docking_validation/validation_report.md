# JNK Docking Selectivity — Validation Report

## 1. Cognate redocking (RMSD)

| PDB_ID   | target   | ligand_name   |   glide_score |   rmsd_A | pass(<2A)   | pass_redock   |
|:---------|:---------|:--------------|--------------:|---------:|:------------|:--------------|
| 3ELJ     | JNK1     | GS7_3ELJ      |      -12.7939 |   0.6634 | True        | True          |
| 4L7F     | JNK1     | AX13587_4L7F  |      -12.2065 |   0.9181 | True        | True          |
| 3E7O     | JNK2     | 35F-3E7O      |      -11.266  |   0.2596 | True        | True          |
| 3TTI     | JNK3     | CC-930_3TTI   |      -12.8575 |   1.501  | True        | True          |
| 4WHZ     | JNK3     | 3NL_4WHZ      |      -10.0919 |   1.8828 | True        | True          |

**Pass rate:** 5/5 structures below 2.0 Å.

All five prepared structures reproduce cognate poses adequately; receptor grids are trustworthy for pose comparison.

## 2. Isoform score aggregation & Δsel (docking)

Formula: `score_isoform = mean(PDB scores)` for JNK1 (3ELJ, 4L7F) and JNK3 (3TTI, 4WHZ); single structure for JNK2 (3E7O).

`Δsel_dock = min(score_JNK2, score_JNK3) − score_JNK1` — **Δsel > 0 → computational JNK1 preference**.

| name       | expected_profile    |   score_JNK1 |   score_JNK2 |   score_JNK3 |   delta_sel_dock |
|:-----------|:--------------------|-------------:|-------------:|-------------:|-----------------:|
| AS602801   | pan-JNK             |      -5.9745 |       -5.233 |      -9.072  |          -3.0975 |
| CC-401     | unknown-isoform     |      -6.365  |      -11.267 |      -5.941  |          -4.902  |
| SP600125   | pan-JNK             |      -8.4395 |      -11.013 |      -7.5265 |          -2.5735 |
| TCS JNK 6O | JNK1-preferring     |      -6.2235 |       -6.626 |      -7.4045 |          -1.181  |
| CC-930     | JNK2/JNK3-biased    |      -7.155  |       -8.372 |     -12.051  |          -4.896  |
| JNK-IN-8   | JNK3-preferring     |      -6.035  |       -4.879 |      -5.6225 |           0.4125 |
| CC-90001   | pan-JNK             |      -6.65   |       -7.575 |      -6.8935 |          -0.925  |
| Q63        | JNK1/JNK3-over-JNK2 |      -7.8395 |       -8.872 |      -8.147  |          -1.0325 |
| E1         | JNK1-preferring     |      -9.3895 |       -5.926 |      -6.3445 |           3.045  |

**MM-GBSA benchmark validation:** Not available (no Prime MM-GBSA on the 9 benchmark compounds). Δsel_mmgbsa for benchmarks is **N/A**.

## 3. Direction reproduction vs experimental IC50

- Spearman(Δsel_dock, −ΔpIC50_sel) = **0.786** (p=0.0208) — sign-aligned so positive correlation = correct direction.
- Spearman(Δsel_dock, raw ΔpIC50_sel) = -0.786 (expected negative if directions consistent).
- Direction accuracy (all compounds with IC50): **22.2%**
- Direction accuracy (JNK1 / JNK23 / pan labels only): **28.6%**

### Key control compounds

| name       | expected_profile   |   delta_sel_dock |   delta_pIC50_sel | exp_dir_pIC50   | pred_dir_dock   | direction_match   |
|:-----------|:-------------------|-----------------:|------------------:|:----------------|:----------------|:------------------|
| SP600125   | pan-JNK            |          -2.5735 |         -0.352183 | JNK1            | JNK23           | False             |
| TCS JNK 6O | JNK1-preferring    |          -1.181  |         -0.550907 | JNK1            | JNK23           | False             |
| CC-930     | JNK2/JNK3-biased   |          -4.896  |          0.940232 | JNK23           | JNK23           | True              |
| E1         | JNK1-preferring    |           3.045  |         -0.84739  | JNK1            | JNK1            | True              |

## 4. Within-isoform activity ranking (docking vs pIC50)

- **JNK1:** Spearman ρ = -0.429
- **JNK2:** Spearman ρ = 0.190
- **JNK3:** Spearman ρ = 0.371

## 5. Limitations (do not ignore)

- JNK1/JNK2/JNK3 ATP pockets are **highly conserved**; Glide score differences are often **1–3 kcal/mol**, near the noise floor.
- XP docking scores rank **binding affinity**, not isoform selectivity; spurious Δsel can arise from protein prep differences across PDBs.
- MM-GBSA was **not** benchmarked here; VSW MM-GBSA selectivity gates are **uncalibrated** for isoform direction.
- Covalent inhibitor (JNK-IN-8) and multi-target profiles violate simple Δsel logic.

## 6. Verdict

**Docking selectivity direction is NOT reliable for isoform ranking** (|Spearman|=0.79, direction accuracy=29%; thresholds: |ρ|≥0.35, accuracy≥55%). **Recommendation:** treat VSW hits as **pan-JNK family** actives; selectivity must come from orthogonal assays (FP, CETSA, panel IC50).