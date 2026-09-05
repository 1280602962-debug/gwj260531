# C5 W1 gate audit — benzbromarone@9DKA (zero new docking)

## Prep QC

| Check | Result |
|---|---|
| Box `--center_*` | **`[107.167, 105.296, 107.628]`** 22³ — matches lock (all 3 seeds) |
| Receptor leftover R75 | **None** (0 R75 / 0 HETATM / 0 Br in `9DKA_receptor.pdbqt`) |
| Ligand microspecies | PDBQT REMARK **phenolate** `[O-]`; sum_q≈−1.0 |
| Crystal R75 centroid (gemmi / SDF) | **`[107.167, 105.297, 107.628]`** — matches box |
| RMSD reference | 9DKA chain A R75, bond-ordered template; symmetry via `GetBestRMS` |

## 9-mode RMSD summary

| seed | CNNscore Top-1 mode | Top-1 `pose_rmsd` | Top-1 **GetBestRMS** | Top-3 min GetBestRMS | best-of-9 GetBestRMS (mode) | CNN Top-1 ≤2? (GetBestRMS) |
|---:|---:|---:|---:|---:|---:|---|
| 42 | 1 | 3.585 | **1.652** | 0.484 | **0.484** (m3) | YES |
| 43 | 1 | 3.603 | **1.645** | 0.486 | **0.482** (m4) | YES |
| 44 | 1 | 3.595 | **1.649** | 0.477 | **0.477** (m3) | YES |

## Per-mode detail (GetBestRMS)

### seed 42

| mode | CNNscore | CNNaffinity | pose_rmsd | GetBestRMS |
|---:|---:|---:|---:|---:|
| 1 | 0.9661 | 7.672 | 3.585 | 1.652 |
| 2 | 0.9633 | 7.656 | 3.617 | 1.680 |
| 3 | 0.9495 | 7.603 | 1.110 | 0.484 |
| 4 | 0.9480 | 7.606 | 1.134 | 0.485 |
| 5 | 0.9240 | 7.298 | 6.361 | 0.780 |
| 6 | 0.9208 | 7.276 | 6.382 | 0.764 |
| 7 | 0.8626 | 7.373 | 1.528 | 0.773 |
| 8 | 0.8613 | 7.074 | 4.299 | 0.925 |
| 9 | 0.8377 | 7.071 | 4.291 | 0.938 |

### seed 43

| mode | CNNscore | CNNaffinity | pose_rmsd | GetBestRMS |
|---:|---:|---:|---:|---:|
| 1 | 0.9662 | 7.668 | 3.603 | 1.645 |
| 2 | 0.9635 | 7.661 | 3.612 | 1.677 |
| 3 | 0.9538 | 7.609 | 1.087 | 0.486 |
| 4 | 0.9486 | 7.603 | 1.131 | 0.482 |
| 5 | 0.9220 | 7.291 | 6.362 | 0.784 |
| 6 | 0.9201 | 7.273 | 6.380 | 0.752 |
| 7 | 0.8237 | 6.793 | 5.769 | 0.943 |
| 8 | 0.8152 | 7.079 | 4.245 | 0.980 |
| 9 | 0.7903 | 6.687 | 5.790 | 0.919 |

### seed 44

| mode | CNNscore | CNNaffinity | pose_rmsd | GetBestRMS |
|---:|---:|---:|---:|---:|
| 1 | 0.9659 | 7.672 | 3.595 | 1.649 |
| 2 | 0.9634 | 7.661 | nan | nan |
| 3 | 0.9499 | 7.603 | 1.109 | 0.477 |
| 4 | 0.9475 | 7.602 | 1.130 | 0.479 |
| 5 | 0.9228 | 7.282 | nan | nan |
| 6 | 0.9193 | 7.279 | 6.367 | 0.775 |
| 7 | 0.8943 | 7.483 | 1.681 | 0.852 |
| 8 | 0.8520 | 7.081 | 4.282 | 0.919 |
| 9 | 0.8461 | 7.071 | 4.303 | 0.916 |

## Interpretation (audit only; gate file still `pass=false` under `pose_rmsd`)

- Prep QC does **not** support Fork 1 (wrong box / leftover R75 / wrong protonation).
- Under **symmetry-corrected GetBestRMS**: CNNscore Top-1 ≈ **1.65 Å** (≤2.0) on all three seeds; best-of-9 ≈ **0.48 Å**.
- Under original gate metric **`pose_rmsd`** (substructure match without full symmetry): Top-1 ≈ **3.59 Å** (fail); best-of-9 ≈ **1.09–1.11 Å**.
- The ~3.6 Å figure is consistent with a **wrong atom mapping / missing ring symmetry** in `pose_rmsd`, not with a missing crystal basin.
- Gate runner used `pose_rmsd` → recorded `pass=false`. No Task2/3 started.

Artifacts: `gate_audit_9mode_rmsd_dual.csv`, `gate_audit_summary_dual.csv`.
