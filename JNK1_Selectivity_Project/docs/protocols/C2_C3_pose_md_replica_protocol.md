# C2 / C3 Protocol — Pose consensus & MD replicas (purchased set)

**Scope:** HIT103871685 (**690**), HIT101201113 (**2157**); optional positive controls E1, CC-90001.  
**Receptors:** JNK1 `3ELJ`, JNK2 `3E7O`, JNK3 `3TTI` (project ensemble).  
**Goal:** Show purchase poses are not single-run flukes; quantify MD metric variance.  
**Claim limit:** These runs support **pose credibility**, not isoform selectivity.

---

## Software / license requirement

| Tool | Requirement |
|------|-------------|
| Glide / Maestro / Desmond (if used) | **Valid institutional or academic Schrödinger license** |
| Alternative docking | AutoDock Vina / Gnina **only if re-docked de novo**; do not relabel Glide poses |
| MD | Desmond (licensed) **or** OpenMM/GROMACS with documented force field |

Do **not** publish methods claiming Schrödinger modules without a compliant license. Prefer re-running C2 under a legal license or an open engine with full re-dock.

---

## C2 — Multi-seed / multi-start redock (≥3)

### Inputs
- Ligand SMILES from `data/shortlist/md_shortlist_final.csv`
- Prepared receptor grids used in the archived VS (same prep protocol)

### Procedure (Glide XP, if licensed)
1. For each ligand × each receptor (`3ELJ`, `3E7O`, `3TTI`):
   - Run XP docking with **≥3 independent seeds** (or ≥3 random starts / distinct ligand conformer expansions).
   - Keep top-1 pose per seed.
2. Align poses in the pocket; compute pairwise heavy-atom RMSD of the ligand.
3. Cluster at 2.0 Å; report whether the archived VS pose falls in the dominant cluster.
4. Record hinge H-bond presence (Met108/Met111/Met149 family hinge — use project residue mapping) for each seed pose.

### Outputs (write under `results/pose_consensus/`)
- `c2_pose_rmsd_matrix.csv`
- `c2_interaction_consistency.csv` (hinge HB yes/no per seed)
- `C2_POSE_CONSENSUS_REPORT.md`

### Pass interpretation (pose QC, not selectivity)
- Dominant cluster contains ≥2/3 seeds **and** archived pose RMSD ≤ 2.0 Å to cluster centroid → **pose stable to protocol noise**.
- Failure → do not over-interpret single archived pose figures.

### Open-source fallback (Vina/Gnina)
- Exhaustiveness ≥ 32; ≥3 independent runs; same RMSD clustering.
- State clearly that scores are **not** comparable to archived Glide XP ranks; use only for pose geometry consensus.

---

## C3 — MD replica mini-panel

### Systems
- 690 and 2157 × {JNK1, JNK2, JNK3} = **6 complexes**
- Seeds: **≥2** independent initial velocity seeds per complex (preferably 3)
- Length: **20–50 ns** production each (longer if affordable)
- **No ligand heavy-atom restraints** in production (restraints only for equilibration if needed)

### Metrics (match project MD QC language)
- Ligand RMSD (protein-aligned)
- Hinge hydrogen-bond occupancy (%)
- Optional: pocket RMSF; end-point visual inspection

### Outputs (`results/md_replicas/`)
- Per-replica time series + `c3_replica_summary.csv` with mean±SD hinge occ / RMSD
- `pass_md_overall` recomputed per replica using archived thresholds (RMSD ≤ 3 Å; hinge HB ≥ 30% — confirm exact numbers from MD shortlist README before locking)
- `C3_MD_REPLICA_REPORT.md`

### Interpretation
- If pass/fail flips across seeds → archived single-replica rank is **fragile**; soften pose claims.
- Do **not** treat replica-mean hinge asymmetry as selectivity proof (E1 / SP600125 already show MD ≠ enzyme SI).

---

## Status

| Item | Status |
|------|--------|
| Protocol written | DONE (this file) |
| C2 execution | **PENDING** — requires licensed docking environment |
| C3 execution | **PENDING** — requires MD compute + force-field setup |

When executed, point PaperSpine `source_map.md` to the new `results/pose_consensus/` and `results/md_replicas/` artifacts.
