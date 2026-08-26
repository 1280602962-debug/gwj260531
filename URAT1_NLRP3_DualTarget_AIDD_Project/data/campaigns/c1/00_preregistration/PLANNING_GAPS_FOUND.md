# C1 planning gaps found during local execution (2026-08-26)

## Blockers / environment
1. **No NVIDIA GPU** on this WSL host → L3 (~9,849 × exh32 × 9 modes) is blocked. Only L0–L2 CPU path is valid.
2. **gnina 1.3.2** present; `campaign_c1.yaml` requires **1.3.1**. Recorded; do not treat as silent pass. Prefer matching frozen-campaign binary before L3 Rank claims.
3. PR branch not present in older worktrees; used `git worktree` at `/home/gwj/work/c1-campaign` @ `5995f6fb`.

## Schema / script traps (already warned in LOCAL_C1 §12 — confirmed)
4. `run_gnina_batch.py` still reports **max CNNaffinity (= C1_P2max)**. Rank must use SDF parser (`parse_c1_sdf_readouts.py`).
5. `campaign_c1.yaml` is NOT a `--config` for gnina batch. Engine yaml: `docking_c1.yaml` / CPU override `docking_c1_cpu.yaml`.
6. PR ships C1 docs/configs but **no** `prepare_ligands_c1.py` — implemented locally before L2.

## Chemistry / numbering
7. Prepared `9DKB_receptor.pdbqt` renumbers paper **Arg477 → ARG A 476**. Distances computed from CIF Arg477 guanidinium coords; label kept as Arg477.
8. **GSK-3008348** at pH 7.4 is a **zwitterion** (COO− + piperidinium). Acid gate (COO− present) passes; net charge 0. Sensitivity: neutral amine + COO− microspecies later if Acid-track keeps GSK.
9. **benzbromarone / dotinurad** are phenols, not carboxylates. Forced-recovery list mixes acid chemotypes; Acid-track carboxylate rule must not be applied blindly to phenols.
10. **NP3-146** Dimorphite form deprotonates sulfonylurea NH (N−). Plausible; if self-dock RMSD fails, re-test neutral microspecies before declaring NLRP3 arm exploratory.

## Process
11. Do **not** start L3 without written auth after L2 gate.
12. Do **not** use `run_funnel_p2.sh` or overwrite `data/repurposing/p2/`.

## Critical: Arg477 ≤ 4 Å vs 9DKB crystal (discovered during L2b)

13. **9DKB crystal lesinurad (A1AIL) carboxylate O–Arg477 N min distance = 6.70 Å.**
    The preregistered gate `acid_arg477_le_A: 4.0` is **stricter than the crystal**.
    Crystal-start `--local_only` recovers RMSD 0.75 Å but Arg ≈ 7.4 Å → fails written gate.
14. Therefore free-dock “Arg fail” mixes two issues: (a) CNNscore selects ~14 Å poses;
    (b) the absolute 4 Å cutoff cannot be satisfied even by the crystal reference.
15. **Do not** silently change 4→8 after seeing results. Amend preregistration in writing
    (crystal-relative tolerance) before any Acid-track pass language.
16. NP3-146@7ALV **passes** (CNNscore-selected RMSD 0.67–0.82 Å) after fixing RMSD
    protonation matching — initial null RMSD was a measurement bug.
