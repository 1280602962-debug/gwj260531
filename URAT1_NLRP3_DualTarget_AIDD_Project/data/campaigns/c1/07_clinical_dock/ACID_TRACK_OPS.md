# C1 Acid-track ops (Amendment A1)

## Locked
- Arg477 ≤ **7.7027 Å** (= 9DKB crystal min 6.7027 + 1.0)
- Rank track **closed**; no `unique_docking_pool` L3
- No percentile ranking; claim = **acid-pose dual-node hypotheses**
- Exclude Vecabrutinib from URAT1 acid claim

## Done
- [x] A1 written (`00_preregistration/AMENDMENT_A1_*.yaml`)
- [x] `05_metrics/pass_fail.json` (Rank FAIL / Acid OPEN)
- [x] Acid pool: 303 acid-equivalents; 156 chemistry soft-pass
- [x] PDBQT prep for 156 chemistry soft-pass (pH 7.4)

## Next (needs GPU preferred)
1. Dock `01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv`
   - 9DKB + 7ALV, `config/docking_c1.yaml` (or `_cpu.yaml` only for smoke)
   - `num_modes: 9`, parse SDF (CNNscore-selected pose)
2. Keep if: carboxylate/equivalent · Arg477 ≤ 7.7027 · both pockets OK · chemistry soft
3. Nominate ≤2 primary + ≤3 backup → `08_nomination/`
4. MD only after shortlist freeze

## Do not
- Reopen Rank / L3 without new preregistration
- Overwrite `data/repurposing/p2/`
- Call results dual inhibitors or rank-prioritized hits
