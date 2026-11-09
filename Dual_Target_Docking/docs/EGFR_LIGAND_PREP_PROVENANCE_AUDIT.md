# EGFR/HER2 ligand-preparation provenance audit

Date: 2026-09-18  
Branch: `cursor/scientific-freeze-c7cc`  
Scope: what actually generated the ligand PDBQT used for current EGFR/HER2 production and five-seed docking.

This audit does **not** treat `protocol.yaml` text as proof that LigPrep produced the docked coordinates.

Decision: **CASE 2** (uniform RDKit ETKDGv3 + MMFF + Meeko rebuild of all 110 ligands).

---

## DIRECT_EVIDENCE

1. **A Schrödinger LigPrep job did run on this machine for panel_v0_40.**  
   Local files (not in git):
   - `/mnt/d/CADD paper exercise/dual target docking/Maestro doc/panel_v0_40_ligprep/panel_v0_40_ligprep.log`
   - `panel_v0_40_ligprep.sh`
   - `panel_v0_40_ligprep.inp`  
   The log records `Program: LigPrep`, `Version: 2025-1 build 129`, `StartTime: 2026-07-24-10:59:13`, command `D:\Schrodinger2025\ligprep.exe ...`, 14/14 subjobs succeeded, and `# of processed structures in "panel_v0_40_ligprep-out.maegz" : 40`.  
   This is evidence that LigPrep **executed** on 40 class-titled structures (`dual` / `A_only` / `B_only` / `neither`). It is **not** evidence that those coordinates were converted to PDBQT or used in the current production/five-seed Vina tables.

2. **Frozen EGFR/HER2 receptor PDBQT blobs still exist in git history.**  
   Recovered from `c213c485` (user-cited `dbe7ac1d` is not in this clone; blob IDs match the requested hashes):
   - `3POZ_receptor.pdbqt` blob `2b398c62a14395212913401055baf8b248e75f24`
   - `3RCD_receptor.pdbqt` blob `c338094b67d11a71304aa7f6c9a851f4967aabdd`  
   Restored to `data/egfr_her2_uniform_rdkit_v1/receptors/`.

3. **Track-B ligand prep code is recoverable and is RDKit ETKDGv3 + Meeko.**  
   Historical `prep_track_b_ligands_v1.py` (`c213c485`, originally `31aba656`): `MolFromSmiles` → largest fragment by heavy atoms → `AddHs` → `ETKDGv3` seed `20260727` → `MMFFOptimizeMolecule(maxIters=200)` → `mk_prepare_ligand.py`. Track-B YAML records `no_ligprep: true`. This applies to the five Track-B pairs, **not** to a recovered EGFR EH40 PDBQT set.

4. **Current production EGFR scores exist as numeric tables, not as ligand PDBQT.**  
   `data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv` (110 ligands) and `data/jcim_multiseed_v0/tables/scores_vina_mode1_EGFR_corrected_box_fiveseed.csv` (1100 rows). No `ligands_pdbqt/` directory exists on this branch for panel120 (0 files at `c213c485`; deleted earlier in slim cleanup).

5. **Corrected boxes exist and match the required heavy-atom construction.**  
   `data/egfr_her2_panel120_v0/boxes/3POZ_box_corrected.json` and `3RCD_box_corrected.json`.

---

## INDIRECT_METADATA

1. Historical `data/egfr_her2_panel120_v0/protocol/protocol.yaml` (git `be1fda60` / `c213c485`) claimed `panel40_reused: LigPrep` and `new_ligands: RDKit ETKDG + meeko`. **Not used as LigPrep proof.** Current yaml records deposited scores as `not_recoverable` and points the CASE 2 rebuild at `data/egfr_her2_uniform_rdkit_v1`.

2. Historical `data/egfr_her2_panel120_v0/MANIFEST.md` claimed mixed LigPrep/RDKit. Current MANIFEST states deposited ligand PDBQT are not recoverable.

3. Manuscript / SI methods (out of this audit’s rewrite scope) still describe a single RDKit ETKDGv3 + MMFF + Meeko pipeline for the eight-pair study. That claim is **not** proven for the deposited EGFR production scores.

4. PIK3CA/mTOR `panel48_rdkit_v0` still contains `ligprep_old` score columns and `prep_delta_vs_ligprep.csv`. Those are a **PIK3CA** historical comparison. M4-style renaming of old scores to `ligprep_old` is not independent EGFR LigPrep evidence (user constraint).

5. Local Maestro files also exist for `pik3ca_mtor_panel48_v0_ligprep/` (log/sh/inp). Same class of metadata: a LigPrep job ran for that panel; not EGFR production PDBQT provenance.

6. `Maestro doc/vina_docking/tmp/ligands.mae` and `ligands_fixed.mae` exist locally. They are not mapped to EH40 IDs in this audit and are not git-tracked production inputs.

---

## CONTRADICTORY_RECORD

1. **PR / Track-B text says “No LigPrep”.**  
   `DOCKING_PLAN_V1.md` / `LOCAL_RECOMPUTE_PACK_V1.md` (still in git history, deleted from the slim tree): “Prep is RDKit ETKDGv3 + meeko. No LigPrep.” That contract is for Track-B five ordinary pairs, not a recovered EGFR EH40 PDBQT chain.

2. **Historical EGFR protocol.yaml claimed LigPrep reuse of panel40** while the eight-pair methods claim a uniform RDKit/Meeko ligand prep. Current yaml no longer presents LigPrep or uniform RDKit as the method of the deposited ablation scores.

3. **User confirmation: LigPrep was never used as the intended production ligand prep.** That conflicts with the LigPrep job log (a job did run) and with protocol.yaml, and cannot be resolved by choosing one narrative.

4. **Output `panel_v0_40_ligprep-out.maegz` is named in the log and in historical protocol text, but is absent** from the Maestro job directory (only `.log` / `.sh` / `.inp` remain) and was **never tracked in git** (`git log --all -- '**/*.maegz'` empty).

5. No script in the current tree reads a `.maegz` or calls `ligprep.exe` to produce the ablation scores.

---

## NOT_RECOVERABLE

1. `panel_v0_40_ligprep-out.maegz` (40 LigPrep output structures).
2. The ligand PDBQT files actually docked for current EGFR production and corrected-box five-seed Vina.
3. A file-level map from those PDBQT files to either LigPrep output or RDKit/Meeko output.
4. `dbe7ac1d` as a commit in this clone (receptor blobs recovered from `c213c485` instead).
5. Original `mk_prepare_ligand.py` logs for EH40 / EH120 production.
6. Proof that the 70 “new” EH120 ligands and the 40 EH40 ligands currently in `ablation_ligand_scores.csv` share one preparation pipeline.

---

## What CASE 1 would have required

Direct, traceable evidence that EH40 production PDBQT were RDKit ETKDGv3 + Meeko **and** that current production/five-seed used those files.

That evidence is absent. The LigPrep log proves a job ran; it does not prove current docking inputs. Protocol.yaml is metadata only.

## Decision

**CASE 2.** Do not guess. Rebuild all 110 ligands from `panel_v0_120.csv` canonical SMILES with the Track-B RDKit/Meeko protocol, restore the frozen receptor PDBQT, keep the corrected boxes, and redock five Vina seeds plus EGFR independent GNINA.

Historical M4 “LigPrep vs RDKit” comparisons are not current scientific evidence for EGFR/HER2.
