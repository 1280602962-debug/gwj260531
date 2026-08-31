# Tooling notes — exhaustiveness_sensitivity_v1

## RTMScore
**RTM unavailable for this sensitivity rerun.**

Attempted to rescore Experiment A poses with:
- env: `conda activate rtmscore`
- script: `/home/gwj/software/RTMScore/example/rtmscore.py`
- model: `trained_models/rtmscore_model1.pth`
- pocket: `receptors/{3POZ,3RCD}_pocket_10.0.pdb`

Failure (first job `E8 3POZ`):
```
ValueError: not enough values to unpack (expected 2, got 0)
```
in `VSDataset` after ligand graph construction — all converted poses filtered as `None` (likely pdbqt→SDF conversion / sanitization mismatch for the batch SDF produced in this rerun).

Therefore:
- No `scores_rtm_experimentA.csv` was written.
- Primary conclusions below rely on **Vina sampling metrics** (RMSD / mode stability / Vina ranks).
- Do **not** invent RTM ranks for this sensitivity package.
- Historical panel40 RTM ranks remain archived under `ROOT/tables/scores_rtm*.csv` and are not re-interpreted under the new fixed seed.

## Vina
- Engine: `/home/gwj/miniconda3/bin/vina` 1.2.7
- Inputs locked via `tables/ligand_input_manifest.csv` (original LigPrep→meeko pdbqt; no re-LigPrep)
