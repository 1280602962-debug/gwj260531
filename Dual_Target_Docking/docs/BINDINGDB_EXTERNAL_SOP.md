# BindingDB external-validation SOP (local; only if time-split fails)

Do not start this until `data/jcim_novelty_v0/analysis/TIME_SPLIT_VERDICT.md` says the primary cutoff cannot be packaged as external validation.

## Freeze before seeing AUROC

1. Target accessions (UniProt / BindingDB): EGFR P00533, HER2 P04626, AChE P22303, BChE P06276, PIK3CA P42336, PIK3CB P42338, mTOR P42345.
2. Keep exact-relation quantitative Ki, Kd, and IC50 separately; do not mix them into one pChEMBL-like scale until the rule is written down.
3. Deduplicate against ChEMBL panel structures (salt, stereo, canonical SMILES).
4. Drop BindingDB records whose source document is already a ChEMBL document used in the panels.
5. Independent means: **structure not repeated, literature not repeated, database independent, unused for method selection**.
6. Evaluation contract: same θ = 6.0 four-state labels, same pocket-matched directional AUROCs, same `summary_min`. No cutoff shopping.

## Stop and keep the formulation-audit claim if

- fewer than two pairs remain
- any primary class stays below 10
- most BindingDB records cite the same papers as ChEMBL
- a favorable AUROC appears only after changing the threshold
- the test set influenced receptor or metric choice

New ligands then still need local docking. Cloud work stops at the protocol.
