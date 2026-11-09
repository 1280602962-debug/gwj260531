# REPRODUCIBILITY BOUNDARY

Date: 2026-09-19

## A. Zero-dock analysis reproducibility

**In scope:** rebuild statistics, tables, and (with the figure scripts) plotted-value JSON from frozen score/activity inputs.

Required:

- Python 3.12.3
- `requirements-analysis.txt` pins: numpy 2.5.2, pandas 3.0.5, scipy 1.18.1, scikit-learn 1.9.0, RDKit 2026.3.5
- Frozen `current_score_master` inputs: panel CSVs + deposited score tables + EGFR uniform five-seed CSV
- Scripts: `check_analysis_env.py` → `build_current_score_master.py` → `compute_canonical_results.py` (+ ECFP/descriptor/detectable) → `verify_freeze_rebuild.py`

This audit independently recomputed primary estimands from the deposited master without calling `bootstrap_metrics` / `compute_canonical_results`. 192/192 comparisons MATCH.

A full clean-temp `rebuild_freeze` was **not** re-executed in this read-only pass. The prior writing-freeze recorded PASS. That prior PASS is not re-used as truth; the independent replay is the evidence for Table 2 arithmetic.

Analysis RDKit 2026.03.5 is **not** the historical docking RDKit (ENV_PIN: 2026.3.1 on the docking machine). Do not back-date.

Python 3.12.13 has **no evidence** in this tree.

## B. Complete docking reproducibility

| Pack | SMILES → PDBQT → dock? | Evidence |
|------|------------------------|----------|
| EGFR/HER2 current uniform | **yes, in principle** | `prep_uniform_ligands.py`, 110 PDBQT, receptors, corrected boxes, Vina/GNINA scripts |
| PIK3CA/mTOR | **no bit-identical rerun** | script recovered; ligand PDBQT absent |
| AChE/BChE | **no** | script recovered; ligand PDBQT absent |
| Track-B five pairs | **no from ligand PDBQT** | prep script + status table; PDBQT never in git; **poses** deposited |
| Receptor rebuild | **forbidden / not possible from protocol** | method NOT_RECOVERABLE |

Do **not** claim that all eight pairs are fully reproducible from SMILES.

EGFR uniform is a **new** CASE 2 campaign. It is the current EGFR truth. It is not a reconstruction of the deleted historical LigPrep/ablation PDBQT.

## Environment registers (do not borrow versions)

| Register | What may be cited | What may not |
|----------|-------------------|--------------|
| ANALYSIS ENVIRONMENT | Python 3.12.3 + `requirements-analysis.txt` | docking-machine RDKit 2026.3.1 |
| DOCKING ENVIRONMENT | Vina 1.2.7 (yaml / ENV_PIN); E and seed from protocol | analysis Python |
| LIGAND PREPARATION | EGFR: Meeko 0.7.1 + current prep script. Others: METADATA_ONLY | “all pairs used the analysis RDKit pin” |
| GNINA ENVIRONMENT | 1.3.2, `--no_gpu` (ENV_PIN + independent-dock yaml) | Vina version |
| RTMSCORE ENVIRONMENT | `rtmscore_model1` path in ENV_PIN | analysis sklearn |

## What this repository can defend

- Table 2 / Table 3 / SI numeric tables that come from frozen scores + the analysis freeze
- EGFR uniform Vina/GNINA as a deposited campaign
- Track-B pose archive as deposited poses, not as input PDBQT

## What it cannot defend

- Bit-identical re-preparation of 14 receptors
- Bit-identical re-docking of PIK3CA / AChE / Track-B from SMILES
- Historical EGFR LigPrep as current input
- External docking validation (none was run)
