# ACTIVE SCRIPT AUDIT

Date: 2026-09-19  
Scope: `scripts/` in Dual_Target_Docking. This is a scientific-writer audit, not lint.

## Classification

| Script | Class | Writes current scientific outputs? |
|--------|--------|------------------------------------|
| `analysis/analysis_config.py` | ACTIVE_SCIENCE | constants only |
| `analysis/build_current_score_master.py` | ACTIVE_SCIENCE | `current_score_master.csv` (canonical + processed) |
| `analysis/compute_canonical_results.py` | ACTIVE_SCIENCE | primary/sensitivity canonical CSVs |
| `analysis/fit_ecfp4_models.py` | ACTIVE_SCIENCE | ECFP4 CSVs |
| `analysis/compute_descriptor_baselines.py` | ACTIVE_SCIENCE | descriptor CSVs |
| `analysis/compute_class_chemistry.py` | ACTIVE_SCIENCE | `class_chemistry_summary.csv` |
| `analysis/compute_leave_one_document.py` | ACTIVE_SCIENCE | `leave_one_document_delta.csv` |
| `analysis/compute_detectable_effect.py` | ACTIVE_SCIENCE | `detectable_effect_simulation.csv` |
| `analysis/bootstrap_metrics.py` | ACTIVE_SCIENCE | library; no file writer |
| `analysis/rebuild_freeze.py` | ACTIVE_SCIENCE | writes a **new** freeze directory, not live canonical unless promoted |
| `analysis/promote_freeze_to_canonical.py` | ACTIVE_SCIENCE | can replace canonical |
| `analysis/rebuild_submission_pack.py` | ACTIVE_SCIENCE | submission_pack copy |
| `analysis/close_publication_from_canonical.py` | HELPER | publication close |
| `analysis/patch_publication_text.py` | HELPER | text patcher; must not invent numbers |
| `analysis/adjudicate_activity_records.py` | HELPER / HISTORICAL | adjudication utility |
| `check_analysis_env.py` | HELPER | env gate |
| `check_docking_env.py` | HELPER | docking env |
| `check_local_env.py` | HELPER | local env |
| `qa/verify_freeze_rebuild.py` | HELPER | compare/verify |
| `qa/check_current_chain.py` | HELPER | current-chain gates |
| `qa/audit_current_paths.py` | HELPER | path audit |
| `qa/compare_rebuild_to_canonical.py` | HELPER | compare |
| `qa/build_receptor_input_registry.py` | HELPER | provenance CSV |
| `qa/run_scientific_traceability_audit.py` | HELPER | older QA |
| `qa/run_final_forensic_audit.py` | HELPER | this forensic runner (QA/provenance only) |

No DEAD script in `scripts/analysis/` writes live canonical by default except the named ACTIVE_SCIENCE writers.

## Unique writer rule

`current_score_master.csv` is written by `build_current_score_master.py` to both `results/canonical/` and `data/processed/`. The two copies are currently byte-identical. That is a **duplicated truth source** (P1 hygiene), not two different numbers.

Other canonical CSVs have a single production writer: `compute_canonical_results.py`, plus the specialized ECFP/descriptor/detectable/leave-one-document/class-chemistry scripts.

`promote_freeze_to_canonical.py` can overwrite the same files after a rebuild. That is intended, but it means two promotion paths exist.

## Hard-coded scientific values

`analysis_config.py` holds the real scientific constants (θ, seeds, B, pairs, receptors). That is correct.

`GNINA_SOURCES` still lists the **historical** EGFR GNINA path. Runtime `current_gnina_sources()` overrides when the uniform file exists. A caller that imports `GNINA_SOURCES` directly would use the wrong file. P1.

No analysis script hard-codes 0.3237 / 0.3341 / 0.0112 as results.

## Silent exception handling

`except Exception` then fallback in:

- `compute_canonical_results.py` scaffold helper → `"unparsed"`
- `compute_descriptor_baselines.py` / `fit_ecfp4_models.py` Murcko → `"acyclic"` / `"unparsed"`

These can hide SMILES parse failures by collapsing them to a dummy scaffold. They do not silently drop docking scores. P2 / document.

No `except: pass` that writes primary AUROC was found.

## Other findings

- `comparable_to_current_primary` is still written as 1 whenever protocol matches, including AChE five-seed where membership differs. Separate flags exist; Figure 5C metadata still plots the fuzzy flag. P1.
- PIK3CA recovered prep script contains a historical absolute ROOT. Historical / unavailable as a current instruction.
- ENV_PIN contains personal absolute paths (`/home/gwj/...`). Snapshot only; not a current analysis instruction.
- Order dependence: GroupKFold is documented non-shuffled; bootstrap uses a fixed Generator seed. Primary AUROC is rank-based and order-insensitive.
- Empty cells are parsed with `parse_finite` → None, not 0. EGFR timeout EH40_31 is skipped, not zero-filled.
- Score sign is consistently `score_S = −energy` in master construction and five-seed loaders.

## Verdict

Active current outputs have named writers. The scientific risk is **wrong default paths / fuzzy flags / duplicate master copies**, not silent AUROC invention.
