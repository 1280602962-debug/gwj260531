#!/usr/bin/env python3
"""Copy submission-facing manuscripts, tables, figures, and scripts into submission_pack/."""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACK = ROOT / "submission_pack"

DOCS = [
    "MANUSCRIPT_JCIM_EN.md",
    "MANUSCRIPT_JCIM_ZH.md",
    "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    "STATISTICAL_LOCK_V1.md",
    "FIGURE_PANEL_LOCK_V3.md",
    "SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md",
    "DATA_AND_SOFTWARE_AVAILABILITY_JCIM_EN.md",
    "DATA_AND_SOFTWARE_AVAILABILITY_JCIM_ZH.md",
    "REFERENCES_JCIM.md",
]

TABLES = [
    "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv",
    "data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv",
    "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv",
    "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv",
    "data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv",
    "data/jcim_novelty_v0/tables/incremental_information_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv",
    "data/jcim_novelty_v0/tables/ecfp4_docking_scaler_sensitivity_v1.csv",
    "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv",
    "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv",
    "data/jcim_multiseed_v0/tables/multiseed_auroc_aggregate_v2.csv",
    "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/fiveseed_summary_min_aggregate_v1.csv",
    "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv",
    "data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv",
    "data/jcim_novelty_v0/tables/DUALFOURCLASS_EVALUATION_CONTRACT_v1.json",
    "data/jcim_novelty_v0/tables/REVISION_CHECKSUM_MANIFEST_v1.csv",
    "data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv",
    "data/jcim_j0j1_v0/tables/j0_strict_label_supply.csv",
    "data/jcim_chembl_universe_v0/tables/track_b_panel_summary_v1.csv",
]

SCRIPTS = [
    "docs/assemble_manuscript_en.py",
    "docs/assemble_manuscript_zh.py",
    "scripts/primary/bootstrap_primary.py",
    "scripts/audit/audit_submission_five_rounds_v1.py",
    "scripts/audit/pack_submission_v1.py",
    "data/jcim_novelty_v0/scripts/validate_revision_v1.py",
    "data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py",
    "data/jcim_novelty_v0/scripts/build_master_results_table_v1.py",
    "data/jcim_novelty_v0/scripts/ecfp4_docking_scaler_sensitivity_v1.py",
    "data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py",
    "data/jcim_bench_v0/scripts/jcim_figure_style.py",
    "data/jcim_bench_v0/scripts/plot_jcim_si_composites_v1.py",
]


def copy_file(rel: str, dest_root: Path) -> None:
    src = ROOT / rel
    dst = dest_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    if PACK.exists():
        shutil.rmtree(PACK)
    (PACK / "manuscript").mkdir(parents=True)
    (PACK / "tables").mkdir()
    (PACK / "scripts").mkdir()
    (PACK / "figures").mkdir()

    for name in DOCS:
        src = ROOT / "docs" / name
        if src.exists():
            shutil.copy2(src, PACK / "manuscript" / name)
    for rel in TABLES:
        src = ROOT / rel
        dst = PACK / "tables" / Path(rel).name
        if src.exists():
            shutil.copy2(src, dst)
    for rel in SCRIPTS:
        src = ROOT / rel
        if src.exists():
            shutil.copy2(src, PACK / "scripts" / Path(rel).name)

    fig_src = ROOT / "figures" / "jcim_article"
    allowed_suffix = {".png", ".pdf", ".tif", ".tiff", ".md", ".json"}
    if fig_src.exists():
        for path in sorted(fig_src.iterdir()):
            if path.suffix.lower() not in allowed_suffix:
                continue
            if path.stem.startswith(("Fig", "TOC", "CAPTIONS", "plotted", "FIGURE_AUDIT", "FIGURE_SHA")):
                shutil.copy2(path, PACK / "figures" / path.name)

    readme = """# JCIM submission pack

This directory is the publication-facing slice of DualFourClass-Bench.
It is not a second copy of the docking pose workspaces.

## Use these files for submission

| Path | Role |
|---|---|
| `manuscript/MANUSCRIPT_JCIM_EN.md` | English manuscript (assemble from section drafts; do not hand-edit the assembled file as the source of truth) |
| `manuscript/MANUSCRIPT_JCIM_ZH.md` | Chinese working manuscript |
| `manuscript/SUPPORTING_INFORMATION_JCIM_EN_V1.md` | English SI Tables S1–S13 |
| `manuscript/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` | Chinese SI |
| `manuscript/STATISTICAL_LOCK_V1.md` | Table 2 / Table 3 estimand lock |
| `manuscript/FIGURE_PANEL_LOCK_V3.md` | Figure-to-CSV map |
| `manuscript/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` | Five-round numeric audit |
| `tables/` | Frozen CSVs cited by Tables 1–3 and S2–S11 |
| `figures/` | Regenerated main and SI figures from `figures/jcim_article/scripts/update_figures_pr32.py` |
| `scripts/` | Assemble, validate, bootstrap-lock reader, figure v3, audit, pack |

## Do not submit as primary evidence

- `plot_jcim_article_figures_v1.py` / `v2.py` (withdrawn-pair leftovers)
- `plot_jcim_si_composites_v1.py` S1–S3 / S9 / S10 (original-set archive; may still tick PIK3CA/PIK3CB)
- `data/pik3ca_pik3cb_panel_v0/` (withdrawn pair archive)
- `external_slice_summary_202608_contract_v1.csv` (legacy three-pair BindingDB snapshot)
- Pose workspaces and multi-GB score dumps (indexed in the repository, not copied here)
- A minted Zenodo DOI (not created)

## Rebuild

```bash
python3 docs/assemble_manuscript_en.py
python3 docs/assemble_manuscript_zh.py
python3 scripts/primary/bootstrap_primary.py
python3 data/jcim_novelty_v0/scripts/validate_revision_v1.py
python3 data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py --check
python3 scripts/audit/audit_submission_five_rounds_v1.py
python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root .
python3 figures/jcim_article/scripts/audit_figures_pr32.py
python3 scripts/audit/pack_submission_v1.py
```
"""
    (PACK / "README.md").write_text(readme, encoding="utf-8")
    inventory = ["# Inventory", ""]
    for folder in ("manuscript", "tables", "figures", "scripts"):
        inventory.append(f"## {folder}")
        for path in sorted((PACK / folder).iterdir()):
            inventory.append(f"- `{path.name}` ({path.stat().st_size} bytes)")
        inventory.append("")
    (PACK / "INVENTORY.md").write_text("\n".join(inventory), encoding="utf-8")
    print("wrote", PACK)


if __name__ == "__main__":
    main()
