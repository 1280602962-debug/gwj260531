#!/usr/bin/env python3
"""QA-only scientific traceability tables. Zero-dock. No protocol retuning."""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    EGFR_UNIFORM_GNINA_REL,
    EGFR_UNIFORM_VINA_REL,
    FIVE_SEEDS,
    PRIMARY_PAIRS,
    RECEPTORS,
    current_egfr_score_source,
    egfr_uniform_ready,
    is_primary_row,
)
from analysis.compute_canonical_results import (  # noqa: E402
    _load_five_seed_scores,
    _smin_point,
    read_csv,
)

QA = ROOT / "results" / "qa"
FIELDS_MATRIX = [
    "pair",
    "channel",
    "role",
    "ligand_panel_source",
    "activity_source",
    "activity_source_kind",
    "label_rule",
    "ligand_smiles_source",
    "ligand_prep_method",
    "ligand_prep_version",
    "ligand_prep_seed",
    "ligand_prep_evidence",
    "receptor_pdb",
    "receptor_chain_or_construct",
    "receptor_input_file",
    "receptor_prep_method",
    "receptor_prep_evidence",
    "cognate_ligand",
    "box_file",
    "box_definition",
    "vina_version",
    "dock_seed",
    "exhaustiveness",
    "num_modes",
    "energy_range",
    "score_definition",
    "score_source_file",
    "sample_membership_rule",
    "analysis_script",
    "canonical_output",
    "figure_table_destination",
    "provenance_status",
    "notes",
]


def exists(rel: str) -> bool:
    return (ROOT / rel).is_file()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def five_seed_fixed_membership() -> list[dict]:
    master = read_csv(ROOT / "results/canonical/current_score_master.csv")
    current = {r["pair"]: r for r in read_csv(ROOT / "results/canonical/five_seed_summary_min.csv")}
    wide = _load_five_seed_scores()
    out = []
    for pair in PRIMARY_PAIRS:
        primary = [
            r
            for r in master
            if r["pair"] == pair and r.get("analysis_set") == "main" and is_primary_row(r)
        ]
        lig_ok = []
        for rec in primary:
            lig = rec["ligand_id"]
            if all("A" in wide.get((pair, seed, lig), {}) and "B" in wide.get((pair, seed, lig), {}) for seed in FIVE_SEEDS):
                lig_ok.append(rec)
        counts = defaultdict(int)
        for rec in lig_ok:
            counts[rec["primary_class_theta6"]] += 1
        n_inter = len(lig_ok)
        for seed in FIVE_SEEDS:
            recs = []
            for rec in lig_ok:
                sc = wide[(pair, seed, rec["ligand_id"])]
                recs.append({"cls": rec["primary_class_theta6"], "score_A": sc["A"], "score_B": sc["B"]})
            da, db, sm, nd, na, nb, n = _smin_point(recs)
            key = (pair, str(seed))
            cur = current.get(pair)
            # five_seed_summary_min has one row per pair×seed
            cur_row = next(
                r
                for r in read_csv(ROOT / "results/canonical/five_seed_summary_min.csv")
                if r["pair"] == pair and str(r["seed"]) == str(seed)
            )
            cur_sm = float(cur_row["summary_min"])
            out.append(
                {
                    "pair": pair,
                    "n_intersection": n_inter,
                    "n_dual": counts["dual"],
                    "n_A_only": counts["A_only"],
                    "n_B_only": counts["B_only"],
                    "n_neither": counts["neither"],
                    "n_primary_activity_eligible_complete_case": len(primary),
                    "n_dropped_vs_primary": len(primary) - n_inter,
                    "seed": seed,
                    "summary_min": None if n_inter == 0 else round(sm, 4),
                    "auroc_D_vs_A_pocketB": None if n_inter == 0 else round(da, 4),
                    "auroc_D_vs_B_pocketA": None if n_inter == 0 else round(db, 4),
                    "current_per_seed_summary_min": round(cur_sm, 4),
                    "current_per_seed_n": cur_row["n"],
                    "delta": None if n_inter == 0 else round(sm - cur_sm, 4),
                    "qualitative_change": "none_if_abs_delta_lt_0.02_else_document",
                    "source_kind": "qa_fixed_membership_sensitivity",
                    "note": "intersection of ligands with both pocket scores on all five seeds; current labels; no redock",
                }
            )
    return out


def path_audit() -> list[dict]:
    docs = [
        "docs/WRITING_INDEX_FREEZE.md",
        "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
        "docs/PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md",
        "docs/PR39_SCIENTIFIC_DATA_AUDIT.md",
        "data/jcim_chembl_universe_v0/local_track_b_v0/README.md",
    ]
    candidates = [
        "results/canonical/primary_summary_min.csv",
        "results/canonical/primary_directional_auroc.csv",
        "results/canonical/class_counts.csv",
        "results/canonical/two_pocket_mean_ranking.csv",
        "results/canonical/top10_operating_points.csv",
        "results/canonical/fixed_score_negative_class_delta.csv",
        "results/canonical/ecfp4_incremental_information.csv",
        "results/canonical/matched_minus_mismatched.csv",
        "results/canonical/holdout_metrics.csv",
        "results/canonical/descriptor_baselines.csv",
        "results/canonical/computational_robustness.csv",
        "results/canonical/receptor_substitution.csv",
        "results/canonical/cognate_rmsd.csv",
        "results/canonical/detectable_effect_simulation.csv",
        "results/canonical/max_vs_median_sensitivity.csv",
        "results/canonical/five_seed_summary_min.csv",
        "results/canonical/current_score_master.csv",
        "figures/jcim_article/plotted_values_postfix.json",
        "docs/HISTORICAL_LEFTOVER_FILES.md",
        "results/freeze_rebuild",
        "figures/jcim_article/plotted_values.json",
        "data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml",
        "data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv",
        "analysis/DOCKING_PLAN_V1.md",
        "analysis/LOCAL_RECOMPUTE_PACK_V1.md",
        "TIER1_DOCKING_ROSTER_V1.md",
        "FEASIBLE_PAIR_LADDER_V1.md",
        "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "data/jcim_multiseed_v0/tables/scores_vina_mode1_EGFR_corrected_box_fiveseed.csv",
        "data/egfr_her2_uniform_rdkit_v1/protocol/protocol.yaml",
        "results/canonical/ecfp4_scaler_sensitivity.csv",
        "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv",
        "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv",
        "data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv",
        "submission_pack/manuscript/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
        "scripts/analysis/analysis_config.py",
        "scripts/analysis/rebuild_freeze.py",
        "figures/jcim_article/scripts/update_figures_pr32.py",
    ]
    cited = []
    for doc in docs:
        text = (ROOT / doc).read_text(encoding="utf-8") if (ROOT / doc).is_file() else ""
        for rel in candidates:
            if rel in text or rel.replace("\\", "/") in text:
                cited.append((doc, rel))
    extra_hist = [
        "formulation_equal_score_negative_v1.csv",
        "unified_threshold_sensitivity_v2.csv",
        "table2_comparable_theta6_v1.csv",
        "eight_pair_ranking_operating_point_v1.csv",
        "wrong_pocket_paired_delta_bootstrap_v1.csv",
        "holdout_pocket_matched_v1.csv",
        "assembled_AChE_BChE.csv",
        "pm110_vs_pm48_pocket_matched_v1.csv",
        "equal_score_cluster_bootstrap_v1.csv",
        "primary_directional_intervals_review_v1.csv",
        "formulation_conventional_vs_directional_v1.csv",
        "ecfp4_docking_scaler_sensitivity_v1.csv",
        "HISTORICAL_LEFTOVER_FILES.md",
    ]
    lock = ROOT / "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md"
    lock_text = lock.read_text(encoding="utf-8") if lock.is_file() else ""
    rows = []
    seen = set()
    for doc, rel in cited:
        p = ROOT / rel
        status = "exists" if p.exists() else "missing"
        role = "current" if rel.startswith("results/canonical/") or rel.endswith("plotted_values_postfix.json") else "check"
        if rel in {"docs/HISTORICAL_LEFTOVER_FILES.md", "results/freeze_rebuild", "figures/jcim_article/plotted_values.json"}:
            role = "must_be_absent"
            verdict = "ok_absent" if not p.exists() else "stale_error_present"
        elif status == "missing":
            verdict = "stale_error" if "current" in doc.lower() or doc.startswith("docs/") else "archive_ok"
            if rel in {"analysis/DOCKING_PLAN_V1.md", "analysis/LOCAL_RECOMPUTE_PACK_V1.md", "TIER1_DOCKING_ROSTER_V1.md", "FEASIBLE_PAIR_LADDER_V1.md"}:
                verdict = "historical_citation_missing_file"
        else:
            verdict = "ok"
        key = (doc, rel)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "citing_file": doc,
                "referenced_path": rel,
                "exists": int(p.exists()),
                "path_kind": "file" if p.is_file() else ("dir" if p.is_dir() else "missing"),
                "verdict": verdict,
                "role": role,
            }
        )
    for name in extra_hist:
        if name in lock_text:
            hits = list(ROOT.rglob(name))
            rows.append(
                {
                    "citing_file": "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
                    "referenced_path": name,
                    "exists": int(bool(hits)),
                    "path_kind": "basename_search",
                    "verdict": "historical_alias_in_lock",
                    "role": "superseded_if_canonical_replacement_exists",
                }
            )
    return rows


def receptor_rows() -> list[dict]:
    rmsd = {r["pdb"]: r for r in read_csv(ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv")}
    slots = [
        ("EGFR/HER2", "A", "3POZ", "A", "03P", "data/egfr_her2_uniform_rdkit_v1/receptors/3POZ_receptor.pdbqt", "data/egfr_her2_panel120_v0/boxes/3POZ_box_corrected.json", "cognate_heavy_atom_AABB+5A_min20", "git blob 2b398c62 from c213c485"),
        ("EGFR/HER2", "B", "3RCD", "A", "03P", "data/egfr_her2_uniform_rdkit_v1/receptors/3RCD_receptor.pdbqt", "data/egfr_her2_panel120_v0/boxes/3RCD_box_corrected.json", "cognate_heavy_atom_AABB+5A_min20", "git blob c338094b from c213c485"),
        ("PIK3CA/mTOR", "A", "4L23", "A", "X6K", "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_receptor.pdbqt", "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4L23_box.json", "cognate_heavy_atom_box_file", "deposited pdbqt; prep script not in current tree"),
        ("PIK3CA/mTOR", "B", "4JT6", "A", "X6K", "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_receptor.pdbqt", "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4JT6_box.json", "cognate_heavy_atom_box_file", "deposited pdbqt; prep script not in current tree"),
        ("AChE/BChE", "A", "4EY7", "AB", "E20", "data/ache_bche_panel_v0/receptors/ACHE_receptor.pdbqt", "data/ache_bche_panel_v0/boxes/4EY7_box.json", "cognate_heavy_atom_box_file", "deposited pdbqt (alias 4EY7_receptor.pdbqt); dimer AB; prep script not in current tree"),
        ("AChE/BChE", "B", "4BDS", "A", "THA", "data/ache_bche_panel_v0/receptors/BCHE_receptor.pdbqt", "data/ache_bche_panel_v0/boxes/4BDS_box.json", "cognate_heavy_atom_box_file", "deposited pdbqt (alias 4BDS_receptor.pdbqt); prep script not in current tree"),
        ("F2/F10", "A", "4UDW", "HIL", "N6L", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/4UDW_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/4UDW_box.json", "cognate_heavy_atom_AABB+5A_min20", "deposited pdbqt; chains H/I/L; Track-B yaml site thrombin S1 chain H"),
        ("F2/F10", "B", "2JKH", "AL", "BI7", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/2JKH_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/2JKH_box.json", "cognate_heavy_atom_AABB+5A_min20", "deposited pdbqt; SI notes OpenBabel SDF issue then CCD for BI7"),
        ("JAK1/JAK2", "A", "6N7A", "AB", "KEV", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6N7A_box.json", "cognate_heavy_atom_AABB+5A_min20", "shared JAK1 receptor with JAK1/TYK2"),
        ("JAK1/JAK2", "B", "8BXH", "A", "C87", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/8BXH_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/8BXH_box.json", "cognate_heavy_atom_AABB+5A_min20", "deposited pdbqt"),
        ("JAK1/TYK2", "A", "6N7A", "AB", "KEV", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6N7A_box.json", "cognate_heavy_atom_AABB+5A_min20", "same file as JAK1/JAK2 pocket A"),
        ("JAK1/TYK2", "B", "3LXP", "A", "IZA", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/3LXP_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/3LXP_box.json", "cognate_heavy_atom_AABB+5A_min20", "JH1 ATP not JH2 (yaml)"),
        ("PPARG/PPARA", "A", "9V8H", "AB", "BRL", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/9V8H_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/9V8H_box.json", "cognate_heavy_atom_AABB+5A_min20", "chain B is 11-residue peptide (yaml: keep PG08-NL peptide); confirmed in pdbqt"),
        ("PPARG/PPARA", "B", "6LXA", "A", "EPA", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6LXA_box.json", "cognate_heavy_atom_AABB+5A_min20", "shared PPARA receptor with PPARA/PPARD"),
        ("PPARA/PPARD", "A", "6LXA", "A", "EPA", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6LXA_box.json", "cognate_heavy_atom_AABB+5A_min20", "same file as PPARG/PPARA pocket B"),
        ("PPARA/PPARD", "B", "5U3Q", "A", "7UJ", "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/5U3Q_receptor.pdbqt", "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/5U3Q_box.json", "cognate_heavy_atom_AABB+5A_min20", "yaml: LBD agonist 1 not PEG"),
    ]
    rows = []
    seen_pdb = set()
    for pair, pocket, pdb, chain, cognate, rec, box, boxdef, note in slots:
        unique = pdb not in seen_pdb
        seen_pdb.add(pdb)
        r = rmsd.get(pdb, {})
        rec_ok = exists(rec)
        box_ok = exists(box)
        # Unique 14 PDBs: file confirmed; prep method not reconstructed.
        file_status = "confirmed" if rec_ok else "not_recoverable"
        prep_status = "not_recoverable"
        rows.append(
            {
                "pair": pair,
                "pocket": pocket,
                "receptor_pdb": pdb,
                "unique_14slot": int(unique),
                "receptor_chain_or_construct": chain,
                "cognate_ligand": cognate,
                "receptor_input_file": rec,
                "receptor_file_exists": int(rec_ok),
                "box_file": box,
                "box_file_exists": int(box_ok),
                "box_definition": boxdef,
                "waters": "none_in_pdbqt",
                "cofactor_peptide": note,
                "protonation": "not_recoverable",
                "missing_loops": "not_recoverable",
                "calcrrms_top1_A": r.get("calcrrms_top1_A", ""),
                "calcrrms_best_A": r.get("calcrrms_best_A", ""),
                "receptor_file_status": file_status,
                "receptor_prep_status": prep_status,
                "notes": note,
            }
        )
    return rows


def matrix_rows() -> list[dict]:
    adj = "data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv"
    rows = []

    def add(**kwargs):
        rec = {k: "" for k in FIELDS_MATRIX}
        rec.update(kwargs)
        rows.append(rec)

    # EGFR CASE 2 uniform RDKit/Meeko is the only current production truth.
    add(
        pair="EGFR/HER2",
        channel="production_vina",
        role="current_primary",
        ligand_panel_source="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass; EH120_059 unresolved_missing_arm excluded; EH40_31 timeout_skipped",
        ligand_smiles_source="panel_v0_120.csv smiles",
        ligand_prep_method="rdkit_etkdgv3_mmff_meeko",
        ligand_prep_version="meeko 0.7.1",
        ligand_prep_seed="20260727",
        ligand_prep_evidence="data/egfr_her2_uniform_rdkit_v1/scripts/prep_uniform_ligands.py + ligands_pdbqt 110/110 + ligand_prep_status.csv",
        receptor_pdb="3POZ+3RCD",
        receptor_chain_or_construct="A / A",
        receptor_input_file="data/egfr_her2_uniform_rdkit_v1/receptors/{3POZ,3RCD}_receptor.pdbqt",
        receptor_prep_method="not_recoverable",
        receptor_prep_evidence="frozen pdbqt blobs restored from c213c485; no PPW/pH/loop record",
        cognate_ligand="03P / 03P",
        box_file="data/egfr_her2_panel120_v0/boxes/{3POZ,3RCD}_box_corrected.json",
        box_definition="cognate_heavy_atom_AABB+5A_min20",
        vina_version="1.2.7",
        dock_seed="20260727",
        exhaustiveness="8",
        num_modes="9",
        energy_range="3",
        score_definition="score_S=-vina_mode1; timeout_skipped at 600s",
        score_source_file=current_egfr_score_source(),
        sample_membership_rule="activity_eligible AND complete_case AND analysis_set=main",
        analysis_script="scripts/analysis/build_current_score_master.py; compute_canonical_results.py",
        canonical_output="results/canonical/primary_summary_min.csv",
        figure_table_destination="Table 2; Fig 2",
        provenance_status="confirmed" if egfr_uniform_ready() else "not_recoverable",
        notes="CASE 2 uniform RDKit/Meeko is the only current EGFR production truth",
    )
    add(
        pair="EGFR/HER2",
        channel="five_seed_vina",
        role="current_robustness",
        ligand_panel_source="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass",
        ligand_smiles_source="panel_v0_120.csv smiles",
        ligand_prep_method="rdkit_etkdgv3_mmff_meeko",
        ligand_prep_evidence="same uniform ligands as production",
        receptor_pdb="3POZ+3RCD",
        receptor_chain_or_construct="A / A",
        receptor_input_file="data/egfr_her2_uniform_rdkit_v1/receptors/{3POZ,3RCD}_receptor.pdbqt",
        receptor_prep_method="not_recoverable",
        receptor_prep_evidence="frozen pdbqt blobs",
        cognate_ligand="03P / 03P",
        box_file="data/egfr_her2_panel120_v0/boxes/{3POZ,3RCD}_box_corrected.json",
        box_definition="cognate_heavy_atom_AABB+5A_min20",
        vina_version="1.2.7",
        dock_seed="20260727,20260811-14",
        exhaustiveness="8",
        num_modes="9",
        energy_range="3",
        score_definition="score_S=-vina_mode1; timeout_skipped at 600s",
        score_source_file=EGFR_UNIFORM_VINA_REL,
        sample_membership_rule="per-seed complete-case on current eligible panel",
        analysis_script="scripts/analysis/compute_canonical_results.py:compute_five_seed",
        canonical_output="results/canonical/five_seed_summary_min.csv",
        figure_table_destination="Fig 5C; Table S7",
        provenance_status="confirmed" if egfr_uniform_ready() else "not_recoverable",
        notes="uniform five-seed is current; historical corrected-box table is not current truth",
    )
    add(
        pair="EGFR/HER2",
        channel="gnina_independent",
        role="current_robustness",
        ligand_panel_source="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass",
        ligand_prep_method="rdkit_etkdgv3_mmff_meeko",
        ligand_prep_evidence="same uniform ligands as production Vina",
        receptor_pdb="3POZ+3RCD",
        cognate_ligand="03P / 03P",
        box_file="data/egfr_her2_panel120_v0/boxes/{3POZ,3RCD}_box_corrected.json",
        vina_version="GNINA 1.3.2 independent search (not Vina rescore)",
        dock_seed="20260727",
        exhaustiveness="8",
        num_modes="9",
        score_definition="score_S=-gnina_mode1; timeout_skipped at 600s",
        score_source_file=EGFR_UNIFORM_GNINA_REL if (ROOT / EGFR_UNIFORM_GNINA_REL).is_file() else "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv",
        sample_membership_rule="activity_eligible complete-case intersect GNINA both pockets",
        analysis_script="scripts/analysis/compute_canonical_results.py",
        canonical_output="results/canonical/computational_robustness.csv",
        figure_table_destination="Fig 5A; Table S7",
        provenance_status="confirmed" if (ROOT / EGFR_UNIFORM_GNINA_REL).is_file() else "not_recoverable",
        notes="independent pose generation on uniform ligands; EH120_109 fail (invalid PDBQT atom type CG0); timeouts skipped",
    )
    add(
        pair="EGFR/HER2",
        channel="historical_deposited_ablation",
        role="historical_not_current",
        ligand_panel_source="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass",
        ligand_smiles_source="panel_v0_120.csv smiles",
        ligand_prep_method="not_recoverable",
        ligand_prep_evidence="deleted production PDBQT; not used after CASE 2 promotion",
        receptor_pdb="3POZ+3RCD",
        score_source_file="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        canonical_output="not current",
        figure_table_destination="none",
        provenance_status="not_recoverable",
        notes="historical deposited scores; not Table 2 after uniform promotion",
    )

    # PIK3CA
    add(
        pair="PIK3CA/mTOR",
        channel="production_vina",
        role="current_primary",
        ligand_panel_source="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass",
        ligand_smiles_source="panel_v0_48.csv smiles",
        ligand_prep_method="rdkit_etkdg_meeko",
        ligand_prep_version="protocol.yaml ligand_prep field",
        ligand_prep_seed="20260727",
        ligand_prep_evidence="protocol.yaml + primary score table; ligand PDBQT not in current tree; local LigPrep job is predecessor only",
        receptor_pdb="4L23+4JT6",
        receptor_chain_or_construct="A / A",
        receptor_input_file="data/pik3ca_mtor_panel48_rdkit_v0/receptors/{4L23,4JT6}_receptor.pdbqt",
        receptor_prep_method="not_recoverable",
        receptor_prep_evidence="deposited pdbqt only",
        cognate_ligand="X6K / X6K",
        box_file="data/pik3ca_mtor_panel48_rdkit_v0/boxes/{4L23,4JT6}_box.json",
        box_definition="deposited box json",
        vina_version="1.2.7",
        dock_seed="20260727",
        exhaustiveness="16",
        num_modes="9",
        energy_range="",
        score_definition="score_S=-vina_mode1",
        score_source_file="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        sample_membership_rule="activity_eligible AND complete_case AND analysis_set=main",
        analysis_script="scripts/analysis/build_current_score_master.py; compute_canonical_results.py",
        canonical_output="results/canonical/primary_summary_min.csv",
        figure_table_destination="Table 2; Fig 2",
        provenance_status="confirmed",
        notes="scores+receptors+boxes+panel exist; ligand 3D inputs not deposited; LigPrep is not current primary",
    )
    add(
        pair="PIK3CA/mTOR",
        channel="five_seed_vina",
        role="current_robustness",
        ligand_panel_source="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated",
        label_rule="theta_6.0_fourclass",
        ligand_prep_method="rdkit_etkdg_meeko",
        ligand_prep_evidence="same panel as production",
        receptor_pdb="4L23+4JT6",
        receptor_input_file="data/pik3ca_mtor_panel48_rdkit_v0/receptors/{4L23,4JT6}_receptor.pdbqt",
        cognate_ligand="X6K / X6K",
        box_file="data/pik3ca_mtor_panel48_rdkit_v0/boxes/{4L23,4JT6}_box.json",
        vina_version="1.2.7",
        dock_seed="20260727,20260811-14",
        exhaustiveness="16",
        num_modes="9",
        score_definition="score_S=-vina_mode1",
        score_source_file="data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv",
        sample_membership_rule="per-seed complete-case on current eligible panel",
        analysis_script="scripts/analysis/compute_canonical_results.py:compute_five_seed",
        canonical_output="results/canonical/five_seed_summary_min.csv",
        figure_table_destination="Fig 5C; Table S7",
        provenance_status="confirmed",
        notes="deposited long scores exist",
    )
    add(
        pair="PIK3CA/mTOR",
        channel="gnina_independent",
        role="current_robustness",
        activity_source_kind="chembl_assay_adjudicated",
        ligand_prep_method="rdkit_etkdg_meeko",
        receptor_pdb="4L23+4JT6",
        vina_version="GNINA independent search",
        score_definition="score_S=-gnina_mode1",
        score_source_file="data/jcim_independent_dock_v0/tables/gnina_dock_scores_PIK3CA_mTOR.csv",
        analysis_script="scripts/analysis/compute_canonical_results.py",
        canonical_output="results/canonical/computational_robustness.csv",
        figure_table_destination="Fig 5A; Table S7",
        provenance_status="confirmed",
        notes="independent pose generation, not CNN rescore of Vina",
    )
    add(
        pair="PIK3CA/mTOR",
        channel="receptor_substitution",
        role="sensitivity",
        receptor_pdb="4JPS / 5DXT / 4JSX",
        score_source_file="data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt*.csv",
        analysis_script="scripts/analysis/compute_canonical_results.py",
        canonical_output="results/canonical/receptor_substitution.csv",
        figure_table_destination="Fig 5B; Table S7",
        provenance_status="confirmed",
        notes="sensitivity only; primary remains 4L23/4JT6",
    )
    add(
        pair="PIK3CA/mTOR",
        channel="protocol_E8_and_PM110",
        role="sensitivity",
        score_source_file="data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv; data/pik3ca_mtor_panel110_rdkit_v0/tables/ablation_ligand_scores.csv",
        canonical_output="results/canonical/protocol_sensitivity.csv",
        figure_table_destination="Fig S2; Table S7",
        provenance_status="confirmed",
        notes="not current primary; E=16 PM48 remains primary",
    )
    add(
        pair="PIK3CA/mTOR",
        channel="historical_ligprep_scores",
        role="historical/superseded",
        ligand_prep_method="LigPrep (predecessor campaign)",
        ligand_prep_evidence="local Maestro log/sh/inp; maegz not in git; prep_delta_vs_ligprep.csv is comparison table",
        score_source_file="data/pik3ca_mtor_panel48_rdkit_v0/tables/prep_delta_vs_ligprep.csv",
        provenance_status="superseded",
        notes="must not be described as current primary",
    )

    # AChE
    add(
        pair="AChE/BChE",
        channel="production_vina",
        role="current_primary",
        ligand_panel_source="data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        activity_source=adj,
        activity_source_kind="chembl_assay_adjudicated + 5 panel_pchembl_no_audit_rows (only AB_056 in primary)",
        label_rule="panel constructed strict_6.5_5.5; primary evaluation theta_6.0; AB_087 unresolved excluded",
        ligand_smiles_source="panel_v0_strict.csv smiles",
        ligand_prep_method="rdkit_etkdg_meeko",
        ligand_prep_seed="20260727",
        ligand_prep_evidence="protocol.yaml; ligand PDBQT not in current tree",
        receptor_pdb="4EY7+4BDS",
        receptor_chain_or_construct="AB / A",
        receptor_input_file="data/ache_bche_panel_v0/receptors/{ACHE,BCHE}_receptor.pdbqt",
        receptor_prep_method="not_recoverable",
        receptor_prep_evidence="deposited pdbqt; alt receptors 5DYW/6QAA/6ZWI sensitivity only",
        cognate_ligand="E20 / THA",
        box_file="data/ache_bche_panel_v0/boxes/{4EY7,4BDS}_box.json",
        vina_version="1.2.7",
        dock_seed="20260727",
        exhaustiveness="8",
        num_modes="9",
        score_definition="score_S=-vina_mode1",
        score_source_file="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        sample_membership_rule="activity_eligible AND complete_case AND analysis_set=main",
        analysis_script="scripts/analysis/build_current_score_master.py; compute_canonical_results.py",
        canonical_output="results/canonical/primary_summary_min.csv",
        figure_table_destination="Table 2; Fig 2; Fig 4",
        provenance_status="confirmed",
        notes="panel label_rule vs theta=6.0 evaluation is documented dual-layer, not two current truths",
    )
    add(
        pair="AChE/BChE",
        channel="five_seed_vina",
        role="current_robustness",
        ligand_panel_source="data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        activity_source_kind="chembl_assay_adjudicated",
        ligand_prep_method="rdkit_etkdg_meeko",
        receptor_pdb="4EY7+4BDS",
        vina_version="1.2.7",
        dock_seed="20260727,20260811-14",
        exhaustiveness="8",
        num_modes="9",
        score_source_file="data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv",
        canonical_output="results/canonical/five_seed_summary_min.csv",
        figure_table_destination="Fig 5C",
        provenance_status="confirmed",
        notes="",
    )
    add(
        pair="AChE/BChE",
        channel="gnina_panel_scores",
        role="historical/superseded_as_independent_gnina",
        score_source_file="data/ache_bche_panel_v0/tables/scores_gnina_best.csv",
        canonical_output="not in GNINA_SOURCES",
        figure_table_destination="none",
        provenance_status="superseded",
        notes="not wired into computational_robustness GNINA channel",
    )

    track_b = [
        ("JAK1/JAK2", "6N7A+8BXH", "KEV / C87", False),
        ("JAK1/TYK2", "6N7A+3LXP", "KEV / IZA", True),
        ("F2/F10", "4UDW+2JKH", "N6L / BI7", False),
        ("PPARG/PPARA", "9V8H+6LXA", "BRL / EPA", False),
        ("PPARA/PPARD", "6LXA+5U3Q", "EPA / 7UJ", False),
    ]
    for pair, pdbs, cognate, has_gnina in track_b:
        panel = {
            "JAK1/JAK2": "panel_JAK1_JAK2_v1.csv",
            "JAK1/TYK2": "panel_JAK1_TYK2_v1.csv",
            "F2/F10": "panel_F2_F10_v1.csv",
            "PPARG/PPARA": "panel_PPARG_PPARA_v1.csv",
            "PPARA/PPARD": "panel_PPARA_PPARD_v1.csv",
        }[pair]
        add(
            pair=pair,
            channel="production_vina",
            role="current_primary",
            ligand_panel_source=f"data/jcim_chembl_universe_v0/tables/track_b_panels/{panel}",
            activity_source=f"data/jcim_chembl_universe_v0/tables/track_b_panels/{panel}",
            activity_source_kind="chembl37_dump_panel",
            label_rule="panel constructed strict_6.5_5.5; primary evaluation theta_6.0 from dump pChEMBL; not assay-adjudicated",
            ligand_smiles_source=f"{panel} canonical_smiles",
            ligand_prep_method="rdkit_etkdgv3_mmff_meeko",
            ligand_prep_version="meeko 0.7.1",
            ligand_prep_seed="20260727",
            ligand_prep_evidence="track_b_local_run_v1.yaml no_ligprep; historical prep_track_b_ligands_v1.py in git; ligands_pdbqt dir absent; production poses deposited",
            receptor_pdb=pdbs,
            receptor_input_file="data/jcim_chembl_universe_v0/local_track_b_v0/receptors/*_receptor.pdbqt",
            receptor_prep_method="not_recoverable",
            receptor_prep_evidence="deposited pdbqt; yaml records site/cognate; no PPW/pH/loop protocol",
            cognate_ligand=cognate,
            box_file="data/jcim_chembl_universe_v0/local_track_b_v0/boxes/{pdb}_box.json",
            box_definition="cognate heavy-atom AABB + 5 A, min edge 20 A",
            vina_version="1.2.7",
            dock_seed="20260727",
            exhaustiveness="8",
            num_modes="9",
            energy_range="3",
            score_definition="score_S=-vina_mode1",
            score_source_file="data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv",
            sample_membership_rule="activity_eligible default 1 AND complete_case AND analysis_set=main",
            analysis_script="scripts/analysis/build_current_score_master.py; compute_canonical_results.py",
            canonical_output="results/canonical/primary_summary_min.csv",
            figure_table_destination="Table 2; Fig 2",
            provenance_status="confirmed",
            notes="Track-B is dump-gated labels; do not call assay-adjudicated",
        )
        add(
            pair=pair,
            channel="five_seed_vina",
            role="current_robustness",
            activity_source_kind="chembl37_dump_panel",
            ligand_prep_method="rdkit_etkdgv3_mmff_meeko",
            receptor_pdb=pdbs,
            vina_version="1.2.7",
            dock_seed="20260727,20260811-14",
            exhaustiveness="8",
            num_modes="9",
            energy_range="3",
            score_source_file="data/jcim_chembl_universe_v0/local_track_b_v0/tables/multiseed/scores_vina_mode1_seed{seed}.csv",
            canonical_output="results/canonical/five_seed_summary_min.csv",
            figure_table_destination="Fig 5C",
            provenance_status="confirmed",
            notes="",
        )
        add(
            pair=pair,
            channel="rtm_best9_rescore",
            role="sensitivity",
            score_source_file="data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_rtm_best9_v1.csv",
            canonical_output="results/canonical/computational_robustness.csv",
            figure_table_destination="Table S7 not Fig 5",
            provenance_status="confirmed",
            notes="rescore of Vina poses, not independent docking",
        )
        add(
            pair=pair,
            channel="gnina_cnn_best9_rescore",
            role="sensitivity",
            score_source_file="data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_cnn_best9_v1.csv",
            canonical_output="results/canonical/computational_robustness.csv",
            figure_table_destination="Table S7 not Fig 5",
            provenance_status="confirmed",
            notes="CNN rescore of Vina poses; distinct from independent GNINA search",
        )
        if has_gnina:
            add(
                pair=pair,
                channel="gnina_independent",
                role="current_robustness",
                activity_source_kind="chembl37_dump_panel",
                receptor_pdb=pdbs,
                vina_version="GNINA independent search",
                score_source_file="data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_independent_jak1_tyk2_v1.csv",
                canonical_output="results/canonical/computational_robustness.csv",
                figure_table_destination="Fig 5A; Table S7",
                provenance_status="confirmed",
                notes="independent search on JAK1/TYK2 only among Track-B pairs",
            )
    return rows


def figure_table_map() -> list[dict]:
    return [
        {"item": "Table 1", "displayed": "pair quotas / PDB / exhaustiveness / n_scored", "canonical_or_frozen_source": "results/canonical/class_counts.csv + panel CSVs + analysis_config.RECEPTORS", "row_or_key": "eight PRIMARY_PAIRS", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists; lock previously named unspecified panel tables"},
        {"item": "Table 2", "displayed": "directional AUROC + summary_min + n_scored", "canonical_or_frozen_source": "results/canonical/primary_directional_auroc.csv; primary_summary_min.csv; class_counts.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists; historical unified_threshold_sensitivity_v2.csv superseded"},
        {"item": "Table 3", "displayed": "two-pocket mean D-vs-neither; Top10%; EF", "canonical_or_frozen_source": "results/canonical/two_pocket_mean_ranking.csv; top10_operating_points.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Table S1", "displayed": "computational settings", "canonical_or_frozen_source": "scripts/analysis/analysis_config.py; requirements-analysis.txt; data/jcim_strengthen_t0t1_v0/ENV_PIN.md", "row_or_key": "constants", "generating_script": "n/a (typeset from config)", "status": "config_exists; analysis freeze Python 3.12.3"},
        {"item": "Table S2", "displayed": "receptors/boxes/RMSD", "canonical_or_frozen_source": "results/canonical/cognate_rmsd.csv copied from all14_cognate_rmsd_calcrrms_v1.csv; pair box json", "row_or_key": "14 pdb slots", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists; receptor prep method not in S2"},
        {"item": "Table S3", "displayed": "threshold / max-vs-median", "canonical_or_frozen_source": "results/canonical/label_aggregation_sensitivity.csv; max_vs_median_sensitivity.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Table S4", "displayed": "fixed-score Δ; cluster bootstrap; detectable-effect", "canonical_or_frozen_source": "results/canonical/fixed_score_negative_class_delta.csv; cluster_bootstrap_sensitivity.csv; detectable_effect_simulation.csv", "row_or_key": "pair/contrast", "generating_script": "compute_canonical_results.py; compute_detectable_effect.py", "status": "canonical_exists"},
        {"item": "Table S5", "displayed": "descriptor / ECFP / incremental", "canonical_or_frozen_source": "results/canonical/descriptor_baselines.csv; descriptor_nested_scaffold_cv.csv; ecfp4_incremental_information.csv; ecfp4_scaler_sensitivity.csv", "row_or_key": "pair/arm/model", "generating_script": "compute_descriptor_baselines.py; fit_ecfp4_models.py", "status": "canonical_exists"},
        {"item": "Table S6", "displayed": "matched vs mismatched main+holdout", "canonical_or_frozen_source": "results/canonical/matched_mismatched_pocket.csv; matched_minus_mismatched.csv; holdout_metrics.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Table S7", "displayed": "GNINA / receptor swap / RTM/CNN / five-seed", "canonical_or_frozen_source": "results/canonical/computational_robustness.csv; receptor_substitution.csv; five_seed_summary_min.csv", "row_or_key": "pair/channel", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Table S8", "displayed": "external eligibility", "canonical_or_frozen_source": "results/canonical/external_eligibility.csv from external_slice_summary_v1.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Table S9", "displayed": "pair eligibility", "canonical_or_frozen_source": "data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv", "row_or_key": "pair", "generating_script": "frozen audit table", "status": "exists; TIER1/FEASIBLE md cited but missing"},
        {"item": "Table S10", "displayed": "ranking operating points + AND filter", "canonical_or_frozen_source": "results/canonical/top10_operating_points.csv; and_filter_operating_points.csv", "row_or_key": "pair", "generating_script": "scripts/analysis/compute_canonical_results.py", "status": "canonical_exists"},
        {"item": "Figure 1A", "displayed": "four-state schematic", "canonical_or_frozen_source": "schematic", "row_or_key": "n/a", "generating_script": "figures/jcim_article/scripts/update_figures_pr32.py", "status": "schematic"},
        {"item": "Figure 1B", "displayed": "directional schematic", "canonical_or_frozen_source": "schematic", "row_or_key": "n/a", "generating_script": "update_figures_pr32.py", "status": "schematic"},
        {"item": "Figure 1C", "displayed": "universe census counts", "canonical_or_frozen_source": "data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv", "row_or_key": "census rows", "generating_script": "update_figures_pr32.py", "status": "exists"},
        {"item": "Figure 2A", "displayed": "fixed-score ΔAUROC", "canonical_or_frozen_source": "results/canonical/fixed_score_negative_class_delta.csv", "row_or_key": "pair × direction", "generating_script": "update_figures_pr32.py", "status": "canonical; lock historical name superseded"},
        {"item": "Figure 2B", "displayed": "two directional AUROCs", "canonical_or_frozen_source": "results/canonical/primary_directional_auroc.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 2C", "displayed": "summary_min vs two-pocket mean", "canonical_or_frozen_source": "results/canonical/primary_summary_min.csv; two_pocket_mean_ranking.csv; class_counts.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 2D", "displayed": "top-10% stacked bars", "canonical_or_frozen_source": "results/canonical/top10_operating_points.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 3A", "displayed": "Vina vs ECFP4 OOF", "canonical_or_frozen_source": "results/canonical/ecfp4_incremental_information.csv", "row_or_key": "pair/arm", "generating_script": "fit_ecfp4_models.py; update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 3B", "displayed": "ΔAUROC ECFP4+docking − ECFP4", "canonical_or_frozen_source": "results/canonical/ecfp4_incremental_information.csv", "row_or_key": "max |Δ|=0.0234 PPARA/PPARD D vs B", "generating_script": "fit_ecfp4_models.py; update_figures_pr32.py", "status": "canonical; 0.0112 superseded"},
        {"item": "Figure 4A", "displayed": "main matched−mismatched", "canonical_or_frozen_source": "results/canonical/matched_mismatched_pocket.csv", "row_or_key": "pair main_panel", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 4B", "displayed": "holdout matched−mismatched", "canonical_or_frozen_source": "results/canonical/holdout_metrics.csv", "row_or_key": "pair holdout", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 5A", "displayed": "independent GNINA", "canonical_or_frozen_source": "results/canonical/computational_robustness.csv", "row_or_key": "EGFR/PIK3CA/JAK1-TYK2", "generating_script": "update_figures_pr32.py", "status": "canonical; EGFR ligand-prep chain broken"},
        {"item": "Figure 5B", "displayed": "PIK3CA receptor substitution", "canonical_or_frozen_source": "results/canonical/receptor_substitution.csv", "row_or_key": "alt pdb", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure 5C", "displayed": "five-seed summary_min range", "canonical_or_frozen_source": "results/canonical/five_seed_summary_min.csv", "row_or_key": "pair × seed", "generating_script": "update_figures_pr32.py", "status": "canonical; EGFR inputs not_recoverable"},
        {"item": "Figure S1A", "displayed": "best descriptor vs Vina", "canonical_or_frozen_source": "results/canonical/descriptor_baselines.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure S1B", "displayed": "AChE TPSA jitter", "canonical_or_frozen_source": "results/canonical/current_score_master.csv", "row_or_key": "AChE/BChE ligands", "generating_script": "update_figures_pr32.py", "status": "canonical; lock assembled_AChE_BChE.csv superseded"},
        {"item": "Figure S2A", "displayed": "PM48 vs PM110", "canonical_or_frozen_source": "results/canonical/protocol_sensitivity.csv", "row_or_key": "PIK3CA/mTOR", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure S2B", "displayed": "E16 vs E8", "canonical_or_frozen_source": "results/canonical/protocol_sensitivity.csv", "row_or_key": "PIK3CA/mTOR", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure S3", "displayed": "14-slot cognate RMSD", "canonical_or_frozen_source": "results/canonical/cognate_rmsd.csv", "row_or_key": "pdb", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure S4A", "displayed": "theta sensitivity", "canonical_or_frozen_source": "results/canonical/label_aggregation_sensitivity.csv", "row_or_key": "pair/theta", "generating_script": "update_figures_pr32.py", "status": "canonical; lock historical name superseded"},
        {"item": "Figure S4B", "displayed": "cluster bootstrap", "canonical_or_frozen_source": "results/canonical/cluster_bootstrap_sensitivity.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
        {"item": "Figure S5A/B", "displayed": "external eligibility", "canonical_or_frozen_source": "results/canonical/external_eligibility.csv", "row_or_key": "pair", "generating_script": "update_figures_pr32.py", "status": "canonical"},
    ]


def main() -> int:
    QA.mkdir(parents=True, exist_ok=True)
    seed_rows = five_seed_fixed_membership()
    for r in seed_rows:
        if r["delta"] is None:
            r["qualitative_change"] = "empty_intersection"
        else:
            r["qualitative_change"] = "no" if abs(float(r["delta"])) < 0.02 else "document"
    write_csv(QA / "five_seed_fixed_membership_sensitivity.csv", seed_rows, list(seed_rows[0].keys()))
    rec_rows = receptor_rows()
    write_csv(QA / "receptor_provenance_14slot.csv", rec_rows, list(rec_rows[0].keys()))
    mat = matrix_rows()
    write_csv(QA / "SCIENTIFIC_PROVENANCE_MATRIX.csv", mat, FIELDS_MATRIX)
    write_csv(ROOT / "docs/SCIENTIFIC_PROVENANCE_MATRIX.csv", mat, FIELDS_MATRIX)
    paths = path_audit()
    write_csv(QA / "current_truth_path_existence.csv", paths, list(paths[0].keys()) if paths else ["citing_file"])
    fig = figure_table_map()
    write_csv(QA / "figure_table_source_audit.csv", fig, list(fig[0].keys()))
    print(f"wrote {QA} rows: matrix={len(mat)} receptors={len(rec_rows)} five_seed={len(seed_rows)} paths={len(paths)} fig={len(fig)}")
    # summary for the report
    for pair in PRIMARY_PAIRS:
        sub = [r for r in seed_rows if r["pair"] == pair]
        deltas = [r["delta"] for r in sub]
        print(f"FIXED {pair:14} n_inter={sub[0]['n_intersection']} dropped={sub[0]['n_dropped_vs_primary']} deltas={deltas}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
