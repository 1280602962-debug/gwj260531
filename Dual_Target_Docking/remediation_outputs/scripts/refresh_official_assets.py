#!/usr/bin/env python3
"""Refresh official score tables that figures and SI still read."""
from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("updated", path)


def update_fiveseed() -> None:
    dest = ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv"
    rows = list(csv.DictReader(dest.open()))
    egfr = {
        r["seed"]: r
        for r in csv.DictReader((ROOT / "remediation_outputs/phase_fiveseed/multiseed_auroc_by_seed_corrected.csv").open())
    }
    ache = {
        r["seed"]: r
        for r in csv.DictReader((ROOT / "remediation_outputs/phase_ache_fiveseed/multiseed_auroc_by_seed_ACHE_corrected.csv").open())
    }
    out = []
    for r in rows:
        src = None
        if r["pair"] == "EGFR/HER2":
            src = egfr.get(r["seed"])
        elif r["pair"] == "AChE/BChE":
            src = ache.get(r["seed"])
        if src:
            r["n_complete"] = src["n_complete"]
            r["n_dual"] = src["n_dual"]
            r["n_A_only"] = src["n_A_only"]
            r["n_B_only"] = src["n_B_only"]
            r["n_neither"] = src["n_neither"]
            r["auroc_dual_vs_A_only"] = src["auroc_dual_vs_A_only"]
            r["auroc_dual_vs_B_only"] = src["auroc_dual_vs_B_only"]
            r["summary_min"] = src["summary_min"]
            r["auroc_dual_vs_neither_vina_mean"] = src["auroc_dual_vs_neither_vina_mean"]
            r["formulation_gap_neither_minus_summary_min"] = src["formulation_gap_neither_minus_summary_min"]
            if "mean_marginal_pocket_auroc_D_vs_neither" in r:
                r["mean_marginal_pocket_auroc_D_vs_neither"] = src["auroc_dual_vs_neither_vina_mean"]
        out.append(r)
    write_csv(dest, out)


def update_wrong_pocket() -> None:
    dest = ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv"
    rows = list(csv.DictReader(dest.open()))
    for r in rows:
        if r["pair"] == "EGFR/HER2" and r["set"] == "main_panel":
            r["matched_D_vs_A"] = "0.6607"
            r["matched_D_vs_B"] = "0.3237"
            r["matched_summary_min"] = "0.3237"
            r["matched_ci_lo"] = "0.1953"
            r["matched_ci_hi"] = "0.4710"
            r["wrong_summary_min"] = "0.2679"
            r["delta_matched_minus_wrong"] = "0.0558"
            r["delta_boot_mean"] = "0.0558"
            r["delta_ci_lo"] = "-0.0440"
            r["delta_ci_hi"] = "0.1600"
            r["ci_excludes_zero"] = "False"
            r["point_matched_gt_wrong"] = "True"
        if r["pair"] == "AChE/BChE" and r["set"] == "main_panel":
            r["n_A_only"] = "26"
            r["matched_D_vs_A"] = "0.6524"
            r["matched_D_vs_B"] = "0.6058"
            r["matched_summary_min"] = "0.6058"
            r["delta_matched_minus_wrong"] = "0.1770"
            r["delta_boot_mean"] = "0.1770"
            r["delta_ci_lo"] = "0.0498"
            r["delta_ci_hi"] = "0.2973"
            r["ci_excludes_zero"] = "True"
        if r["pair"] == "AChE/BChE" and r["set"] == "unused_pool_holdout":
            r["matched_D_vs_A"] = "0.615"
            r["matched_summary_min"] = "0.615"
            r["delta_matched_minus_wrong"] = "0.025"
            r["delta_boot_mean"] = "0.025"
            r["delta_ci_lo"] = "-0.0900"
            r["delta_ci_hi"] = "0.1075"
            r["ci_excludes_zero"] = "False"
            r["point_matched_gt_wrong"] = "True"
    write_csv(dest, rows)


def update_ml_and_holdout() -> None:
    dest = ROOT / "data/jcim_strengthen_t0t1_v0/tables/ligand_ml_baseline_scaffold_cv_v1.csv"
    rows = list(csv.DictReader(dest.open()))
    dock = {
        ("EGFR/HER2", "D_vs_A"): 0.6607,
        ("EGFR/HER2", "D_vs_B"): 0.3237,
        ("AChE/BChE", "D_vs_A"): 0.6524,
        ("AChE/BChE", "D_vs_B"): 0.6058,
    }
    for r in rows:
        key = (r["pair"], r["contrast"])
        if key in dock:
            r["auroc_dock_pocket_matched"] = f"{dock[key]:.4f}"
            r["delta_ml_minus_dock"] = f"{float(r['auroc_ml']) - dock[key]:.4f}"
    write_csv(dest, rows)

    hold = ROOT / "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv"
    rows = list(csv.DictReader(hold.open()))
    for r in rows:
        if r["prefix"] == "HOAB" and r["variant"] == "pocket_matched_vina":
            r["auroc_D_vs_A"] = "0.615"
            r["summary_min"] = "0.615"
            r["delta_vs_main_panel"] = f"{0.615 - 0.6058:.4f}"
    write_csv(hold, rows)


def append_ache_assembled() -> None:
    dest = ROOT / "data/jcim_bench_v0/tables/assembled_AChE_BChE.csv"
    rows = list(csv.DictReader(dest.open()))
    if any(r["ligand"] == "AB_056" for r in rows):
        print("assembled already has AB_056")
        return
    smiles = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3cccc(Cl)c3)CC1)C2"
    mol = Chem.MolFromSmiles(smiles)
    rows.append(
        {
            "pair": "AChE/BChE",
            "ligand": "AB_056",
            "cls": "A_only",
            "smiles": smiles,
            "label_rule": "strict_6.5_5.5",
            "prep": "rdkit_meeko",
            "pA": "6.94",
            "pB": "4.3",
            "min_pchembl": "4.3",
            "heavy": str(float(mol.GetNumHeavyAtoms())),
            "mw": str(float(Descriptors.MolWt(mol))),
            "clogp": str(float(Crippen.MolLogP(mol))),
            "tpsa": str(float(Descriptors.TPSA(mol))),
            "vina_mean": "10.9135",
            "vina_min": "9.05",
            "rtm_min": "",
            "rtm_mean": "",
            "rtm_min_z": "",
            "gnina_cnn_min": "",
        }
    )
    write_csv(dest, rows)


def update_cognate() -> None:
    src = list(csv.DictReader((ROOT / "remediation_outputs/phase1_cognate_redock/cognate_rank_rmsd_corrected_box.csv").open()))
    dest = ROOT / "data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv"
    old = list(csv.DictReader(dest.open()))
    keep = [r for r in old if r["pdb"] not in {"3POZ", "3RCD"}]
    for r in src:
        keep.append(
            {
                "target": r["target"],
                "pdb": r["pdb"],
                "exhaustiveness": r["exhaustiveness"],
                "n_modes_deposited": r["n_modes_deposited"],
                "pose_rank": r["pose_rank"],
                "rmsd_A": r["rmsd_A"],
                "best_top1_A": r["best_top1_A"],
                "best_top3_A": r["best_top3_A"],
                "best_all_deposited_A": r["best_all_deposited_A"],
                "pass_top1_lt2": r["pass_top1_lt2"],
                "pass_top3_lt2": r["pass_top3_lt2"],
                "pass_all_deposited_lt2": r["pass_all_deposited_lt2"],
                "mapping_method": r["mapping_method"],
                "mapping_max_distance_A": r["mapping_max_distance_A"],
                "rmsd_method": r["rmsd_method"],
            }
        )
    write_csv(dest, keep)

    summ = ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv"
    rows = list(csv.DictReader(summ.open()))
    by = {}
    for r in src:
        if int(r["pose_rank"]) == 1:
            by[r["pdb"]] = r
    for r in rows:
        if r["pdb"] in by:
            s = by[r["pdb"]]
            r["pose_status"] = "canonical_heavy_atom_redock"
            r["calcrrms_top1_A"] = f"{float(s['best_top1_A']):.3f}"
            r["calcrrms_top3_A"] = f"{float(s['best_top3_A']):.3f}"
            r["calcrrms_best_A"] = f"{float(s['best_all_deposited_A']):.3f}"
            r["calcrrms_best_mode"] = "1"
            r["pass_best_lt2"] = s["pass_all_deposited_lt2"]
            r["pass_top1_lt2"] = s["pass_top1_lt2"]
            r["pass_top3_lt2"] = s["pass_top3_lt2"]
            r["mapping_method"] = s["mapping_method"]
    write_csv(summ, rows)


def update_gnina_summary() -> None:
    dest = ROOT / "data/jcim_independent_dock_v0/tables/independent_dock_summary_v1.csv"
    rows = list(csv.DictReader(dest.open()))
    for r in rows:
        if r["pair"] == "EGFR/HER2":
            r["gnina_summary_min"] = "0.2645"
            r["gnina_D_vs_neither_mean"] = "0.7370"
            r["vina_summary_min_ref"] = "0.3237"
            r["vina_D_vs_neither_ref"] = "0.7589"
            r["delta_summary_min_vs_vina"] = f"{0.2645 - 0.3237:.4f}"
            r["verdict"] = "gap_remains"
    write_csv(dest, rows)


def main() -> int:
    update_fiveseed()
    update_wrong_pocket()
    update_ml_and_holdout()
    append_ache_assembled()
    update_cognate()
    update_gnina_summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
