#!/usr/bin/env python3
"""Score stage 3B poses with the locked A2 and NLRP3 structural gates."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/gwj/work/stage3b-refdock/URAT1_NLRP3_DualTarget_AIDD_Project")
sys.path.insert(0, str(ROOT / "scripts"))
from c1_acid_pose_selection import evaluate_urat1_acid_sdf  # noqa: E402
from c1_nlrp3_pose_metrics import evaluate_nlrp3_structural  # noqa: E402
from parse_c1_sdf_readouts import load_poses, load_ref_centroid  # noqa: E402

REFS = ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
ARG = REFS / "arg477_coords.json"
REF_U = load_ref_centroid(REFS / "lesinurad_crystal_ref.sdf")
BASE = ROOT / "results/nlrp3_model_v2/stage3b_local"
NEW = BASE / "new32"
GAP = BASE / "gap_REP_07837"
HIST = ROOT / "data/campaigns/c1/07_clinical_dock"
OUT = BASE / "analysis"
SEEDS = (42, 43, 44)


def hist_sdf(ligand: str, target: str, seed: int) -> tuple[Path, str]:
    if ligand == "REP_07837" and target == "nlrp3_7alv" and seed in (43, 44):
        return GAP / "rerun" / f"seed{seed}" / "poses_out.sdf", "stage3b_gap_fill"
    folder = "acid_dual" if seed == 42 else "acid_dual_a2"
    return HIST / folder / target / f"seed{seed}" / f"{ligand}_out.sdf", "historical_deposited"


def new_sdf(ligand: str, target: str, seed: int) -> Path:
    return NEW / "poses" / target / f"seed{seed}" / ligand / "poses_out.sdf"


def dock_state(sdf: Path, status_path: Path | None) -> tuple[str, int]:
    if status_path is not None and status_path.exists():
        rec = json.loads(status_path.read_text())
        state = rec.get("state", "failed")
        n = int(rec.get("n_poses") or 0)
        if state == "completed" and sdf.exists() and sdf.stat().st_size > 0:
            return "dock_completed", n
        if state == "timeout":
            return "timeout", n
        return "dock_failed", n
    if not sdf.exists() or sdf.stat().st_size == 0:
        return "dock_failed", 0
    n = len(load_poses(sdf))
    if n < 1:
        return "dock_failed", 0
    return "dock_completed", n


def eval_pair(ligand: str, seed: int, u_sdf: Path, n_sdf: Path, u_status: Path | None, n_status: Path | None, source: str) -> dict:
    u_state, u_n = dock_state(u_sdf, u_status)
    n_state, n_n = dock_state(n_sdf, n_status)
    row = {
        "ligand_id": ligand,
        "seed": seed,
        "pose_source": source,
        "urat1_dock_state": u_state,
        "urat1_n_poses": u_n,
        "nlrp3_dock_state": n_state,
        "nlrp3_n_poses": n_n,
        "urat1_geometry_state": None,
        "nlrp3_geometry_state": None,
        "urat1_pass": False,
        "nlrp3_structural_pass": False,
        "dual_seed_pass": False,
    }
    if u_state == "dock_completed":
        ev = evaluate_urat1_acid_sdf(u_sdf, ARG, REF_U, ligand, seed, rule="a2")
        row.update({
            "urat1_geometry_state": ev.get("pose_selection_status"),
            "urat1_selected_mode": ev.get("selected_mode"),
            "urat1_n_geometry_pass_modes": ev.get("n_geometry_pass_modes"),
            "urat1_CNNscore": ev.get("CNNscore"),
            "urat1_CNNaffinity": ev.get("CNNaffinity"),
            "urat1_acid_arg477_min_A": ev.get("acid_arg477_min_A"),
            "urat1_centroid_A": ev.get("centroid_to_crystal_lig_A"),
            "urat1_pass": bool(ev.get("keep_urat1_acid")),
        })
        if not row["urat1_pass"]:
            row["urat1_geometry_state"] = ev.get("pose_selection_status") or "no_geometry_pass"
    else:
        row["urat1_geometry_state"] = u_state
    if n_state == "dock_completed":
        evn = evaluate_nlrp3_structural(n_sdf, ligand, seed)
        row.update({
            "nlrp3_selected_mode": evn.get("selected_mode"),
            "nlrp3_CNNscore": evn.get("CNNscore"),
            "nlrp3_CNNaffinity": evn.get("CNNaffinity"),
            "nlrp3_centroid_A": evn.get("centroid_to_crystal_lig_A"),
            "nlrp3_overlap": evn.get("pocket_overlap_frac"),
            "nlrp3_ifp": evn.get("ifp_jaccard_vs_np3146"),
            "nlrp3_n_key_contacts": evn.get("n_key_contacts"),
            "nlrp3_clash": "CLASH:ANY" in (evn.get("ifp_bits") or ""),
            "nlrp3_structural_pass": bool(evn.get("keep_nlrp3_structural")),
            "nlrp3_loose_pass": bool(evn.get("keep_nlrp3_pose")),
        })
        row["nlrp3_geometry_state"] = "structural_pass" if row["nlrp3_structural_pass"] else "no_structural_pass"
    else:
        row["nlrp3_geometry_state"] = n_state
    row["dual_seed_pass"] = bool(row["urat1_pass"] and row["nlrp3_structural_pass"])
    return row


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    decisions = pd.read_csv(ROOT / "results/nlrp3_model_v2/stage3a/new_molecule_docking_decisions.csv")
    lib = pd.read_csv(ROOT / "results/nlrp3_model_v2/stage2/clinical_library_scores.csv")
    dock_ids = decisions.loc[decisions.docking_decision == "dock", "repurposing_id"].tolist()
    exclude = decisions.loc[decisions.docking_decision != "dock"].copy()
    prep = pd.read_csv(NEW / "ligand_prep" / "ligand_prep_report.csv")

    # job execution table
    jobs = json.loads((NEW / "job_list.json").read_text())
    job_rows = []
    for job in jobs:
        status = Path(job["out"]) / "status.json"
        rec = json.loads(status.read_text()) if status.exists() else {"state": "not_started"}
        job_rows.append({
            "kind": job["kind"], "ligand": job["ligand"], "target": job["target"], "seed": job["seed"],
            "state": rec.get("state"), "n_poses": rec.get("n_poses"), "scores_readable": rec.get("scores_readable"),
            "exit_code": rec.get("exit_code"), "elapsed_s": rec.get("elapsed_s"), "error": rec.get("error"),
            "skipped": rec.get("skipped", False),
        })
    job_df = pd.DataFrame(job_rows)
    job_df.to_csv(OUT / "task_execution.csv", index=False)

    seed_rows = []
    for lig in dock_ids:
        for seed in SEEDS:
            u = new_sdf(lig, "urat1_9dkb", seed)
            n = new_sdf(lig, "nlrp3_7alv", seed)
            seed_rows.append(eval_pair(
                lig, seed, u, n,
                u.parent / "status.json", n.parent / "status.json",
                "stage3b_new32",
            ))
    seed_df = pd.DataFrame(seed_rows)
    seed_df.to_csv(OUT / "new32_per_seed.csv", index=False)

    chem = decisions.merge(
        lib[["repurposing_id", "chembl_id", "max_phase", "canonical_smiles"]],
        on="repurposing_id", how="left",
    )
    summaries = []
    for lig in dock_ids:
        sub = seed_df[seed_df.ligand_id == lig]
        meta = chem[chem.repurposing_id == lig].iloc[0]
        pr = prep[prep.repurposing_id == lig].iloc[0]
        n_dual = int(sub.dual_seed_pass.sum())
        structural = n_dual >= 2
        beta = bool(meta.beta_lactam)
        clean = bool(meta.alert_clean_pains_brenk)
        chem_ok = structural and (not beta) and clean
        if structural and not chem_ok:
            chem_state = "structural_pass_chemistry_fail"
        elif structural and chem_ok:
            chem_state = "structural_and_chemistry_pass"
        else:
            chem_state = "not_structural_pass"
        summaries.append({
            "ligand_id": lig,
            "name": meta["name"],
            "chembl_id": meta.chembl_id,
            "max_phase": meta.max_phase,
            "cohort": "stage3b_new32",
            "n_dual_seeds": n_dual,
            "structural_pass": structural,
            "beta_lactam": beta,
            "alert_clean_pains_brenk": clean,
            "chemistry_nomination_pass": chem_ok,
            "nomination_state": chem_state,
            "formal_charge": pr.formal_charge,
            "selected_smiles": pr.selected_smiles,
            "n_unassigned_stereo": pr.selected_n_unassigned_stereo,
            "prep_status": pr.status,
        })
    sum_df = pd.DataFrame(summaries)
    sum_df.to_csv(OUT / "new32_summary.csv", index=False)

    # historical deposited cohort, same functions; gap fills marked
    deposited = sorted(p.name.replace("_out.sdf", "") for p in (HIST / "acid_dual/urat1_9dkb/seed42").glob("*_out.sdf"))
    hist_rows = []
    for lig in deposited:
        for seed in SEEDS:
            u, su = hist_sdf(lig, "urat1_9dkb", seed)
            n, sn = hist_sdf(lig, "nlrp3_7alv", seed)
            source = "stage3b_gap_fill" if "gap" in (su, sn) else "historical_deposited"
            # gap URAT1 stays historical; only NLRP3 43/44 is new
            row = eval_pair(lig, seed, u, n, None, None, source)
            if lig == "REP_07837" and seed in (43, 44):
                row["nlrp3_pose_origin"] = "stage3b_gap_fill"
                row["urat1_pose_origin"] = "historical_deposited"
            else:
                row["nlrp3_pose_origin"] = "historical_deposited"
                row["urat1_pose_origin"] = "historical_deposited"
            hist_rows.append(row)
    hist_seed = pd.DataFrame(hist_rows)
    hist_seed.to_csv(OUT / "historical_per_seed_rescored.csv", index=False)

    frozen = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/reportable_candidates.csv")
    primary = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/primary_candidates.csv")
    primary_ids = set(primary.ligand_id)
    hist_sum = []
    for lig, sub in hist_seed.groupby("ligand_id"):
        n_dual = int(sub.dual_seed_pass.sum())
        fr = frozen[frozen.ligand_id == lig]
        name = fr.iloc[0]["name"] if len(fr) else ""
        pains = bool(fr.iloc[0].pains_any) if len(fr) else None
        brenk = bool(fr.iloc[0].brenk) if len(fr) else None
        beta = bool(fr.iloc[0].beta_lactam_flag) if len(fr) else None
        structural = n_dual >= 2
        if pains is None:
            chem_ok = None
            chem_state = "structural_pass_chemistry_not_in_frozen_table" if structural else "not_structural_pass"
        else:
            chem_ok = bool(structural and (not beta) and (not pains) and (not brenk))
            if structural and not chem_ok:
                chem_state = "structural_pass_chemistry_fail"
            elif structural and chem_ok:
                chem_state = "structural_and_chemistry_pass"
            else:
                chem_state = "not_structural_pass"
        hist_sum.append({
            "ligand_id": lig,
            "name": name,
            "cohort": "historical_deposited",
            "in_frozen_primary12": lig in primary_ids,
            "in_frozen_reportable": lig in set(frozen.ligand_id),
            "n_dual_seeds": n_dual,
            "structural_pass": structural,
            "beta_lactam": beta,
            "pains_any": pains,
            "brenk": brenk,
            "chemistry_nomination_pass": chem_ok,
            "nomination_state": chem_state,
            "frozen_reportable_flag": bool(fr.iloc[0].reportable) if len(fr) else None,
        })
    hist_sum_df = pd.DataFrame(hist_sum)
    hist_sum_df.to_csv(OUT / "historical_summary_rescored.csv", index=False)

    # PF specific
    pf = hist_seed[hist_seed.ligand_id == "REP_07580"].copy()
    pf.to_csv(OUT / "PF03882845_REP_07580_per_seed.csv", index=False)

    exclude_rows = exclude.merge(lib[["repurposing_id", "name", "chembl_id"]], on="repurposing_id", how="left", suffixes=("", "_lib"))
    exclude_rows = exclude_rows.assign(
        cohort="stage3a_exclude",
        disposition="rule_not_applicable_noncarboxylate",
        note="acyl-sulfonamide or tetrazole; carboxylate-Arg477 gate was not applied and is not evidence of absent URAT1 activity",
    )
    exclude_rows.to_csv(OUT / "excluded_noncarboxylate.csv", index=False)

    unified = pd.concat([
        sum_df.assign(table_arm="new32")[[
            "ligand_id", "name", "cohort", "n_dual_seeds", "structural_pass",
            "chemistry_nomination_pass", "nomination_state",
        ]],
        hist_sum_df.loc[hist_sum_df.structural_pass | hist_sum_df.in_frozen_reportable, [
            "ligand_id", "name", "cohort", "n_dual_seeds", "structural_pass",
            "chemistry_nomination_pass", "nomination_state",
        ]],
    ], ignore_index=True)
    unified.to_csv(OUT / "unified_structural_candidates.csv", index=False)

    counts = {
        "new32_prepared": int((prep.status == "prepared").sum()),
        "new32_prep_failed": int((prep.status != "prepared").sum()),
        "jobs_total": int(len(job_df)),
        "jobs_by_state": job_df.state.value_counts(dropna=False).to_dict(),
        "jobs_with_fewer_than_9_poses": int(((job_df.state == "completed") & (job_df.n_poses < 9)).sum()),
        "new32_structural_pass": int(sum_df.structural_pass.sum()),
        "new32_chemistry_pass": int(sum_df.chemistry_nomination_pass.sum()),
        "new32_structural_but_chemistry_fail": int((sum_df.nomination_state == "structural_pass_chemistry_fail").sum()),
        "historical_structural_pass": int(hist_sum_df.structural_pass.sum()),
        "historical_chemistry_pass": int(hist_sum_df.chemistry_nomination_pass.fillna(False).sum()),
        "pf_n_dual_seeds": int(hist_sum_df.loc[hist_sum_df.ligand_id == "REP_07580", "n_dual_seeds"].iloc[0]),
        "pf_structural_pass": bool(hist_sum_df.loc[hist_sum_df.ligand_id == "REP_07580", "structural_pass"].iloc[0]),
        "exclude_noncarboxylate": int(len(exclude)),
        "rep_07837_gap_states": job_df.loc[job_df.kind == "gap_REP_07837", ["seed", "state", "n_poses"]].to_dict("records"),
    }
    (OUT / "counts.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
