#!/usr/bin/env python3
"""Analyze exhaustiveness sensitivity (Exp A/B). Must run in conda env cadd_tools."""
from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0")
OUT = ROOT / "analysis" / "exhaustiveness_sensitivity_v1"
TABLES = OUT / "tables"
POSES = OUT / "poses"
LIG_SDF = Path("/mnt/d/CADD paper exercise/dual target docking/Maestro doc/vina_docking/ligands_sdf/EH40_01.sdf")
FIXED_SEED = 20260727
TARGETS = ["3POZ", "3RCD"]
TMPL = Chem.RemoveHs(Chem.SDMolSupplier(str(LIG_SDF), removeHs=False)[0])


def parse_vina_score(path: Path) -> float | None:
    if not path.exists():
        return None
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def crystal_mol(pdb_path: Path) -> Chem.Mol:
    lines = [
        line
        for line in pdb_path.read_text(errors="ignore").splitlines()
        if line.startswith(("HETATM", "ATOM  ")) and line[17:20].strip() == "03P"
    ]
    if not lines:
        lines = [
            line
            for line in pdb_path.read_text(errors="ignore").splitlines()
            if line.startswith(("HETATM", "ATOM  "))
        ]
    mol = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=True, sanitize=False)
    mol = AllChem.AssignBondOrdersFromTemplate(TMPL, mol)
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def pose_mol(pdbqt_path: Path) -> Chem.Mol:
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tf:
        out = tf.name
    subprocess.check_call(
        ["obabel", "-ipdbqt", str(pdbqt_path), "-osdf", "-O", out],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    mol = Chem.SDMolSupplier(out, removeHs=False, sanitize=False)[0]
    mol = AllChem.AssignBondOrdersFromTemplate(TMPL, Chem.RemoveHs(mol))
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def symmetry_corrected_rmsd(ref: Chem.Mol, pose: Chem.Mol) -> float:
    """Heavy-atom RMSD with template-constrained symmetry matches.

    Matches historical panel40 redock numbers (as-run 3POZ mode1 ~9.51, mode2 ~1.02).
    Uses min CalcRMS over RDKit substructure matches to EH40_01 template.
    """
    matches_ref = ref.GetSubstructMatches(TMPL, uniquify=False)
    matches_pose = pose.GetSubstructMatches(TMPL, uniquify=False)
    if not matches_ref or not matches_pose:
        raise RuntimeError("no template substructure match for RMSD")
    best = None
    for a in matches_ref:
        for b in matches_pose:
            rms = float(rdMolAlign.CalcRMS(pose, ref, map=[list(zip(b, a))]))
            if best is None or rms < best:
                best = rms
    return float(best)


REF = {
    "3POZ": crystal_mol(TABLES / "3POZ_cocrystal_03P.pdb"),
    "3RCD": crystal_mol(TABLES / "3RCD_cocrystal_03P.pdb"),
}


def pose_path(exhaustiveness: int, seed: int, target: str, lig: str, mode: int) -> Path:
    return POSES / f"E{exhaustiveness}_seed{seed}" / target / lig / f"mode_{mode:02d}.pdbqt"


def subset() -> list[str]:
    with (TABLES / "subset_ligands.csv").open() as fh:
        return [r["ligand_id"] for r in csv.DictReader(fh)]


def build_vina_long() -> pd.DataFrame:
    rows = []
    for exhaustiveness in [8, 16, 32]:
        for target in TARGETS:
            for lig in subset():
                for mode in range(1, 10):
                    p = pose_path(exhaustiveness, FIXED_SEED, target, lig, mode)
                    rows.append(
                        {
                            "experiment": "A",
                            "target": target,
                            "ligand_id": lig,
                            "exhaustiveness": exhaustiveness,
                            "seed": FIXED_SEED,
                            "mode": mode,
                            "vina_score": parse_vina_score(p),
                            "pose_path": str(p),
                        }
                    )
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "scores_vina_experimentA_long.csv", index=False)
    return df


def eh40_01_rmsd_table() -> pd.DataFrame:
    rtm = None
    rtm_path = TABLES / "scores_rtm_experimentA.csv"
    if rtm_path.exists():
        rtm = pd.read_csv(rtm_path)
    rows = []
    for exhaustiveness in [8, 16, 32]:
        for target in TARGETS:
            rmsds = []
            for mode in range(1, 10):
                pose = pose_mol(pose_path(exhaustiveness, FIXED_SEED, target, "EH40_01", mode))
                rms = symmetry_corrected_rmsd(REF[target], pose)
                rmsds.append((mode, rms))
            mode1 = next(r for m, r in rmsds if m == 1)
            best_mode, best_rmsd = min(rmsds, key=lambda x: x[1])
            rec = {
                "target": target,
                "exhaustiveness": exhaustiveness,
                "seed": FIXED_SEED,
                "rmsd_mode1": round(mode1, 3),
                "rmsd_best_of_9": round(best_rmsd, 3),
                "best_of_9_mode": best_mode,
                "pass_mode1_lt2": mode1 < 2.0,
                "pass_best_of_9_lt2": best_rmsd < 2.0,
                "rmsd_rtm_best": None,
                "rtm_best_mode": None,
            }
            if rtm is not None:
                hit = rtm[
                    (rtm["target"] == target)
                    & (rtm["ligand_id"] == "EH40_01")
                    & (rtm["exhaustiveness"] == exhaustiveness)
                ]
                if not hit.empty:
                    rtm_mode = int(hit.iloc[0]["best_rtm_mode"])
                    rtm_pose = pose_mol(pose_path(exhaustiveness, FIXED_SEED, target, "EH40_01", rtm_mode))
                    rec["rtm_best_mode"] = rtm_mode
                    rec["rmsd_rtm_best"] = round(symmetry_corrected_rmsd(REF[target], rtm_pose), 3)
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "eh40_01_rmsd_by_E.csv", index=False)
    return df


def ranks_experimentA(vina_long: pd.DataFrame) -> pd.DataFrame:
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_40.csv").rename(columns={"panel_id": "ligand_id"})
    rows = []
    for exhaustiveness in [8, 16, 32]:
        mode1 = vina_long[(vina_long["exhaustiveness"] == exhaustiveness) & (vina_long["mode"] == 1)]
        wide = (
            mode1.pivot(index="ligand_id", columns="target", values="vina_score")
            .rename(columns={"3POZ": "vina_3POZ_mode1", "3RCD": "vina_3RCD_mode1"})
            .reset_index()
        )
        wide["vina_mean_mode1_hb"] = -(wide["vina_3POZ_mode1"] + wide["vina_3RCD_mode1"]) / 2.0
        wide["exhaustiveness"] = exhaustiveness
        rows.append(wide)
    ranks = pd.concat(rows, ignore_index=True).merge(
        panel[["ligand_id", "class", "pref_name"]], on="ligand_id", how="left"
    )
    out = []
    for _, sub in ranks.groupby("exhaustiveness"):
        order = sub.sort_values("vina_mean_mode1_hb", ascending=False).reset_index(drop=True)
        order["vina_mode1_rank"] = np.arange(1, len(order) + 1)
        out.append(order)
    ranks = pd.concat(out, ignore_index=True)
    rtm_path = TABLES / "scores_rtm_experimentA.csv"
    if rtm_path.exists():
        rtm = pd.read_csv(rtm_path)
        rtm_rows = []
        for exhaustiveness in [8, 16, 32]:
            sub = rtm[rtm["exhaustiveness"] == exhaustiveness]
            wide = (
                sub.pivot(index="ligand_id", columns="target", values="rtmscore")
                .rename(columns={"3POZ": "rtm_3POZ", "3RCD": "rtm_3RCD"})
                .reset_index()
            )
            wide["rtm_min"] = wide[["rtm_3POZ", "rtm_3RCD"]].min(axis=1)
            wide["exhaustiveness"] = exhaustiveness
            rtm_rows.append(wide)
        ranks = ranks.merge(pd.concat(rtm_rows, ignore_index=True), on=["ligand_id", "exhaustiveness"], how="left")
        out = []
        for _, sub in ranks.groupby("exhaustiveness"):
            order = sub.sort_values("rtm_min", ascending=False).reset_index(drop=True)
            order["rtm_best_rank"] = np.arange(1, len(order) + 1)
            out.append(order)
        ranks = pd.concat(out, ignore_index=True)
    ranks = ranks.sort_values(["exhaustiveness", "vina_mode1_rank", "ligand_id"])
    ranks.to_csv(TABLES / "ranks_experimentA.csv", index=False)
    return ranks


def seed_noise() -> pd.DataFrame:
    rows = []
    for exhaustiveness in [8, 16]:
        for target in TARGETS:
            vals = {}
            for seed in [20260727, 7, 42]:
                best = None
                for mode in range(1, 10):
                    pose = pose_mol(pose_path(exhaustiveness, seed, target, "EH40_01", mode))
                    rms = symmetry_corrected_rmsd(REF[target], pose)
                    if best is None or rms < best:
                        best = rms
                vals[seed] = round(best, 3)
            arr = list(vals.values())
            rows.append(
                {
                    "metric": "EH40_01_best_of_9_rmsd",
                    "target": target,
                    "exhaustiveness": exhaustiveness,
                    "min_value": min(arr),
                    "max_value": max(arr),
                    "range_value": round(max(arr) - min(arr), 3),
                    "values_by_seed": json.dumps(vals),
                }
            )
        # relative ranks among {01,18,23} by vina mode1 mean
        for lig in ["EH40_18", "EH40_23"]:
            ranks_by_seed = {}
            for seed in [20260727, 7, 42]:
                scores = []
                for lid in ["EH40_01", "EH40_18", "EH40_23"]:
                    poz = parse_vina_score(pose_path(exhaustiveness, seed, "3POZ", lid, 1))
                    rcd = parse_vina_score(pose_path(exhaustiveness, seed, "3RCD", lid, 1))
                    scores.append((lid, -((poz + rcd) / 2.0)))
                order = sorted(scores, key=lambda x: x[1], reverse=True)
                ranks_by_seed[seed] = {lid: i + 1 for i, (lid, _) in enumerate(order)}[lig]
            arr = list(ranks_by_seed.values())
            rows.append(
                {
                    "metric": f"{lig}_rank_among_01_18_23_vina_mode1",
                    "target": "dual_aggregate_key3",
                    "exhaustiveness": exhaustiveness,
                    "min_value": min(arr),
                    "max_value": max(arr),
                    "range_value": max(arr) - min(arr),
                    "values_by_seed": json.dumps(ranks_by_seed),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "seed_noise_E8_E16.csv", index=False)
    return df


def asrun_compare(ranks: pd.DataFrame, rmsd_df: pd.DataFrame) -> pd.DataFrame:
    vina = pd.read_csv(ROOT / "tables" / "scores_vina.csv").rename(columns={"ligand": "ligand_id"})
    sub = subset()
    asrun = vina[vina["ligand_id"].isin(sub)][["ligand_id", "3POZ_affinity", "3RCD_affinity"]].copy()
    asrun["vina_mean_mode1_hb_asrun"] = -(asrun["3POZ_affinity"] + asrun["3RCD_affinity"]) / 2.0
    asrun = asrun.sort_values("vina_mean_mode1_hb_asrun", ascending=False).reset_index(drop=True)
    asrun["vina_mode1_rank_asrun_subset8"] = np.arange(1, len(asrun) + 1)
    fixed = ranks[ranks["exhaustiveness"] == 8][["ligand_id", "vina_mode1_rank"]].rename(
        columns={"vina_mode1_rank": "vina_mode1_rank_fixedseed_E8_subset8"}
    )
    out = asrun.merge(fixed, on="ligand_id", how="left")
    out.to_csv(TABLES / "asrun_vs_fixedseed_E8.csv", index=False)

    # EH40_01 RMSD side table: as-run vs fixed
    rows = []
    for target in TARGETS:
        asrun_rms = []
        for mode in range(1, 10):
            pose = pose_mol(ROOT / "poses" / target / "EH40_01" / f"mode_{mode:02d}.pdbqt")
            asrun_rms.append((mode, symmetry_corrected_rmsd(REF[target], pose)))
        mode1 = next(r for m, r in asrun_rms if m == 1)
        best = min(asrun_rms, key=lambda x: x[1])
        fix = rmsd_df[(rmsd_df["target"] == target) & (rmsd_df["exhaustiveness"] == 8)].iloc[0]
        rows.append(
            {
                "target": target,
                "asrun_rmsd_mode1": round(mode1, 3),
                "asrun_rmsd_best_of_9": round(best[1], 3),
                "asrun_best_mode": best[0],
                "fixed_E8_rmsd_mode1": fix["rmsd_mode1"],
                "fixed_E8_rmsd_best_of_9": fix["rmsd_best_of_9"],
                "fixed_E8_best_mode": fix["best_of_9_mode"],
                "note": "as-run uses per-job random seeds; fixed uses seed=20260727",
            }
        )
    pd.DataFrame(rows).to_csv(TABLES / "asrun_vs_fixedseed_E8_eh40_01_rmsd.csv", index=False)
    return out


def main() -> None:
    print("Building Exp A vina long...", flush=True)
    vina_long = build_vina_long()
    print("Computing EH40_01 RMSD by E...", flush=True)
    rmsd_df = eh40_01_rmsd_table()
    print(rmsd_df.to_string(index=False), flush=True)
    print("Building ranks...", flush=True)
    ranks = ranks_experimentA(vina_long)
    print("Seed noise...", flush=True)
    noise = seed_noise()
    print(noise.to_string(index=False), flush=True)
    print("As-run compare...", flush=True)
    asrun_compare(ranks, rmsd_df)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
